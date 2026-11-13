/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

    This class is modified from OpenFOAM's source code
    applications/solvers/compressible/rhoSimpleFoam

    OpenFOAM: The Open Source CFD Toolbox

    Copyright (C): 2011-2016 OpenFOAM Foundation

    OpenFOAM License:

        OpenFOAM is free software: you can redistribute it and/or modify it
        under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
    
        OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
        ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
        FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
        for more details.
    
        You should have received a copy of the GNU General Public License
        along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DAHisaFoam.H"
#include "wallDist.H"
#include "characteristicBase.H"
#include "fluidThermo.H"
#include "transform.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAHisaFoam, 0);
addToRunTimeSelectionTable(DASolver, DAHisaFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAHisaFoam::DAHisaFoam(
    char* argsAll,
    PyObject* pyOptions)
    : DASolver(argsAll, pyOptions),
      hisaSolverPtr_(nullptr)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAHisaFoam::initSolver()
{
    /*
    Description:
        Initialize variables for DASolver
    */

    Info << "Initializing fields for DAHisaFoam" << endl;
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();

    dictionary dict;
    dict.add("solver", "hisa");
    hisaSolverPtr_.reset(
        solverModule::New(Foam::polyMesh::defaultRegion, runTime, mesh, dict));

    solverModule& solver = hisaSolverPtr_();

    solver.initialise();

    // v2506 wallDist owns an unregistered yWall field. Expose that same
    // live field to the DAFoam turbulence adapter, without copying/fixing it.
    // Its mesh-motion updates and active coordinate derivatives are retained.
    if (!mesh.foundObject<volScalarField>("yWall"))
    {
        volScalarField& distance =
            const_cast<volScalarField&>(wallDist::New(mesh).y());
        if (!distance.checkIn())
            FatalErrorInFunction << "Cannot register wallDist field for DAHisaFoam"
                                 << exit(FatalError);
    }


    // read the RAS model from constant/turbulenceProperties
    const word turbModelName(
        IOdictionary(
            IOobject(
                "turbulenceProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false))
            .subDict("RAS")
            .lookup("RASModel"));
    daTurbulenceModelPtr_.reset(DATurbulenceModel::New(turbModelName, mesh, daOptionPtr_()));

#include "createAdjoint.H"
    haveFreestream_ = readFreestream(previousFreestream_, previousFreestreamP_, previousFreestreamT_);
}

// Sample actual restart boundary references, not the uniform 0/ initial field.
// All ranks participate, including ranks with no farfield faces.
bool DAHisaFoam::readFreestream(vector& velocity, scalar& pressure, scalar& temperature)
{
    const auto& U = meshPtr_->lookupObject<volVectorField>("U");
    vector sumU = vector::zero;
    scalar sumP = 0, sumT = 0, count = 0;
    forAll(U.boundaryField(), patchI)
    {
        const auto& patch = U.boundaryField()[patchI];
        if (patch.type() != "characteristicFarfieldVelocity" || !patch.size()) continue;
        const auto& ref = refCast<const characteristicBase>(patch);
        scalar n = patch.size();
        sumU += n * ref.URef(); sumP += n * ref.pRef(); sumT += n * ref.TRef(); count += n;
    }
    reduce(sumU, sumOp<vector>()); reduce(sumP, sumOp<scalar>());
    reduce(sumT, sumOp<scalar>()); reduce(count, sumOp<scalar>());
    if (!count) return false;
    velocity = sumU/count; pressure = sumP/count; temperature = sumT/count;
    scalar spread = 0;
    forAll(U.boundaryField(), patchI)
    {
        const auto& patch = U.boundaryField()[patchI];
        if (patch.type() != "characteristicFarfieldVelocity" || !patch.size()) continue;
        const auto& ref = refCast<const characteristicBase>(patch);
        spread = max(spread, mag(ref.URef()-velocity));
    }
    reduce(spread, maxOp<scalar>());
    if (spread > 1e-10*max(mag(velocity),scalar(1)))
        FatalErrorInFunction << "Freestream correction requires a uniform farfield velocity" << exit(FatalError);
    return true;
}

void DAHisaFoam::correctFreestreamInitialGuess()
{
    freestreamCorrected_ = false;
#if !defined(CODI_ADR) && !defined(CODI_ADF)
    const auto& options = daOptionPtr_->getAllOptions();
    const bool enabled = options.get<bool>("infChangeCorrection");
    const scalar tolerance = options.get<scalar>("infChangeCorrectionTol");
    const word method = options.get<word>("infChangeCorrectionType");
    if (!std::isfinite(tolerance) || tolerance < 0 || (method != "offset" && method != "rotate"))
        FatalErrorInFunction << "Invalid infChangeCorrectionTol/Type (offset or rotate)" << exit(FatalError);
    vector next; scalar p, T;
    if (!readFreestream(next,p,T)) return;
    if (enabled && haveFreestream_)
    {
        // Current support is patchVelocity (speed/AoA) at fixed p,T.
        // Do not silently claim a general thermodynamic-state correction.
        if (mag(p-previousFreestreamP_) > 1e-10*max(p,scalar(1)) ||
            mag(T-previousFreestreamT_) > 1e-10*max(T,scalar(1)))
            FatalErrorInFunction << "infChangeCorrection currently supports velocity changes at fixed farfield p,T. Disable it for thermodynamic changes." << exit(FatalError);
        auto& thermo = meshPtr_->lookupObject<fluidThermo>("thermophysicalProperties");
        // ADflow uRef=sqrt(pInf/rhoInf)=sqrt(R*TInf), for the ideal gas.
        tmp<volScalarField> gas = thermo.Cp()-thermo.Cv();
        scalar gasR = gAverage(gas().primitiveField());
        scalar scale2 = gasR*T;
        if (scale2 <= 0) FatalErrorInFunction << "Invalid freestream velocity scale" << exit(FatalError);
        vector delta = next-previousFreestream_;
        scalar deltaEnergy = 0.5*(magSqr(next)-magSqr(previousFreestream_))/scale2;
        scalar change = sqrt(magSqr(delta)/scale2 + sqr(deltaEnergy));
        if (change > 0 && change >= tolerance)
        {
            auto& U = meshPtr_->lookupObjectRef<volVectorField>("U");
            const vectorField oldVelocity(U.primitiveField());
            if (method == "offset") U.primitiveFieldRef() += delta;
            else
            {
                scalar oldSpeed=mag(previousFreestream_), newSpeed=mag(next);
                if (oldSpeed <= SMALL || newSpeed <= SMALL)
                    FatalErrorInFunction << "rotate requires nonzero old/new speed; use offset" << exit(FatalError);
                tensor rotation = rotationTensor(previousFreestream_/oldSpeed,next/newSpeed);
                forAll(U,cellI) U[cellI] = (newSpeed/oldSpeed)*(rotation & U[cellI]);
            }
            // ADflow stores rho,u,v,w,rhoE. Preserve density and add the
            // freestream total-energy increment; convert back to p,T.
            auto& temperature = meshPtr_->lookupObjectRef<volScalarField>("T");
            auto& pressure = meshPtr_->lookupObjectRef<volScalarField>("p");
            tmp<volScalarField> cv = thermo.Cv();
            forAll(U,cellI)
            {
                scalar oldT = temperature[cellI];
                scalar rhoInfOverRho = (p/T)*(oldT/pressure[cellI]);
                scalar dEnergy = rhoInfOverRho*deltaEnergy*scale2
                    - 0.5*(magSqr(U[cellI])-magSqr(oldVelocity[cellI]));
                scalar newT = oldT + dEnergy/cv()[cellI];
                if (newT <= 0 || !std::isfinite(newT))
                    FatalErrorInFunction << "Freestream correction gives nonpositive T; disable it or reduce the condition change" << exit(FatalError);
                temperature[cellI] = newT;
                pressure[cellI] *= newT/oldT;
            }
            U.correctBoundaryConditions();
            temperature.correctBoundaryConditions(); pressure.correctBoundaryConditions();
            freestreamCorrected_ = true;
            Info << "Freestream initial-guess correction: " << method << " deltaU " << delta
                 << " nondimensionalChange " << change << " tolerance " << tolerance << endl;
        }
    }
    previousFreestream_=next; previousFreestreamP_=p; previousFreestreamT_=T; haveFreestream_=true;
#endif
}

label DAHisaFoam::solvePrimal()
{
    /*
    Description:
        Call the primal solver to get converged state variables
    */

    correctFreestreamInitialGuess();

    if (meshPtr_->solutionDict().subOrEmptyDict("hisaHybrid").lookupOrDefault<bool>("enabled", false))
        return solveHybrid();

    Info << "Running HiSA " << endl;

    solverModule& solver = hisaSolverPtr_();
    Time& runTime = runTimePtr_();
    runTime.setTime(0.0, 0);

    bool steadyState = solver.steadyState();
    // Opt-in, HiSA-only diagnostics. Shared turbulence models are unchanged.
    const dictionary diagnostics = meshPtr_->solutionDict().subOrEmptyDict("hisaDiagnostics");
    const bool freezeSA = diagnostics.lookupOrDefault<bool>("freezeSA", false);
    const bool auditSA = diagnostics.lookupOrDefault<bool>("auditSA", false);
    if ((freezeSA || auditSA) && (!steadyState || !meshPtr_->foundObject<volScalarField>("nuTilda")))
        FatalErrorInFunction << "SA diagnostics require steady HiSA with nuTilda" << exit(FatalError);


    if (!steadyState)
    {
#include "createTimeControls.H"

        scalar scale = solver.timeStepScaling(maxCo);
        //Generate a pretend Courant number so we can use setInitialDeltaT.H unmodified
        scalar CoNum = maxCo / scale;
#include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (steadyState)
    {
        Info << "\nStarting pseudotime iteration loop\n"
             << endl;
    }
    else
    {
        Info << "\nStarting time loop\n"
             << endl;
    }

    //scalar prevCpuTime = runTime.elapsedCpuTime();
    while (runTime.run())
    {

        if (!steadyState)
        {
#include "createTimeControls.H"

            // Scale deltaT
            scalar scale = solver.timeStepScaling(maxCo);
            if (adjustTimeStep)
            {
                scalar deltaT = scale * runTime.deltaTValue();
                deltaT = min(
                    min(deltaT, runTime.deltaTValue() + 0.1 * deltaT),
                    1.2 * runTime.deltaTValue());
                runTime.setDeltaT(min(deltaT, maxDeltaT));

                Info << "deltaT = " << runTime.deltaTValue() << endl;
            }
        }

        solver.preTimeIncrement();

        runTime++;

        if (steadyState)
        {
            printToScreen_ = this->isPrintTime(runTime, printInterval_);
        }
        else
        {
            printToScreen_ = 1;
        }

        if (printToScreen_)
        {
            Info << "Time = " << runTime.timeName() << nl << endl;
            // this->calcPrimalResidualStatistics("print");
        }

        solver.beginTimeStep();

        // --- Outer corrector loop
        bool notFinished;

        while ((notFinished = solver.solnControl().loop()))
        {
            solver.outerIteration();

            // NOTE: we removed the turbulence correct in hisaModule, so we need to call it here
            // NOTE: this trick works for ESI version of OpenFOAM only!
            if (!freezeSA) daTurbulenceModelPtr_->correct(printToScreen_);
            if (auditSA && runTime.value() >= runTime.endTime().value())
            {
                dictionary opts;
                opts.set("isPC", 0);
                daTurbulenceModelPtr_->calcResiduals(opts);
                const scalarField lagged(meshPtr_->lookupObject<volScalarField>("nuTildaRes").primitiveField());
                const scalarField oldPhi(meshPtr_->lookupObject<surfaceScalarField>("phi").primitiveField());
                // Recompute flux at the updated flow state, without a new SA solve.
                daResidualPtr_->calcResiduals(opts);
                daTurbulenceModelPtr_->calcResiduals(opts);
                const scalarField refreshed(meshPtr_->lookupObject<volScalarField>("nuTildaRes").primitiveField());
                const scalarField newPhi(meshPtr_->lookupObject<surfaceScalarField>("phi").primitiveField());
                this->updateStateBoundaryConditions();
                this->calcResiduals();
                const scalarField full(meshPtr_->lookupObject<volScalarField>("nuTildaRes").primitiveField());
                Info << "SAAudit: laggedL2 " << sqrt(gSum(sqr(lagged)))
                     << " refreshedL2 " << sqrt(gSum(sqr(refreshed)))
                     << " fullL2 " << sqrt(gSum(sqr(full)))
                     << " refreshDifferenceL2 " << sqrt(gSum(sqr(refreshed-lagged)))
                     << " boundaryDifferenceL2 " << sqrt(gSum(sqr(full-refreshed)))
                     << " phiDifferenceL2 " << sqrt(gSum(sqr(newPhi-oldPhi))) << endl;
            }

            if (steadyState)
            {
                break;
            }
        }

        // calculate all functions
        this->calcAllFunctions(printToScreen_);
        // calculate yPlus
        daTurbulenceModelPtr_->printYPlus(printToScreen_);

        // print run time
        this->printElapsedTime(runTime, printToScreen_);

        // if we want to use the function std as the convergence criteria, we need to compute the std
        if (primalFuncStdNames_.size() > 0)
        {
            this->calcFuncStd();
            // also compute the function slope over the same window (diagnostic only)
            this->calcFuncSlope();
        }
        if (this->isFuncStdConverged() && runTime.timeIndex() > primalMinIters_)
        {

            Info << "All functions satisfied their prescribed std and slope tolerances" << endl
                 << endl;
            notFinished = false;
        }

        if (steadyState && !notFinished)
        {
            runTime.writeNow();
            break;
        }
        else
        {
            runTime.write();
        }
    }

    // write the mesh to files
    meshPtr_->write();

    primalFinalTimeIndex_ = runTime.timeIndex();

    Info << "End\n"
         << endl;

    return 0;
}

} // End namespace Foam

// ************************************************************************* //
