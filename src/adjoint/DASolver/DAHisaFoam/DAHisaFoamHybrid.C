#include "DAHisaFoam.H"
#include "HisaPetscNK.H"
#include "IOdictionary.H"
#include "uniformDimensionedFields.H"
#include <iostream>
#include <streambuf>
#include <memory>
namespace {
struct NullOutput : std::streambuf {int overflow(int c) override {return c;}};
struct Quiet {
    NullOutput sink; std::streambuf* old;
    explicit Quiet(bool quiet):old(quiet?std::cout.rdbuf(&sink):nullptr){}
    ~Quiet(){if(old)std::cout.rdbuf(old);}
};
}
namespace Foam {
label DAHisaFoam::solveHybrid()
{
#ifdef CODI_ADR
    FatalErrorInFunction<<"Hybrid primal must run through the native solver, not the reverse AD instance"<<exit(FatalError);
    return 1;
#else
    const dictionary cfg=meshPtr_->solutionDict().subDict("hisaHybrid");
    solverModule& solver=hisaSolverPtr_(); Time& time=runTimePtr_();
    if(!solver.steadyState() || daOptionPtr_->getOption<word>("adjStateOrdering")!="state")
        FatalErrorInFunction<<"Hybrid requires steady state and adjStateOrdering=state"<<exit(FatalError);
    const label nc=getNLocalCells(),ns=getNLocalAdjointStates();
    if(ns!=6*nc || !meshPtr_->foundObject<volScalarField>("nuTilda"))
        FatalErrorInFunction<<"Hybrid currently supports U,p,T,nuTilda RANS only"<<exit(FatalError);
    const label maxOuter=cfg.lookupOrDefault<label>("ncycles",8000);
    const label ptcMin=cfg.lookupOrDefault<label>("ptcMinIterations",10);
    const label interval=cfg.lookupOrDefault<label>("monitorInterval",10);
    const label recoveryMin=cfg.lookupOrDefault<label>("ptcRecoverySteps",10);
    const double switchTol=cfg.lookupOrDefault<scalar>("nkSwitchTol",1e-6);
    const double finalTol=cfg.lookupOrDefault<scalar>("steadyRelTol",1e-8);
    const double absTol=cfg.lookupOrDefault<scalar>("steadyAbsTol",1e-14);
    wordList monitors;
    if(cfg.found("monitors")) cfg.lookup("monitors")>>monitors;
    else {monitors.setSize(6);monitors[0]="relative";monitors[1]="absolute";monitors[2]="ksp";monitors[3]="linear";monitors[4]="alpha";monitors[5]="cfl";}
    auto show=[&](const word& key){return monitors.found(key);};
    forAll(monitors,i) if(monitors[i]!="relative"&&monitors[i]!="absolute"&&monitors[i]!="blocks"&&monitors[i]!="combined"&&monitors[i]!="ksp"&&monitors[i]!="linear"&&monitors[i]!="alpha"&&monitors[i]!="cfl") FatalErrorInFunction<<"Unknown hybrid monitor "<<monitors[i]<<exit(FatalError);
    const bool detailed=cfg.lookupOrDefault<bool>("detailed",false);
    if(maxOuter<1||ptcMin<0||interval<1||recoveryMin<0||finalTol<=0||switchTol<finalTol)
        FatalErrorInFunction<<"Invalid hisaHybrid thresholds or iteration limits"<<exit(FatalError);
    hisaNK::Options opts;opts.subspace=cfg.lookupOrDefault<label>("nkSubspace",60);opts.maxLinear=cfg.lookupOrDefault<label>("nkMaxLinearIterations",60);
    opts.overlap=cfg.lookupOrDefault<label>("nkASMOverlap",1);opts.fill=cfg.lookupOrDefault<label>("nkILUFill",2);opts.lag=cfg.lookupOrDefault<label>("nkJacobianLag",20);opts.rtol=cfg.lookupOrDefault<scalar>("nkLinearTol",.3);opts.ew=cfg.lookupOrDefault<bool>("nkUseEW",true);
    const dictionary scaling=daOptionPtr_->getAllOptions().subDict("normalizeStates");
    scalarList scales(4);wordList names(4);names[0]="U";names[1]="p";names[2]="T";names[3]="nuTilda";
    forAll(scales,i) {scales[i]=scaling.get<scalar>(names[i]);if(scales[i]<=0)FatalErrorInFunction<<"Invalid state scale"<<exit(FatalError);}
    auto block=[&](label i)->label{return i<3*nc?0:(i<4*nc?1:(i<5*nc?2:3));};
    auto state=[&](){scalarField raw(ns);getOFFields(raw.data());hisaNK::Array z(ns);forAll(raw,i){double v;assignValueCheckAD(v,raw[i]);z[i]=v/scales[block(i)];}return z;};
    auto residual=[&](const hisaNK::Array&z){Quiet q(!detailed);scalarField w(ns);forAll(w,i)w[i]=z[i]*scales[block(i)];updateOFFields(w.data());hisaNK::Array r(ns);getResiduals(r.data());return r;};
    auto admissible=[&](const hisaNK::Array&z){int ok=1,global=0;for(label i=0;i<ns;++i)if(!std::isfinite(z[i])||(i>=3*nc&&i<5*nc&&z[i]<=0)||(i>=5*nc&&z[i]<0))ok=0;MPI_Allreduce(&ok,&global,1,MPI_INT,MPI_MIN,PETSC_COMM_WORLD);return global!=0;};
    auto norms=[&](const hisaNK::Array&r){double sum[4]={0,0,0,0},total[4];for(label i=0;i<ns;++i)sum[block(i)]+=r[i]*r[i];MPI_Allreduce(sum,total,4,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);scalarList result(4);label globalN=returnReduce(nc,sumOp<label>());forAll(result,i)result[i]=std::sqrt(total[i]/(globalN*(i==0?3:1)));return result;};
    hisaNK::Array z=state(),r=residual(z);scalarList initial=norms(r);
    label iter=0,nkCount=0,recovery=recoveryMin;word stage="PTC";
    IOdictionary saved(IOobject("hisaHybridState",time.timeName(),"uniform",meshPtr_(),IOobject::READ_IF_PRESENT,IOobject::NO_WRITE,false));
    bool reset=cfg.lookupOrDefault<bool>("resetBaseline",false);
    if(saved.found("initialResidual")&&!reset){saved.lookup("initialResidual")>>initial;saved.lookup("stage")>>stage;iter=saved.get<label>("outerIteration");nkCount=saved.get<label>("nkIteration");recovery=saved.get<label>("recoverySteps");scalarList oldScales;saved.lookup("stateScales")>>oldScales;
        if(oldScales!=scales||initial.size()!=4||saved.get<label>("globalCells")!=returnReduce(nc,sumOp<label>()))FatalErrorInFunction<<"Hybrid checkpoint incompatible with state scales or mesh size"<<exit(FatalError);
    }else if(time.value()!=0&&!reset)FatalErrorInFunction<<"Restart missing hisaHybridState. Supply a compatible checkpoint or explicitly resetBaseline."<<exit(FatalError);
    if(stage!="PTC"&&stage!="NK")FatalErrorInFunction<<"Invalid saved stage"<<exit(FatalError);
    forAll(initial,i)if(!std::isfinite(initial[i])||initial[i]<0)FatalErrorInFunction<<"Invalid initial residual"<<exit(FatalError);
    hisaNK::Array weights(ns);for(label i=0;i<ns;++i)weights[i]=1/std::max(double(initial[block(i)]),1e-30);
    const label globalCells=returnReduce(nc,sumOp<label>());
    double initialSquared=0;
    forAll(initial,i) {double a=initial[i]/std::max(double(initial[i]),1e-30);initialSquared+=globalCells*(i==0?3:1)*a*a;}
    double initialTotal=std::sqrt(initialSquared);
    if(!reset && saved.found("initialTotalResidual")) initialTotal=saved.get<scalar>("initialTotalResidual");
    if(!std::isfinite(initialTotal)||initialTotal<0) FatalErrorInFunction<<"Invalid total residual baseline"<<exit(FatalError);
    label ptcCount=iter-nkCount;
    if(ptcCount<0) FatalErrorInFunction<<"Invalid checkpoint iteration counts"<<exit(FatalError);
    bool colored=false;
    auto pcMatrix=[&](){Quiet q(!detailed);residual(z);if(!colored){runColoring();colored=true;}Mat m;hisaNK::check(MatCreate(PETSC_COMM_WORLD,&m));calcdRdWT(1,m);hisaNK::check(MatTranspose(m,MAT_INPLACE_MATRIX,&m));return m;};
    std::unique_ptr<hisaNK::Solver> nk(new hisaNK::Solver(opts,weights,residual,admissible,pcMatrix));
    auto persist=[&](const word&status){IOdictionary out(IOobject("hisaHybridState",time.timeName(),"uniform",meshPtr_(),IOobject::NO_READ,IOobject::NO_WRITE,false));out.add("initialResidual",initial);out.add("initialTotalResidual",initialTotal);out.add("ptcIteration",ptcCount);out.add("stateScales",scales);out.add("globalCells",returnReduce(nc,sumOp<label>()));out.add("outerIteration",iter);out.add("nkIteration",nkCount);out.add("stage",stage);out.add("recoverySteps",recovery);out.add("status",status);out.regIOobject::write();};
    Info<<"Hybrid initialResidual "<<initial<<" initialTotalResidual "<<initialTotal<<" cumulative iteration "<<iter<<" stage "<<stage<<endl;
    word status="maxIterations";hisaNK::Result linear;
    while(true){scalarList values=norms(r);
        hisaNK::Array scaled=r;
        for(label i=0;i<ns;++i) scaled[i]*=weights[i];
        const double ab=hisaNK::norm(scaled);
        const double rel=ab/std::max(initialTotal,1e-30);
        if(!std::isfinite(ab)) FatalErrorInFunction<<"Nonfinite hybrid residual"<<exit(FatalError);
        const bool converged=(ab<=absTol || rel<=finalTol);
        if(iter%interval==0||converged||iter>=maxOuter) {
            Info<<"Hybrid "<<iter<<' '<<stage;
            if(show("relative")) Info<<" relative "<<rel;
            if(show("absolute")) Info<<" absolute "<<ab;
            if(show("blocks")) Info<<" blockRMS "<<values;
            if(show("combined")) Info<<" combined "<<rel;
            if(show("cfl")) {
                if(stage=="PTC") {
                    if(meshPtr_->foundObject<uniformDimensionedScalarField>("pseudoCoNum"))
                        Info<<" CFL "<<meshPtr_->lookupObject<uniformDimensionedScalarField>("pseudoCoNum").value();
                    if(meshPtr_->foundObject<volScalarField>("pseudoCoField")) {
                        const scalarField& co=meshPtr_->lookupObject<volScalarField>("pseudoCoField").primitiveField();
                        Info<<" localCFL[min,max] "<<gMin(co)<<' '<<gMax(co);
                    }
                } else Info<<" CFL n/a";
            }
            if(show("ksp")) Info<<" KSP "<<linear.iterations;
            if(show("linear")) Info<<" linRel "<<linear.linear;
            if(show("alpha")) Info<<" alpha "<<linear.alpha;
            Info<<endl;
        }
        if(converged){status="converged";break;}if(iter>=maxOuter){status="maxIterations";break;}bool running;{Quiet q(!detailed);running=time.run();}if(!running){status="endTime";break;}
        if(stage=="NK"&&(rel>switchTol||ptcCount<ptcMin)){Info<<"Hybrid NK -> PTC: threshold or minimum PTC requirement"<<endl;stage="PTC";recovery=0;nk->restart();}
        if(stage=="PTC"&&rel<=switchTol&&ptcCount>=ptcMin&&recovery>=recoveryMin){Info<<"Hybrid PTC -> NK: switch threshold reached"<<endl;stage="NK";nk->restart();}
        {Quiet q(!detailed);solver.preTimeIncrement();time++;}
        ++iter;
        if(stage=="NK"){
            ++nkCount;linear=nk->step(z);r=residual(z);
            if(!linear.accepted){Info<<"Hybrid NK -> PTC: rejected update, KSP reason "<<linear.reason<<endl;stage="PTC";recovery=0;nk->restart();}
        }else{
            {Quiet q(!detailed);solver.beginTimeStep();if(!solver.solnControl().loop())FatalErrorInFunction<<"PTC stopped before hybrid controller"<<exit(FatalError);solver.outerIteration();daTurbulenceModelPtr_->correct(detailed);}
            z=state();r=residual(z);++recovery;++ptcCount;
        }
        {Quiet q(!detailed);calcAllFunctions(detailed);}
        if(time.writeTime()){Quiet q(!detailed);time.write();persist("running");}
    }
    residual(z);{Quiet q(!detailed);time.writeNow();persist(status);}primalFinalTimeIndex_=time.timeIndex();
    Info<<"Hybrid finished: "<<status<<" iteration "<<iter<<" NK iterations "<<nkCount<<endl;
    return status=="converged"?0:1;
#endif
}
}
