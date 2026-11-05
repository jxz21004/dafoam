#include "fvCFD.H"
#include "fluidThermo.H"
#include "fluxScheme.H"
int main(int argc,char *argv[])
{
 timeSelector::addOptions();
 #include "setRootCase.H"
 #include "createTime.H"
 instantList times=timeSelector::select0(runTime,args);runTime.setTime(times.last(),0);
 #include "createMesh.H"
 autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));
 volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
 volScalarField rho("rho",thermo->rho());
 volVectorField rhoU("rhoU",rho*U);
 volScalarField rhoE("rhoE",rho*(thermo->he()+0.5*magSqr(U)));
 surfaceScalarField phi("auditPhi",linearInterpolate(rhoU)&mesh.Sf());
 surfaceVectorField phiUp("auditPhiUp",linearInterpolate(rhoU*U+thermo->p()*tensor::I)&mesh.Sf());
 surfaceScalarField phiEp("auditPhiEp",linearInterpolate((rhoE+thermo->p())*U)&mesh.Sf());
 surfaceVectorField Up("auditUp",linearInterpolate(U)*mesh.magSf());
 autoPtr<fluxScheme> scheme(fluxScheme::New(mesh.schemesDict(),thermo(),rho,U,rhoU,rhoE));
 scheme->calcFlux(phi,phiUp,phiEp,Up);
 label patch=mesh.boundaryMesh().findPatchID("wing");
 const auto& sf=mesh.Sf().boundaryField()[patch];const auto& p=thermo->p().boundaryField()[patch];
 scalar maxMass=0,maxEnergy=0,maxMismatch=0,maxPdiff=0;
 OFstream out(runTime.path()/"wall_flux_audit.csv");out.precision(17);
 out<<"face,x,y,z,pwall,powner,Towner,rhoowner,Ux,Uy,Uz,normalVelocityOwner,massFlux,energyFlux,momentumMismatch,pFlux"<<nl;
 forAll(sf,i)
 {
  label cell=mesh.boundary()[patch].faceCells()[i];vector normal=sf[i]/mag(sf[i]);vector centre=mesh.Cf().boundaryField()[patch][i];
  vector mismatch=phiUp.boundaryField()[patch][i]-p[i]*sf[i];
  scalar pf=(phiUp.boundaryField()[patch][i]&normal)/mag(sf[i]);
  maxMass=max(maxMass,mag(phi.boundaryField()[patch][i]));maxEnergy=max(maxEnergy,mag(phiEp.boundaryField()[patch][i]));maxMismatch=max(maxMismatch,mag(mismatch));maxPdiff=max(maxPdiff,mag(pf-p[i]));
  out<<i<<','<<centre.x()<<','<<centre.y()<<','<<centre.z()<<','<<p[i]<<','<<thermo->p()[cell]<<','<<thermo->T()[cell]<<','<<rho[cell]<<','<<U[cell].x()<<','<<U[cell].y()<<','<<U[cell].z()<<','<<(U[cell]&normal)<<','<<phi.boundaryField()[patch][i]<<','<<phiEp.boundaryField()[patch][i]<<','<<mag(mismatch)<<','<<pf<<nl;
 }
 Info().precision(17);Info<<"WALL_AUDIT maxMass "<<maxMass<<" maxEnergy "<<maxEnergy<<" maxMomentumMismatch "<<maxMismatch<<" maxPressureDifference "<<maxPdiff<<nl;
 return 0;
}
