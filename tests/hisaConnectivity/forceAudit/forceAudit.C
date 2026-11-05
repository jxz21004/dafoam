#include "fvCFD.H"
#include "fluidThermo.H"
int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList times = timeSelector::select0(runTime, args);
    runTime.setTime(times.last(), times.size()-1);
    #include "createMesh.H"
    autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    volScalarField nut(IOobject("nut",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    volScalarField rho("auditRho",thermo->rho());
    volScalarField mu("auditMu",thermo->mu());
    volSymmTensorField tau("auditTau",-(mu+rho*nut)*dev(twoSymm(fvc::grad(U))));
    label patch = mesh.boundaryMesh().findPatchID("wing");
    if (patch<0) FatalErrorInFunction << "wing patch missing" << exit(FatalError);
    const auto& sf = mesh.Sf().boundaryField()[patch];
    const auto& p = thermo->p().boundaryField()[patch];
    vector fp=vector::zero, fv=vector::zero, area=vector::zero, fpGauge=vector::zero;
    scalar pinf=26499.89786492552;
    OFstream faces(runTime.path()/"force_audit_faces.csv");
    faces.precision(17);
    faces << "x,y,z,sx,sy,sz,p,mu,nut,rho,fpx,fpy,fpz,fvx,fvy,fvz" << nl;
    forAll(sf,i)
    {
        vector pn=sf[i]*p[i]; vector tv=sf[i]&tau.boundaryField()[patch][i];
        fp+=pn; fv+=tv; area+=sf[i]; fpGauge+=sf[i]*(p[i]-pinf);
        vector c=mesh.Cf().boundaryField()[patch][i];
        faces<<c.x()<<','<<c.y()<<','<<c.z()<<','<<sf[i].x()<<','<<sf[i].y()<<','<<sf[i].z()<<','<<p[i]<<','<<mu.boundaryField()[patch][i]<<','<<nut.boundaryField()[patch][i]<<','<<rho.boundaryField()[patch][i]<<','<<pn.x()<<','<<pn.y()<<','<<pn.z()<<','<<tv.x()<<','<<tv.y()<<','<<tv.z()<<nl;
    }
    Info().precision(17);
    Info<<"AUDIT pressure "<<fp<<" viscous "<<fv<<" gaugePressure "<<fpGauge<<" areaVector "<<area<<" faces "<<sf.size()<<nl;
    return 0;
}
