surfaceScalarField rAUrelf("rAUrelf", fvc::interpolate(rAUrel));

volVectorField HbyA("HbyA", Urel);
HbyA = rAUrel*UrelEqn().H(); // H/A = "H operation source" divided by the
                             //       "central coefficient"
// note: H is the RHS momentum terms excluding pressure
//       == S - sum_over_neighbor_cells{A*U}
//       where S is the discretization source term (i.e. RHS terms?)
//             A is the matrix coefficient
//             U is the uncorrected velocity
//
//       The momentum equation == H - grad{p}
//       ref: http://openfoamwiki.net/index.php/OpenFOAM_guide/H_operator

if (pimple.nCorrPISO() <= 1)
{
    UrelEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + rAUrelf*fvc::ddtCorr(Urel, Urelf) //SRFPimpleFoam: ddtCorr(Urel,phi)
);

if (p.needReference()) // not in SRFPimpleFoam
{
    fvc::makeRelative(phiHbyA, Urel); // not in SRFPimpleFoam
    adjustPhi(phiHbyA, Urel, p);
    fvc::makeAbsolute(phiHbyA, Urel); // not in SRFPimpleFoam
}

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        // SRFPimpleFoam uses rAUrel, not interpolated rAUrel
        fvm::laplacian(rAUrelf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux(); // same as SRFPimpleFoam
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
Urel = HbyA - rAUrel*fvc::grad(p);
Urel.correctBoundaryConditions();
fvOptions.correct(Urel);

{ // not in SRFPimpleFoam
    // for non-moving meshes, Urelf should be equal to phi...
    Urelf = fvc::interpolate(Urel);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Urelf += n*(phi/mesh.magSf() - (n & Urelf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, Urel); // not in SRFPimpleFoam

