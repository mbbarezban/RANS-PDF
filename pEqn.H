
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

tUEqn.clear();

surfaceScalarField phiHbyA("phiHbyA",fvc::interpolate(rho)*fvc::flux(HbyA));

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf);

adjustPhi(phiHbyA, U, p);

while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rhorAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U = HbyA - rAU * fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
//U.relax();

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();    
}

rho = thermo.rho(); // rho = p / (R * T)
