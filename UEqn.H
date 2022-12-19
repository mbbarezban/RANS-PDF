// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
          fvm::div(phi, U)
        + turbulence->divDevRhoReff(U)
        ==
          fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (simple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
}
