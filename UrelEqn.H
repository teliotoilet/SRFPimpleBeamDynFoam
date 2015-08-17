// Solve the Momentum equation

tmp<fvVectorMatrix> UrelEqn
(
    fvm::ddt(Urel)
  + fvm::div(phi, Urel)
  + turbulence->divDevReff(Urel)
  + SRF->Su()
 ==
    fvOptions(Urel)
);

UrelEqn().relax();

fvOptions.constrain(UrelEqn());

volScalarField rAUrel(1.0/UrelEqn().A());

if (pimple.momentumPredictor())
{
    solve(UrelEqn() == -fvc::grad(p));

    fvOptions.correct(Urel);
}
