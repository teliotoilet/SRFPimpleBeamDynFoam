if (mesh.changing())
{
    forAll(Urel.boundaryField(), patchI)
    {
        if (Urel.boundaryField()[patchI].fixesValue())
        {
            Urel.boundaryField()[patchI].initEvaluate();
        }
    }

    forAll(Urel.boundaryField(), patchI)
    {
        if (Urel.boundaryField()[patchI].fixesValue())
        {
            Urel.boundaryField()[patchI].evaluate();

            phi.boundaryField()[patchI] =
                Urel.boundaryField()[patchI]
              & mesh.Sf().boundaryField()[patchI];
        }
    }
}

{
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUrelf("rAUrelf", dimTime, 1.0);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUrelf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

    #include "continuityErrs.H"
}
