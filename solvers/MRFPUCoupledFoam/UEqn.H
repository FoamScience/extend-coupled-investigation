{
    // Update boundary velocity for MRF consistency
    mrfZones.correctBoundaryVelocity(U);

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // Handle MRF implicit sources
    tmp<volTensorField> tTU;

    if (addMRF)
    {
        // Create tensor field for MRF
        tTU = tmp<volTensorField>
        (
            new volTensorField
            (
                IOobject
                (
                    "TU",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedTensor("zero", dimless/dimTime, tensor::zero)
            )
        );
        volTensorField& TU = tTU();

        // Add implicit MRF source (Hodge dual of omega)
        TU += *mrfZones.omega();

        // Inverse including MRF contribution
        trTU = inv(TU + tensor(I)*UEqn.A());
        trTU().rename("rAU");
    }
    else
    {
        // No MRF: simple scalar rAU
        trAU = 1.0/UEqn.A();
        trAU().rename("rAU");
    }

    // Under-relax momentum
    UEqn.relax();

    // Insert momentum equation into block system
    UpEqn.insertEquation(0, UEqn);

#   include "addImplicitMRF.H"
#   include "addBlockCoupledBC.H"
}
