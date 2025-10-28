    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    // Block vector field for velocity (first entry) and pressure (second
    // entry).
    Info << "Creating field Up\n" << endl;
    volVector4Field Up
    (
        IOobject
        (
            "Up",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector4("zero", dimless, vector4::zero)
    );

    Info<< "Creating field rAU\n" << endl;
    // Note: depending on the MRF treatment, only one of the rAU fields
    // will be used.
    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    mesh.schemesDict().setFluxRequired(p.name());
    mesh.schemesDict().setFluxRequired("pcorr");
