Info<< "\nReading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

Info<< "\nCreating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "\nReading field U\n" << endl;
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

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, simple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<mcCompTurbulenceModel> turbulence
(
    new mcCompTurbulenceModel(rho, U, phi, thermo)
);

Info<< "\nCreating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo)
);

Info<< "\nCreating Radiation model\n" << endl;
autoPtr<radiationModel> radiation(radiationModel::New(thermo.T()));

Info<< "\nCreating PDF Cloud\n" << endl;
reactPdfCloud pdfCloud
(
        turbulence(),
        thermo,
        reaction(),
        radiation(),
        U,
        rho
);

//volScalarField Qdot
//(
//    IOobject
//    (
//        "Qdot",
//        runTime.timeName(),
//        mesh,
//        IOobject::READ_IF_PRESENT,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
//);

#include "createFvOptions.H"
