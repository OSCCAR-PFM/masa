    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//------------------Modified---------------------//
Info<< "Reading field T\n" << endl;

	volScalarField T
(
	IOobject
	(
		"T",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
//---------------------End----------------------//

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



    Info<< "Reading field U\n" << endl;
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

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi);

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

////////////////////////////////////////////////////////////////////////////////
const dimensionedScalar& cp1 = twoPhaseProperties.cp1();
const dimensionedScalar& cp2 = twoPhaseProperties.cp2();
////////////////////////////////////////////////////////////////////////////////



///////////

/* Extract the dictionary from the database
const dictionary& transportProperties = db().lookupObject<IOdictionary>
(
 "transportProperties"
);
*/
//////////////////


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();

////////////////////////////////////////
/*Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixture::mu() const
{
const volScalarField limitedAlpha1
(
min(max(alpha1_, scalar(0)), scalar(1))
);
return tmp<volScalarField>

(
	new volScalarField
		(
	 	"mu",
		limitedAlpha1*rho1_*nuModel1_->nu()
		+ (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
		)
	);
}
*/


////////////////////////////////////////



    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );
////////////////////////////////////////////////////////////////////////////



 
//////////////////////////////////////////////////////////////////////////
    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */


//------------------Modified---------------------//
  Info<< "Calculating rho*cp\n" << endl;
 
    volScalarField rhoCp
    (
	    IOobject
	    (
	        "rho*Cp",
	        runTime.timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
        ),
        alpha1*rho1*cp1 + rho2*cp2*(scalar(1) - alpha1),
        alpha1.boundaryField().types()
    );
         rhoCp.oldTime();

    Info<< "Reading / calculating rho*phi*cp\n" << endl;
    surfaceScalarField rhoPhiCpf
    (
	IOobject
	(
		"rho*phi*cpf",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	rhoPhi*cp1
    );

///////////////////////////////////////////////////////

 Info<< "Calculating cpcp\n" << endl;
 
    volScalarField CpCp
    (
	    IOobject
	    (
	        "CpCp",
	        runTime.timeName(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
        ),
        alpha1 * cp1 + cp2*(scalar(1) - alpha1),
        alpha1.boundaryField().types()
    );
         CpCp.oldTime();
//---------------------End-----------------------//

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );



    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
/////////////////////////////////////////////////////////

 //volScalarField myVectorField("myVectorField", myVector & mesh.C()); 
   


/////////////////////////////////////////////////////
  /* Info<< "Reading field criterion\n" << endl;


     volScalarField criterion
    (
        IOobject
        (
            "criterion",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),

mesh
    ); */
/////////////////////////////////////////////////////////


volScalarField cellVolume
(
    IOobject
    (
        "cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
//dimensionedScalar("0", dimensionSet(0, 3, 0, 0, 0), 0)  
    dimensionedScalar("zero",dimVolume,0.0)
);
    cellVolume.internalField() = mesh.V();

////////////////////////////////////////////////

volScalarField TLS
(
    IOobject
    (
        "TLS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
   dimensionedScalar("0", dimensionSet(0, 1, 0, 0, 0), 0)


);
////////////////////////////////////////////////

volScalarField Lzero
(
    IOobject
    (
        "Lzero",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
   dimensionedScalar("0", dimensionSet(0, 1, 0, 0, 0), 0)


);

////////////////////////////////////////////////////

   Info<< "Reading field criterion2\n" << endl;


    volScalarField criterion2
    (
        IOobject
        (
            "criterion2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),

mesh
    );

///////////////////////////////////////////////////
    Info<< "Reading field Ndrop\n" << endl;

volScalarField Ndrop
(
    IOobject
    (
        "Ndrop",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        mesh,
   dimensionedScalar("Ndrop",dimensionSet(0, 0, 0, 0, 0), scalar(0.0))


);  

//////////////////////////////////////
/*
    volScalarField Vdrop
    (
        IOobject
        (
            "Vdrop",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedScalar("Vdrop",dimensionSet(0, 3, 0, 0, 0), 0)
    );
 
 */
 
///////////////////////////////////////////
      volVectorField momSp   // volVectorField
    (
        IOobject
        (
            "momSp",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
  // dimensionedScalar("momSp",dimensionSet(0, 1, -2, 0, 0), scalar(0.0))
    );

/////////////////////////////////////////////////
 /*
     volScalarField Alnegative
    (
        IOobject
        (
            "Alnegative",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedScalar("Alnegative",dimensionSet(0, 0, 0, 0, 0), scalar(0.0))
    );
  */
////////////////////////////////////




    volScalarField criterion2LowerLimit
    (
        IOobject
        (
            "criterion2LowerLimit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "criterion2LowerLimit",
            dimensionSet(1, 2, -2, 0, 0, 0, 0),
            scalar(3e-10)  ///scalar(1.0e-10) 
        )
    );

/////////////////////////////////////////////
volScalarField VOLVOL
    (
        IOobject
        (
            "VOLVOL",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "VOLVOL",
            dimensionSet(0, 3, 0, 0, 0, 0, 0),
            scalar(1.0e-9)
        )
    );
VOLVOL.internalField() = mesh.V();

///////////////////////////////////////////////////


     volScalarField AlSu
    (
        IOobject
        (
            "AlSu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedScalar("AlSu",dimensionSet(0, 0, -1, 0, 0), scalar(0.0))
  // dimensionedScalar("AlSu",dimensionSet(0, 0, 0, 0, 0), scalar(0.0))
    ); 
 //////////////////////////////////////////////

     volScalarField MassAlSu
    (
        IOobject
        (
            "MassAlSu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedScalar("MassAlSu",dimensionSet(1, 0, -1, 0, 0), scalar(0.0))
    ); 
///////////////////////////////////////////////
//
volScalarField frequ
(
    IOobject
    (
        "frequ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
        mesh,
   dimensionedScalar("frequ",dimensionSet(0, 0, -1, 0, 0), scalar(0.0))


); 
///////////////////////////////////////////////


   // Info<< "Reading field U\n" << endl;
    volVectorField Vorticity
    (
        IOobject
        (
            "Vorticity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ); 
///////////////////////////////////////////////


   // Info<< "Reading field U\n" << endl;
    volVectorField fluctuation1
    (
        IOobject
        (
            "fluctuation1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedVector("fluctuation1",dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0,0,0))
    );

    volVectorField fluctuation2
    (
        IOobject
        (
            "fluctuation2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedVector("fluctuation2",dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0,0,0))
    );

   volVectorField grad1Alpha
    (
        IOobject
        (
            "grad1Alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedVector("grad1Alpha",dimensionSet(0, -1, 0, 0, 0, 0, 0), vector(0,0,0))
    );

    volVectorField gradAlphaUSE
    (
        IOobject
        (
            "gradAlphaUSE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedVector("gradAlphaUSE",dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(0,0,0))
    );

    volVectorField UFLUK
    (
        IOobject
        (
            "UFLUK",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
   dimensionedVector("UFLUK",dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0,0,0))
    );


volScalarField CPCPPARTICLE
(
    IOobject
    (
        "CPCPPARTICLE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
        mesh,
   dimensionedScalar("CPCPPARTICLE",dimensionSet(0, 2, -2, -1, 0, 0, 0), scalar(4200))


);

