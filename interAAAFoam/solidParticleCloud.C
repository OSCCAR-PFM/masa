/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidParticleCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "twoPhaseMixture.H"
//#include "incompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticleCloud::solidParticleCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    bool readFields
)
:
    Cloud<solidParticle>(mesh, cloudName, false),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            "particleProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    rhop_(dimensionedScalar(particleProperties_.lookup("rhop")).value()),
    e_(dimensionedScalar(particleProperties_.lookup("e")).value()),
    mu_(dimensionedScalar(particleProperties_.lookup("mu")).value()),

/////////////////////////////////////////edit////////////////////////////////////
    mudrop_(dimensionedScalar(particleProperties_.lookup("mudrop")).value()),
    sigmad_(dimensionedScalar(particleProperties_.lookup("sigmad")).value()),
    PrPr_(dimensionedScalar(particleProperties_.lookup("PrPr")).value()),
    kapkap_(dimensionedScalar(particleProperties_.lookup("kapkap")).value()),
    CpCppart_(dimensionedScalar(particleProperties_.lookup("CpCppart")).value()),


    posP1_(dimensionedVector(particleProperties_.lookup("posP1")).value()),
    dP1_(dimensionedScalar(particleProperties_.lookup("dP1")).value()),
    UP1_(dimensionedVector(particleProperties_.lookup("UP1")).value()),
    TP1_(dimensionedScalar(particleProperties_.lookup("TP1")).value()),
    posP2_(dimensionedVector(particleProperties_.lookup("posP2")).value()),
    dP2_(dimensionedScalar(particleProperties_.lookup("dP2")).value()),
    TP2_(dimensionedScalar(particleProperties_.lookup("TP2")).value()),
    UP2_(dimensionedVector(particleProperties_.lookup("UP2")).value()),

///2way coupling/////

    smom_(mesh_.nCells(), vector::zero),
    sTemp_(mesh_.nCells(), 0),

    correctalpha1_(mesh_.nCells(), 0),
    correctU_(mesh_.nCells(), vector::zero),

    sourcecorrectT_(mesh_.nCells(), 0), //Temp
    sourcecorrectalpha1_(mesh_.nCells(), 0),
    sourcecorrectU_(mesh_.nCells(), vector::zero),
////////////////////

    tInjStart_(dimensionedScalar(particleProperties_.lookup("tInjStart")).value()),
    tInjEnd_(dimensionedScalar(particleProperties_.lookup("tInjEnd")).value())

/////////////////////////////////////////////////////////////////////////////////

{
    if (readFields)
    {
        solidParticle::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
///////////////////////////////////////////////////////////////////////////////////////
void Foam::solidParticleCloud::inject(Foam::solidParticle::trackingData &td) 
{
    Info<< "*** ONE" << endl;

    //const volScalarField& alpha1 = mesh_.lookupObject<const volScalarField>("alpha1");
    //const volScalarField& alpha1 = twoPhaseProperties.alpha1();
    const volScalarField& alpha1 = mesh_.lookupObject<const volScalarField>("alpha.phase1");

    Info<< "*** TWO" << endl;

    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volVectorField& fluctuation1 = mesh_.lookupObject<const volVectorField>("fluctuation1");
    const volVectorField& fluctuation2 = mesh_.lookupObject<const volVectorField>("fluctuation2");
    const volVectorField& UFLUK = mesh_.lookupObject<const volVectorField>("UFLUK");
    const volScalarField& criterion2 = mesh_.lookupObject<const volScalarField>("criterion2");
    const volScalarField& k = mesh_.lookupObject<const volScalarField>("k");
    const volScalarField& TLS = mesh_.lookupObject<const volScalarField>("TLS");
    const volScalarField& Ndrop = mesh_.lookupObject<const volScalarField>("Ndrop");
    const volScalarField& Lzero = mesh_.lookupObject<const volScalarField>("Lzero");
    const volScalarField& T = mesh_.lookupObject<const volScalarField>("T");

    interpolationCellPoint<scalar> alpha1Interp(alpha1);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<vector> fluctuation1Interp(fluctuation1);
    interpolationCellPoint<vector> fluctuation2Interp(fluctuation2);
    interpolationCellPoint<vector> UFLUKInterp(UFLUK);
    interpolationCellPoint<scalar> kInterp(k);
    interpolationCellPoint<scalar> criterion2Interp(criterion2);
    interpolationCellPoint<scalar> TLSInterp(TLS);
    interpolationCellPoint<scalar> LzeroInterp(Lzero);
    interpolationCellPoint<scalar> NdropInterp(Ndrop);
    interpolationCellPoint<scalar> TInterp(T);

    Info<< "*** THREE" << endl;

//////
    forAll(mesh_.cells(), cellI)
    {
        if (alpha1[cellI]  < scalar(0.9) && alpha1[cellI]  > scalar(0.1)  )
        {
           	if ( TLS[cellI] <  pow(mesh_.V()[cellI],(1.0/3.0)) && criterion2[cellI] > scalar(0) && Ndrop[cellI] > scalar(1.0) )   
            {
                label tetFaceI = 1; //-1
                label tetPtI = 1; //-1
                const point& pt=mesh_.C()[cellI];
                mesh_.findTetFacePt(cellI,pt,tetFaceI,tetPtI);

                solidParticle* ptr1= new solidParticle(td.cloud().mesh_,mesh_.C()[cellI],cellI,tetFaceI,tetPtI, TLS[cellI], UFLUK[cellI]); 

                Cloud<solidParticle>::addParticle(ptr1);

             /*   if (Ndrop[cellI] > scalar(1.99))
                {
                    solidParticle* ptr2= new solidParticle(td.cloud().mesh_,mesh_.C()[cellI],cellI,tetFaceI,tetPtI, TLS[cellI], UFLUK[cellI]); 
                    solidParticle* ptr3= new solidParticle(td.cloud().mesh_,mesh_.C()[cellI],cellI,tetFaceI,tetPtI, TLS[cellI], UFLUK[cellI]); 
                    Cloud<solidParticle>::addParticle(ptr2);
                    Cloud<solidParticle>::addParticle(ptr3);
                }*/

                sourcecorrectalpha1_=scalar(1.0)*.52*pow(TLS[cellI],3)/(mesh_.time().deltaT().value()*mesh_.V()[cellI]); 
                sourcecorrectU_=  (scalar(1.0) * .52 *pow(TLS[cellI],3)/(mesh_.V()[cellI]))* U[cellI]; 
                sourcecorrectT_=( scalar(1.0) * .52 *pow(TLS[cellI],3)/(mesh_.V()[cellI]))*T[cellI];
            }                     
         }
    }
}

////////////////////////////////////////////////////////////////////////////////////////
bool Foam::solidParticleCloud::hasWallImpactDistance() const
{
    return true;
}


void Foam::solidParticleCloud::move(const dimensionedVector& g)
{
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");

 // const volScalarField& alpha1 = mesh_.lookupObject<const volScalarField>("alpha1");
    const volScalarField& alpha1 = mesh_.lookupObject<const volScalarField>("alpha.phase1");

    const volScalarField& T = mesh_.lookupObject<const volScalarField>("T");
    const volScalarField& criterion2 = mesh_.lookupObject<const volScalarField>("criterion2");
    const volScalarField& k = mesh_.lookupObject<const volScalarField>("k");
    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> nuInterp(nu);
    interpolationCellPoint<scalar> kInterp(k);
    interpolationCellPoint<scalar> criterion2Interp(criterion2);
    interpolationCell<scalar> alpha1Interp(alpha1);
    interpolationCellPoint<scalar> TInterp(T);  //Temp

    ///////////////////////////////////////////////////
        smom_=vector::zero;
        sTemp_=0;
        correctalpha1_=0;
        correctU_=vector::zero;
        sourcecorrectalpha1_=0;
        sourcecorrectU_=vector::zero;
        sourcecorrectT_=0;
       //correctT_=0;
    //////////////////////////////////////////////////



    solidParticle::trackingData td(*this, rhoInterp, UInterp, nuInterp, alpha1Interp, TInterp, g.value());  ///

    Cloud<solidParticle>::move(td, mesh_.time().deltaTValue());

    if(mesh_.time().value() > td.cloud().tInjStart_ && mesh_.time().value() < td.cloud().tInjEnd_)      
    {
        this->inject(td);
    }
}


// ************************************************************************* //
