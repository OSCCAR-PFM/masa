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

Class
    Foam::solidParticleCloud

Description
    A Cloud of solid particles

SourceFiles
    solidParticleCloudI.H
    solidParticleCloud.C
    solidParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticleCloud_H
#define solidParticleCloud_H

#include "Cloud.H"
#include "solidParticle.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class solidParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class solidParticleCloud
:
    public Cloud<solidParticle>
{
    // Private data

        const fvMesh& mesh_;

        IOdictionary particleProperties_;

        scalar rhop_;
        scalar e_;
        scalar mu_;
/////////////////////////
vectorField smom_;
scalarField sTemp_;
                scalarField correctalpha1_;
                scalarField sourcecorrectT_;
                vectorField correctU_;
                scalarField sourcecorrectalpha1_;//
                vectorField sourcecorrectU_;//
        scalar mudrop_;
        scalar sigmad_;
        scalar PrPr_;
        scalar kapkap_;
        scalar CpCppart_;

scalar TP2_;
scalar TP1_;
vector posP1_;
scalar dP1_;
vector UP1_;
vector posP2_;
scalar dP2_;
vector UP2_;
scalar tInjStart_;
scalar tInjEnd_;
////////////////////



    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidParticleCloud(const solidParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const solidParticleCloud&);


public:

    // Constructors

        //- Construct given mesh
        solidParticleCloud
        (
            const fvMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );


    // Member Functions

        // Access

            virtual bool hasWallImpactDistance() const;

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mu() const;
            inline scalar mudrop() const;
            inline scalar sigmad() const;
            inline scalar PrPr() const;
            inline scalar kapkap() const;
            inline scalar CpCppart() const;
                                       // inline scalar T() const;
////////////////////////////////////////////////2way coupling//////////////////



inline  vectorField& smom();
inline  const vectorField& smom() const;
inline  tmp<volVectorField> momentumSource() const;


inline  scalarField& sTemp();
inline  const scalarField& sTemp() const;
inline  tmp<volScalarField> sTempsTemp() const;



                inline  scalarField& correctalpha1();
                inline  const scalarField& correctalpha1() const;
                inline  tmp<volScalarField> AddTOAlpha() const;
                inline  vectorField& correctU();
                inline  const vectorField& correctU() const;
                inline  tmp<volVectorField> AddToU() const;
                inline  scalarField& sourcecorrectT();  //Temp

                inline  const scalarField& sourcecorrectT() const; //Temp
                inline  tmp<volScalarField> sourceAddTOT() const; //Temp
               //inline  scalarField& sourcecorrectT();  //Temp
               // inline  const scalarField& correctT() const; //Temp
               // inline  tmp<volScalarField> AddTOT() const; //Temp

                inline  scalarField& sourcecorrectalpha1();
                inline  const scalarField& sourcecorrectalpha1() const;
                inline  tmp<volScalarField> sourceAddTOAlpha() const;
                inline  vectorField& sourcecorrectU();
                inline  const vectorField& sourcecorrectU() const;
                inline  tmp<volVectorField> sourceAddToU() const;





//////////////////////////////////////////////

        // Edit

            //- Move the particles under the influence of the given
            //  gravitational acceleration
            void move(const dimensionedVector& g);
void move2(const dimensionedVector& g);
//////////////////////////////////////////////////////////////
//- Inject particles according to the dictionnary particleProperties
void inject(solidParticle::trackingData &td);
//void removeAL(solidParticle::trackingData &td);
///////////////////////////////////////////////////////////////////////////////
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
