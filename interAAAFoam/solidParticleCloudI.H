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
namespace Foam
{
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::solidParticleCloud::mesh() const
{
    return mesh_;
}


inline Foam::scalar Foam::solidParticleCloud::rhop() const
{
    return rhop_;
}


inline Foam::scalar Foam::solidParticleCloud::e() const
{
    return e_;
}


inline Foam::scalar Foam::solidParticleCloud::mu() const
{
    return mu_;
}

inline Foam::scalar Foam::solidParticleCloud::mudrop() const
{
    return mudrop_;
}

inline Foam::scalar Foam::solidParticleCloud::sigmad() const
{
    return sigmad_;
}

inline Foam::scalar Foam::solidParticleCloud::PrPr() const
{
    return PrPr_;
}

inline Foam::scalar Foam::solidParticleCloud::kapkap() const
{
    return kapkap_;
}

inline Foam::scalar Foam::solidParticleCloud::CpCppart() const
{
    return CpCppart_;
}


 ///////////////////////////////// injaaa

inline tmp<volVectorField> solidParticleCloud::momentumSource() const
{
	tmp<volVectorField> tsource
	(
		new volVectorField
		(
			IOobject
			(
				"smom",
				mesh_.time().timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh_,
			dimensionedVector
			(
				"zero",
				dimensionSet(1, -2, -2, 0, 0),
				vector::zero
			)
		)
	
	);
	tsource().internalField() = (- smom_/(mesh_.time().deltaT().value()*mesh_.V())) ;//-( rhop_ * sourcecorrectU_);    //sourcecorrectalpha1_*mesh_.time().deltaTValue()*mesh_.V()
	return tsource;
}


inline Foam::vectorField& Foam::solidParticleCloud::smom()
{
  	return smom_;
}

inline const Foam::vectorField& Foam::solidParticleCloud::smom() const
{
	return smom_;
}
///////////////////////////////////TEMP_MAI_2014/////////////////////////
inline tmp<volScalarField> solidParticleCloud::sTempsTemp() const
{
    tmp<volScalarField> ttsource
    (
        new volScalarField
        (       
            IOobject
            (
                "sTemp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(1, -3, -1, 1, 0),
                0
            )
        )
    );
    ttsource().internalField() =   (- sTemp_/(mesh_.time().deltaT().value()*mesh_.V())); 
                                        //kg.k   t                           m3  
    return ttsource;
}


inline Foam::scalarField& Foam::solidParticleCloud::sTemp()
{
  	return sTemp_;
}

inline const Foam::scalarField& Foam::solidParticleCloud::sTemp() const
{
	return sTemp_;
}


///////////////////////////////////////
inline tmp<volScalarField> solidParticleCloud::sourceAddTOAlpha() const
{
    tmp<volScalarField> sourcealphasource
    (
        new volScalarField
        (       
            IOobject
            (
                "sourcecorrectalpha1",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 0, -1, 0, 0),
                0
            )
        )
    );
   // sourcealphasource().internalField() = ( (-1)* sourcecorrectalpha1_)/(mesh_.time().deltaT().value()*mesh_.V()); //passed
   // sourcealphasource().internalField() = ( (-1)* sourcecorrectalpha1_)/(mesh_.time().deltaTValue()*mesh_.V()); 
      sourcealphasource().internalField() =  (-1)* sourcecorrectalpha1_; 
    return sourcealphasource;
}
/////////////////////////////////////////
inline tmp<volVectorField> solidParticleCloud::sourceAddToU() const
{
    tmp<volVectorField> sourceUsource
    (
        new volVectorField
        (
            IOobject
            (
                "sourcecorrectU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "sourcecorrectU",
                dimensionSet(0, 1, -2, 0, 0),   //dimensionSet(0, 1, -1, 0, 0),
                vector(0,0,0)
            )
        )
    );

     // sourceUsource().internalField() = ((-1)* sourcecorrectU_)/(mesh_.time().deltaT().value()*mesh_.V()); //passed
      //sourceUsource().internalField() = ((-1)* sourcecorrectU_)/(mesh_.time().deltaTValue()*mesh_.V());
       sourceUsource().internalField() = (-1)* sourcecorrectU_;
//runTime_.deltaT().value()*
    return sourceUsource;
}
///////////////////////////////////////////
inline tmp<volScalarField> solidParticleCloud::sourceAddTOT() const
{
    tmp<volScalarField> sourceTsource
    (
        new volScalarField
        (
            IOobject
            (
                "sourcecorrectT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0, 0, -1, 1, 0),
                0
            )
        )
    );
//sourceTsource().internalField() = (-1) * (sourcecorrectT_) / (mesh_.time().deltaT().value()*mesh_.V());//passed
//sourceTsource().internalField() = (-1) * (sourcecorrectT_) / (mesh_.time().deltaTValue()*mesh_.V());
        sourceTsource().internalField() = (-1.0) * (sourcecorrectT_/mesh_.time().deltaTValue()) ;
        return sourceTsource;
}


///////////////////////////////////////////LPT---->VOF/////////////////////////
inline tmp<volScalarField> solidParticleCloud::AddTOAlpha() const
{
    tmp<volScalarField> alphasource
    (
        new volScalarField
        (       
            IOobject
            (
                "correctalpha1",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar

            (
                "zero",
                dimensionSet(0, 0, 0, 0, 0),
                0
            )
        )
    );

 
    alphasource().internalField() =   correctalpha1_/mesh_.V();   
    return alphasource;
}
inline tmp<volVectorField> solidParticleCloud::AddToU() const
{
    tmp<volVectorField> Usource
    (
        new volVectorField
        (
            IOobject
            (
                "correctU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "correctU",
                dimensionSet(0, 1, -1, 0, 0),   //dimensionSet(0, 1, -1, 0, 0),
                vector(0,0,0)
            )
        )
    );
   // Usource().internalField() = - correctU_;    //- smom_/(mesh_.time().deltaT().value()*mesh_.V());
     Usource().internalField() =  correctU_  -  sourcecorrectU_; //-  sourcecorrectU_;///(mesh_.time().deltaT().value()*mesh_.V());
    return Usource;
}
///////////////////////////////////////////
//Convert particle to VOF
inline Foam::scalarField& Foam::solidParticleCloud::sourcecorrectalpha1()
{
return sourcecorrectalpha1_;
}
inline const Foam::scalarField& Foam::solidParticleCloud::sourcecorrectalpha1() const
{
return sourcecorrectalpha1_;
}
inline Foam::vectorField& Foam::solidParticleCloud::sourcecorrectU()
{
return sourcecorrectU_;
}
inline const Foam::vectorField& Foam::solidParticleCloud::sourcecorrectU() const
{
return sourcecorrectU_;
}

inline Foam::scalarField& Foam::solidParticleCloud::sourcecorrectT()
{
return sourcecorrectT_;
}
inline const Foam::scalarField& Foam::solidParticleCloud::sourcecorrectT() const
{
return sourcecorrectT_;
}

inline Foam::scalarField& Foam::solidParticleCloud::correctalpha1()
{
return  correctalpha1_;
}
inline const Foam::scalarField& Foam::solidParticleCloud::correctalpha1() const
{
return  correctalpha1_;
}
inline Foam::vectorField& Foam::solidParticleCloud::correctU()
{
return  correctU_;
}
inline const Foam::vectorField& Foam::solidParticleCloud::correctU() const
{
return  correctU_;
}

//////////////////////////////////////////

}//end namespace Foam








// ************************************************************************* //
