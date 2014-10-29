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
#include "mytwoPhaseMixture.H"
//#include "twoPhaseMixture.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<solidParticle>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        if (debug)
        {
            Info<< "Time = " << mesh_.time().timeName()
                << " trackTime = " << trackTime
                << " tEnd = " << tEnd
                << " steptFraction() = " << stepFraction() << endl;
        }

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        label cellI = cell();

        dt *= trackToFace(position() + dt*U_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

        cellPointWeight cpw(mesh_, position(), cellI, face());
        scalar rhoc = td.rhoInterp().interpolate(cpw);
        vector Uc = td.UInterp().interpolate(cpw);
        scalar nuc = td.nuInterp().interpolate(cpw);
      //  scalar CpCpc = td.CpCpInterp().interpolate(cpw);

       // scalar kappafc = td.kappafInterp().interpolate(cpw);
//kappaf
                                    scalar Tc = td.TInterp().interpolate(cpw);
                              // scalar Td = td.TInterp().interpolate(cpw); ///
        scalar rhop = td.cloud().rhop();
        scalar magUr = mag(Uc - U_);
                           scalar magT = Tc - T_;

        scalar ReFunc = 1.0;
        scalar Re = magUr*d_/nuc;

///////////////////////////////////////////////////////////////////////////////////////

     /*   IOdictionary particleProperties
        (
            IOobject
            (
                "particleProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
*/
             //Info<< "*** " << td.cloud().PrPr()<< endl;
             //  Info<< "*** " << particleProperties.lookup("KaKa") << endl;
 
             scalar Nuss = 2.0 + 0.6 * pow((td.cloud().PrPr()),(1.0/3.0)) * pow(Re,(1.0/2.0));   //?
             //scalar Sdh0= (6 * Nuss *  td.cloud().kapkap() ) / (td.cloud().rhop()* d_ * d_ * td.cloud().CpCppart()) ;  //maghale openfoam 2009
             scalar Sdh0= (6 * Nuss *  td.cloud().kapkap() ) / (td.cloud().rhop()* d_ * d_ * td.cloud().CpCppart()) ;    //maghale ETH shokolat dagh  
  //Info<< "Sdh0 = " << Sdh0 << endl;
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////-------new drag coefficient----////////////////////////////////
       scalar aaaRe = 0.0;
      if (Re > 0.1)
        {
          ReFunc += 18.5 * pow(Re, 3.6) + pow(.5*Re, 11) ;
          aaaRe += .44 *  pow(Re, 0.8)  / (330 + pow(Re, .8));
        }
         
         scalar Dc = (24.0*nuc/d_)*(pow(ReFunc,(.033))+aaaRe)*(3.0/4.0)*(rhoc/(d_*rhop));
 
//////////////////////////////////////////////////////////
 
///////////////////////////////////
                scalar m=rhop*4/3*(3.1415)*pow(d_/2,3);

                vector oldMom=U_*m;  
///////////////////////////////////////////////////////...Dc= 1/tau!!!
                 scalar oldTemp= T_*m;
/////////////////////////////////////
                U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc);
                T_ = (T_ + dt*(Sdh0*Tc));
                //T_ = (T_ + dt*(Sdh0*Tc))/(1.0 + dt*Sdh0);   //ba - bad nabood


////////////////////////////////////
                vector newMom=U_*m;//!?  
                scalar newTemp=T_*m; 
                td.cloud().smom()[cellI] += newMom-oldMom;
                td.cloud().sTemp()[cellI] += newTemp-oldTemp;
 
/////////////////////////////////////

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::solidParticle::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::solidParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::solidParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    vector nw = tetIs.faceTri(mesh_).normal();
    nw /= mag(nw);

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;
//////////////////////////////////////// here I can modify ////////
    for(solidParticleCloud::iterator elmnt = td.cloud().begin();elmnt != td.cloud().end();++elmnt)
        {   
            //Info << "d=" << elmnt().d() << endl;
         //   Info << "U=" << elmnt().U() << endl;  
           // Info << "T=" << elmnt().T() << endl; 
                     label cellI = cell();
                     scalar ReNum= (td.cloud().rhop() * mag(elmnt().U()) * elmnt().d())/td.cloud().mudrop(); 
                     scalar WeNum= (td.cloud().rhop() * mag(elmnt().U()) * mag(elmnt().U()) * elmnt().d())/td.cloud().sigmad();
 

//////////  These two criteria are the same. the only difference is that it is empowered by 1.6  ~~~  K = pow (k', 1.6) !!
                       scalar KNum= pow(ReNum,.4)  *  pow(WeNum,.8);

                    if(KNum < 657 ) //|| KNum <  57.7 )     ///  K number is computed here
                     //if(KNum > 0.001)
                            {
                       // Info<< "particle deposition..." << nl << endl;
                       td.keepParticle = false;

                       td.cloud().correctalpha1()[cellI] = 10 * .52 * pow(elmnt().d(),3);
                       td.cloud().correctU()[cellI] = 0.0 * elmnt().U();

                            } 


        } 
 
///////////////////////////////////////////////////////////7

 


////////////////////////////////////////////////////////////////7       
 
///////////////////////////////////////////////////////////////////
    if (Un > 0)
    {
        U_ -= (1.0 + td.cloud().e())*Un*nw;
    }

    U_ -= td.cloud().mu()*Ut;
}
///////////////////////////////////////////////////////////////////////////////

void Foam::solidParticle::hitPatch
(
    const polyPatch&,
    trackingData& td
)
{
    td.keepParticle = false;
}


void Foam::solidParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}


void Foam::solidParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


Foam::scalar Foam::solidParticle::wallImpactDistance(const vector&) const
{
    return 0.5*d_;
}


// ************************************************************************* //
