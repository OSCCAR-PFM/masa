label nAlphaCorr(readLabel(pimple.dict().lookup("nAlphaCorr")));

label nAlphaSubCycles(readLabel(pimple.dict().lookup("nAlphaSubCycles")));




if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum(0.0*rhoPhi);

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {

	/*	if (ALFA  = 1)// && alpha1  > 0.05 && criterion2 > 0)
		{
           		
   				#include "alphaEqnA.H"
     			  	rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
		}

		else
		{ */

    	  	  #include "alphaEqn.H"
      		  rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
	//	}




    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

interface.correct();

rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;

//-------------------Modified--------------------//
rhoCp == alpha1*rho1*cp1 + (scalar(1) - alpha1)*rho2*cp2;
 
//----------------------End-----------------------//
