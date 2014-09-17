 Lzero=(1.22*.09*turbulence->k()*sqrt(turbulence->k())/turbulence->epsilon());

 	  criterion = (.66 * rho1 *  turbulence->k()  + 0.5 * rho2* pow(mag(U),2))*pow((1.22*.09*turbulence->k()*sqrt(turbulence->k())/turbulence->epsilon()),3) -interface.sigma()*pow((1.22*.09 *turbulence->k()*sqrt(turbulence->k())/turbulence->epsilon()),2);

//
volSymmTensorField rey(turbulence->R());
surfaceScalarField QEQ(interface.nHatf());

//dynamic_cast <Foam::incompressible::RASModels::myrealizableKE*> (turbulence.ptr())->r2Cmu(fvc::grad(U));

 TLS= sqrt(           0.5*   10 * twoPhaseProperties.nu() * tr(rey) / turbulence->epsilon()           );

	//criterion2 = (rho1 *   rey.component(3) + 2.5 * rho2* pow(mag(U),2))* pow(TLS,3)  - .5 * interface.sigma()*pow(TLS,2);
         // criterion2 = (rho1 *   (rey.component(0)+rey.component(5)) + 3.7 * rho2* pow(mag(U),2))* pow(TLS,3)  - 1.6 * interface.sigma()*pow(TLS,2);

              criterion2 = pos(alpha1 - 0.11)*neg(alpha1 - .85)*(rho1 *   (rey.component(0)+rey.component(5)) + 3.7 * rho2* pow(mag(U),2))* pow(TLS,3)  - 1.6 * interface.sigma()*pow(TLS,2);

