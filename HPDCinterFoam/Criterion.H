 //=(1.22*.09*turbulence->k()*sqrt(turbulence->k())/turbulence->epsilon());


//
volSymmTensorField rey(turbulence->R());
surfaceScalarField QEQ(interface.nHatf());




//volVectorField grad1Alpha(fvc::grad(alpha1));
//volVectorField gradAlphaUSE(grad1Alpha/(mag(grad1Alpha)+interface.deltaN()));
 

grad1Alpha=(fvc::grad(alpha1));
gradAlphaUSE=(-grad1Alpha/(mag(grad1Alpha)+interface.deltaN()));
 

//fluctuation1 = rey & (gradAlpha/(mag(gradAlpha)+interface.deltaN()));

volScalarField AAAup (sqrt(rey.component(0)));
volScalarField AAAvp (sqrt(rey.component(3)));
volScalarField AAAwp (sqrt(rey.component(5)));

tensor tx1(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
tensor tx2(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
tensor tx3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

vector vx1(1.0, 0.0, 0.0);
vector vx2(0.0, 1.0, 0.0);
vector vx3(0.0, 0.0, 1.0);

volTensorField fluctensor = AAAup*tx1 + AAAvp*tx2 + AAAwp*tx3 ;

//UFLUK = U + fluctuation1;

fluctuation1 = fluctensor & gradAlphaUSE;     // simon suggests
fluctuation2 = AAAup*vx1 + AAAvp*vx2 + AAAwp*vx3 ;  //simply use

UFLUK =  (U ^ vx1) + 3 *fluctuation1;
 //UFLUK =  U  + 2.5 *fluctuation1;

 


  TLS=  sqrt(           0.5*   10 * twoPhaseProperties.nu() * tr(rey) / turbulence->epsilon()           );   //teylor micro scale
  Lzero=  (1.22*.09*tr(rey)*sqrt(tr(rey))/turbulence->epsilon());    ///integral length sclae


 

 
criterion2 = ((rho1 *   (rey.component(0)+rey.component(3)+rey.component(5)) + 3.5 * rho2* pow(mag(U),2))* pow(TLS,3)  - 1 * interface.sigma()*pow(TLS,2));   //3 , 1  baraye maghale
 
 
 

Vorticity= fvc::curl(U);
