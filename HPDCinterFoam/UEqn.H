    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)  //temporal term
      + fvm::div(rhoPhi, U)   //convection
      - fvm::laplacian(muEff, U)  //diffusion term
      - (fvc::grad(U) & fvc::grad(muEff))  //change of mu
     // + (rho1-rho2) * particles.AddToU()
      // +fvm::Sp((rho1-rho2)*fvc::ddt(alpha1), U)   //source term
       // +
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)  //surface tension
                  - ghf*fvc::snGrad(rho) // body force
                  - fvc::snGrad(p_rgh)  //pressure grad
                ) * mesh.magSf()
            )
//momentum source term
            + particles.momentumSource()
            + (rho1-rho2) * particles.sourceAddToU()  //1-2
           //+ momSp
            //  +  (rho1-rho2) *AlSu*U
           // - fvm::Sp((-rho1+rho2)*AlSu,U)     // source term "value" should always be added here without considering sign. negative or positive effects (adding or removing momentum) come from the factors. here a negative factor multiplies by AlSu to make it sink term!  
        );
    }
