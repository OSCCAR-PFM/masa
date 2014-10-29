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
 
            + particles.momentumSource()
            + (rho1-rho2) * particles.sourceAddToU()  //1-2
 
        );
    }
