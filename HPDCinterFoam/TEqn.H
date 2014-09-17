
surfaceScalarField kappaf = twoPhaseProperties.kappaf();

//volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + kappaf);

//+ fvc::interpolate(rho*turbulence->nut())

        fvScalarMatrix TEqn
        (
          fvm::ddt(rhoCp, T)
        + fvm::div(rhoPhiCpf, T)
        - fvm::laplacian(kappaf, T) // - fvm::laplacian(kappaEff, T)
        - (rho1-rho2) * CPCPPARTICLE  * particles.sourceAddTOT()
        -  CPCPPARTICLE * particles.sTempsTemp()
        );
//TEqn.relax();
TEqn.solve();

                //Info<< "*** CpCp " <<  CpCp << endl;
