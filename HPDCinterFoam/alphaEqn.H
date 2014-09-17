{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*interface.nHatf());
///////////////////////////////////////////////////////////////////////////////////

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
             
        );

          //   MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);
//////////////////////////////////////////////////////////////////////////////////////
      volScalarField divU = fvc::div(phi);
        volVectorField GRAD = fvc::grad(alpha1);
        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
           (divU*alpha1)-(divU*alpha1)

            //U & GRAD
          //+ vDotcAlphal
        );



//////////////////////////////////////////////////////////////////////////////////////
              MULES::implicitSolve
                        (
                            geometricOneField(),
                            alpha1,
                            phi,
                            phiAlpha,
                            AlSu, //AlSp  implicit (when source term is positive)  //mass transfer
                            Su,// AlSp,  //AlSu  explicit(when source tesm is neg) // divergence term
                            1,
                            0
                        );


////////////////////////////////////////////////////////////////////////////////////

              rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2 ;
//-------------------Modified--------------------//
              rhoPhiCpf = phiAlpha*(rho1*cp1 - rho2*cp2) + phi*rho2*cp2;
//----------------------End-----------------------//

//////////////////////////////////////////////////////
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
