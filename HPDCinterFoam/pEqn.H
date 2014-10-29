{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    U = rAU*UEqn.H();
    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    adjustPhi(phiU, U, p_rgh);

    phi = phiU +
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
      - ghf*fvc::snGrad(rho)
    )*rAUf*mesh.magSf();


///////////////////////////////////////

///////////////////////////////////////////

    while (pimple.correctNonOrthogonal())
///////////////////////////////////
    {
 

    		     fvScalarMatrix p_rghEqn
     			   (
       			      fvm::laplacian(rAUf, p_rgh) == fvc::div(phi) //+ ((rho2-rho1)/rho) * AlSu //((rho2-rho1)/rho2) * (AlSu)//source term "value" should always be added here without considering sign. negative or positive effects (adding or removing momentum) come from the factors. here a negative factor multiplies by AlSu to make it sink term!  

       			   );
		 
////////////////////////////////////
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
