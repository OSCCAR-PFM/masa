

    ////////////////////////
 

       frequ= pos(alpha1-0.01) * neg(alpha1-.99) * sqrt(rey.component(0))/Lzero;
 
       Ndrop= pos(alpha1-0.01) * neg(alpha1-.99) * (sqrt(rey.component(0)) / Lzero) * runTime.deltaT() * pow((pow(VOLVOL,1.0/3.0) / TLS),2);
 

        //AlSu = fvc::div(U);
        AlSu = particles.sourceAddTOAlpha(); 
      //  momSp = particles.sourceAddToU();

 

 















