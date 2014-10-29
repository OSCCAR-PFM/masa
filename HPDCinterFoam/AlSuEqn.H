

    ////////////////////////
 


       frequ= pos(alpha1-0.01) * neg(alpha1-.99) * sqrt(rey.component(0))/Lzero;
 
       Ndrop= pos(alpha1-0.01) * neg(alpha1-.99) * (sqrt(rey.component(0)) / Lzero) * runTime.deltaT() * pow((pow(VOLVOL,1.0/3.0) / TLS),2);

 

        AlSu =  particles.sourceAddTOAlpha();///runTime.deltaT();  //+ particles.AddTOAlpha2()/runTime.deltaT();
        momSp = particles.sourceAddToU();

 
 

        MassAlSu= rho1 * VOLVOL * AlSu;
















