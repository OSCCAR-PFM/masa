

    ////////////////////////
       // Alnegative= 
    //    neg(TLS-(1.91)*pow(VOLVOL,1.0/3.0)) * ((0.52 * pow(TLS,3.0)) / (VOLVOL))
    //  + pos(TLS-(1.91)*pow(VOLVOL,1.0/3.0)) * (0.0);

        // Ndrop= (0.81 * sqrt(turbulence->k()) / TLS) * runTime.deltaT() * ((1.91)*pow(VOLVOL,1.0/3.0) / TLS);


       frequ= pos(alpha1-0.01) * neg(alpha1-.99) * sqrt(rey.component(0))/Lzero;
      // Ndrop= pos(alpha1-0.01) * neg(alpha1-.99) * scalar(100.0) * (sqrt(rey.component(0)) / Lzero) * runTime.deltaT() * (pow(VOLVOL,1.0/3.0) / TLS);
       Ndrop= pos(alpha1-0.01) * neg(alpha1-.99) * (sqrt(rey.component(0)) / Lzero) * runTime.deltaT() * pow((pow(VOLVOL,1.0/3.0) / TLS),2);

     //  fluct44= pos(alpha1-0.01) * neg(alpha1-.99) 

    //    Ndrop= neg(TLS-(1.91)*pow(VOLVOL,1.0/3.0)) * (0.81 * sqrt(rey.component(0)+rey.component(3)+rey.component(5)) / TLS) * runTime.deltaT() * ((1.91)*pow(VOLVOL,1.0/3.0) / TLS);
///////////////////////////////.................type_1...............///////////////////////////////

     //   AlSu = pos(criterion2 - criterion2LowerLimit) * pos(alpha1 - 0.05) * neg(alpha1 - .50) * (Alnegative/runTime.deltaT());

        AlSu =  particles.sourceAddTOAlpha();///runTime.deltaT();  //+ particles.AddTOAlpha2()/runTime.deltaT();
        momSp = particles.sourceAddToU();

///////////////////////////////.................smoothing...............////////////////////////////
 /*AlSp = 
        pos(criterion2 - criterion2LowerLimit)
      * pos(alpha1 - 0.10) 
      * neg(alpha1 - .75) 
      * (Ndrop*(0+Alnegative)/runTime.deltaT());

solve 
( 
            AlSu
          - fvm::laplacian(Lzero*Lzero, AlSu)
==
             AlSp
); 

Info<< "Calculating AlSu field  \n" << endl;
 
*/
/////////////////////////////...............s....................////////////////////////////////////
   

MassAlSu= rho1 * VOLVOL * AlSu;
















