///////////////////////////////////////////////////////////////////////////////
for(solidParticleCloud::iterator elmnt = particles.begin();elmnt != particles.end();++elmnt)
{
vector pp= elmnt.position()-wallDist.position();






}




volScalarField yy(wallDist(mesh).y());
vector p= particles.position() 

sacalar distance= mag(p)- yy;
scalar RR= .5 * particles.d();

        if (RR > distance)
         {
          impact=true;
          Info<<"wall impact occurs"<<endl;
         }

////////////////////////////////////////////////////////////////////////////////
    if (impact)
    {
       
  scalar ReNum = rho1*mag(particle.U())*particle.d() / mu1;
  scalar WeNum = rho1*mag(particle.U())*mag(particle.U())*particle.d() / sigma;
  scalar KaNum = (WeNum^(.25))*sqrt(ReNum);
            if (KaNum < 57.7)
            {
              deposition1=true;
            }
            
            else {rebound1=true;}

    } 

