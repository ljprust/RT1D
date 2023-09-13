
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double rho, P, v, X;

   double Eej = 1.0e51; // 1.0e51;
   double Mej = 2.5*2.0e33; // 2.0e33;
   double t0  = 4.32e6; // 50 days
   //double rhoISM = 5.0e-25; // 1.6e-24;
   double vmax = 1.72e9;

   double v0 = sqrt(4.0/3.0*Eej/Mej);
   double r0 = vmax*t0;
   double vr = vmax*r/r0;
   double rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);

   double npower = 10.0;
   double deltapower = 1.1;
   double K = (npower-3.0)*(3.0-deltapower)/4.0/3.14159/(npower-deltapower);
   double vt = sqrt((npower-5.0)*(5.0-deltapower)/(npower-3.0)/(3.0-deltapower)*2.0*Eej/Mej);
   double rt = vt*t0;
   double rhoprefactor = K*Mej/rt/rt/rt;
   double rhoOut = rhoprefactor*pow(r/rt,-npower);
   double rhoIn  = rhoprefactor*pow(r/rt,-deltapower);

   double vwind = 10.0e5;
   double Mdot = 4.0e-5*2.0e33/365.25/24.0/3600.0;
   double rhoISM = Mdot/4.0/3.14159/r/r/vwind;
/*
   if( r < r0 ){
      rho = rhoSunny;
      v = vr;
      X = 1.0;
   } else {
      rho = rhoISM;
      v = 0.0;
      X = 0.0;
   }
*/
   if( r < rt ){
      rho = rhoIn;
      v = vr;
      X = 1.0;
   } else if(r < r0){
      rho = rhoOut;
      v = vr;
      X = 1.0;
   } else {
      rho = rhoISM;
      v = 0.0;
      X = 0.0;
   }
   
   P = 1.0e-5*rho*vmax*vmax;
 
   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
