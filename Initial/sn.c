
#include "../paul.h"
#include <stdbool.h>

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double rho, P, v, X;
   double Eej, Mej, t0, vmax, rhoISM;
   double v0, r0, vr, rhoSunny;
   double npower, deltapower, K, vt, rt;
   double rhoprefactor, rhoOut, rhoIn, Mdot, vwind;
   double Msun = 2.0e33;
   double yr = 365.25*24.0*3600.0; // sec
   double day = 24.0*3600.0;
   bool wind, powerlaw;

   wind = true;
   powerlaw = true;

   Eej    = 1.0e51; // 1.0e51;
   Mej    = 2.5*Msun; // 2.0e33;
   t0     = 50.0*day;
   vmax   = 1.72e9;
   vwind  = 10.0e5;
   Mdot   = 4.0e-5*Msun/yr;
   rhoISM = 5.0e-25; // 1.6e-24;

   v0 = sqrt(4.0/3.0*Eej/Mej);
   r0 = vmax*t0;
   vr = vmax*r/r0;
   rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);

   npower = 10.0;
   deltapower = 1.1;
   K = (npower-3.0)*(3.0-deltapower)/4.0/3.14159/(npower-deltapower);
   vt = sqrt((npower-5.0)*(5.0-deltapower)/(npower-3.0)/(3.0-deltapower)*2.0*Eej/Mej);
   rt = vt*t0;
   rhoprefactor = K*Mej/rt/rt/rt;
   rhoOut = rhoprefactor*pow(r/rt,-npower);
   rhoIn  = rhoprefactor*pow(r/rt,-deltapower);

   if (wind) {
      rhoISM = Mdot/4.0/3.14159/r/r/vwind;
   }

   if(powerlaw){
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
   } else {
      if( r < r0 ){
         rho = rhoSunny;
         v = vr;
         X = 1.0;
      } else {
         rho = rhoISM;
         v = 0.0;
         X = 0.0;
      }
   }
   
   P = 1.0e-5*rho*vmax*vmax;
 
   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
