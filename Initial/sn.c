
#include "../paul.h"
#include <stdbool.h>

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r , double densRead, double vrRead ){

   double rho, P, v, X;
   double Eej, Mej, t0, vmax, rhoISM;
   double v0, r0, vr, rhoSunny;
   double npower, deltapower, K, vt, rt;
   double rhoprefactor, rhoOut, rhoIn, Mdot, vwind;
   double fh, mpower, thetah, thetap, kasenA, theta, kasenFactor;
   double Msun = 2.0e33;
   double yr = 365.25*24.0*3600.0; // sec
   double day = 24.0*3600.0;
   bool wind, powerlaw, readrho, readvr, kasen, ejecta;
   double Rgas = 8.314e7; // cgs
   double molarMass = 0.6504; // 63% H, 37% He
   double constTemp = 100.0; // K

   wind     = true;
   powerlaw = true;
   readrho  = false;
   readvr   = false;
   kasen    = false;

   Eej    = 1.0e51; // 0.97e51; // 1.31e51;
   Mej    = 2.5*Msun; // 1.789623e33;
   t0     = 1.0*day; // 26.614*yr;
   vmax   = 1.72e9; // 2.53e9;
   vwind  = 10.0e5;
   Mdot   = 4.0e-5*Msun/yr;
   rhoISM = 6.31e-25; // 5.0e-25; // 1.7e-24;

   // Kasen fit parameters
   fh     = 0.1;
   mpower = 8.0;
   thetah = 30.0;
   thetap = 15.0;
   kasenA = 1.8;
   theta  = 90.0; // 36.7567567568;

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

   kasenFactor = fh+(1.0-fh)*pow(theta/thetah,mpower)/(1.0+pow(theta/thetah,mpower)) * (1.0+kasenA*exp(-pow(theta/thetah-1.0,2.0)/pow(thetap/thetah,2.0)));

   if (wind) {
      rhoISM = Mdot/4.0/3.14159/r/r/vwind;
   }

   if (readrho) {
      if (densRead > 0.0 && r < 2.0*r0) {
         ejecta = true;
      } else {
         ejecta = false;
      }
   } else if (r < r0) {
      ejecta = true;
   } else {
      ejecta = false;
   }

   if (ejecta) {
      X = 1.0;
      if (vrRead>0.0) { // read in density & velocity profile
         rho = densRead;
         v   = vrRead;
      } else if(densRead>0.0) { // read in only density profile
         rho = densRead;
         v   = vr;
      } else if(powerlaw) { // ----- tony broken power law -----
         if( r < rt ){
            rho = rhoIn;
            v = vr;
         } else {
            rho = rhoOut;
            v = vr;
         }
      } else { // ------ sunny gaussian -------
         rho = rhoSunny;
         v = vr;
      }
      if( kasen ) {
         rho *= kasenFactor;
      }
   } else {
      rho = rhoISM;
      v = 0.0;
      X = 0.0;
   }
   
   P = 1.0e-5*rho*vmax*vmax;
   //P = 1.0e4*Rgas/molarMass*constTemp*rho;
 
   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
