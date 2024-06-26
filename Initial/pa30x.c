
#include "../paul.h"
#include <stdbool.h>

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r , double densRead, double vrRead ){

   double rho, P, v, X;
   double Eej, Mej, t0, rhoISM;
   double vr, vmax;
   double npower, delta;
   double Mdot, vwind, rhoWind;
   double muISM, nISM, amu, A_Chevalier, rsh, Pratio;
   double numerator, denominator, gn, g, Rc;
   double rbr, rho_br, rho_br_in, rho_br_out;
   double Msun = 2.0e33;
   double yr = 365.25*24.0*3600.0; // sec
   double day = 24.0*3600.0;
   bool sunny;
   double prefactor, rhoSunny, v0sq, r0;

   Eej         = 1.0e49; // 0.131e51; // 1.0e48;
   Mej         = 0.35*Msun; // 0.353*Msun; // 0.18*Msun;
   t0          = 10.0*yr;
   vwind       = 1.5e9;
   Mdot        = 1.0e-6*Msun/yr;
   npower      = 6.0;
   delta       = 1.5;
   muISM       = 0.615;
   nISM        = 0.1; // 1.0;
   amu         = 1.66e-24; // g
   A_Chevalier = 2.4;
   rsh         = -1.0; // outer radius of wind
   Pratio      = 1.0e-5;
   sunny       = true;

   rhoISM = muISM*amu*nISM;
   rhoWind = Mdot/4.0/3.14159/r/r/vwind;

   vr = r/t0;

   numerator = 2.0*(5.0-delta)*(npower-5.0)*Eej;
   denominator = (3.0-delta)*(npower-3.0)*Mej;
   gn = 1.0/4.0/3.14159/(npower-delta)
      * pow(numerator,(npower-3.0)/2.0)
      / pow(denominator,(npower-5.0)/2.0);
   g = pow(gn,1.0/npower);

   Rc = pow(A_Chevalier,1.0/6.0) * g / pow(rhoISM,1.0/6.0) * pow(t0,0.5);

   rbr = pow(numerator/denominator,0.5)*t0;
   rho_br = pow(rbr/t0/g,-npower)/t0/t0/t0;
   rho_br_in  = rho_br * pow( r/rbr, -delta  );
   rho_br_out = rho_br * pow( r/rbr, -npower );

   /////// GAUSSIAN STUFF //////////

   prefactor = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5);
   v0sq = 4.0/3.0*Eej/Mej;
   vmax = sqrt( -1.0*log(rhoISM/prefactor*t0*t0*t0)*v0sq );
   r0 = vmax*t0;
   rhoSunny = prefactor /t0/t0/t0 * exp(-vr*vr/v0sq);

   if (sunny) {
      if (r > r0) { // ISM
         X = 0.0;
         v = 0.0;
         rho = rhoISM;
      } else { // gaussian
         X = 1.0;
         v = vr;
         rho = rhoSunny;
      }
   } else {
      if (r > Rc) { // ISM
         X = 0.0;
         v = 0.0;
         rho = rhoISM;
      } else if (r > rbr) { // n power law
         X = 1.0;
         v = vr;
         rho = rho_br_out;
      } else if (r > rsh) { // delta power law
         X = 1.0;
         v = vr;
         rho = rho_br_in;
      } else { // wind
         X = 1.0;
         v = vr;
         rho = rhoWind;
      }
   }

   if (!sunny) vmax = Rc/t0;
   P = Pratio*rho*vmax*vmax;

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
