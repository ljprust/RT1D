
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double rho, P, v, X;

   double Eej = 1.3e51; // 1.0e51;
   double Mej = 1.4*2.0e33; // 2.0e33;
   double t0  = 3.15576e8; // 10 yr // 7.8894e8; // 25 yr
   double rhoISM = 5.0e-25; // 1.6e-24;
   double vmax = 2.0e9;

   double v0 = sqrt(4.0/3.0*Eej/Mej);
   double r0 = vmax*t0;
   double vr = vmax*r/r0;
   double rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);

   if( r < r0 ){
      rho = rhoSunny;
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
