
#include "../paul.h"
#include <stdbool.h>

static double Cfit       = 0.0;
static double Rfit       = 0.0;
static double pfit       = 0.0;
static double qfit       = 0.0;
static double sfit       = 0.0;
static double tfit       = 0.0;
static double Eej_SN     = 0.0;
static double Mej_SN     = 0.0;
static double vmax_SN    = 0.0;
static double vwindinner = 0.0;
static double vwindouter = 0.0;
static double Mdotinner  = 0.0;
static double Mdotouter  = 0.0;
static double rhoISM     = 0.0;
static double t0         = 0.0;
static double pratio     = 0.0;

void setICparams( struct domain * theDomain ){
   Cfit       = theDomain->theParList.C_fit;
   Rfit       = theDomain->theParList.R_fit;
   pfit       = theDomain->theParList.p_fit;
   qfit       = theDomain->theParList.q_fit;
   sfit       = theDomain->theParList.s_fit;
   tfit       = theDomain->theParList.t_fit;
   Eej_SN     = theDomain->theParList.SN_Eej;
   Mej_SN     = theDomain->theParList.SN_Mej;
   vmax_SN    = theDomain->theParList.SN_vmax;
   Mdotinner  = theDomain->theParList.Mdot_inner;
   Mdotouter  = theDomain->theParList.Mdot_outer;
   vwindinner = theDomain->theParList.v_wind_inner;
   //vwindouter = theDomain->theParList.v_wind_outer;
   vwindouter = 1000.0*7.0e10/(3000.0*24.0*3600.0);
   rhoISM     = theDomain->theParList.rho_ISM;
   t0         = theDomain->theParList.T_Start;
   pratio     = theDomain->theParList.Pressure_Ratio;
}

void initial( double * prim , double r , double densRead, double vrRead ){

   double rho, P, v, X;
   double v0, r0, vr, rhoSunny;
   double npower, deltapower, K, vt, rt;
   double rhoprefactor, rhoOut, rhoIn;
   double fh, mpower, thetah, thetap, kasenA, theta, kasenFactor;
   bool readrho, readvr, kasen, ejecta;
   double Rsun = 7.0e10;
   //double Msun = 2.0e33;
   //double yr = 365.25*24.0*3600.0; // sec
   //double day = 24.0*3600.0;
   //double Rgas = 8.314e7; // cgs
   //double molarMass = 0.6504; // 63% H, 37% He
   //double constTemp = 100.0; // K

   //v0 = sqrt(4.0/3.0*Eej_SN/Mej_SN);
   //vr = vmax*r/r0;
   //rhoSunny = pow(3.0/4.0/3.14159, 1.5) * pow(Mej, 2.5)/pow(Eej, 1.5) /t0/t0/t0 * exp(-vr*vr/v0/v0);

   double rho_wind_inner, rho_wind_outer, rho_CEE, rho_SN;
   double r_max;
   bool isInner=false;
   bool isOuter=false;
   bool isSN=false;
   bool isCEE=false;
   bool isISM=false;

   rho_wind_inner = Mdotinner/4.0/3.14159/r/r/vwindinner;
   rho_wind_outer = Mdotouter/4.0/3.14159/r/r/vwindouter;
   rho_CEE = Cfit*pow(r/10.0/Rsun,pfit)*pow(tfit/t0,3.0)*pow(1.0+pow(Rfit/r,qfit),sfit);
  
   r_max = Rfit*pow(-pfit/(pfit-qfit*sfit),-1.0/qfit);

   if( r > r_max ) {
      if(rho_wind_outer>rho_CEE && rho_wind_outer>rhoISM) isOuter=true;
      if(rho_CEE>rho_wind_outer && rho_CEE>rhoISM) isCEE=true;
      if(rhoISM>rho_CEE && rhoISM>rho_wind_outer) isISM=true;
   } else {
      if(rho_wind_inner>rho_CEE) isInner=true;
      if(rho_CEE>rho_wind_inner) isCEE=true;
   }

   if(isCEE) {
      rho = rhoCEE;
      v = r*t0;
      X = 1.0;
   } else if(isInner) {
      rho = rho_wind_inner;
      v = vwindinner;
      X = 0.0;
   } else if(isOuter) {
      rho = rho_wind_outer;
      v = vwindouter;
      X = 0.0;
   } else {
      rho = rhoISM;
      v = 0.0;
      X = 0.0;
   }

   P = pratio*0.5*rho*vmax*vmax;
   //P = 1.0e4*Rgas/molarMass*constTemp*rho;
 
   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
