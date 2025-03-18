
#include "../paul.h"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;
static double grav_G = 0.0;
static int USE_RT = 1;
static double rt_A = 0.0;
static double rt_B = 0.0;
static double rt_C = 0.0;
static double rt_D = 0.0;
static double ar = 0.0;
static double mu = 0.0;
static double mProton = 0.0;
static double kB = 0.0;

void setHydroParams( struct domain * theDomain ){
   //GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   USE_RT = theDomain->theParList.rt_flag;
   grav_G = theDomain->theParList.grav_G;
   rt_A = theDomain->theParList.rt_A;
   rt_B = theDomain->theParList.rt_B;
   rt_C = theDomain->theParList.rt_C;
   rt_D = theDomain->theParList.rt_D;

   // constants
   ar      = 7.5646e-15;
   mu      = 2.0;
   mProton = 1.6726e-24;
   kB      = 1.3807e-16;
}

//takes in temperature and density and outputs pressure
double pressureEq(double temp, double rho) {
  return 1.0/3.0*ar*temp*temp*temp*temp + rho*kB*temp/(mu*mProton);
}

//takes in temperature and density and outputs internal energy
double energyEq(double temp, double rho) {
  return 1.5*kB*temp*rho/mu/mProton + ar*temp*temp*temp*temp;
}

//Solves the cubic x^4+4Bx-A^2=0 to get a root of our repressed cubic, formula off of wolfram
double findRootCubic(double A, double B) {
  double z3 = pow(81.0*pow(A,4.0)+768.00*pow(B,3.0),0.5)+9.0*A*A;
  double numerator = pow(2.0,1.0/3.0)*pow(z3, 2.0/3.0)-8.0*pow(3.0,1.0/3.0)*B;
  double denominator = pow(6.0, 2.0/3.0) * pow(z3, 1.0/3.0);
  return numerator/denominator;
}

//takes in pressure and density and solves the quartic analytically to get you temperature
Real calcTemperaturePressure(Real rho, Real pres) {
  Real A, B, temp, y;

  A = 3.0*kB*rho/(ar*mu*mProton);
  B = 3.0*pres/ar;
  y = findRootCubic(A,B);
  temp = pow(y,0.5)*(pow(2.0*A/pow(y*y*y,0.5)-1.0,0.5)-1.0)/2.0;
  return temp;       
} 

//takes in energy and density and solves the quartic analytically to get you temperature
double calcTemperatureEnergy(Real rho, Real energy) {
  double A, B, temp, y;

  A = 3.0*kB*rho/(2.0*ar*mu*mProton);
  B = energy/(ar);
  y = findRootCubic(A,B);
  temp = pow(y,0.5)*(pow(2.0*A/pow(y*y*y,0.5)-1.0,0.5)-1.0)/2.0;
  return temp;
}

//takes in pressure and density and returns gamma
double calcGamma(Real rho, Real pres) {
  double gasPres, beta, temp;

  temp = calcTemperaturePressure(rho, pres);
  gasPres = rho*kB*temp/(mu*mProton);
  beta = gasPres/pres;
  return (32.0-24.0*beta-3.0*beta*beta)/(24.0-21.0*beta);
}

double get_vr( double * prim ){
   return( prim[VRR] );
}

void prim2cons( double * prim , double * cons , double GMr , double dV ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   //double gam = GAMMA_LAW;

   double temp = calcTemperaturePressure(rho, Pp);
   double rhoe = energyEq(temp, rho);

   //double rhoe = Pp/(gam-1.);

   double egrav = -rho*GMr;

   cons[DDD] = rho*dV;
   cons[SRR] = rho*vr*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + egrav)*dV;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      cons[q] = cons[DDD]*prim[q];
   }
}

void cons2prim( double * cons , double * prim , double GMr , double dV ){

   double rho = cons[DDD]/dV;
   double Sr  = cons[SRR]/dV;
   double E   = cons[TAU]/dV;

   double egrav = -rho*GMr;

   double vr = Sr/rho;
   double v2 = vr*vr;
   double rhoe = E - .5*rho*v2 - egrav;
   //double gam = GAMMA_LAW;

   double temp = calcTemperatureEnergy(rho, rhoe);
   double Pp = pressureEq(temp, rho);

   //double Pp = (gam-1.)*rhoe;

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = vr;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void get_cs( double rho , double P ) {
   double gamma1 = calcGamma(rho,P);
   return sqrt(gamma1*P/rho);
}

void getUstar( double * prim , double * Ustar , double Sk , double Ss ){

   double rho = prim[RHO];
   double vr  = prim[VRR];
   double Pp  = prim[PPP];
   double v2  = vr*vr;

   //double gam = GAMMA_LAW;

   double temp = calcTemperaturePressure(rho, Pp);
   double rhoe = energyEq(temp, rho);

   //double rhoe = Pp/(gam-1.);

   double rhostar = rho*(Sk - vr)/(Sk - Ss);
   double Pstar = Pp*(Ss - vr)/(Sk - Ss);
   double Us = rhoe*(Sk - vr)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = rhostar*( Ss );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vr) + Pstar;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void flux( double * prim , double * flux ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double v2  = vr*vr;
   //double gam = GAMMA_LAW;

   double temp = calcTemperaturePressure(rho, Pp);
   double rhoe = energyEq(temp, rho);

   //double rhoe = Pp/(gam-1.);
 
   flux[DDD] = rho*vr;
   flux[SRR] = rho*vr*vr + Pp;
   flux[TAU] = (.5*rho*v2 + rhoe + Pp )*vr;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      flux[q] = flux[DDD]*prim[q];
   }

}

void source( double * prim , double * cons , double rp , double rm , double dVdt ){
   double Pp  = prim[PPP];
   double r  = .5*(rp+rm);
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   cons[SRR] += 2.*Pp*(r/r2)*dVdt;
}

void source_alpha( double * prim , double * cons , double * grad_prim , double r , double dVdt ){

   double A = rt_A;//2e-5/1.7; //2e-5;//1e-4;
   double B = rt_B;//1.2;//0.9;
   double D = rt_D;//0.0;

   //double gam = GAMMA_LAW;
   double alpha = prim[AAA];

   double Pp = prim[PPP];
   double rho = prim[RHO];
   double P1 = grad_prim[PPP];
   double rho1 = grad_prim[RHO];
 
   double g2 = -P1*rho1;
   if( g2 < 0.0 ) g2 = 0.0;
   //double cs = sqrt(gam*fabs(Pp/rho));
   double cs = get_cs(rho, Pp);

   cons[AAA] += ( (A+B*alpha)*sqrt(g2) - D*rho*alpha*cs/r )*dVdt;
   if( cons[AAA] < 0. ) cons[AAA] = 0.;

}

double get_eta( double * prim , double * grad_prim , double r ){

   double C = rt_C;//0.06*1.7;//0.03;
   //double gam = GAMMA_LAW;

   //double cs = sqrt( gam*fabs(prim[PPP]/prim[RHO]) );
   double cs = get_cs(prim[RHO], prim[PPP]);

   double alpha = prim[AAA];
   if( alpha < 0.0 ) alpha = 0.0;

   double u_eddy = cs*sqrt( alpha );
   double lambda = r*sqrt(alpha);
   
   double eta = C*u_eddy*lambda;

   return( eta );

}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss ){
   
   //double gam = GAMMA_LAW;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[VRR];

   //double cs1 = sqrt(fabs(gam*P1/rho1));
   double cs1 = get_cs(rho1, P1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[VRR];

   //double cs2 = sqrt(fabs(gam*P2/rho2));
   double cs2 = get_cs(rho2, P2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;
   
}

double mindt( double * prim , double w , double r , double g , double dr ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   //double gam = GAMMA_LAW;

//   Pp += g*g/8./M_PI/grav_G;

   //double cs = sqrt(fabs(gam*Pp/rho));
   double cs = get_cs(rho, Pp);
   double eta = get_eta( prim , NULL , r );

   double maxvr = cs + fabs( vr - w );
   double dt = dr/maxvr;
   double dt_eta = dr*dr/eta;
   if( dt > dt_eta && USE_RT ) dt = dt_eta;

   return( dt );

}

