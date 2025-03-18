
#include "../paul.h"
#include <stdbool.h>
#include <stdio.h>

#define NUMTOREAD 90

static double t0     = 0.0;
static double vmax   = 0.0;
static double rhoISM = 0.0;
static double Rmax   = 0.0;
static double rhoRead[90];
static double vrRead[90];
static double ejectaRead[90];

void setICparams( struct domain * theDomain ){
   t0       = theDomain->theParList.t_initial;
   Rmax     = theDomain->theParList.rmax;
   vmax     = theDomain->theParList.v_max;
   rhoISM   = theDomain->theParList.rho_ISM;

   FILE *infile_rho;
   FILE *infile_vr;
   FILE *infile_ejecta;

   char filename_rho[256];
   char filename_vr[256];
   char filename_ejecta[256];

   sprintf(filename_rho,"rt1dinput_rho.txt");
   sprintf(filename_vr, "rt1dinput_vx.txt" );
   sprintf(filename_ejecta,"rt1dinput_ejecta.txt");

   infile_rho = fopen(filename_rho,"r");
   infile_vr  = fopen(filename_vr,"r");
   infile_ejecta = fopen(filename_ejecta,"r");

   printf("Opened %s\n", filename_rho);
   printf("Opened %s\n", filename_vr);
   printf("Opened %s\n", filename_ejecta);

   int i;
   for( i=0 ; i<NUMTOREAD ; ++i ){
      fscanf(infile_rho,"%lf",&rhoRead[i]);
      fscanf(infile_vr, "%lf",&vrRead[i] );
      fscanf(infile_ejecta,"%lf",&ejectaRead[i]);
   }

   fclose(infile_rho);
   fclose(infile_vr);
   fclose(infile_ejecta);
   printf("Closed input files\n");
}

void initial( double * prim , double r , int debug ){

   double rho, P, v, X, Y, Z;
   double vr;
   int index, i;
   double dist, minDist;

   vr = r/t0;

   index = -1;
   minDist = 2.0*Rmax/t0;
   for( i=0 ; i<NUMTOREAD ; ++i ) {
      dist = vr - vrRead[i];
      if( dist < 0.0 ) dist = -dist;
      if( dist < minDist ) {
         minDist = dist;
         index = i;
      }
   }

   if(debug==1) printf("vr vr_in index: %5.3e %5.3e %d\n",vr,vrRead[index],index);
   if(index<0) printf("NEIGHBOR SEARCH FAILED!!!\n");

   bool ejecta = false;
   if( ejectaRead[index] > 0.5 ) ejecta = true;

   if (ejecta) {
      X   = 1.0;
      Y   = 0.75;
      Z   = 0.66667;
      rho = rhoRead[index];
      v   = vr;
   } else {
      X   = 0.0;
      Y   = 0.25;
      Z   = 0.33333;
      rho = rhoISM;
      v   = 0.0;
   }
   
   P = 1.0e-5*rho*vmax*vmax;
   //P = 1.0e4*Rgas/molarMass*constTemp*rho;
 
   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[YYY] = Y;
   prim[ZZZ] = Z;
   prim[AAA] = 0.0;

}
