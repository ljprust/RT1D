
#include "paul.h"

double get_moment_arm( double , double );
void initial( double * , double ); 
double get_g( struct cell * );
double get_GMr( struct cell * );
double get_dV( double , double );
void prim2cons( double * , double * , double , double );

void boundary( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int gE = theDomain->theParList.grav_e_mode;
   
   if( rank == size-1 ){
      struct cell * cB = theCells+Nr-1;
      double rp = cB->riph;
      double rm = rp-cB->dr;
      double r = get_moment_arm(rp,rm);
      initial( cB->prim , r );
      double dV = get_dV( rp , rm );

      double GMr = 0.0;
      if( gE == 3 ) GMr = get_GMr( cB );
      prim2cons( cB->prim , cB->cons , GMr , dV ); 

   }
/*
   if( rank == 0 ){
      struct cell * cB = theCells+0;
      struct cell * cP = theCells+1;
      cB->prim[RHO] = cP->prim[RHO];
      cB->prim[PPP] = cP->prim[PPP];
   }
*/

   int WIND_R0 = theDomain->theParList.Wind_BC;
   if( WIND_R0 && rank==0 ){
      double Mdotinner  = theDomain->theParList.Mdot_inner;
      double vwindinner = theDomain->theParList.v_wind_inner;
      struct cell * c4 = theCells;
      double r_min = theDomain->theParList.rmin;
      double rho_wind_inner = Mdotinner/4.0/3.14159/r_min/r_min/vwindinner;
      int q;
      c4->prim[RHO] = rho_wind_inner;
      c4->prim[VRR] = vwindinner;
      c4->prim[PPP] = 1.0e-5*0.5*rho_wind_inner*vwindinner*vwindinner;
      for( q=3 ; q<NUM_Q ; ++q ){
         c4->prim[q] = 0.0;
      }
   }

   int ABSORB_R0 = theDomain->theParList.Absorb_BC;
   if( ABSORB_R0 && rank==0 ){
      struct cell * c3 = theCells+1;
      struct cell * c4 = theCells;
 
      //struct cell * cB = theCells+0;
      //struct cell * cP = theCells+1;
      //if( cP->prim[VRR] < 0.0 ){
         //cB->prim[RHO] = cP->prim[RHO];
         //cB->prim[PPP] = cP->prim[PPP];
         //cB->prim[VRR] = cP->prim[VRR];
      //}
int q;
      //if( c3->prim[VRR] < 0.0 ){
         for( q=0 ; q<NUM_Q ; ++q ){
            c4->prim[q] = c3->prim[q];
         }    
      //c4->prim[AAA] = 0.0;
      //if( c4->prim[VRR] > 0.0 ) c4->prim[VRR] *= .5;
      //c4->prim[VRR] = 0.0;
      //}
   }    

}
