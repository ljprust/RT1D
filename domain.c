
#include "paul.h"
#include <stdbool.h>
#include <stdio.h>

double get_moment_arm( double , double );
double get_dV( double , double );

int num_diagnostics( void );

void setICparams( struct domain * );
void setHydroParams( struct domain * );
void setRiemannParams( struct domain * );
void setGravityParams( struct domain * );

void setupDomain( struct domain * theDomain ){

   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->point_mass = theDomain->theParList.grav_pointmass;

   theDomain->count_steps = 0;
   theDomain->final_step = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   setICparams( theDomain );
   setHydroParams( theDomain );
   setRiemannParams( theDomain );
   setGravityParams( theDomain );
 
}

void initial( double * , double , double , double ); 
void prim2cons( double * , double * , double , double );
void cons2prim( double * , double * , double , double );
void restart( struct domain * ); 
void calculate_mass( struct domain * );
void calculate_pot( struct domain * );
double get_g( struct cell * );
double get_GMr( struct cell * );
void calc_dp( struct domain * );
void set_wcell( struct domain * );

void setupCells( struct domain * theDomain ){

   int i;
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   // read in density profile
   bool readRho = false;
   bool readVr =  false;
   int startDay = 100;
   int numToRead = 1752;
   double densIn[numToRead];
   double vrIn[numToRead];
   FILE *infile_rho;
   FILE *infile_vr;
   if (readRho) {
      char filename_rho[256];
      //sprintf(filename,"kundu_day%i.txt",startDay);
      sprintf(filename_rho,"midscaled_rho.txt");
      infile_rho = fopen(filename_rho,"r");
      printf("opened %s\n", filename_rho);
      for( i=0 ; i<numToRead ; ++i ){
         fscanf(infile_rho,"%lf",&densIn[i]);
         printf("read density %5.3e\n",densIn[i]);
      }
   }
   if (readVr) {
      char filename_vr[256];
      //sprintf(filename,"kundu_day%i.txt",startDay);
      sprintf(filename_vr,"midscaled_vr.txt");
      infile_vr = fopen(filename_vr,"r");
      printf("opened %s\n", filename_vr);
      for( i=0 ; i<numToRead ; ++i ){
         fscanf(infile_vr,"%lf",&vrIn[i]);
         printf("read velocity %5.3e\n",vrIn[i]);
      }
   }

   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = &(theCells[i]);
      double rp = c->riph;
      double rm = rp - c->dr;
      c->wiph = 0.0; 
      double r = get_moment_arm( rp , rm );
      double dV = get_dV( rp , rm );
      if(readVr && i<numToRead) {
         initial( c->prim , r , densIn[i], vrIn[i] ); 
      } else if(readRho && i<numToRead) {
         initial( c->prim , r , densIn[i], 0.0 );
      } else {
         initial( c->prim , r , 0.0 , 0.0 );
      }
      prim2cons( c->prim , c->cons , 0.0 , dV );
      cons2prim( c->cons , c->prim , 0.0 , dV );
   }

   if(readRho) fclose(infile_rho);
   if(readVr)  fclose(infile_vr);

   int gE = theDomain->theParList.grav_e_mode;
   if( gE ){
      calculate_mass( theDomain );
      if( gE == 2 ) calculate_pot( theDomain );

      for( i=0 ; i<Nr ; ++i ){
         struct cell * c = &(theCells[i]);
         double rp = c->riph;
         double rm = rp - c->dr;
         double dV = get_dV( rp , rm );

         double GMr = 0.0;
         if( gE == 3 ) GMr = get_GMr( c );
         prim2cons( c->prim , c->cons , GMr , dV );
         cons2prim( c->cons , c->prim , GMr , dV );
      }
   }
}

void freeDomain( struct domain * theDomain ){
   free( theDomain->theCells );
}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
   }

   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );
   if( final ) theDomain->final_step = 1;

}

void report( struct domain * );
void snapshot( struct domain * , char * );
void output( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      report( theDomain );
      if( theDomain->rank==0 ) printf("t = %.3e\n",t);
   }

   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
         sprintf(filename,"checkpoint_%04d",n0);
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }

}


