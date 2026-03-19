enum{VAR_INT,VAR_DOUB,VAR_STR};

#include "paul.h"
#include <string.h>

int readvar( char * filename , char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512];
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 ) return(1);

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( ptr , stringval );
   }

   return(0);
}

int read_par_file( struct domain * theDomain ){

   MPI_Comm_rank(MPI_COMM_WORLD,&(theDomain->rank));
   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));

   int rank = theDomain->rank;
   int size = theDomain->size;

   struct param_list * theList = &( theDomain->theParList );

   char pfile[] = "in.par";

   int err=0;  

   int nrank;
   for( nrank=0 ; nrank<size ; ++nrank ){
      if( rank==nrank ){
         err += readvar( pfile , "Num_R"             , VAR_INT  , &(theList->Num_R)           );
         err += readvar( pfile , "Num_Reports"       , VAR_INT  , &(theList->NumRepts)        );
         err += readvar( pfile , "Num_Snapshots"     , VAR_INT  , &(theList->NumSnaps)        );
         err += readvar( pfile , "Num_Checkpoints"   , VAR_INT  , &(theList->NumChecks)       );
         err += readvar( pfile , "T_Start"           , VAR_DOUB , &(theList->t_min)           );
         err += readvar( pfile , "T_End"             , VAR_DOUB , &(theList->t_max)           );
         err += readvar( pfile , "R_Min"             , VAR_DOUB , &(theList->rmin)            );
         err += readvar( pfile , "R_Max"             , VAR_DOUB , &(theList->rmax)            );
         err += readvar( pfile , "Use_Logtime"       , VAR_INT  , &(theList->Out_LogTime)     );
         err += readvar( pfile , "Log_Zoning"        , VAR_INT  , &(theList->LogZoning)       );
         err += readvar( pfile , "Log_Radius"        , VAR_DOUB , &(theList->LogRadius)       );
         err += readvar( pfile , "CFL"               , VAR_DOUB , &(theList->CFL)             );
         err += readvar( pfile , "PLM"               , VAR_DOUB , &(theList->PLM)             );
         err += readvar( pfile , "Adiabatic_Index"   , VAR_DOUB , &(theList->Adiabatic_Index) );
         err += readvar( pfile , "Density_Floor"     , VAR_DOUB , &(theList->Density_Floor)   );
         err += readvar( pfile , "Pressure_Floor"    , VAR_DOUB , &(theList->Pressure_Floor)  );
         err += readvar( pfile , "Mesh_Motion"       , VAR_INT  , &(theList->Mesh_Motion)     );
         err += readvar( pfile , "Riemann_Solver"    , VAR_INT  , &(theList->Riemann_Solver)  );
         err += readvar( pfile , "Use_RT"            , VAR_INT  , &(theList->rt_flag)    );
         err += readvar( pfile , "Measure_Growth"    , VAR_INT  , &(theList->grow_flag)    );
         err += readvar( pfile , "Use_Logtime"       , VAR_INT  , &(theList->Out_LogTime)     );
         err += readvar( pfile , "Max_Aspect_Short"  , VAR_DOUB , &(theList->MaxShort)        );
         err += readvar( pfile , "Max_Aspect_Long"   , VAR_DOUB , &(theList->MaxLong)         );
         err += readvar( pfile , "RT_A"              , VAR_DOUB , &(theList->rt_A)            );
         err += readvar( pfile , "t_initial"         , VAR_DOUB , &(theList->t_initial)       );
         err += readvar( pfile , "E_ejecta"          , VAR_DOUB , &(theList->E_ejecta)        );
         err += readvar( pfile , "M_ejecta"          , VAR_DOUB , &(theList->M_ejecta)        );
         err += readvar( pfile , "v_max"             , VAR_DOUB , &(theList->v_max)           );
         err += readvar( pfile , "Mdot_wind"         , VAR_DOUB , &(theList->Mdot_wind)       );
         err += readvar( pfile , "v_wind"            , VAR_DOUB , &(theList->v_wind)          );
         err += readvar( pfile , "rho_ISM"           , VAR_DOUB , &(theList->rho_ISM)         );
         err += readvar( pfile , "Use_PowerLaw"      , VAR_INT  , &(theList->Use_PowerLaw)    );
         err += readvar( pfile , "Use_Wind"          , VAR_INT  , &(theList->Use_Wind)        );
         err += readvar( pfile , "C_fit"             , VAR_DOUB , &(theList->C_fit)           );
         err += readvar( pfile , "R_fit"             , VAR_DOUB , &(theList->R_fit)           );
         err += readvar( pfile , "p_fit"             , VAR_DOUB , &(theList->p_fit)           );
         err += readvar( pfile , "q_fit"             , VAR_DOUB , &(theList->q_fit)           );
         err += readvar( pfile , "s_fit"             , VAR_DOUB , &(theList->s_fit)           );
         err += readvar( pfile , "t_fit"             , VAR_DOUB , &(theList->t_fit)           );
         err += readvar( pfile , "SN_Eej"            , VAR_DOUB , &(theList->SN_Eej)          );
         err += readvar( pfile , "SN_Mej"            , VAR_DOUB , &(theList->SN_Mej)          );
         err += readvar( pfile , "SN_vmax"           , VAR_DOUB , &(theList->SN_vmax)         );
         err += readvar( pfile , "Mdot_inner"        , VAR_DOUB , &(theList->Mdot_inner)      );
         err += readvar( pfile , "Mdot_outer"        , VAR_DOUB , &(theList->Mdot_outer)      );
         err += readvar( pfile , "v_wind_inner"      , VAR_DOUB , &(theList->v_wind_inner)    );
         err += readvar( pfile , "v_wind_outer"      , VAR_DOUB , &(theList->v_wind_outer)    );
         err += readvar( pfile , "Wind_BC"           , VAR_INT  , &(theList->Wind_BC)         );
         err += readvar( pfile , "RT_B"              , VAR_DOUB , &(theList->rt_B)            );
         err += readvar( pfile , "RT_C"              , VAR_DOUB , &(theList->rt_C)            );
         err += readvar( pfile , "RT_D"              , VAR_DOUB , &(theList->rt_D)            );
         err += readvar( pfile , "Absorb_Inner"      , VAR_INT  , &(theList->Absorb_BC)       );
         err += readvar( pfile , "Gravity_Switch"    , VAR_INT  , &(theList->grav_flag)       );
         err += readvar( pfile , "Gravity_G"         , VAR_DOUB , &(theList->grav_G)          );
         err += readvar( pfile , "Gravity_Pointmass" , VAR_DOUB , &(theList->grav_pointmass)  );
         err += readvar( pfile , "Gravity_E_Mode"    , VAR_INT  , &(theList->grav_e_mode)     );
         err += readvar( pfile , "Gravity_Balanced"  , VAR_INT  , &(theList->grav_bal)        );
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   int errtot;
   MPI_Allreduce( &err , &errtot , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );

   if( errtot > 0 ){
      printf("Read Failed\n");
      return(1);
   }

   return(0);

}


