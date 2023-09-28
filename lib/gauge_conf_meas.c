﻿#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/function_pointers.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/geometry.h"


#include<time.h> // DEBUG

// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix);
   }


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  long r,
                                  int i,
                                  int j)
   {
   GAUGE_GROUP matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix)+I*imtr(&matrix);
   }


// computation of the plaquette (matrix) in position r and positive directions i,j
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix)
   {
   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ j
//       |   (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> i
//       r   (1)
//

   equal(matrix, &(GC->lattice[r][i]));
   times_equal(matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(matrix, &(GC->lattice[r][j]));
   }


// compute the four-leaf clover in position r, in the plane i,j and save it in M
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M)
   {
   GAUGE_GROUP aux;
   long k, p;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(i >= STDIM || j >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param;
   #endif

   zero(M);

//
//                   i ^
//                     |
//             (14)    |     (3)
//         +-----<-----++-----<-----+
//         |           ||           |
//         |           ||           |
//   (15)  V      (13) ^V (4)       ^ (2)
//         |           ||           |
//         |   (16)    || r   (1)   |
//    p    +----->-----++----->-----+------>   j
//         +-----<-----++-----<-----+
//         |    (9)    ||   (8)     |
//         |           ||           |
//    (10) V      (12) ^V (5)       ^ (7)
//         |           ||           |
//         |           ||           |
//         +------>----++----->-----+
//              (11)   k      (6)
//
   // avanti-avanti
   equal(&aux, &(GC->lattice[r][i]) );                           // 1
   times_equal(&aux, &(GC->lattice[nnp(geo, r, i)][j]) );        // 2
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, r, j)][i]) );    // 3
   times_equal_dag(&aux, &(GC->lattice[r][j]) );                 // 4
   plus_equal(M, &aux);

   k=nnm(geo, r, j);

   // avanti-indietro
   equal_dag(&aux, &(GC->lattice[k][j]) );                       // 5
   times_equal(&aux, &(GC->lattice[k][i]) );                     // 6
   times_equal(&aux, &(GC->lattice[nnp(geo, k, i)][j]) );        // 7
   times_equal_dag(&aux, &(GC->lattice[r][i]) );                 // 8
   plus_equal(M, &aux);

   p=nnm(geo, r, i);

   // indietro-indietro
   equal_dag(&aux, &(GC->lattice[p][i]) );                       // 9
   times_equal_dag(&aux, &(GC->lattice[nnm(geo, k, i)][j]) );    // 10
   times_equal(&aux, &(GC->lattice[nnm(geo, k, i)][i]) );        // 11
   times_equal(&aux, &(GC->lattice[k][j]) );                     // 12
   plus_equal(M, &aux);

   // indietro-avanti
   equal(&aux, &(GC->lattice[r][j]) );                            // 13
   times_equal_dag(&aux, &(GC->lattice[nnp(geo, p, j)][i]) );     // 14
   times_equal_dag(&aux, &(GC->lattice[p][j]) );                  // 15
   times_equal(&aux, &(GC->lattice[p][i]) );                      // 16
   plus_equal(M, &aux);
   }


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt)
   {
   long r;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i, j;
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, param, r, i, j);
         }
     
      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, param, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
     }
   else
     {
     ps=0.0;
     }

   pt*=param->d_inv_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }


// compute the clover discretization of
// sum_{\mu\nu}  Tr(F_{\mu\nu}F_{\mu\nu})/2
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy)
  {
  long r;
  double ris;

  ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<param->d_volume; r++)
     {
     int i, j;
     GAUGE_GROUP aux1, aux2;

     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           clover(GC, geo, param, r, i, j, &aux1);

           ta(&aux1);
           equal(&aux2, &aux1);
           times_equal(&aux1, &aux2);
           ris+=-NCOLOR*retr(&aux1)/16.0;
           }
        }
     }

  *energy=ris*param->d_inv_vol;
  }


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly)
   {
   long rsp;
   double rep, imp;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      rep+=retr(&matrix);
      imp+=imtr(&matrix);
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }


// compute the mean Polyakov loop in the adjoint representation (the trace of)
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double *repoly,
                  double *impoly)
   {
   long rsp;
   double rep, imp;
   double complex tr;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      GAUGE_GROUP matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }
      tr=NCOLOR*retr(&matrix)+NCOLOR*imtr(&matrix)*I;

      #if NCOLOR==1
        rep+=0.0;
      #else
        rep+=(cabs(tr)*cabs(tr)-1)/(NCOLOR*NCOLOR-1);
      #endif

      imp+=0.0;
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }



// compute the mean Polyakov loop and its powers (trace of) in the presence of trace deformation
void polyakov_with_tracedef(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double *repoly,
                            double *impoly)
   {
   long rsp;
   double **rep, **imp;
   int j, err;
   long i;

   for(j=0;j<(int)floor(NCOLOR/2);j++)
      {
      repoly[j]=0.0;
      impoly[j]=0.0;
      }

   err=posix_memalign((void**)&rep, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   err=posix_memalign((void**)&imp, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double*));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(i=0; i<param->d_space_vol; i++)
      {
      err=posix_memalign((void**)&(rep[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      err=posix_memalign((void**)&(imp[i]), (size_t)DOUBLE_ALIGN, (size_t) (int)floor(NCOLOR/2) * sizeof(double));
      if(err!=0)
        {
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      for(j=0; j<(int)floor(NCOLOR/2); j++)
         {
         rep[i][j] = 0.0;
         imp[i][j] = 0.0;
         }
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int k;
      GAUGE_GROUP matrix, matrix2;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(k=0; k<param->d_size[0]; k++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

       rep[rsp][0] = retr(&matrix);
       imp[rsp][0] = imtr(&matrix);

       equal(&matrix2, &matrix);

      for(k=1; k<(int)floor(NCOLOR/2.0); k++)
         {
         times_equal(&matrix2, &matrix);
         rep[rsp][k] = retr(&matrix2);
         imp[rsp][k] = imtr(&matrix2);
         }
      }

    for(j=0; j<(int)floor(NCOLOR/2); j++)
       {
       for(i=0; i<param->d_space_vol; i++)
          {
          repoly[j] += rep[i][j];
          impoly[j] += imp[i][j];
          }
       }

   for(j=0; j<(int)floor(NCOLOR/2.0); j++)
      {
      repoly[j] *= param->d_inv_space_vol;
      impoly[j] *= param->d_inv_space_vol;
      }

   for(i=0; i<param->d_space_vol; i++)
      {
      free(rep[i]);
      free(imp[i]);
      }
   free(rep);
   free(imp);
   }


// compute the local topological charge at point r
// see readme for more details
double loc_topcharge(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     long r)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     (void) GC;
     (void) geo;
     (void) param;
     (void) r;
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;

   #if (STDIM==4 && NCOLOR>1)
     GAUGE_GROUP aux1, aux2, aux3;
     double real1, real2, loc_charge;
     const double chnorm=1.0/(128.0*PI*PI);
     int i, dir[4][3], sign;

     dir[0][0] = 0;
     dir[0][1] = 0;
     dir[0][2] = 0;

     dir[1][0] = 1;
     dir[1][1] = 2;
     dir[1][2] = 3;

     dir[2][0] = 2;
     dir[2][1] = 1;
     dir[2][2] = 1;

     dir[3][0] = 3;
     dir[3][1] = 3;
     dir[3][2] = 2;

     sign=-1;
     loc_charge=0.0;

     for(i=0; i<3; i++)
        {
        clover(GC, geo, param, r, dir[0][i], dir[1][i], &aux1);
        clover(GC, geo, param, r, dir[2][i], dir[3][i], &aux2);

        times_dag2(&aux3, &aux2, &aux1); // aux3=aux2*(aux1^{dag})
        real1=retr(&aux3)*NCOLOR;

        times(&aux3, &aux2, &aux1); // aux3=aux2*aux1
        real2=retr(&aux3)*NCOLOR;

        loc_charge+=((double) sign*(real1-real2));
        sign=-sign;
        }
     ris=(loc_charge*chnorm);
   #endif

   #if (STDIM==2 && NCOLOR==1)
     GAUGE_GROUP u1matrix;
     double angle;

     plaquettep_matrix(GC, geo, param, r, 0, 1, &u1matrix);
     angle=atan2(cimag(u1matrix.comp), creal(u1matrix.comp))/PI2;

     ris=angle;
   #endif

   return ris;
   }


// compute the topological charge
// see readme for more details
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double ris;
   long r;

   ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      ris+=loc_topcharge(GC, geo, param, r);
      }

   return ris;
   }

// chi^\prime = (1/8) int d^4x |x|^2 <q(x)q(0)> = < (1/8) int d^4x |x|^2 q(x) q(0) > = < G2 >
// This function computes the quantity (q(0)/8) sum_{x} d(x,0)^2 q(x) = a^2 G2, whose mean over the ensamble is a^2 chi^\prime
// d(x,y) = lattice distance between sites x and y keeping periodic boundary conditions into account (i.e., the shortest distance between x and y)
double topo_chi_prime(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	double ris=0.0, factor=0.125; // factor = 1/(2D) = 1/8
	long r;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris)
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		double d2 = square_distance(r, 0, param); // d(r,0)^2
		ris += d2 * loc_topcharge(GC, geo, param, r);
	}
	ris *=  loc_topcharge(GC, geo, param, 0) * factor; // ris *= q(0) / 8
	
	return ris;
}

void topcharge_timeslices(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, double *ris, int ncool, FILE *topchar_tprof_filep){
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	long r;
	int N_t = param->d_size[0];

	for (int i=0; i<N_t; i++) ris[i]=0.0;
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris[:N_t])
	#endif
	for(r=0; r<(param->d_volume); r++)
	{
		int t = geo->d_timeslice[r];
		ris[t] += loc_topcharge(GC, geo, param, r);
	}

	fprintf(topchar_tprof_filep, "%ld %d ", GC->update_index, ncool);
	for (int i=0; i<param->d_size[0]; i++) fprintf(topchar_tprof_filep, " %.12g", ris[i]);
	fprintf(topchar_tprof_filep, "\n");
}

void FT_topcharge_timeslices(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, double complex *ris, int ncool, FILE *topchar_MOM_tprof_filep){
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
   {
      fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   long r;
   int N_t = param->d_size[0];
   int N_s = param->d_size[1]; 
   int coordinate[STDIM];     
   int n_i;                    
   double inv_N_s = 1./N_s;
   double fase;

   for(int kappa = 0; kappa<1; kappa++){       
      for (int i=0; i<N_t; i++) {ris[i]=0.0+I*0.0;}
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris[:N_t])
      #endif
      for(r=0; r<(param->d_volume); r++)
      {
         int t = geo->d_timeslice[r];
         si_to_cart(coordinate,r,param);
         n_i = coordinate[param->d_topcharge_MOM_tprof_dir];

         fase = (double) PI2 * n_i * inv_N_s * kappa ;
         ris[t] += loc_topcharge(GC, geo, param, r)*cexp(I*fase);
      }

      fprintf(topchar_MOM_tprof_filep, "%ld %d %d ", GC->update_index, ncool, kappa);
      for (int i=0; i<param->d_size[0]; i++) fprintf(topchar_MOM_tprof_filep, " %.12g %.12g", creal(ris[i]), cimag(ris[i]));
      fprintf(topchar_MOM_tprof_filep, "\n");
   }

   double complex **zp; //dummy array for storing measurements
   int err=posix_memalign((void**)&zp, (size_t)2*DOUBLE_ALIGN, (size_t)N_s/2 * sizeof(complex*));
   if (err!=0){
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
   for(int i = 0; i<N_s/2;i++){
      err=posix_memalign((void**)&zp[i], (size_t)2*DOUBLE_ALIGN, (size_t)N_t*sizeof(complex));
      if (err!=0){
        fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
   }


   for(int kappa = 1; kappa<=(N_s)/2; kappa++){       
      for (int i=0; i<N_t; i++) {ris[i]=0.0+I*0.0;}
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+:ris[:N_t])
      #endif
      for(r=0; r<(param->d_volume); r++)
      {
         int t = geo->d_timeslice[r];
         si_to_cart(coordinate,r,param);
         n_i = coordinate[param->d_topcharge_MOM_tprof_dir];

         fase = (double) PI2 * n_i * inv_N_s * kappa ;
         ris[t] += loc_topcharge(GC, geo, param, r)*cexp(I*fase);
      }

      fprintf(topchar_MOM_tprof_filep, "%ld %d %d ", GC->update_index, ncool, kappa);
      for (int i=0; i<param->d_size[0]; i++){
         fprintf(topchar_MOM_tprof_filep, " %.12g %.12g", creal(ris[i]), cimag(ris[i]));
         zp[kappa-1][i] = ris[i];
      }
      fprintf(topchar_MOM_tprof_filep, "\n");
   }
   for(int kappa = (N_s)/2+1; kappa<N_s; kappa++){       
      

      fprintf(topchar_MOM_tprof_filep, "%ld %d %d ", GC->update_index, ncool, kappa);
      for (int i=0; i<param->d_size[0]; i++) {
         fprintf(topchar_MOM_tprof_filep, " %.12g %.12g", creal(zp[N_s-kappa-1][i]), -cimag(zp[N_s-kappa-1][i]));

      }
      fprintf(topchar_MOM_tprof_filep, "\n");
   }
   free(zp);
}

// compute topological observables (Q, chi_prime, Q(t)) after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topo_obs_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
											 double *chi_prime,
                       double *meanplaq,	double *sum_q_timeslices, FILE *topchar_tprof_filep,
                                          double complex *sum_q_MOM_timeslices, FILE *topchar_MOM_tprof_filep)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

				if (param->d_chi_prime_meas == 1)
				{
					ris=topo_chi_prime(&helperconf, geo, param);
					chi_prime[iter]=ris;
				}
				else
				{
					(void) chi_prime;
				}
			
				if (param->d_topcharge_tprof_meas == 1)	topcharge_timeslices(&helperconf, geo, param, sum_q_timeslices, (iter+1)*param->d_coolsteps, topchar_tprof_filep);
				else
				{
					(void) sum_q_timeslices;
				}

            if (param->d_topcharge_MOM_tprof_meas == 1)   FT_topcharge_timeslices(&helperconf, geo, param, sum_q_MOM_timeslices, (iter+1)*param->d_coolsteps, topchar_MOM_tprof_filep);
            else
            {
               (void) sum_q_MOM_timeslices;
            }

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

		 if (param->d_topcharge_tprof_meas == 1) fflush(topchar_tprof_filep);
		 else (void) topchar_tprof_filep;
       if (param->d_topcharge_MOM_tprof_meas == 1) fflush(topchar_MOM_tprof_filep);
       else (void) topchar_MOM_tprof_filep;
     free_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, ris2, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
		 /*if (param->d_chi_prime_meas == 1) ris2=topo_chi_prime(GC, geo, param);
		 else (void) ris2;*/
     ris2=topo_chi_prime(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
				chi_prime[iter]=ris2;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }

/*---------------------------------------------*/
// OBSERVABLE NEEDED JUST TO CHECK HOW COOLING DESTROYS TOPOLOGICAL CORRELATIONS
void check_correlation_decay_cooling(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, double *ratio)
{
	if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
	{
		fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if(param->d_coolsteps>0)  // if using cooling
	{  
		Gauge_Conf helperconf; 
		double Q, satd;
		init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
		for(int iter=0; iter<(param->d_coolrepeat); iter++)
		{
			cooling(&helperconf, geo, param, param->d_coolsteps);
			Q = fabs(topcharge(&helperconf, geo, param));
			satd = sum_abs_topcharge_dens(&helperconf, geo, param);
			ratio[iter] = (satd-Q)/satd;
		}
		free_gauge_conf(&helperconf, param); 
	}
}

double sum_abs_topcharge_dens(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param)
{
	double sum=0.0;
	for (long r=0; r<(param->d_volume); r++)
	{
		sum += fabs(loc_topcharge(GC, geo, param, r));
	}
	return sum;
}

/*---------------------------------------------*/

// compute the topological charge after some cooling
// in the cooling procedure the action at theta=0 is minimized
void topcharge_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   if(param->d_coolsteps>0)  // if using cooling
     {  
     Gauge_Conf helperconf; 
     double ris, plaqs, plaqt;
     int iter;

     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);
     // helperconf is a copy of the configuration
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        cooling(&helperconf, geo, param, param->d_coolsteps);

        ris=topcharge(&helperconf, geo, param);
        charge[iter]=ris;

        plaquette(&helperconf, geo, param, &plaqs, &plaqt);
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }

     free_gauge_conf(&helperconf, param); 
     }
   else   // no cooling
     {
     double ris, plaqs, plaqt; 
     int iter;

     ris=topcharge(GC, geo, param);
     plaquette(GC, geo, param, &plaqs, &plaqt);
  
     for(iter=0; iter<(param->d_coolrepeat); iter++)
        {
        charge[iter]=ris;
        #if(STDIM==4)
          meanplaq[iter]=0.5*(plaqs+plaqt);
        #else
          meanplaq[iter]=plaqt;
        #endif
        }
     } 
   }


// compute the correlator of the local topological charge
// after "ncool" cooling steps up to spatial distance "dist"
void loc_topcharge_corr(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    int ncool,
                    int dist,
                    double *ris)
   {
   if(!(STDIM==4 && NCOLOR>1) && !(STDIM==2 && NCOLOR==1) )
     {
     fprintf(stderr, "Wrong number of dimensions or number of colors! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   double *topch;
   long r;
   int err, i;

   err=posix_memalign((void**) &(topch), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating memory! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   // compute the local topological charge
   if(ncool>0)
     {
     Gauge_Conf helperconf;

     // helperconf is a copy of GC
     init_gauge_conf_from_gauge_conf(&helperconf, GC, param);

     // cool helperconf
     cooling(&helperconf, geo, param, ncool);

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        topch[r]=loc_topcharge(&helperconf, geo, param, r);
        }

     // free helperconf
     free_gauge_conf(&helperconf, param);
     }
   else
     {
     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r)
     #endif
     for(r=0; r<param->d_volume; r++)
        {
        topch[r]=loc_topcharge(GC, geo, param, r);
        }
     }

   // compute correlators
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(i)
   #endif
   for(i=0; i<dist; i++)
      {
      double aux;
      long r1, r2;
      int j, dir;

      ris[i]=0.0;

      for(r1=0; r1<param->d_volume; r1++)
         {
         aux=0.0;

         for(dir=1; dir<STDIM; dir++)
            {
            r2=r1;
            for(j=0; j<i; j++)
               {
               r2=nnp(geo, r2, dir);
               }
            aux+=topch[r2];
            }
         aux/=(double)(STDIM-1);

         ris[i]+=aux*topch[r1];
         }
      ris[i]*=param->d_inv_vol;
      }

   // free memory
   free(topch);
   }


void perform_measures_localobs(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep, FILE *chiprimefilep, FILE *topchar_tprof_filep, FILE *topchar_MOM_tprof_filep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )
     int i, err;
     double plaqs, plaqt, polyre, polyim, *charge, *chi_prime, *meanplaq, charge_nocooling, chi_prime_nocooling;
		 double *sum_q_timeslices;
       double complex *sum_q_MOM_timeslices;

		 // measure plaquette and polyakov loop
     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);
		 
		 // measure topological charge, chi_prime and topological charge profile at 0 cooling
     charge_nocooling=topcharge(GC, geo, param);
		 if (param->d_chi_prime_meas == 1 ) chi_prime_nocooling=topo_chi_prime(GC, geo, param);

		 // refresh topological charge of periodic replica (only for multicanonic)
		 GC->stored_topo_charge = charge_nocooling;

     fprintf(datafilep, "%ld %.12g %.12g %.12g %.12g %.12g ", GC->update_index, plaqs, plaqt, polyre, polyim, charge_nocooling);
		 if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld 0 %.12lg\n", GC->update_index, chi_prime_nocooling);
	
		 if (param->d_topcharge_tprof_meas == 1)
		 {
			int err=posix_memalign((void**) &(sum_q_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double));
			if(err!=0)
			{
				fprintf(stderr, "Problems in allocating the aux vector for topcharge tprof meas! (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			topcharge_timeslices(GC, geo, param, sum_q_timeslices, 0, topchar_tprof_filep);
  		fflush(topchar_tprof_filep);
		 }
		 else { (void) sum_q_timeslices; }

      if (param->d_topcharge_MOM_tprof_meas == 1)
       {
         int err=posix_memalign((void**) &(sum_q_MOM_timeslices), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0]*sizeof(double complex));
         if(err!=0)
         {
            fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr MOM meas! (%s, %d)\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
         }
         //topcharge_timeslices(GC, geo, param, Resum_q_timeslices, 0, topchar_tprof_filep); 
         FT_topcharge_timeslices(GC, geo, param, sum_q_MOM_timeslices, 0, topchar_MOM_tprof_filep);
      fflush(topchar_MOM_tprof_filep);
       }
       else { (void) sum_q_MOM_timeslices;  }

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_chi_prime_meas == 1)
		 {
     	err=posix_memalign((void**)&chi_prime, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     	if(err!=0)
       	{
       	fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       	exit(EXIT_FAILURE);
       	}
		 }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

		 if (param->d_chi_prime_meas == 1 
       || param->d_topcharge_tprof_meas == 1
       || param->d_topcharge_MOM_tprof_meas == 1) topo_obs_cooling(GC, geo, param, charge, chi_prime, meanplaq, sum_q_timeslices, topchar_tprof_filep, sum_q_MOM_timeslices, topchar_MOM_tprof_filep);
		 else topcharge_cooling(GC, geo, param, charge, meanplaq);

     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
				if (param->d_chi_prime_meas == 1 ) fprintf(chiprimefilep, "%ld %d %.12lg\n", GC->update_index, (i+1)*param->d_coolsteps, chi_prime[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);
		 if (param->d_chi_prime_meas == 1 ) fflush(chiprimefilep);

     free(charge);
		 if (param->d_chi_prime_meas == 1 ) free(chi_prime);
     else 
			{
				(void) chiprimefilep;
			}
		 free(meanplaq);

   #else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);
   #endif
   }


// perform local observables in the case of trace deformation, it computes all the order parameters
void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep)
   {
   #if( (STDIM==4 && NCOLOR>1) || (STDIM==2 && NCOLOR==1) )

     int i, err;
     double plaqs, plaqt, polyre[NCOLOR/2+1], polyim[NCOLOR/2+1]; // +1 just to avoid warning if NCOLOR=1
     double *charge, *meanplaq, charge_nocooling;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov_with_tracedef(GC, geo, param, polyre, polyim);
     charge_nocooling=topcharge(GC, geo, param);

     fprintf(datafilep, "%.12g %.12g ", plaqs, plaqt);

     for(i=0; i<(int)floor(NCOLOR/2); i++)
        {
        fprintf(datafilep, "%.12g %.12g ", polyre[i], polyim[i]);
        }
     fprintf(datafilep, "%.12g ", charge_nocooling);

     err=posix_memalign((void**)&charge, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     err=posix_memalign((void**)&meanplaq, (size_t)DOUBLE_ALIGN, (size_t) param->d_coolrepeat * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

     topcharge_cooling(GC, geo, param, charge, meanplaq);
     for(i=0; i<param->d_coolrepeat; i++)
        {
        fprintf(datafilep, "%.12g %.12g ", charge[i], meanplaq[i]);
        }
     fprintf(datafilep, "\n");

     fflush(datafilep);

     free(charge);
     free(meanplaq);

   #else

     double plaqs, plaqt, polyre, polyim;

     plaquette(GC, geo, param, &plaqs, &plaqt);
     polyakov(GC, geo, param, &polyre, &polyim);

     fprintf(datafilep, "%.12g %.12g %.12g %.12g ", plaqs, plaqt, polyre, polyim);
     fprintf(datafilep, "\n");
     fflush(datafilep);

   #endif
   }


// to optimize the number of hits to be used in multilevel
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double complex poly_corr;
  double poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }
        poly_array[r]=retr(&matrix)+I*imtr(&matrix);
        }

     // average correlator computation
     poly_corr=0.0+I*0.0;
     poly_corr_abs=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr += poly_array[r]*conj(poly_array[r2]);
        poly_corr_abs += cabs(poly_array[r]*conj(poly_array[r2]));
        }
     poly_corr*=param->d_inv_space_vol;
     poly_corr_abs*=param->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1
        poly_corr_fluct+=cabs( poly_array[r]*conj(poly_array[r2]) - poly_corr );
        }
     poly_corr_fluct*=param->d_inv_space_vol;


     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel
void optimize_multilevel_polycorr(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   int i, err;
   long r;
   double complex poly_corr;
   double poly_corr_abs, poly_corr_fluct;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorr(GC,
                       geo,
                       param,
                       param->d_size[0]);
   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct += cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*= sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*= sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g ", poly_corr_abs);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator with the multilevel algorithm
void perform_measures_polycorr(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorr(GC,
                geo,
                param,
                param->d_size[0]);

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorr(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorr(GC, geo, param, datafilep);
   #endif
   }


// to optimize the number of hits to be used in multilevel for the adjoint representation
void optimize_multihit_polycorradj(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep)
  {
  const int max_hit=50;
  const int dir=1;

  int i, mh, t_tmp, err;
  long r, r1, r2;
  double poly_corr, poly_corr_abs, poly_corr_fluct, diff_sec;
  double complex tr;
  double *poly_array;
  time_t time1, time2;
  GAUGE_GROUP matrix, tmp;

  err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef THETA_MODE
   compute_clovers(GC, geo, param, 0);
  #endif

  fprintf(datafilep, "Multihit optimization: \n");
  fprintf(datafilep, "the smaller the value the better the multihit\n");

  for(mh=1; mh<max_hit; mh++)
     {
     time(&time1);

     // polyakov loop in the adjoint representation computation
     for(r=0; r<param->d_space_vol; r++)
        {
        one(&matrix);
        for(i=0; i<param->d_size[0]; i++)
           {
           multihit(GC,
                    geo,
                    param,
                    sisp_and_t_to_si(geo, r, i),
                    0,
                    mh,
                    &tmp);
           times_equal(&matrix, &tmp);
           }

        //trace of the matix in the fundamental representation
        tr=NCOLOR*(retr(&matrix)+I*imtr(&matrix));

        //trace of the matrix in adjoint representation
        poly_array[r]=cabs(tr*conj(tr))-1.0;
        #if NCOLOR != 1
          poly_array[r]/=(NCOLOR*NCOLOR-1);
        #endif
        }

     // average correlator computation
     poly_corr=0.0;
     poly_corr_abs=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr+= poly_array[r]*poly_array[r2];
        poly_corr_abs+=fabs(poly_array[r]*poly_array[r2]);
        }
     poly_corr*=param->d_inv_space_vol;
     poly_corr_abs*=param->d_inv_space_vol;

     // fluctuation of the average correlator computation
     poly_corr_fluct=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        r1=sisp_and_t_to_si(geo, r, 0);
        for(i=0; i<param->d_dist_poly; i++)
           {
           r1=nnp(geo, r1, dir);
           }
        si_to_sisp_and_t(&r2, &t_tmp, geo, r1); // r2 is the spatial value of r1

        poly_corr_fluct+=fabs(poly_array[r]*poly_array[r2]-poly_corr);
        }
     poly_corr_fluct*=param->d_inv_space_vol;

     time(&time2);
     diff_sec = difftime(time2, time1);

     fprintf(datafilep, "%d  %.12g  %.12g (time:%g)\n", mh, poly_corr_abs*sqrt(mh), poly_corr_fluct*sqrt(mh), diff_sec);

     fflush(datafilep);
     }

  free(poly_array);
  }


// to optimize the multilevel (adjoint representation)
void optimize_multilevel_polycorradj(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr, poly_corr_abs, poly_corr_fluct;
   double *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   multilevel_polycorradj(GC,
                          geo,
                          param,
                          param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
         }
      }

   // averages
   poly_corr=0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=fabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuations
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct+=fabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalizations
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "\n");
   fflush(datafilep);

   free(poly_array);
   }


// perform the computation of the polyakov loop correlator in the adjoint representation with the multilevel algorithm
void perform_measures_polycorradj(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep)
   {
   #ifndef OPT_MULTIHIT
   #ifndef OPT_MULTILEVEL
     double ris;
     long r;
     int i;

     multilevel_polycorradj(GC,
                            geo,
                            param,
                            param->d_size[0]);

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProdAdj(&(GC->ml_polycorradj[0][0][r]), &(GC->ml_polycorradj[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProdAdj(&(GC->ml_polycorradj[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   #endif

   #ifdef OPT_MULTIHIT
     optimize_multihit_polycorradj(GC, geo, param, datafilep);
   #endif

   #ifdef OPT_MULTILEVEL
     optimize_multilevel_polycorradj(GC, geo, param, datafilep);
   #endif
   }


// to optimize the multilevel
void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       GParam const * const param,
                                       FILE *datafilep)
   {
   int i, err;
   long r;
   double poly_corr_abs, poly_corr_fluct;
   double complex poly_corr;
   double complex *poly_array;

   err=posix_memalign((void**)&poly_array, (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(double complex));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   fprintf(datafilep, "Multilevel optimization: ");
   fprintf(datafilep, "the smaller the value the better the update\n");

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   // average
   poly_corr=0.0+I*0.0;
   poly_corr_abs=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_array[r]=retr_TensProd(&(GC->ml_polycorr[0][0][r]))+I*imtr_TensProd(&(GC->ml_polycorr[0][0][r]));

      poly_corr+=poly_array[r];
      poly_corr_abs+=cabs(poly_array[r]);
      }
   poly_corr*=param->d_inv_space_vol;
   poly_corr_abs*=param->d_inv_space_vol;

   // fluctuation
   poly_corr_fluct=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      poly_corr_fluct+=cabs(poly_array[r]-poly_corr);
      }
   poly_corr_fluct*=param->d_inv_space_vol;

   // normalization
   for(i=0; i<NLEVELS; i++)
      {
      poly_corr_abs*=sqrt(param->d_ml_upd[i]);
      poly_corr_fluct*=sqrt(param->d_ml_upd[i]);
      }
   poly_corr_abs*=sqrt(param->d_ml_level0_repeat);
   poly_corr_fluct*=sqrt(param->d_ml_level0_repeat);

   poly_corr_abs*=sqrt(param->d_multihit);
   poly_corr_fluct*=sqrt(param->d_multihit);

   fprintf(datafilep, "%.12g %.12g ", poly_corr_abs, poly_corr_fluct);
   for(i=0; i<NLEVELS; i++)
      {
      fprintf(datafilep, "(%d, %d) ", param->d_ml_step[i], param->d_ml_upd[i]);
      }
   fprintf(datafilep, "(1, %d) ", param->d_multihit);
   fprintf(datafilep, "(%d) ", param->d_ml_level0_repeat);
   fprintf(datafilep, "\n");

   fflush(datafilep);

   free(poly_array);
   }


// print the value of the polyakov loop correlator that has been computed by multilevel
void perform_measures_polycorr_long(Gauge_Conf *GC,
                                    GParam const * const param,
                                    FILE *datafilep)
   {
   #ifdef OPT_MULTILEVEL
      optimize_multilevel_polycorr_long(GC, param, datafilep);
   #else
     double ris;
     long r;
     int i;

     for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
        {
        #ifdef OPENMP_MODE
        #pragma omp parallel for num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<param->d_space_vol; r++)
           {
           times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
           }
        }

     ris=0.0;
     for(r=0; r<param->d_space_vol; r++)
        {
        ris+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
        }
     ris*=param->d_inv_space_vol;

     fprintf(datafilep, "%.12g\n", ris);
     fflush(datafilep);
   #endif
   }


// perform the computation of the string width with the
// disconnected correlator using the multilevel algorithm
void perform_measures_tube_disc(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   multilevel_tube_disc(GC,
                        geo,
                        param,
                        param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// perform the computation of the string width with the
// connected correlator using the multilevel algorithm
void perform_measures_tube_conn(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   multilevel_tube_conn(GC,
                        geo,
                        param,
                        param->d_size[0]);

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }


// print the value of the the string width with the
// connected correlator that has been computed by multilevel
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep)
   {
   double risr, risi;
   long r;
   int i;

   for(i=1; i<param->d_size[0]/param->d_ml_step[0]; i++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<param->d_space_vol; r++)
         {
         times_equal_TensProd(&(GC->ml_polycorr[0][0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaq[0][r]), &(GC->ml_polycorr[0][i][r]) );
         times_equal_TensProd(&(GC->ml_polyplaqconn[0][r]), &(GC->ml_polycorr[0][i][r]) );
         }
      }

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polycorr[0][0][r]));
      risi+=imtr_TensProd(&(GC->ml_polycorr[0][0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaq[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaq[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   risr=0.0;
   risi=0.0;
   for(r=0; r<param->d_space_vol; r++)
      {
      risr+=retr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      risi+=imtr_TensProd(&(GC->ml_polyplaqconn[0][r]));
      }
   risr*=param->d_inv_space_vol;
   risi*=param->d_inv_space_vol;
   fprintf(datafilep, "%.12g %.12g ", risr, risi);

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

#endif
