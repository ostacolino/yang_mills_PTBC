#ifndef SUN_H
#define SUN_H

#include"macro.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>

typedef struct SuN {
   double complex comp[NCOLOR*NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
} SuN;
//
//  the element [i][j] can be obtained by matrix.comp[m(i,j)] with m(i,j) defined in macro.h
//

typedef struct SuNAdj {
   #if NCOLOR!=1
     double comp[(NCOLOR*NCOLOR-1)*(NCOLOR*NCOLOR-1)] __attribute__((aligned(DOUBLE_ALIGN)));
   #else // this will never be used, is defined just to avoid warnings
     double comp[1] __attribute__((aligned(DOUBLE_ALIGN)));
   #endif

} SuNAdj;
//
//  the element [i][j] can be obtained by matrix.comp[madj(i,j)] with madj(i,j) defined in macro.h
//


// A=1
inline void one_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }

  for(i=0; i<NCOLOR; i++)
     {
     A->comp[m(i,i)]=1.0+0.0*I;
     }
  }


// A=0
inline void zero_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }
  }


// A=B
inline void equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A=B^{dag}
inline void equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=conj(B->comp[m(j,i)]);
        }
     }
  }


// A+=B
inline void plus_equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A+=B^{dag}
inline void plus_equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]+=conj(B->comp[m(j,i)]);
        }
     }
  }


// A-=B
inline void minus_equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=(r*B)
inline void minus_equal_times_real_SuN(SuN * restrict A, SuN const * const restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=(r*B->comp[i]);
     }
  }


// A-=B^{dag}
inline void minus_equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]-=conj(B->comp[m(j,i)]);
        }
     }
  }


// A=b*B+c*C
inline void lin_comb_SuN(SuN * restrict A,
                  double b, SuN const * const restrict B,
                  double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=b*(B->comp[i])+c*(C->comp[i]);
     }
  }


// A=b*B^{dag}+c*C
inline void lin_comb_dag1_SuN(SuN * restrict A,
                       double b, SuN const * const restrict B,
                       double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*conj(B->comp[m(j,i)])+c*(C->comp[m(i,j)]);
        }
     }
  }


// A=b*B+c*C^{dag}
inline void lin_comb_dag2_SuN(SuN * restrict A,
                       double b, SuN const * const restrict B,
                       double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*(B->comp[m(i,j)])+c*conj(C->comp[m(j,i)]);
        }
     }
  }


// A=b*B^{dag}+c*C^{dag}
inline void lin_comb_dag12_SuN(SuN * restrict A,
                        double b, SuN const * const restrict B,
                        double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*conj(B->comp[m(j,i)])+c*conj(C->comp[m(j,i)]);
        }
     }
  }


// A*=r
inline void times_equal_real_SuN(SuN * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=r
inline void times_equal_complex_SuN(SuN * restrict A, double complex r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=B
inline void times_equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*(B->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A*=B^{dag}
inline void times_equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*conj(B->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C
inline void times_SuN(SuN * restrict A,
                      SuN const * const restrict B,
                      SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C
inline void times_dag1_SuN(SuN * restrict A,
                    SuN const * const restrict B,
                    SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[m(k,i)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C^{dag}
inline void times_dag2_SuN(SuN * restrict A,
                    SuN const * const restrict B,
                    SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*conj(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C^{dag}
inline void times_dag12_SuN(SuN * restrict A,
                     SuN const * const restrict B,
                     SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[m(k,i)])*conj(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// SU(N) random matrix
// generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices
void rand_matrix_SuN(SuN *A);


// generate a matrix in the algebra of SuN with gaussian
// random components in the base T_i such that Tr(T_iT_j)=delta_{ij}
void rand_algebra_gauss_matrix_SuN(SuN *A);


// l2 norm of the matrix
inline double norm_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     aux=cabs(A->comp[i]);
     ris+=aux*aux;
     }
  return sqrt(ris);
  }


// real part of the trace /N
inline double retr_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=creal(tr)/(double)NCOLOR;
  return ris;
  }


// imaginary part of the trace /N
inline double imtr_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=cimag(tr)/(double)NCOLOR;
  return ris;
  }


// LU decomposition with partial pivoting
void LU_SuN(SuN const * const A, SuN *ris, int *sign);


// determinant
inline complex double det_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  #if NCOLOR==3
    complex double ris=0.0+0.0*I;

    ris+=(A->comp[m(0,0)])*(A->comp[m(1,1)])*(A->comp[m(2,2)]);
    ris+=(A->comp[m(1,0)])*(A->comp[m(2,1)])*(A->comp[m(0,2)]);
    ris+=(A->comp[m(2,0)])*(A->comp[m(0,1)])*(A->comp[m(1,2)]);
    ris-=(A->comp[m(2,0)])*(A->comp[m(1,1)])*(A->comp[m(0,2)]);
    ris-=(A->comp[m(1,0)])*(A->comp[m(0,1)])*(A->comp[m(2,2)]);
    ris-=(A->comp[m(0,0)])*(A->comp[m(2,1)])*(A->comp[m(1,2)]);

    return ris;
  #else
    int i;
    double complex ris;
    SuN lu;

    LU_SuN(A, &lu, &i);

    if(i>0)
      {
      ris=1.0+0.0*I;
      }
    else
     {
     ris=-1.0+0.0*I;
     }

    for(i=0; i<NCOLOR; i++)
       {
       ris*=(lu.comp[m(i,i)]);
       }

    return ris;
  #endif
  }


// gives 0 if the matrix is in SU(N) and 1 otherwise
int scheck_SuN(SuN const * const A);


// sunitarize
void unitarize_SuN(SuN *A);


// takes the traceless antihermitian part
inline void ta_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  SuN aux, aux1;
  double complex trace;
  int i;

  equal_SuN(&aux, A);
  equal_dag_SuN(&aux1, A);
  minus_equal_SuN(&aux, &aux1);
  times_equal_real_SuN(&aux, 0.5); // now aux=(A-A^{dag})/2

  trace=aux.comp[m(0,0)];
  for(i=1; i<NCOLOR; i++)
     {
     trace+=aux.comp[m(i,i)];
     }
  trace/=(double)NCOLOR;

  for(i=0; i<NCOLOR; i++)
     {
     aux.comp[m(i,i)]-=trace;
     }

  equal_SuN(A, &aux);
  }


// eponential of the traceless antihermitian part
inline void taexp_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  SuN aux, uno, ris;

  equal_SuN(&aux, A);
  ta_SuN(&aux);

  one_SuN(&uno);

  // now aux is the traceless antihermitian part of the initial matrix
  // and we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.2);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.25);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.33333333333333333333);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.5);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  plus_equal_SuN(&ris, &uno);

  unitarize_SuN(&ris);
  equal_SuN(A, &ris);
  }



// return 0 if matrix is traceless antihermitian, 1 otherwise
inline int ta_check_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double complex aux;
  int i, j, ris;

  ris=0;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux+=(A->comp[m(i,j)]+conj(A->comp[m(j,i)]));
        }
     }
  if(cabs(aux)>MIN_VALUE) ris=1;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     aux+=A->comp[m(i,i)];
     }
  if(cabs(aux)>MIN_VALUE) ris=1;

  return ris;
  }


// exponential of a TA matrix
inline void exp_of_ta_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  // we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  #ifdef DEBUG
  if(ta_check_SuN(A)!=0)
    {
    fprintf(stderr, "Trying to exp. a non TA matrix! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  SuN aux, uno;

  one_SuN(&uno);
  equal_SuN(&aux, A); // in aux the initial matrix is stored

  times_equal_real_SuN(A, 1.0/5.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/4.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/3.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/2.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  plus_equal_SuN(A, &uno);

  unitarize_SuN(A);
  }


// print on screen
void print_on_screen_SuN(SuN const * const A);


// print on file
void print_on_file_SuN(FILE *fp, SuN const * const A);


// print on binary file without changing endiannes
void print_on_binary_file_noswap_SuN(FILE *fp, SuN const * const A);


// print on binary file changing endiannes
void print_on_binary_file_swap_SuN(FILE *fp, SuN const * const A);


// print on binary file in bigendian
void print_on_binary_file_bigen_SuN(FILE *fp, SuN const * const A);


// read from file
void read_from_file_SuN(FILE *fp, SuN *A);


// read from binary file without changing endiannes
void read_from_binary_file_noswap_SuN(FILE *fp, SuN *A);


// read from binary file changing endianness
void read_from_binary_file_swap_SuN(FILE *fp, SuN *A);


// read from binary file written in bigendian
void read_from_binary_file_bigen_SuN(FILE *fp, SuN *A);


// initialize tensor product
inline void TensProd_init_SuN(TensProd * restrict TP, SuN const * const restrict A1, SuN const * const restrict A2)
  {
  #ifdef DEBUG
  if(A1==A2)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(A2->comp), DOUBLE_ALIGN);
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              TP->comp[i][j][k][l]=conj(A1->comp[m(i,j)])*A2->comp[m(k,l)];
              }
           }
        }
     }
  }


// convert the fundamental representation matrix B to the adjoint representation matrix A
inline void fund_to_adj_SuN(SuNAdj * restrict A, SuN const * const restrict B)
  {
  (void) A;
  (void) B;

  fprintf(stderr, "The function fund_to_adj_SuN still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
inline void TensProdAdj_init_SuN(TensProdAdj * restrict TP, SuN const * const restrict A1, SuN const * const restrict A2)
  {
  (void) TP;
  (void) A1;
  (void) A2;

  fprintf(stderr, "The function TensProd_adj_init_SuN still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


#endif
