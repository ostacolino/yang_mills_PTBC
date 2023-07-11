#ifndef SU2_H
#define SU2_H

#include"macro.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>

//
// An Su2 matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
//
typedef struct Su2 {
     double comp[4] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2;


typedef struct Su2Adj {
     double comp[9] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2Adj;


inline void init_Su2(Su2 * restrict A, double vec[4])
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=vec[0];
  A->comp[1]=vec[1];
  A->comp[2]=vec[2];
  A->comp[3]=vec[3];
  }


// A=1
inline void one_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=1.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=0
inline void zero_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=0.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=B
inline void equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  A->comp[0]=B->comp[0];
  A->comp[1]=B->comp[1];
  A->comp[2]=B->comp[2];
  A->comp[3]=B->comp[3];
  }


// A=B^{dag}
inline void equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  A->comp[0]= B->comp[0];
  A->comp[1]=-B->comp[1];
  A->comp[2]=-B->comp[2];
  A->comp[3]=-B->comp[3];
  }


// A+=B
inline void plus_equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  A->comp[0]+=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A+=B^{dag}
inline void plus_equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  A->comp[0]+=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=B
inline void minus_equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  A->comp[0]-=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=(r*B)
inline void minus_equal_times_real_Su2(Su2 * restrict A, Su2 const * const restrict B, double r)
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

  A->comp[0]-=(r*B->comp[0]);
  A->comp[1]-=(r*B->comp[1]);
  A->comp[2]-=(r*B->comp[2]);
  A->comp[3]-=(r*B->comp[3]);
  }


// A-=B^{dag}
inline void minus_equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  A->comp[0]-=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A=b*B+c*C
inline void lin_comb_Su2(Su2 * restrict A,
                  double b, Su2 const * const restrict B,
                  double c, Su2 const * const restrict C)
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

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] + c*C->comp[1];
  A->comp[2]= b*B->comp[2] + c*C->comp[2];
  A->comp[3]= b*B->comp[3] + c*C->comp[3];
  }


// A=b*B^{dag}+c*C
inline void lin_comb_dag1_Su2(Su2 * restrict A,
                       double b, Su2 const * const restrict B,
                       double c, Su2 const * const restrict C)
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

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] + c*C->comp[1];
  A->comp[2]= -b*B->comp[2] + c*C->comp[2];
  A->comp[3]= -b*B->comp[3] + c*C->comp[3];
  }


// A=b*B+c*C^{dag}
inline void lin_comb_dag2_Su2(Su2 * restrict A,
                       double b, Su2 const * const restrict B,
                       double c, Su2 const * const restrict C)
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

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] - c*C->comp[1];
  A->comp[2]= b*B->comp[2] - c*C->comp[2];
  A->comp[3]= b*B->comp[3] - c*C->comp[3];
  }


// A=b*B^{dag}+c*C^{dag}
inline void lin_comb_dag12_Su2(Su2 * restrict A,
                        double b, Su2 const * const restrict B,
                        double c, Su2 const * const restrict C)
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

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] - c*C->comp[1];
  A->comp[2]= -b*B->comp[2] - c*C->comp[2];
  A->comp[3]= -b*B->comp[3] - c*C->comp[3];
  }


// A*=r
inline void times_equal_real_Su2(Su2 * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]*=r;
  A->comp[1]*=r;
  A->comp[2]*=r;
  A->comp[3]*=r;
  }


// A*=r
inline void times_equal_complex_Su2(Su2 * restrict A, double complex r)
  {
  #ifdef DEBUG
  if(fabs(cimag(r))>MIN_VALUE)
    {
    fprintf(stderr, "Trying to multiply SU(2) matrix by a non-real number (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]*=creal(r);
  A->comp[1]*=creal(r);
  A->comp[2]*=creal(r);
  A->comp[3]*=creal(r);
  }


// A*=B
inline void times_equal_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  register double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];

  A->comp[0]= a0*B->comp[0] - a1*B->comp[1] - a2*B->comp[2] - a3*B->comp[3];
  A->comp[1]= a0*B->comp[1] + a1*B->comp[0] - a2*B->comp[3] + a3*B->comp[2];
  A->comp[2]= a0*B->comp[2] + a2*B->comp[0] + a1*B->comp[3] - a3*B->comp[1];
  A->comp[3]= a0*B->comp[3] + a3*B->comp[0] - a1*B->comp[2] + a2*B->comp[1];
  }


// A*=B^{dag}
inline void times_equal_dag_Su2(Su2 * restrict A, Su2 const * const restrict B)
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

  register double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];

  A->comp[0]=  a0*B->comp[0] + a1*B->comp[1] + a2*B->comp[2] + a3*B->comp[3];
  A->comp[1]= -a0*B->comp[1] + a1*B->comp[0] + a2*B->comp[3] - a3*B->comp[2];
  A->comp[2]= -a0*B->comp[2] + a2*B->comp[0] - a1*B->comp[3] + a3*B->comp[1];
  A->comp[3]= -a0*B->comp[3] + a3*B->comp[0] + a1*B->comp[2] - a2*B->comp[1];
  }


// A=B*C
inline void times_Su2(Su2 * restrict A,
               Su2 const * const restrict B,
               Su2 const * const restrict C)
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

  A->comp[0]= B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C
inline void times_dag1_Su2(Su2 * restrict A,
                    Su2 const * const restrict B,
                    Su2 const * const restrict C)
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

  A->comp[0]= B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B*C^{dag}
inline void times_dag2_Su2(Su2 * restrict A,
                    Su2 const * const restrict B,
                    Su2 const * const restrict C)
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

  A->comp[0]=  B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C^{dag}
inline void times_dag12_Su2(Su2 * restrict A,
                     Su2 const * const restrict B,
                     Su2 const * const restrict C)
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

  A->comp[0]=  B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// random SU(2) matrix
void rand_matrix_Su2(Su2 *A);


// random SU(2) matrix with p0 given (used in the update)
void rand_matrix_p0_Su2(double p0, Su2 *A);


// sqrt of the determinant
inline double sqrtdet_Su2(Su2 const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  return sqrt(A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1]
             + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3]);
  }


// l2 norm of the matrix
inline double norm_Su2(Su2 const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  return sqrtdet_Su2(A);
  }


// real part of the trace /2
inline double retr_Su2(Su2 const * const restrict A)
  {
  return A->comp[0];
  }


// imaginary part of the trace /2
inline double imtr_Su2(Su2 const * const restrict A)
  {
  (void) A; // to suppress compilation warning of unused variable
  return 0.0;
  }


// unitarize the matrix
inline void unitarize_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double p;

  p=A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1] + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3];
  p=1.0/sqrt(p);

  A->comp[0]*=p;
  A->comp[1]*=p;
  A->comp[2]*=p;
  A->comp[3]*=p;
  }


// traceless antihermitian part
inline void ta_Su2(Su2 * restrict A)
  {
  A->comp[0]=0;
  }


// exponential of the traceless antihermitian part
inline void taexp_Su2(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double v1, v2, v3, norm, s;

  // comp[0] is neglected since we consider the ta part
  norm=(A->comp[1]*A->comp[1]);
  norm+=(A->comp[2]*A->comp[2]);
  norm+=(A->comp[3]*A->comp[3]);
  norm=sqrt(norm);

  v1=A->comp[1]/norm;
  v2=A->comp[2]/norm;
  v3=A->comp[3]/norm;

  s=sin(norm);

  A->comp[0]=cos(norm);
  A->comp[1]=v1*s;
  A->comp[2]=v2*s;
  A->comp[3]=v3*s;
  }


// print on screen
void print_on_screen_Su2(Su2 const * const A);


// print on file
void print_on_file_Su2(FILE *fp, Su2 const * const A);


// print on binary file without changing endiannes
void print_on_binary_file_noswap_Su2(FILE *fp, Su2 const * const A);


// print on binary file changing endiannes
void print_on_binary_file_swap_Su2(FILE *fp, Su2 const * const A);


// print on binary file in big endian format
void print_on_binary_file_bigen_Su2(FILE *fp, Su2 const * const A);


// read from file
void read_from_file_Su2(FILE *fp, Su2 *A);


// read from binary file without changing endiannes
void read_from_binary_file_noswap_Su2(FILE *fp, Su2 *A);


// read from binary file changing endiannes
void read_from_binary_file_swap_Su2(FILE *fp, Su2 *A);


// read from binary file written in big endian
void read_from_binary_file_bigen_Su2(FILE *fp, Su2 *A);


// initialize tensor product
inline void TensProd_init_Su2(TensProd * restrict TP, Su2 const * const restrict A1, Su2 const * const restrict A2)
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
  double complex aux1[4] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex aux2[4] __attribute__((aligned(DOUBLE_ALIGN)));

  #define m2(X,Y) ((X)*2 + (Y))

  // reconstruct the complex form of the matrices
  aux1[m2(0,0)] = A1->comp[0]+I*A1->comp[3];
  aux1[m2(0,1)] = A1->comp[2]+I*A1->comp[1];
  aux1[m2(1,0)] =-A1->comp[2]+I*A1->comp[1];
  aux1[m2(1,1)] = A1->comp[0]-I*A1->comp[3];

  aux2[m2(0,0)] = A2->comp[0]+I*A2->comp[3];
  aux2[m2(0,1)] = A2->comp[2]+I*A2->comp[1];
  aux2[m2(1,0)] =-A2->comp[2]+I*A2->comp[1];
  aux2[m2(1,1)] = A2->comp[0]-I*A2->comp[3];

  for(i=0; i<2; i++)
     {
     for(j=0; j<2; j++)
        {
        for(k=0; k<2; k++)
           {
           for(l=0; l<2; l++)
              {
              TP->comp[i][j][k][l]=conj(aux1[m2(i,j)])*aux2[m2(k,l)];
              }
           }
        }
     }

  #undef m2
  }


// convert the fundamental representation matrix B to the adjoint representation matrix A
inline void fund_to_adj_Su2(Su2Adj * restrict A, Su2 const * const restrict B)
  {
  A->comp[0]= (pow(B->comp[0],2)+ pow(B->comp[1],2)-pow(B->comp[2],2)-pow(B->comp[3],2));
  A->comp[1]= 2*(B->comp[0]*B->comp[3]+B->comp[1]*B->comp[2]);
  A->comp[2]= 2*(-B->comp[0]*B->comp[2]+B->comp[1]*B->comp[3]);
  A->comp[3]= 2*(-B->comp[0]*B->comp[3]+B->comp[1]*B->comp[2]);
  A->comp[4]= (pow(B->comp[0],2)- pow(B->comp[1],2)+pow(B->comp[2],2)-pow(B->comp[3],2));
  A->comp[5]= 2*(B->comp[0]*B->comp[1]+B->comp[2]*B->comp[3]);
  A->comp[6]= 2*(B->comp[0]*B->comp[2]+B->comp[1]*B->comp[3]);
  A->comp[7]= 2*(-B->comp[0]*B->comp[1]+B->comp[2]*B->comp[3]);
  A->comp[8]= (pow(B->comp[0],2)- pow(B->comp[1],2)-pow(B->comp[2],2)+pow(B->comp[3],2));
  }


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
inline void TensProdAdj_init_Su2(TensProdAdj * restrict TP, Su2 const * const restrict A1, Su2 const * const restrict A2)
  {
  #ifdef DEBUG
  if(A1==A2)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, k, l;

  #define m2adj(X,Y) ((X)*3 + (Y))

  Su2Adj A1adj  __attribute__((aligned(DOUBLE_ALIGN)));
  Su2Adj A2adj  __attribute__((aligned(DOUBLE_ALIGN)));

  fund_to_adj_Su2(&A1adj, A1);
  fund_to_adj_Su2(&A2adj, A2);

  for(i=0; i<3; i++)
     {
     for(j=0; j<3; j++)
        {
        for(k=0; k<3; k++)
           {
           for(l=0; l<3; l++)
              {
              TP->comp[i][j][k][l]=(A1adj.comp[m2adj(i,j)])*(A2adj.comp[m2adj(k,l)]);
              }
           }
        }
     }

  #undef m2adj
  }

#endif
