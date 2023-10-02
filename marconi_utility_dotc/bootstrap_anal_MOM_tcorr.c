/*------------------------------------------------------*/
//                                                      //
//  This code generates a bootstrap sample of the       //
//  time correlator of any observable O read from       //
//  input. If input file is formatted as O(0), O(1),    //
//  O(2), ..., O(T-1), where i=0, ..., N_t-1 is the     //
//  time in lattice units and O(t)=sum_x O(x,t), then   //
//  this program computes G(t) = <O(t1)O(t2)>/N_s^3     //
//  for every time separation t = t1-t2 <= N_t/2.       //
//  Bootstrap is binned with given block size.          //
//                                                      //
/*------------------------------------------------------*/

// Author: Claudio Bonanno (claudiobonanno93@gmail.com)

/*-----------------------------------------------------------*/

#define _POSIX_C_SOURCE 200809L // to activate posix_memalign correctly define this macro BEFORE including stdlib.h

#include <stdio.h>   // standard io functions
#include <stdlib.h>  // define EXIT_FAILURE and EXIT_SUCCESS macros
#include <math.h>    // sqrt function
#include <string.h>  // string manipulation
#include <time.h>    // compute elapsed time
#include <float.h>   // define DBL_EPSILON macro

/*------------ useful macros ------------*/
#define MAX_STRING_LENGTH 500   // max length for unknown string
#define DOUBLE_ALIGN 32         // double type memory alignement for posix_memalign function
#define INT_ALIGN 16            // int type memory alignement for posix_memalign function

/*------------ ran2 macros ------------*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

// rng ran2 state struct
typedef struct RNG_State {
        long idum;
        long idum2;
        long iy;
        long iv[NTAB];
} RNG_State;

/*------------ subroutine prototypes ------------*/
void ran2_std_init_rng_state(RNG_State *, long const);
double random(RNG_State *);
void resampling(double const * const, double *, long const, long const, RNG_State *);
double sum_vec_elements(double const * const, long const, long const);
void jack_knife_mean(double const * const, long const, long const, double *);

/*------------ main ------------*/
int main(int argc, char **argv)
{
        if (argc!=8)
        {
                fprintf(stdout, "This program computes G(t) = <O(t1)O(t2)> / Ns^3 with t=t1-t2 <= N_t/2 and O(t)=sum_x O(x,t) through standard binned bootstrap analisys\n");
                fprintf(stdout, "Compiled from file %s\n", __FILE__);
                fprintf(stdout, "Author: Claudio Bonanno    (claudiobonanno93@gmail.com)\n");
                #ifdef __INTEL_COMPILER
                fprintf(stdout, "Compiled with icc\n");
                #elif defined( __GNUC__ )
                fprintf(stdout, "Compiled with gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
                #endif
    fprintf(stdout, "Error! Program usage: %s input_file_name output_file_name block_size num_resampling random_seed Nt Ns\n", argv[0]);
    exit(EXIT_FAILURE);
        }
        else if (strlen(argv[1]) >= MAX_STRING_LENGTH)
        {
                fprintf(stderr, "Input file name too long. Increase MAX_STRING_LENGTH macro in %s\n", __FILE__);
                exit(EXIT_FAILURE);
        }
        else if (strlen(argv[2]) >= MAX_STRING_LENGTH)
        {
                fprintf(stderr, "Output file name too long. Increase MAX_STRING_LENGTH macro in %s\n", __FILE__);
                exit(EXIT_FAILURE);
        }
        else
        {
                FILE *fp;
                int L, T, err, num_resamplings;
                long seed, num_data=0, B;
                char input_file_name[MAX_STRING_LENGTH], output_file_name[MAX_STRING_LENGTH];
                double **Re_sample, **Im_sample, **sample_resampled;
                double *single_line, *sample_resampled_1, *sample_resampled_2;
                int *counter;
                double res[2];
                RNG_State rng_state;

                strcpy(input_file_name, argv[1]);
                strcpy(output_file_name, argv[2]);
                B=(long)(atoi(argv[3]));
                num_resamplings=atoi(argv[4]);
                seed=(long)(atoi(argv[5]));
                T=atoi(argv[6]);
                L=atoi(argv[7]);
                if (seed>0) seed=-seed;
                ran2_std_init_rng_state(&rng_state, seed);

                // count data in input file
                fp = fopen(input_file_name, "r");
                if (fp == NULL)
                {
                        fprintf(stderr, "Can't open input file %s\n", input_file_name);
                        return(EXIT_FAILURE);
                }

                // it is very important that c is int and not char, otherwise this loop my give wrong results on some machines
                for (int c = getc(fp); c != EOF; c = getc(fp)) // slide the file
                {
                        if (c == '\n') num_data++; // increase data counter when newline is found
                }
                fclose(fp);

                // allocate memory to store data read from a single line
                err=posix_memalign((void**) &(single_line), (size_t) DOUBLE_ALIGN, (size_t) 2*T*sizeof(double));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store single line data \n");
                        exit(EXIT_FAILURE);
                }

                // allocate memory to store data of the whole Re_sample
                // T/2 + 1 different values for the time correlator
                err=posix_memalign((void**) &(Re_sample), (size_t) DOUBLE_ALIGN, (size_t) (T/2+1)*sizeof(double));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store data Re_sample\n");
                        exit(EXIT_FAILURE);
                }
                // for delta_t = 0,...,T/2-1, there are T*num_data data, for delta_t = T/2 there are (T/2)*num_data data
                for(int delta_t=0; delta_t<(T/2); delta_t++)
                {
                        err=posix_memalign((void**) &(Re_sample[delta_t]), (size_t) DOUBLE_ALIGN, (size_t) T*num_data*sizeof(double));
                        if(err!=0)
                        {
                                fprintf(stderr, "Problems in allocating memory to store data Re_sample\n");
                                exit(EXIT_FAILURE);
                        }
                }
                err=posix_memalign((void**) &(Re_sample[(T/2)]), (size_t) DOUBLE_ALIGN, (size_t) (T/2)*num_data*sizeof(double));
                if(err!=0)
                {
        fprintf(stderr, "Problems in allocating memory to store data Re_sample\n");
        exit(EXIT_FAILURE);
                }

                // allocate memory to store counters to count data for each value of delta_t
                err=posix_memalign((void**) &(counter), (size_t) INT_ALIGN, (size_t) (T/2+1)*sizeof(int));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store counters\n");
                        exit(EXIT_FAILURE);
                }
                for (int delta_t=0; delta_t<(T/2+1); delta_t++) counter[delta_t]=0;

                // allocate memory to store data of the whole Im_sample
                // T/2 + 1 different values for the time correlator
                err=posix_memalign((void**) &(Im_sample), (size_t) DOUBLE_ALIGN, (size_t) (T/2+1)*sizeof(double));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store data Im_sample\n");
                        exit(EXIT_FAILURE);
                }
                // for delta_t = 0,...,T/2-1, there are T*num_data data, for delta_t = T/2 there are (T/2)*num_data data
                for(int delta_t=0; delta_t<(T/2); delta_t++)
                {
                        err=posix_memalign((void**) &(Im_sample[delta_t]), (size_t) DOUBLE_ALIGN, (size_t) T*num_data*sizeof(double));
                        if(err!=0)
                        {
                                fprintf(stderr, "Problems in allocating memory to store data Im_sample\n");
                                exit(EXIT_FAILURE);
                        }
                }
                err=posix_memalign((void**) &(Im_sample[(T/2)]), (size_t) DOUBLE_ALIGN, (size_t) (T/2)*num_data*sizeof(double));
                if(err!=0)
                {
        fprintf(stderr, "Problems in allocating memory to store data Im_sample\n");
        exit(EXIT_FAILURE);
                }

                // allocate memory to store counters to count data for each value of delta_t
                err=posix_memalign((void**) &(counter), (size_t) INT_ALIGN, (size_t) (T/2+1)*sizeof(int));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store counters\n");
                        exit(EXIT_FAILURE);
                }
                for (int delta_t=0; delta_t<(T/2+1); delta_t++) counter[delta_t]=0;

                // allocate memory to store resamplings
                err=posix_memalign((void**) &(sample_resampled_1), (size_t) DOUBLE_ALIGN, (size_t) T*num_data*sizeof(double));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store resamplings\n");
                        exit(EXIT_FAILURE);
                }
                for (int k=0; k<(T*num_data); k++) sample_resampled_1[k]=0.0;

                err=posix_memalign((void**) &(sample_resampled_2), (size_t) DOUBLE_ALIGN, (size_t) (T/2)*num_data*sizeof(double));
                if(err!=0)
                {
                        fprintf(stderr, "Problems in allocating memory to store resamplings\n");
                        exit(EXIT_FAILURE);
                }
                for (int k=0; k<((T/2)*num_data); k++) sample_resampled_2[k]=0.0;

                // read data from file
                fp = fopen(input_file_name, "r");
                if (fp == NULL)
                {
                        fprintf(stderr, "Can't open input file %s\n", input_file_name);
                        return(EXIT_FAILURE);
                }

                printf("Reading %ld data from %s\n", num_data, input_file_name);
                for(int i=0; i<num_data; i++) // slide the file
                {
                        err=0;
                        for (int t=0; t<2*T; t++) err+=fscanf(fp, "%lf", &(single_line[t]));
                        if(err!=2*T)
                        {
                                fprintf(stderr, "Problems while reading input file %s: expecting %d field at line %d but %d found\n", input_file_name, T, i+1, err);
                                exit(EXIT_FAILURE);
                        }
                        for(int t1=0; t1<T; t1++)
                        {
                                for(int t2=0; t2<=t1; t2++)
                                {
                                        int delta_t = t1-t2;
                                        if (delta_t > T/2) delta_t = T - delta_t;
                                        Re_sample[delta_t][counter[delta_t]] = (single_line[2*t1]*single_line[2*t2]+single_line[2*t1+1]*single_line[2*t2+1])/((double)(L*L*L));
                                        Im_sample[delta_t][counter[delta_t]] = (single_line[2*t2]*single_line[2*t1+1]-single_line[2*t1]*single_line[2*t2+1])/((double)(L*L*L));
                                        counter[delta_t]++;
                                }
                        }
                }
                fclose(fp);
                free(single_line);

                // open output file
                fp = fopen(output_file_name, "w");
                if (fp == NULL)
                {
                        fprintf(stderr, "Can't open output file %s\n", output_file_name);
                        exit(EXIT_FAILURE);
                }
                fprintf(fp, "# bootstrap_index  t/a   a^5G(t) err_a^5G(t)     (block_size = %ld)\n", B);

                // Real part bootstrap begin
                clock_t start_time=clock();
                for(int i_res=0; i_res<num_resamplings; i_res++)
                {
                        // ANALYSIS (binned jack-knife with blocking B)
                        // NOTE: effective block is B*T for delta_t < T/2 and B*(T/2) for delta_t = T/2
                        // delta_t from 0 to (T/2)-1
                        for (int delta_t=0; delta_t<(T/2); delta_t++)
                        {
                                resampling(Re_sample[delta_t], sample_resampled_1, T*num_data, T*B, &rng_state); // generate binned bootstrap resampling of input data
                                jack_knife_mean(sample_resampled_1, T*num_data, T*B, res);
                                fprintf(fp, "%d   %d   %.15lg %.15lg\n", i_res+1, delta_t, res[0], res[1]);
                        }
                        // delta_t = T/2
                        resampling(Re_sample[(T/2)], sample_resampled_2, (T/2)*num_data, (T/2)*B, &rng_state); // generate binned bootstrap resampling of input data
                        jack_knife_mean(sample_resampled_2, (T/2)*num_data, (T/2)*B, res);
                        fprintf(fp, "%d   %d   %.15lg %.15lg\n", i_res+1, (T/2), res[0], res[1]);
                }
                //fprintf(stdout, "#Analysis time: %.10lf s\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
                // bootstrap end

                // Imaginary bootstrap begin
                //clock_t start_time=clock();
                for(int i_res=0; i_res<num_resamplings; i_res++)
                {
                        // ANALYSIS (binned jack-knife with blocking B)
                        // NOTE: effective block is B*T for delta_t < T/2 and B*(T/2) for delta_t = T/2
                        // delta_t from 0 to (T/2)-1
                        for (int delta_t=0; delta_t<(T/2); delta_t++)
                        {
                                resampling(Im_sample[delta_t], sample_resampled_1, T*num_data, T*B, &rng_state); // generate binned bootstrap resampling of input data
                                jack_knife_mean(sample_resampled_1, T*num_data, T*B, res);
                                fprintf(fp, "%d   %d   %.15lg %.15lg\n", i_res+1, delta_t, res[0], res[1]);
                        }
                        // delta_t = T/2
                        resampling(Im_sample[(T/2)], sample_resampled_2, (T/2)*num_data, (T/2)*B, &rng_state); // generate binned bootstrap resampling of input data
                        jack_knife_mean(sample_resampled_2, (T/2)*num_data, (T/2)*B, res);
                        fprintf(fp, "%d   %d   %.15lg %.15lg\n", i_res+1, (T/2), res[0], res[1]);
                }
                fprintf(stdout, "#Analysis time: %.10lf s\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
                // bootstrap end

                // close output file
                fclose(fp);

                // free allocated memory
                for(int delta_t=0; delta_t<(T/2+1); delta_t++) free(Re_sample[delta_t]);
                for(int delta_t=0; delta_t<(T/2+1); delta_t++) free(Im_sample[delta_t]);
                free(Re_sample);
                free(Im_sample);
                free(counter);
                free(sample_resampled_1);
                free(sample_resampled_2);

                // close program
                fprintf(stdout, "Ciao!\n");
                return(EXIT_SUCCESS);
        }
}

/*---------------------------------------------------------------------------------------------------------------------------------------------------*/

// standard initialization of ran2 RNG state
void ran2_std_init_rng_state(RNG_State *rng_state, long const seed)
{
        long k;
        rng_state->idum=seed;
        if (rng_state->idum > 0) (rng_state->idum) *= -1; // if seed is positive, make it negative
        rng_state->idum2=123456789;
        rng_state->iy=0;
        for (int i=0; i<NTAB; i++) rng_state->iv[i]=0;

        if (rng_state->idum <= 0)
        {
        if (-(rng_state->idum) < 1)
                rng_state->idum=1;
        else
                rng_state->idum = -(rng_state->idum);
        rng_state->idum2=(rng_state->idum);
        for (int j=NTAB+7; j>=0; j--)
        {
                k=(rng_state->idum)/IQ1;
                rng_state->idum=IA1*(rng_state->idum-k*IQ1)-k*IR1;
                if (rng_state->idum < 0)
                        rng_state->idum += IM1;
                if (j < NTAB)
                        rng_state->iv[j] = rng_state->idum;
        }
        rng_state->iy=rng_state->iv[0];
        }
}

double random(RNG_State *rng_state)
{
        int j;
        long k;
        double temp;

        k=(rng_state->idum)/IQ1;
        rng_state->idum=IA1*(rng_state->idum-k*IQ1)-k*IR1;
        if (rng_state->idum < 0)
                rng_state->idum += IM1;
        k=rng_state->idum2/IQ2;
        rng_state->idum2=IA2*(rng_state->idum2-k*IQ2)-k*IR2;
        if (rng_state->idum2 < 0)
                rng_state->idum2 += IM2;
        j=(int)(rng_state->iy/NDIV);
        rng_state->iy=rng_state->iv[j]-rng_state->idum2;
        rng_state->iv[j] = rng_state->idum;
        if (rng_state->iy < 1)
                rng_state->iy += IMM1;

        if ( ( temp = AM * ((double)rng_state->iy) ) > RNMX)
                return RNMX;
        else
                return temp;
}

// generete binned bootstrap resampling of x and store it in x_resampled
void resampling(double const * const x, double *x_resampled, long const num_data, long const B, RNG_State *rng_state)
{
        // num_data is the total number of data
        // B is the block size
        // k is the number of blocks
        // N = k*B is the total size of the sample
        // if B doesn't divide num_data, then start_index = (num_data - N) data are discarded to make B a divisor of the sample size

        long k, start_index, i_resampled=0;
        if (num_data % B == 0)
        {
                k = num_data/B;
                start_index = 0;
        }
        else
        {
                k = ( (int) floor( ( (double) num_data) / ( (double) B) ) );
                start_index = num_data - k*B;
        }
        // start to fill x_resampled from i = start_index
        i_resampled = start_index;

        for (long j=0; j<k; j++) // extract k times block indices with repetitions
        {
                double rand_num = random(rng_state); // random double in (0, 1)
                rand_num *= ((double)(k)); // random double in (0, k)
                long rand_block_index = ( (long) floor(rand_num) ); // random long in [0, k-1]

                // include in the resampling all entries within a block of size B
                for(long i=(start_index+rand_block_index*B); i<(start_index+(rand_block_index+1)*B); i++)
                {
                        x_resampled[i_resampled] = x[i];
                        i_resampled++;
                }
        }
}

// sum of vector elements from index i_begin to index i_end
double sum_vec_elements(double const * const x, long const i_end, long const i_begin)
{
        double sum = 0.0;
        for (long i=i_begin; i<i_end; i++)
                sum += x[i];
        return sum;
}

// compute jack-knife mean of vector x with its related jack-knife error at block size B (mean and error stored in res)
void jack_knife_mean(double const * const x, long const num_data, long const B, double *res)
{
        // num_data is the total number of data
        // B is the block length
        // k is the number of blocks
        // N is the effective number of data used for the analysis

        int err;
        long k, N, i_begin;
        double sum, aux, m, err_m;
        double *E;

        if (num_data % B==0) // block size divides sample size
        {
                k = num_data/B;
                i_begin = 0;
                N = num_data;
        }
        else // discard first few data to make block size divide sample size
        {
                k = (long) (floor( ((double) num_data) / ((double) B) ));
                i_begin = num_data-k*B;
                N = num_data-i_begin;
        }

        // allocate aux vector
        err=posix_memalign((void**) &(E), (size_t) DOUBLE_ALIGN, (size_t) k*sizeof(double));
        if(err!=0)
        {
                fprintf(stderr, "Problems in allocating aux memory for jack-knife\n");
                exit(EXIT_FAILURE);
        }

        sum = sum_vec_elements(x, num_data, i_begin); // sum of all elements of x
        for (long j=0; j<k; j++) // slide the blocks
        {
                aux = sum_vec_elements(x, i_begin+(j+1)*B, i_begin+j*B); // sum of elements of x into the j-th block
                E[j] = (sum - aux) / ((double)(N-B)); // <x> on whole sample minus the j-th block
        }
        m = sum_vec_elements(E,k,0) / ((double) k); // jack-knife mean
        res[0] = m;

        err_m = 0.0;
        for (long j=0; j<k; j++)
                err_m += (E[j]-m)*(E[j]-m); // sum of square deviations from m
        err_m *= ((double)(k-1))/((double)k);
        err_m  = sqrt(err_m); // jack-knife error
        res[1] = err_m;

        // free aux vector
        free(E);
}
