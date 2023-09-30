#define _POSIX_C_SOURCE 200809L // to activate posix_memalign correctly define this macro BEFORE including any library
#define MAX_STRING_LENGTH 250   
#define DOUBLE_ALIGN 32         
#define INT_ALIGN 16            

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// subroutine
double sum_vec_elements(double *, int, int);
void jack_knife_mean(double *, int, int, double *);

int main(int argc, char** argv){

	if (argc!=7){

		printf("This program computes a^5_G(t/a) = (1/L^3)sum_x_space_y_space <q(t/a,x) q(0,y)> through standard binned jack-knife analysis\n\n");
                printf("Compiled from file %s\n", __FILE__);
                printf("Author: Claudio Bonanno    (claudiobonanno93@gmail.com)\n");
                #ifdef __INTEL_COMPILER
                printf("Compiled with icc\n");
                #elif defined( __GNUC__ )
                printf("Compiled with gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
                #endif
                printf("Usage: %s input_file_name output_file_name block_max block_spacing L T\n", argv[0]);
                exit(EXIT_FAILURE);
	}
	else if (strlen(argv[1]) >= MAX_STRING_LENGTH){
		fprintf(stderr, "Input file too long. Increase MAX_STRING_LENGTH macro in %s\n", __FILE__);
		exit(EXIT_FAILURE);
	}
	else if (strlen(argv[2]) >= MAX_STRING_LENGTH){
		fprintf(stderr, "Input file too long. Increase MAX_STRING_LENGTH macro in %s\n", __FILE__);
		exit(EXIT_FAILURE);
	}
	else{

		int L, T, err, B, B_max, B_spacing, num_data=0;
		int *counter;
		char input_file_name[MAX_STRING_LENGTH], output_file_name[MAX_STRING_LENGTH];
		FILE *fp;
		double res[2];
		double *single_line; //singola riga letta da input
		double **Re_sample; //matrice bidimensoionale dell'osservabile 
		double **Im_sample;
		double naive_error, tau, err_tau;

		strcpy(input_file_name, argv[1]);
		strcpy(output_file_name, argv[2]);
		B_max = atoi(argv[3]);
		B_spacing = atoi(argv[4]);
		L = atoi(argv[5]);
		T = atoi(argv[6]);

		fp = fopen(input_file_name, "r");
		if (fp == NULL){

			fprintf(stderr, "Can't open input file %s\n", input_file_name);
			return(EXIT_FAILURE);
		}

		for (int c = getc(fp); c!= EOF; c=getc(fp)){
			if (c == '\n') num_data++;
		}
		fclose(fp);

		// saturate B_max to num_data/2 if B_max > num_data/2
        if ( B_max> ((int)(floor((double)(num_data)/2.0))) ) B_max = 1 + ((int)(floor((double)(num_data)/2.0)));
        else B_max++;


        //allocate memory to store data read from a single line 
        err=posix_memalign((void**) &(single_line), (size_t) DOUBLE_ALIGN, (size_t) 2*T*sizeof(double));
        if(err!=0){
        	fprintf(stderr, "Problems in allocating memory to store single line data\n");
        	exit(EXIT_FAILURE);
        }

        //allocate memory to store counters to count data for each value of delta_t
        err=posix_memalign((void**) &(counter), (size_t) INT_ALIGN, (size_t) (T/2+1)*sizeof(int));
        if(err!=0){
        	fprintf(stderr, "Problems in allocating memory to store counters\n");
        	exit(EXIT_FAILURE);
        } 
        for (int delta_t = 0; delta_t < (T/2+1); delta_t++) counter[delta_t]=0;

        //allocate memory to store data of the WHOLE sample
        // T/2 +1 complex value for each values
        err=posix_memalign((void**) &Re_sample, (size_t) DOUBLE_ALIGN, (size_t)(T/2+1)*sizeof(double));
    	if(err!=0){

    		fprintf(stderr, "Problems in allocating memory to store Re_data sample\n");
            exit(EXIT_FAILURE);
    	}

    	err=posix_memalign((void**) &Im_sample, (size_t) DOUBLE_ALIGN, (size_t)(T/2+1)*sizeof(double));
    	if(err!=0){

    		fprintf(stderr, "Problems in allocating memory to store Im_data sample\n");
            exit(EXIT_FAILURE);
    	}

    	//for delta_t = 0, ... , T/2-1, there are T*num_data complex value,
    	//for delta_t = T/2 there are (T/2)*num_data complex value

    	for(int delta_t=0; delta_t<(T/2);delta_t++){

    		err=posix_memalign((void**) &Re_sample[delta_t], (size_t) DOUBLE_ALIGN, (size_t) T*num_data*sizeof(double));
    		if (err!=0){

    			fprintf(stderr, "Problems in allocating memory to store data simple\n");
    			exit(EXIT_FAILURE);
    		}

    		err=posix_memalign((void**) &Im_sample[delta_t], (size_t) DOUBLE_ALIGN, (size_t) T*num_data*sizeof(double));
    		if (err!=0){

    			fprintf(stderr, "Problems in allocating memory to store data simple\n");
    			exit(EXIT_FAILURE);
    		}
    	}
    	err=posix_memalign((void**) &(Re_sample[(T/2)]), (size_t) DOUBLE_ALIGN, (size_t) (T/2)*num_data*sizeof(double));
        if(err!=0){
        
                fprintf(stderr, "Problems in allocating memory to store data sample\n");
                exit(EXIT_FAILURE);
        	}
        err=posix_memalign((void**) &(Im_sample[(T/2)]), (size_t) DOUBLE_ALIGN, (size_t) (T/2)*num_data*sizeof(double));
        if(err!=0){
        
        		fprintf(stderr, "Problems in allocating memory to store data sample\n");
                exit(EXIT_FAILURE);
            }

        //read data from file
        fp = fopen(input_file_name, "r");
        if (fp== NULL){

        	fprintf(stderr, "Can't open input file %s\n", input_file_name);
        	return(EXIT_FAILURE);
        }

        printf("Reading %d data from %s\n", num_data, input_file_name);
        for (int i = 0; i< num_data; i++){

        	err = 0;
        	for (int t = 0 ; t < 2*T; t++) err+=fscanf(fp, "%lf", &(single_line[t]));
        		if(err!=T){

        			fprintf(stderr, "Problems while reading input file %s: expecting %d field at line %d but %d found\n", input_file_name, 2*T, i+1, err);
        			exit(EXIT_FAILURE);
        		}
        	for (int t1 = 0; t1<T; t1++){
        		for(int t2 = 0; t2<=t1; t2++){

        			int delta_t = t1-t2;
        			if(delta_t > T/2) delta_t = T - delta_t;
        			//Calcolo del correlatore ma con Re_sample e Im_sample
        			//G(t1-t2,k)=<Q(t1,k)Q*(t2,k)>
        			Re_sample[delta_t][counter[delta_t]] = (single_line[2*t1]*single_line[2*t2]+single_line[2*t1+1]*single_line[2*t2+1])/((double)(L*L*L));
        			Im_sample[delta_t][counter[delta_t]] = (single_line[2*t2]*single_line[2*t1+1]-single_line[2*t1]*single_line[2*t2+1])/((double)(L*L*L));
        			counter[delta_t]++;
        		}
        	}
        }
        fclose(fp);
        free(single_line);
        free(counter);

        //compute mean value of sample with binned jack-knife
        //open output file
        fp = fopen(output_file_name, "w");
                if (fp == NULL){
                        
                        fprintf(stderr, "Can't open output file %s\n", output_file_name);
                        return(EXIT_FAILURE);
                }
        

        //ANALYSIS (binned Jack-knife with blocking B)
        //NOTE: effective block is B*T for delta_t < T/2 and B*(T/2) for delta_t = T/2

        //PRIMA LA PARTE REALE

        fprintf(fp, "#B t/a a^5_Re{G(at)} err_a^5_Re{G(at)} tau(G) err_tau(G)\n");


        // delta_t from 0 to (T/2)-1
        for (int delta_t=0; delta_t<(T/2); delta_t++){
        	// compute naive error (B=1)
            B=1;
            jack_knife_mean(Re_sample[delta_t], T*num_data, T*B, res);
            naive_error=res[1];
            fprintf(fp, "%d   %d   %.15lg %.15lg\n", B, delta_t, res[0], res[1]);
            // compute error for B>1
            for (B=2; B<B_max; B++){
                if (B % B_spacing == 0){
                                        
                    jack_knife_mean(Re_sample[delta_t], T*num_data, T*B, res);
                    tau = ( (res[1]/naive_error)*(res[1]/naive_error) - 1.0)/2.0; // tau = {[err(B)/err(B=1)]^2 - 1}/2
                    err_tau = tau * sqrt(((double) 2*B) / ((double) num_data)); // err_tau / tau = sqrt(2B/num_data)
                    fprintf(fp, "%d   %d   %.15lg %.15lg   %.15lg %.15lg\n", B, delta_t, res[0], res[1], tau, err_tau);
                    }
                }
            }

            // delta_t = T/2
            // compute naive error (B=1)
            B=1;
            jack_knife_mean(Re_sample[(T/2)], (T/2)*num_data, (T/2)*B, res);
            naive_error=res[1];
            fprintf(fp, "%d   %d   %.15lg %.15lg\n", B, (T/2), res[0], res[1]);
            // compute error for B>1
            for (B=2; B<B_max; B++){
                if (B % B_spacing == 0){
                                
                    jack_knife_mean(Re_sample[(T/2)], (T/2)*num_data, (T/2)*B, res);
                    tau = ( (res[1]/naive_error)*(res[1]/naive_error) - 1.0)/2.0;     // tau = {[err(B)/err(B=1)]^2 - 1}/2
                    err_tau = tau * sqrt(((double) 2*B) / ((double) num_data)); // err_tau / tau = sqrt(2B/num_data)
                    fprintf(fp, "%d   %d   %.15lg %.15lg   %.15lg %.15lg\n", B, (T/2), res[0], res[1], tau, err_tau);
                }
            }


            //Poi la PARTE IMMAGINARIA

        fprintf(fp, "#B t/a a^5_Im{G(at)} err_a^5_Im{G(at)} tau(G) err_tau(G)\n");


        // delta_t from 0 to (T/2)-1
        for (int delta_t=0; delta_t<(T/2); delta_t++){
        	// compute naive error (B=1)
            B=1;
            jack_knife_mean(Im_sample[delta_t], T*num_data, T*B, res);
            naive_error=res[1];
            fprintf(fp, "%d   %d   %.15lg %.15lg\n", B, delta_t, res[0], res[1]);
            // compute error for B>1
            for (B=2; B<B_max; B++){
                if (B % B_spacing == 0){
                                        
                    jack_knife_mean(Im_sample[delta_t], T*num_data, T*B, res);
                    tau = ( (res[1]/naive_error)*(res[1]/naive_error) - 1.0)/2.0; // tau = {[err(B)/err(B=1)]^2 - 1}/2
                    err_tau = tau * sqrt(((double) 2*B) / ((double) num_data)); // err_tau / tau = sqrt(2B/num_data)
                    fprintf(fp, "%d   %d   %.15lg %.15lg   %.15lg %.15lg\n", B, delta_t, res[0], res[1], tau, err_tau);
                    }
                }
            }

            // delta_t = T/2
            // compute naive error (B=1)
            B=1;
            jack_knife_mean(Im_sample[(T/2)], (T/2)*num_data, (T/2)*B, res);
            naive_error=res[1];
            fprintf(fp, "%d   %d   %.15lg %.15lg\n", B, (T/2), res[0], res[1]);
            // compute error for B>1
            for (B=2; B<B_max; B++){
                if (B % B_spacing == 0){
                                
                    jack_knife_mean(Im_sample[(T/2)], (T/2)*num_data, (T/2)*B, res);
                    tau = ( (res[1]/naive_error)*(res[1]/naive_error) - 1.0)/2.0;     // tau = {[err(B)/err(B=1)]^2 - 1}/2
                    err_tau = tau * sqrt(((double) 2*B) / ((double) num_data)); // err_tau / tau = sqrt(2B/num_data)
                    fprintf(fp, "%d   %d   %.15lg %.15lg   %.15lg %.15lg\n", B, (T/2), res[0], res[1], tau, err_tau);
                }
            }

            fclose(fp);

            //free allocated memory
            for(int delta_t=0; delta_t<(T/2+1); delta_t++){
            	free(Re_sample[delta_t]);
            	free(Im_sample[delta_t]);
            }
            free(Re_sample);
            free(Im_sample);

            //close program
            return(EXIT_SUCCESS);

    }
}

double sum_vec_elements(double *x, int i_end, int i_begin){

        double sum = 0.;
        for (int i = i_begin; i<i_end;i++) sum += x[i];

        return sum;
}

// compute jack-knife mean of vector x with its related jack-knife error at block size B (mean and error stored in res)
void jack_knife_mean(double *x, int num_data, int B, double *res)
{
        // num_data is the total number of data
        // B is the block length
        // k is the number of blocks
        // N is the effective number of data used for the analysis

        int j, k, N, i_begin, err;
        double sum, aux, m, err_m;
        double *E;

        if (num_data%B==0) // block size divides sample size
        {
                k = num_data/B;
                i_begin = 0;
                N = num_data;
        }
        else // discard first few data to make block size divide sample size
        {
                k = (int) floor( ((double) num_data) / ((double) B) );
                i_begin = num_data-k*B;
                N = num_data-i_begin;
        }

        // allocate aux vector
        err=posix_memalign((void**) &(E), (size_t) DOUBLE_ALIGN, (size_t) k*sizeof(double));
        if(err!=0){
                fprintf(stderr, "Problems in allocating aux memory for jack-knife\n");
                exit(EXIT_FAILURE);
        }

        sum = sum_vec_elements(x, num_data, i_begin); // sum of all elements of x
        for (j=0; j<k; j++){
                aux = sum_vec_elements(x, i_begin+(j+1)*B, i_begin+j*B); // sum of elements of x into the j-th block
                E[j] = (sum - aux) / ((double)(N-B)); // <x> on whole sample minus the j-th block
        }
        m = sum_vec_elements(E,k,0) / ((double) k); // jack-knife mean
        res[0] = m;

        err_m = 0.0;
        for (j=0; j<k; j++) err_m += (E[j]-m)*(E[j]-m); // sum of square deviations from m
        err_m *= ((double)(k-1))/((double)k);
        err_m = sqrt(err_m); // jack-knife error
        res[1] = err_m;

        // free aux vector
        free(E);
}

