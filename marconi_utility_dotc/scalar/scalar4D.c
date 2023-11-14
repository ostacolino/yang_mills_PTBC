/*

Programma per l'esame di Metodi Numerici che simula un campo scalare neutro

Update: Calcolo del correlatore wall-wall del campo scalare neutro!!!
Questo correlatore mi darà informazioni sulla legge di disperzione della particella, in particolare \sqrt{m**2+\abs(p)**2}

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N_s 20
#define N_t 20
#define N_max 200

#define nlatt N_s*N_s*N_s*N_t
static const double PI=3.141592653589793238462643383279502884197169399375105820974944;
static const double PI2=6.283185307179586476925286766559005768394338798750211641949889;
static const double HALF_PI=1.570796326794896619231321691639751442098584699687552910487472;

double ran2(long *);
void print_lattice(int, int);
void geometry();
void init_lattice(int);
void update_heatbath();
void update_overrelax();
void energy(double *, double *, double *);
void MOM_q_timeslices(double complex *, int, int);



double mass, mass2, mass2p8, sigma2;

double obs1, obs2, obs3;
double complex *sum_q_MOM_timeslices = NULL;

long seed=42;
int npp[N_max][2], nmm[N_max][2];
int iflag, measures, i_decorrel,dir;
double xene_mass, xene_spat, xene_temp;
FILE *finput, *fseed, *fout;

double field[N_t][N_s][N_s][N_s];

int main(){

	
	//fout=fopen("profilo.dat","w");
	//fseed = fopen("seed.txt","r");
	//fscanf(fseed,"%ld",&seed);
	finput=fopen("input.txt","r");
	fscanf(finput, "%d %d %d %d %lf",&iflag, &measures, &i_decorrel, &dir, &mass);
	mass;
	mass2= mass*mass;
	mass2p8 = mass2 + 8.0;
	sigma2 = 1./mass2p8;

	geometry();
	init_lattice(iflag);
	printf("#N_t=%d N_s=%d m=%.12g dir=%d\n",N_t,N_s, mass, dir);
	//print_lattice(N_t,N_s);

	/*int therm = 10000;
	for(int i = 0; i<therm;i++){
		for (int j = 0; j<i_decorrel; j++){

			update_heatbath();
			update_overrelax();
			update_overrelax();
            update_overrelax();
            update_overrelax();
            update_overrelax();

		}

		//energy(&xene_mass, &xene_spat, &xene_temp);
	}*/

	sum_q_MOM_timeslices = malloc(N_t*sizeof(complex));
         		if(!sum_q_MOM_timeslices)
         		{
            		fprintf(stderr, "Problems in allocating the aux vector for topcharge tcorr MOM meas! (%s, %d)\n", __FILE__, __LINE__);
            		exit(EXIT_FAILURE);
         		}
	//MOM_q_timeslices(sum_q_MOM_timeslices,dir,0);

	for (int i = 1; i <= measures; i++){

			
            for (int j =  0; j< i_decorrel; j++){

                    update_heatbath();
                    update_overrelax();
                    update_overrelax();
                    update_overrelax();
                    update_overrelax();
                    update_overrelax();
                }
                

            //energy(&xene_mass, &xene_spat, &xene_temp);
            //printf("%6d %f %f %f\n",i,xene_mass,xene_spat,xene_temp);
            MOM_q_timeslices(sum_q_MOM_timeslices,dir,i);
         	//MOM_q_timeslices(sum_q_MOM_timeslices,2,i);
         	//MOM_q_timeslices(sum_q_MOM_timeslices,3,i);
            
        }

        free(sum_q_MOM_timeslices);
    

        /*obs1 /= (float)(measures-therm);
        obs2 /= (float)(measures-therm);
        obs3 /= (float)(measures-therm);*/



        
	//print_lattice(N_t,N_s);
	//printf("Hello World\n");


	//fclose(fseed);
	fclose(finput);

	return 0;
}

void energy(double *O1, double *O2, double *O3){

	double phi, force_s, force_t;

	xene_mass = 0. ;
	xene_spat = 0. ;
	xene_temp = 0. ;
	for (int it = 0; it <N_t; it++){
		for(int ix = 0; ix < N_s; ix++){
			for(int iy = 0; iy < N_s; iy++){
				for(int iz = 0; iz<N_s;iz++){


					phi = field[it][ix][iy][iz];
					//printf("%f \n",phi);
					force_t = field[npp[it][0]][ix][iy][iz];
					force_s = field[it][npp[ix][1]][iy][iz]+
							  field[it][ix][npp[iy][1]][iz]+
							  field[it][ix][iy][npp[iz][1]];

					xene_mass += mass2*phi*phi;
					xene_temp += 2.0*phi*phi;
					xene_temp -= 2.0*phi*force_t;
					xene_spat += 6.0*phi*phi;
					xene_spat -= 2.0*phi*force_s;

				}
			}
		}
	}

	*O1 = 1.*xene_mass/(N_s*N_s*N_s*N_t);
	*O2 = 1.*xene_spat/(N_s*N_s*N_s*N_t);
	*O3 = 1.*xene_temp/(N_s*N_s*N_s*N_t);

}

void update_heatbath(){

	double force, phi, aver, x, y;

	for(int it = 0; it < N_t; it++){
		for (int ix=0; ix< N_s; ix++){
			for(int iy=0; iy < N_s; iy++){
				for (int iz=0;iz<N_s;iz++){

					force = 0.;
					phi = field[it][ix][iy][iz];
					force += field[npp[it][0]][ix][iy][iz];
					force += field[nmm[it][0]][ix][iy][iz];
					force += field[it][npp[ix][1]][iy][iz];
					force += field[it][nmm[ix][1]][iy][iz];
					force += field[it][ix][npp[iy][1]][iz];
					force += field[it][ix][nmm[iy][1]][iz];
					force += field[it][ix][iy][npp[iz][1]];
					force += field[it][ix][iy][nmm[iz][1]];
					
					aver = force*sigma2;

					//printf("force:%f\taver:%f\n",force,aver);

					x = sqrt(sigma2)*sqrt(-2.*log(ran2(&seed)));
					y = x*cos(PI2*ran2(&seed)) +aver;
					field[it][ix][iy][iz] = y;
				}
			}
		}
	}

}

void update_overrelax(){

	double phi, force, aver;

	for(int it = 0; it < N_t; it++){
		for (int ix=0; ix< N_s; ix++){
			for(int iy=0; iy < N_s; iy++){
				for (int iz=0;iz<N_s;iz++){

					force = 0.;
					phi = field[it][ix][iy][iz];
					force += field[npp[it][0]][ix][iy][iz];
					force += field[nmm[it][0]][ix][iy][iz];
					force += field[it][npp[ix][1]][iy][iz];
					force += field[it][nmm[ix][1]][iy][iz];
					force += field[it][ix][npp[iy][1]][iz];
					force += field[it][ix][nmm[iy][1]][iz];
					force += field[it][ix][iy][npp[iz][1]];
					force += field[it][ix][iy][nmm[iz][1]];

					aver = force / mass2p8;
					field[it][ix][iy][iz] = 2.*aver - phi;
				}
			}
		}
	}
}

void print_lattice(int nt, int ns){

	for (int it=0; it<nt;it++){
		for (int ix=0; ix<ns;ix++){
			for (int iy=0; iy<ns;iy++){
				for (int iz=0; iz<ns;iz++){
					printf("%f ",field[it][ix][iy][iz]);
				}
			}
		}
		printf("\n======================\n");
	}
}

void geometry(){

	for (int i = 0; i<N_t; i++){
                npp[i][0] = i+1;
                nmm[i][0] = i-1;
        }

        npp[N_t-1][0] = 0;
        nmm[0][0] = N_t - 1;

        for (int i = 0; i<N_s; i++){
                npp[i][1] = i+1;
                nmm[i][1] = i-1;
        }

        npp[N_s-1][1] = 0;
        nmm[0][1] = N_s - 1;

        
}

void init_lattice(int iflag){

        if (iflag ==0){
                return;
        }

        if (iflag == 1){

                for (int it=0; it<N_t;it++){
					for (int ix=0; ix<N_s;ix++){
						for (int iy=0; iy<N_s;iy++){
							for (int iz=0; iz<N_s;iz++){
								field[it][ix][iy][iz] = 1.0 - 2.0 * ran2(&seed);
								
				}
			}
		}
	}

 	}
 		if (iflag == 2){

                for (int it=0; it<N_t;it++){
					for (int ix=0; ix<N_s;ix++){
						for (int iy=0; iy<N_s;iy++){
							for (int iz=0; iz<N_s;iz++){
								field[it][ix][iy][iz] = 1.0 ;
								
				}
			}
		}
	}

 	}

}

 // Costanti indispensabili per ran2
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
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum){

        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;}

void MOM_q_timeslices(double complex *ris, int direction, int MC_time){

	double fase;
	double inv_N_s = 1./N_s;
	int n_i;

	for(int kappa = 0;kappa<(N_s/2+1);kappa++){
		for(int i =0; i<N_t;i++) {ris[i]=0.0+I*0.0;}
		
		for(int it = 0; it<N_t;it++){
			for (int ix=0;ix<N_s;ix++){
				for(int iy=0;iy<N_s;iy++){
					for (int iz=0;iz<N_s;iz++){


						//Modifica qui se vuoi altri vettori sul reticolo!
						if(direction==1)	n_i = ix;
						if(direction==2)	n_i = iy;
						if(direction==3)	n_i = iz;
						if((direction!=1)&&(direction!=2)&&(direction!=3)) {printf("Errore con dir!!\ndir=%d\n",direction); exit(EXIT_FAILURE);}

						fase = (double) PI2* n_i * inv_N_s *kappa;
						ris[it] += field[it][ix][iy][iz]*cexp(I*fase)/sqrt(N_s*N_s*N_s);
					}
				}
			}
		//Questa riga divide tutte le somme di un fattore sqrt(N_s**3) che serve per ovviare alla divergenza, su sphaleron questo è fatto dopo...
		//ris[it]/=sqrt(N_s*N_s*N_s);
		}
		printf("%d %d %d ", MC_time, dir, kappa);
      	for (int i=0; i<N_t; i++) printf(" %.12g %.12g", creal(ris[i]), cimag(ris[i]));
      	printf("\n");
	}
}