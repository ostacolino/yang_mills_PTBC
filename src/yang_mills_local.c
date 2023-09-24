#ifndef YM_LOCAL_C
#define YM_LOCAL_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
{
	Gauge_Conf GC;
	Geometry geo;
	GParam param;

	char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
	int count;
	FILE *datafilep, *chiprimefilep, *topchar_tprof_filep, *topchar_MOM_tprof_filep;
	time_t time1, time2;

	// to disable nested parallelism
	#ifdef OPENMP_MODE
		// omp_set_nested(0); // deprecated
		omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
	#endif

	// read input file
	readinput(in_file, &param);

	// initialize random generator
	initrand(param.d_randseed);

	// open data_file
	init_data_file(&datafilep, &chiprimefilep, &topchar_tprof_filep, &topchar_MOM_tprof_filep, &param);

	// initialize geometry
	init_indexing_lexeo();
	init_geometry(&geo, &param);

	// if measure-only mode is active conf must be read from file => start=2 ignoring the value found in input file
	if (param.d_sample == 0)
	{
		fprintf(stdout, "MEASURE-ONLY MODE: performing measures on configuration read from file %s, no update will be performed\n", param.d_conf_file);
		param.d_start=2;
	}

	// initialize gauge configuration
	init_gauge_conf(&GC, &param);

	// --- TO BE REMOVED ---
/*	double* ratio;
	int errs=posix_memalign((void**)&ratio, (size_t)DOUBLE_ALIGN, (size_t) param.d_coolrepeat * sizeof(double));
 	if(errs!=0)
	{
		fprintf(stderr, "Problems in allocating a vector (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	check_correlation_decay_cooling(&GC, &geo, &param, ratio);
	FILE* ratio_fp;
	ratio_fp=fopen("CHECK_charge_cool.dat", "w");
	if (ratio_fp == NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", "CHECK_charge_cool.dat", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for (int iter=0; iter<param.d_coolrepeat; iter++) fprintf(ratio_fp, "%d  %.12lg\n", (iter+1)*param.d_coolsteps, ratio[iter]);
	free(ratio);
	fclose(ratio_fp);
*/
	//---------------------------

	// Monte Carlo begins
	time(&time1);
	if (param.d_sample == 0) // no update is done, only measures are performed on read configuration
	{
		perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep, topchar_MOM_tprof_filep);
	}
	else
	{
		for(count=0; count < param.d_sample; count++)
		{
			// update conf
			update(&GC, &geo, &param);

			// measure local observables
			if(GC.update_index % param.d_measevery == 0 && GC.update_index >= param.d_thermal)
			{
				perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep, topchar_MOM_tprof_filep);
			}

			// save configuration for backup
			if(param.d_saveconf_back_every!=0)
			{
				if(GC.update_index % param.d_saveconf_back_every == 0 )
				{
					// simple
					write_conf_on_file(&GC, &param);

					// backup copy
					write_conf_on_file_back(&GC, &param);
				}
			}

			// save configuration for offline analysis
			if(param.d_saveconf_analysis_every!=0)
			{
				if(GC.update_index % param.d_saveconf_analysis_every == 0 )
				{
					strcpy(name, param.d_conf_file);
					strcat(name, "_step_");
					sprintf(aux, "%ld", GC.update_index);
					strcat(name, aux);
					write_conf_on_file_with_name(&GC, &param, name);
				}
			}
		}
	}
	time(&time2);
	// Monte Carlo ends

	// close data file
	fclose(datafilep);
	if (param.d_chi_prime_meas==1) fclose(chiprimefilep);
	if (param.d_topcharge_tprof_meas==1) fclose(topchar_tprof_filep);
	if (param.d_topcharge_MOM_tprof_meas==1) fclose(topchar_MOM_tprof_filep);

	// save configuration
	if(param.d_saveconf_back_every!=0)
	{
		write_conf_on_file(&GC, &param);
	}

	// print simulation details
	print_parameters_local(&param, time1, time2);

	// free gauge configuration
	free_gauge_conf(&GC, &param);

	// free geometry
	free_geometry(&geo, &param);
}

void print_template_input(void)
{
	FILE *fp;
	fp=fopen("template_input.in", "w");

	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
  else
	{
		fprintf(fp, "size 4 4 4 4\n");
		fprintf(fp,"\n");
		fprintf(fp, "beta 5.705\n");
		fprintf(fp, "theta 1.5\n");
		fprintf(fp,"\n");
		fprintf(fp, "sample    10\n");
		fprintf(fp, "thermal   0\n");
		fprintf(fp, "overrelax 5\n");
		fprintf(fp, "measevery 1\n");
		fprintf(fp,"\n");
		fprintf(fp, "start                    0  # 0=ordered  1=random  2=from saved configuration\n");
		fprintf(fp, "saveconf_back_every      5  # if 0 does not save, else save backup configurations every ... updates\n");
		fprintf(fp, "saveconf_analysis_every  5  # if 0 does not save, else save configurations for analysis every ... updates\n");
		fprintf(fp, "\n");
		fprintf(fp, "coolsteps                  3  # number of cooling steps to be used\n");
		fprintf(fp, "coolrepeat                 5  # number of times 'coolsteps' are repeated\n");
		fprintf(fp, "chi_prime_meas             0  # 1=YES, 0=NO\n");
		fprintf(fp, "topcharge_tprof_meas       0  # 1=YES, 0=NO\n");
		fprintf(fp, "topcharge_MOM_tprof_meas   0  # 1=YES, 0=NO\n");
		fprintf(fp, "topcharge_MOM_tprof_dir    1  # direction of k vector: 1=x, 2=y, 3=z\n");
		fprintf(fp,"\n");
		fprintf(fp, "#output files\n");
		fprintf(fp, "conf_file                 conf.dat\n");
		fprintf(fp, "data_file                 dati.dat\n");
		fprintf(fp, "chiprime_data_file        chiprime_cool.dat\n");
		fprintf(fp, "topcharge_tprof_file      topo_tcorr_cool.dat\n");
		fprintf(fp, "topcharge_MOM_tprof_file  topo_MOM_tcorr_cool.dat\n");
		fprintf(fp, "log_file                  log.dat\n");
		fprintf(fp, "\n");
		fprintf(fp, "randseed 0    #(0=time)\n");
		fclose(fp);
	}
}

int main (int argc, char **argv)
{
	char in_file[STD_STRING_LENGTH];
	if(argc != 2)
	{
		printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
		printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
		printf("Usage: %s input_file\n\n", argv[0]);

		printf("Compilation details:\n");
		printf("\tN_c (number of colors): %d\n", NCOLOR);
		printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
		printf("\tNum_levels (number of levels): %d\n", NLEVELS);
		printf("\n");
		printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
		printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

		#ifdef DEBUG
			printf("\n\tDEBUG mode\n");
		#endif

		#ifdef OPENMP_MODE
			printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
		#endif

		#ifdef THETA_MODE
			printf("\n\tusing imaginary theta\n");
		#endif

		printf("\n");

		#ifdef __INTEL_COMPILER
			printf("\tcompiled with icc\n");
		#elif defined(__clang__)
			printf("\tcompiled with clang\n");
		#elif defined( __GNUC__ )
			printf("\tcompiled with gcc version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		#endif

		print_template_input();

		return EXIT_SUCCESS;
	}
	else
	{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
		{
			fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
		}
		else
		{
			strcpy(in_file, argv[1]);
		}
	}
	real_main(in_file);
	return EXIT_SUCCESS;
}

#endif
