#include <stdio.h>
#include <math.h>
#include "schedule.h"
#include "photoNs.h"

static int NUMBER_OUTPUT_ANALYSIS;
static int *index_output_analysis;


void create_output_list(int step_number, cpuType ai, cpuType af, cpuType aa_i, cpuType aa_f, int analysis_number){
	if (fabs(ai - af) < 1.0e-15)
		return;
	int i,r, num_base;
	cpuType slop = 3.3;
	NUMBER_OUTPUT_ANALYSIS = 7;
	index_output_analysis = (int*)malloc(sizeof(int)*NUMBER_OUTPUT_ANALYSIS);	

	index_output_analysis[0] = 380;
	index_output_analysis[1] = 858;
	index_output_analysis[2] = 1546;
	index_output_analysis[3] = 1622;
	index_output_analysis[4] = 1714;
	index_output_analysis[5] = 1832;
	index_output_analysis[6] = 2000;

	
	
	index_output_analysis[0] = 190;
	index_output_analysis[1] = 429;
	index_output_analysis[2] = 773;
	index_output_analysis[3] = 811;
	index_output_analysis[4] = 857;
	index_output_analysis[5] = 916;
	index_output_analysis[6] = 1000;

	/*
	cpuType lna_tot = log(af) - log(ai);
	cpuType lnaa = log(aa_f) - log(aa_i);
	cpuType lnaaf = log(af) - log(aa_f);
	int num_step_aa = (int) ( lnaa * step_number / lna_tot  ); 
//	printf(" num_step_aa = %d %f %f\n", num_step_aa, lnaa, lna_tot);
	num_base = (int)( num_step_aa / ( slop-1.0) / ( (cpuType) analysis_number ) ); 	
	if (num_base < 1) num_base = 1;

	//	r = analysis_number -1;
	//	index_output_analysis[r] = num_base ;
	index_output_analysis[analysis_number-1] = 0;
	for (i=1; i<analysis_number; i++) {
		r = analysis_number - i -1;
		index_output_analysis[r] = (int) ( num_base * (slop * i / analysis_number + 1.0)) ;
		if (index_output_analysis[r]<1) index_output_analysis[r] = 1;
		//		printf(" output [ %d ] %d \n", r,  index_output_analysis[r]);	
	}

	for (i=1; i<analysis_number; i++) {
		r = analysis_number - i -1;
		index_output_analysis[r] += index_output_analysis[r+1];
	}
	num_step_aa = index_output_analysis[0]; 
	if (num_step_aa > step_number)  { 
		printf("ERROR in schedule \n");
		exit(0);
	}
	//	printf(" tot analysis number = %d\n", index_output_analysis[0]);	
	int idx_aa_f = (int) lnaaf; 
	for (i=0; i<analysis_number; i++) {
		index_output_analysis[i] = step_number  - idx_aa_f - index_output_analysis[i];
		//		printf(" 2 output [ %d ] %d \n", i,  index_output_analysis[i]);	
	}
*/
	//////////////////////

	if ( 0 == PROC_RANK ) {
		int loop_start = 0;
		cpuType dloga =  (log(af) - log(ai) ) / (step_number);
		cpuType loga_i = log( ai);
		cpuType loga_f = dloga + log( ai);	
		int loop;

		char fname[128];
		sprintf(fname, "%s/output_schedue.txt", PathSnapshot);
		FILE *fsch = fopen(fname, "w");
		if (fsch == NULL){ printf(" canot open %s\n", fname); exit(0);}

		fprintf(fsch, "#idx #loop   red-shift scale-factor    time(Gyr)\n");
		fprintf(fsch, "   0     0  %10lf %12lf   %10lf\n", 1.0/ai - 1.0, ai , t_flat_lcdm_a(ai) );
	//	fprintf(fsch, "#snap  Redshift   Time(Gyr)\n");
		for (loop = 0; loop <=step_number; loop++) {
			loga_i = loop * dloga + log( ai);

			cpuType a_i = exp(loga_i);	

			if (analysis_required(loop)) {
				int nn =   analysis_output_index(loop) ;

			//	fprintf(fsch, "%5d  %3.5lf  %2.5lf\n", nn, 1.0/a_f - 1.0, t_flat_lcdm_a(a_f) );
			//	fprintf(fsch, " %03d %5d %e %e %e\n", nn, loop, 1.0/a_f - 1.0, a_f , t_flat_lcdm_a(a_f) );
				fprintf(fsch, " %3d %5d  %10lf %12lf   %10lf\n", nn, loop, 1.0/a_i - 1.0, a_i , t_flat_lcdm_a(a_i) );
			}
		}	
		fclose(fsch);
	}

}





void create_schedule(int step_number, cpuType ai, cpuType af, cpuType aa_i, cpuType aa_f, int analysis_number){
	if (fabs(ai - af) < 1.0e-15)
		return;
	int i,r, num_base;
	cpuType slop = 3.5;
//	cpuType slop = 3.3;
	NUMBER_OUTPUT_ANALYSIS = analysis_number;
	index_output_analysis = (int*)malloc(sizeof(int)*NUMBER_OUTPUT_ANALYSIS);	

	cpuType lna_tot = log(af) - log(ai);
	cpuType lnaa = log(aa_f) - log(aa_i);
	cpuType lnaaf = log(af) - log(aa_f);
	int num_step_aa = (int) ( lnaa * step_number / lna_tot  ); 
//	printf(" num_step_aa = %d %f %f\n", num_step_aa, lnaa, lna_tot);
	num_base = (int)( num_step_aa / ( slop-1.0) / ( (cpuType) analysis_number ) ); 	
	if (num_base < 1) num_base = 1;

	//	r = analysis_number -1;
	//	index_output_analysis[r] = num_base ;
	index_output_analysis[analysis_number-1] = 0;
	for (i=1; i<analysis_number; i++) {
		r = analysis_number - i -1;
		index_output_analysis[r] = (int) ( num_base * (slop * i / analysis_number + 1.0)) ;
		if (index_output_analysis[r]<1) index_output_analysis[r] = 1;
		//		printf(" output [ %d ] %d \n", r,  index_output_analysis[r]);	
	}

	for (i=1; i<analysis_number; i++) {
		r = analysis_number - i -1;
		index_output_analysis[r] += index_output_analysis[r+1];
	}
	num_step_aa = index_output_analysis[0]; 
	if (num_step_aa > step_number)  { 
		printf("ERROR in schedule \n");
		exit(0);
	}
	//	printf(" tot analysis number = %d\n", index_output_analysis[0]);	
	int idx_aa_f = (int) lnaaf; 
	for (i=0; i<analysis_number; i++) {
		index_output_analysis[i] = step_number  - idx_aa_f - index_output_analysis[i];
		//		printf(" 2 output [ %d ] %d \n", i,  index_output_analysis[i]);	
	}

	//////////////////////

	if ( 0 == PROC_RANK ) {
		int loop_start = 0;
		cpuType dloga =  (log(af) - log(ai) ) / (step_number);
		cpuType loga_i = log( ai);
		cpuType loga_f = dloga + log( ai);	
		int loop;

		char fname[128];
		sprintf(fname, "%s/output_schedue.txt", PathSnapshot);
		FILE *fsch = fopen(fname, "w");
		if (fsch == NULL){ printf(" canot open %s\n", fname); exit(0);}

		fprintf(fsch, "#idx #loop   red-shift scale-factor    time(Gyr)\n");
		fprintf(fsch, "   0     0  %10lf %12lf   %10lf\n", 1.0/ai - 1.0, ai , t_flat_lcdm_a(ai) );
	//	fprintf(fsch, "#snap  Redshift   Time(Gyr)\n");
		for (loop = 0; loop <=step_number; loop++) {
			loga_i = loop * dloga + log( ai);

			cpuType a_i = exp(loga_i);	

			if (analysis_required(loop)) {
				int nn =   analysis_output_index(loop) ;

			//	fprintf(fsch, "%5d  %3.5lf  %2.5lf\n", nn, 1.0/a_f - 1.0, t_flat_lcdm_a(a_f) );
			//	fprintf(fsch, " %03d %5d %e %e %e\n", nn, loop, 1.0/a_f - 1.0, a_f , t_flat_lcdm_a(a_f) );
				fprintf(fsch, " %3d %5d  %10lf %12lf   %10lf\n", nn, loop, 1.0/a_i - 1.0, a_i , t_flat_lcdm_a(a_i) );
			}
		}	
		fclose(fsch);
	}

}


int analysis_output_index(int step){
	int i, idx = 0;
	for (i=0; i<NUMBER_OUTPUT_ANALYSIS; i++){
		if ( step == index_output_analysis[i]) {
			idx = NUMBER_OUTPUT_ANALYSIS -  i -1;
			idx = i + 1;
		}
	}
	return idx;
}

int snapshot_required(int step){
	int i, yes = 0;

	for (i=0; i<NUMBER_OUTPUT_ANALYSIS; i++){
		if (step == index_output_analysis[i]) {
			int idx = analysis_output_index(step);
//			if  (  0 == idx ) yes = 1;
//			if  (  4 == idx ) yes = 1;
//			if  ( 12 == idx ) yes = 1;
		
		}
	}
	return yes;
}




int meshout_required(int step){
	int i, yes = 0;
	for (i=0; i<NUMBER_OUTPUT_ANALYSIS; i++){
		if ( i< NUMBER_OUTPUT_ANALYSIS ){
			if (step == index_output_analysis[i]) {
				yes = 1;
			}
		}
	}
	return yes;
}



int analysis_required(int step){
	int i, yes = 0;
	for (i=0; i<NUMBER_OUTPUT_ANALYSIS; i++){
		if ( i>= 0 ){
			if (step == index_output_analysis[i]) {
				yes = 1;
			}
		}
	}
	return yes;
}



void free_schedule(){
	if (index_output_analysis != NULL)
	free(index_output_analysis);
}


