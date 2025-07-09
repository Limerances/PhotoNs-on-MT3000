#include "../inc/typedef.h"
#include "../inc/setparam.h"
#include <stdio.h>
#include <stdlib.h>

typedef char charArr[9];

#define PPN 9

//#define DBG

int main(void) {
	int PROC_SIZE = NSIDE0PROC * NSIDE1PROC;

	printf(" %d %d %d %d %d\n", NDOMinSUDOM, NSIDE0SUDOM, NSIDE1SUDOM, NSIDE0PROC, NSIDE1PROC);
	printf(" PROC_SIZE = %d generating rankfile...\n", PROC_SIZE);	
	int dside0 = NSIDE0PROC/NSIDE0SUDOM;
	int dside1 = NSIDE1PROC/NSIDE1SUDOM;
	int n;
	int new_rank;
	FILE* of = fopen("./rankfile", "w");
	for (n=0; n<PROC_SIZE; n++) {
		int idx = n;
		int is_sudom;

		int py = idx / NSIDE1PROC ; 
		int pz = idx % NSIDE1PROC ; 

		int sy = py / dside0;
		int sz = pz / dside1;

		int sudom_r  = sy * NSIDE1SUDOM + sz;

		int iy = py -  sy*dside0 ;
		int iz = pz -  sz*dside1 ;

		int sudom_size = dside0 * dside1;

		new_rank = sudom_r * sudom_size + iy * dside1 + iz;

		if ( (0 == py % dside0) && (0 == pz % dside1) ) {
			is_sudom = 1;
//			printf("[%d] su %d  size =%d\n", n, sudom_r, sudom_size);
		}
		else {
			is_sudom = 0;
//			printf("[%d] cd %d \n", n, sudom_r);

		}
#ifdef DBG	
		printf(" (%d)  %d >>> %d\n",sudom_r, n, new_rank);
#endif
		
		fprintf(of, "rank %d=+n%d slot=%d\n", n, new_rank/PPN, new_rank%PPN);	
	}

	fclose(of);

}


