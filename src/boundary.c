#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "photoNs.h"

void upload_bnd(cpuType width) ;
void download_bnd(cpuType width);
void global_comm(cpuType width) ;

void boundary() {

	if ( isSUDOM  ) {
		part_bnd = mymalloc(sizeof(BodyBnd)* MAXNPART_BND_SUD, 30);	
		part_com = mymalloc(sizeof(BodyBnd)* MAXNPART_BND_SUD * 3, 31);	
	}


	double t0, t1, t2, t3;


	cpuType width_bnd =  cutoffRadius;


	t0 = dtime();
	upload_bnd( width_bnd );

	t1 = dtime();

	global_comm( width_bnd );

	t2 = dtime();
	download_bnd( width_bnd  );

	t3 = dtime();

//	if ( 0 == PROC_RANK % 8 ) {
//		printf(" up = %lf com = %lf down = %lf, tot=%lf\n", t1-t0, t2-t1, t3-t2, t3- t0);
//	}

	if (isSUDOM ) {
		myfree(part_bnd, 30);
		myfree(part_com, 31);
	}

}

void upload_bnd(cpuType wid_bnd) {
	int n, c;
	int nsend ;
	int nrecv[NDOMinSUDOM];

	MPI_Request req;
	MPI_Status  sta;
	MPI_Status  status[NDOMinSUDOM];


	int box_dom_l[3], box_dom_h[3];

	box_dom_l[0] = (int) (( BOX_DOM_L[0] + wid_bnd ) * POS2INT  );
	box_dom_l[1] = (int) (( BOX_DOM_L[1] + wid_bnd ) * POS2INT  );
	box_dom_l[2] = (int) (( BOX_DOM_L[2] + wid_bnd ) * POS2INT  );

	box_dom_h[0] = (int) (( BOX_DOM_H[0] - wid_bnd ) * POS2INT  );
	box_dom_h[1] = (int) (( BOX_DOM_H[1] - wid_bnd ) * POS2INT  );
	box_dom_h[2] = (int) (( BOX_DOM_H[2] - wid_bnd ) * POS2INT  );

	if ( COM_DOM ) {	
		for (n=0, c=0; n< NPART; n++ ){


#ifndef INTXYZ
			if (
					part[n].pos[0] < BOX_DOM_L[0] + wid_bnd ||
					part[n].pos[1] < BOX_DOM_L[1] + wid_bnd ||
					part[n].pos[2] < BOX_DOM_L[2] + wid_bnd ||
					part[n].pos[0] >=BOX_DOM_H[0] - wid_bnd ||
					part[n].pos[1] >=BOX_DOM_H[1] - wid_bnd ||
					part[n].pos[2] >=BOX_DOM_H[2] - wid_bnd ) 
			{
				part_b[c].pos[0] = part[n].pos[0];
				part_b[c].pos[1] = part[n].pos[1];
				part_b[c].pos[2] = part[n].pos[2];
				c ++;
			}
#else
				if (
					part[n].posi[0] < box_dom_l[0] ||
					part[n].posi[1] < box_dom_l[1] ||
					part[n].posi[2] < box_dom_l[2] ||
					part[n].posi[0] >=box_dom_h[0] ||
					part[n].posi[1] >=box_dom_h[1] ||
					part[n].posi[2] >=box_dom_h[2] ) 
			{
				part_b[c].posi[0] = part[n].posi[0];
				part_b[c].posi[1] = part[n].posi[1];
				part_b[c].posi[2] = part[n].posi[2];
				c ++;
			}





#endif
		}
		nsend = c;

		MPI_Isend(&nsend, 1, MPI_INT, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(&(nrecv[n]), 1, MPI_INT, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
		}
	}

	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);		
	}
	int nrecv_total = 0;
	if ( isSUDOM ) {
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv_total += nrecv[n];
		}
		if (nrecv_total > MAXNPART_BND_SUD) {
			printf(" [%d] bd  error upload_bnd %d %d\n", PROC_RANK, nrecv_total, MAXNPART_BND_SUD);
			exit(0);
		}
		NPART_BND = nrecv_total;
	}

	//////////////////////////////////////
	if ( COM_DOM ) {	
		MPI_Isend(part_b, nsend*sizeof(BodyBnd), MPI_BYTE, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		int ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(part_bnd+ip, nrecv[n]*sizeof(BodyBnd), MPI_BYTE, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}

	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);		
	}

	/////////////////////////////


}




void show_boundary(){
	int n;
	if ( isSUDOM ) {
		char fname[80];
		sprintf(fname ,"bnd_com_%d",SUDOM_RANK); 
		FILE *fd = fopen(fname, "w");
		for (n=0; n<NPART_BND; n++) {
	//	fprintf(fd, "%lf %lf %lf\n", part_bnd[n].pos[0], part_bnd[n].pos[1], part_bnd[n].pos[2]);
		}
		for (n=0; n<NPART_BND_COM; n++) {
//			fprintf(fd, "%lf %lf %lf\n", part_com[n].pos[0], part_com[n].pos[1], part_com[n].pos[2]);
		}
	}
	if ( COM_DOM ) {
		char fname[80];
//		sprintf(fname ,"bnd_%d", PROC_RANK); 
//		FILE *fd = fopen(fname, "w");
//		for (n=0; n<NPART_BND; n++) {
//		if (part_b[n].pos[0] > 900000.0)
//			fprintf(fd, "%lf %lf %lf\n", part_b[n].pos[0], part_b[n].pos[1], part_b[n].pos[2]);
//		}
//		for (n=0; n<NPART; n++) {
//		if (part[n].pos[0] > 900000.0)
//			fprintf(fd, "%lf %lf %lf\n", part[n].pos[0], part[n].pos[1], part[n].pos[2]);
//		}
	}


}


#ifndef INTXYZ

void download_bnd(cpuType wid_bnd) {
	int n, c, ip,m;
	int nrecv;
	int nsend[NDOMinSUDOM];
	int nsend_total;
	cpuType width = BOXSIZE/NSIDEMESH;
	cpuType x0, x1;
	cpuType x2, x3;
	MPI_Status  status;	
	MPI_Status  status1;	

	MPI_Request req[NDOMinSUDOM];
	MPI_Status  sta[NDOMinSUDOM];

	MPI_Request req1[NDOMinSUDOM];
	MPI_Status  sta1[NDOMinSUDOM];
	BodyBnd* sendbuff;

	if ( isSUDOM ) {	
		sendbuff = (BodyBnd*)mymalloc(sizeof(BodyBnd) * MAXNPART_BND_SUD * 3, 32);
		ip = 0;	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			c =0;
			cpuType disp  = 0.0;

			x0 = dom_grp_mesh_start[m] * width - wid_bnd;
			x1 = dom_grp_mesh_start[m] * width;
			x2 = dom_grp_mesh_end[m] * width;
			x3 = dom_grp_mesh_end[m] * width + wid_bnd;
			for (n=0; n< NPART_BND_COM; n++ ){
				if( part_com[n].pos[0] > x0 && part_com[n].pos[0] < x3 ){
					sendbuff[ip].pos[0] = part_com[n].pos[0];
					sendbuff[ip].pos[1] = part_com[n].pos[1];
					sendbuff[ip].pos[2] = part_com[n].pos[2];
					c ++;
					ip ++;
				}
			}




			if (NDOM_HEAD == m) {
				x0 += BOXSIZE;
				x1 += BOXSIZE;
			}

			if ( m == NDOMinSUDOM - 1 ) {
				x2 -= BOXSIZE;
				x3 -= BOXSIZE;
			}

			for (n=0; n< NPART_BND; n++ ){
				if (( part_bnd[n].pos[0]>x0 && part_bnd[n].pos[0]<x1)
				||  ( part_bnd[n].pos[0]>x2 && part_bnd[n].pos[0]<x3) ) {
					sendbuff[ip].pos[0] = part_bnd[n].pos[0];
					sendbuff[ip].pos[1] = part_bnd[n].pos[1];
					sendbuff[ip].pos[2] = part_bnd[n].pos[2];
					c ++;
					ip ++;
				}
			}

			if (NDOM_HEAD==m)
				for (n=0; n< NPART_BND_COM; n++ ){
					if( part_com[n].pos[0] >  BOXSIZE - wid_bnd) {
						sendbuff[ip].pos[0] = part_com[n].pos[0]- BOXSIZE;
						sendbuff[ip].pos[1] = part_com[n].pos[1];
						sendbuff[ip].pos[2] = part_com[n].pos[2];
						c ++;
						ip ++;
					}
				}

			if (m== NDOMinSUDOM - 1)
				for (n=0; n< NPART_BND_COM; n++ ){
					if ( part_com[n].pos[0] < wid_bnd ){
						sendbuff[ip].pos[0] = part_com[n].pos[0] + BOXSIZE;
						sendbuff[ip].pos[1] = part_com[n].pos[1];
						sendbuff[ip].pos[2] = part_com[n].pos[2];
						c ++;
						ip ++;
					}
				}

			nsend[m] = c;

			MPI_Isend(&(nsend[m]), 1, MPI_INT, dom_grp[m], 1, MPI_COMM_WORLD, &(req[m]));
		}
		nsend_total = ip;
	}
	/////////////

	if ( COM_DOM ) {
		MPI_Recv(&nrecv, 1, MPI_INT, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &status);

		if (nrecv > MAXNPART_BND) {
			printf(" [%d] bd  error download_bnd %d %d\n", PROC_RANK, nrecv, MAXNPART_BND);
			exit(0);
		}
		NPART_BND = nrecv;
	}

	if ( isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req[m]), &(sta[m]));
		}
	}

	//////////////////////////////////////
	if (  isSUDOM ) {
		//		for (n=0; n< nsend_total; n++ ){
		//			if ( sendbuff[n].pos[0] < wid_bnd )
		//				sendbuff[n].pos[0] += BOXSIZE;
		//
		//			if ( sendbuff[n].pos[0] > BOXSIZE - wid_bnd )
		//				sendbuff[n].pos[0] -= BOXSIZE;
		//		}	

		ip = 0;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
				MPI_Isend(sendbuff+ip, nsend[m]*sizeof(BodyBnd), MPI_BYTE, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			ip += nsend[m];
		}

	}

	if ( COM_DOM ) {
		MPI_Recv(part_b, nrecv*sizeof(BodyBnd), MPI_BYTE, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);

		cpuType x0 = BOX_DOM_L[0] - wid_bnd;	
		cpuType y0 = BOX_DOM_H[0] + wid_bnd;	
		cpuType x1 = BOX_DOM_L[1] - wid_bnd;	
		cpuType y1 = BOX_DOM_H[1] + wid_bnd;	
		cpuType x2 = BOX_DOM_L[2] - wid_bnd;	
		cpuType y2 = BOX_DOM_H[2] + wid_bnd;	
	
		for (n=0; n<nrecv; n++) { 
			if (part_b[n].pos[0] < x0) part_b[n].pos[0] += BOXSIZE;
			if (part_b[n].pos[1] < x1) part_b[n].pos[1] += BOXSIZE;
			if (part_b[n].pos[2] < x2) part_b[n].pos[2] += BOXSIZE;

			if (part_b[n].pos[0] > y0) part_b[n].pos[0] -= BOXSIZE;
			if (part_b[n].pos[1] > y1) part_b[n].pos[1] -= BOXSIZE;
			if (part_b[n].pos[2] > y2) part_b[n].pos[2] -= BOXSIZE;

		}

	}
	if (  isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req1[m]), &(sta1[m]));
		}
	}




	if (  isSUDOM ) {	
		myfree(sendbuff, 32);
	}

}

void global_comm(cpuType wid_bnd) {
//	MPI_Barrier(MPI_COMM_WORLD);
	if ( isSUDOM ) {
		int thisidx = SUDOM_RANK;		
		int rank_send[8];
		int rank_recv[8];
		int tr, ts;	
		int x_this, y_this;
		x_this = thisidx / NSIDE1SUDOM;
		y_this = thisidx % NSIDE1SUDOM;
		//printf(" <%d> %d %d\n", thisidx, x_this, y_this);
		int n, m, c;
		int np, mp;
		c=0;
		for (n=-1; n<2; n++) {
			for (m=-1; m<2; m++) {
				if (n==0 && m==0)
					continue;

				np = ( x_this + n + NSIDE0SUDOM ) % NSIDE0SUDOM ;
				mp = ( y_this + m + NSIDE1SUDOM ) % NSIDE1SUDOM ;
				ts = np * NSIDE1SUDOM + mp;

				np = ( x_this - n + NSIDE0SUDOM ) % NSIDE0SUDOM ;
				mp = ( y_this - m + NSIDE1SUDOM ) % NSIDE1SUDOM ;
				tr = np * NSIDE1SUDOM + mp;

				rank_send[c] = proc_rank_sudom[ts];
				rank_recv[c] = proc_rank_sudom[tr];

				//printf(" <%d> %d : %d %d :: %d %d \n",  SUDOM_RANK,c, ts, tr, rank_send[c], rank_recv[c]);

				c++;
			}
		}

		for (n=0; n<8; n++) {
			PROC_RANK_SEND_NGB[n] = rank_send[n];
			PROC_RANK_RECV_NGB[n] = rank_recv[n];

		}

		int nrecv[8];
		int nsend[8];
		int i, j;
		m=0;
		cpuType y0, y1, z0, z1;

		y0 = BOX0L +  wid_bnd;
		y1 = BOX0H -  wid_bnd;
		z0 = BOX1L +  wid_bnd;
		z1 = BOX1H -  wid_bnd;

		BodyBnd* sendbuff = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND_SUD*3, 33);
		int ip=0;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[1] < y0 && part_bnd[n].pos[2] < z0 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];

				//				if (x_this == 0) 
				//					sendbuff[ip].pos[1] += BOXSIZE;
				//				if (x_this == NSIDE0SUDOM-1)
				//					sendbuff[ip].pos[1] -= BOXSIZE;
				//
				//				if (y_this == 0) 
				//					sendbuff[ip].pos[1] += BOXSIZE;
				//				if (y_this == NSIDE0SUDOM-1)
				//					sendbuff[ip].pos[1] -= BOXSIZE;

				c ++;
				ip ++;
			}

		}
		nsend[0] = c;
//		printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[1] < y0 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}
		}
		nsend[1] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[1] < y0 && part_bnd[n].pos[2] > z1 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}
		}
		nsend[2] = c;

		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[2] < z0 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}
		}
		nsend[3] = c;

		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[2] > z1 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}
		}
		nsend[4] = c;

		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[1] > y1 && part_bnd[n].pos[2] < z0 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}

		}
		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		nsend[5] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[1] > y1 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}
		}
		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		nsend[6] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].pos[1] > y1 && part_bnd[n].pos[2] > z1 ) {
				sendbuff[ip].pos[0] = part_bnd[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd[n].pos[2];
				c ++;
				ip ++;
			}
		}
		nsend[7] = c;
//		printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);

		int nsend_total = ip;
		//
		//		for (n=0; n< nsend_total; n++ ){
		//			if ( sendbuff[n].pos[1] < wid_bnd )
		//				sendbuff[n].pos[1] += BOXSIZE;
		//
		//			if ( sendbuff[n].pos[1] > BOXSIZE - wid_bnd )
		//				sendbuff[n].pos[1] -= BOXSIZE;
		//
		//			if ( sendbuff[n].pos[2] < wid_bnd )
		//				sendbuff[n].pos[2] += BOXSIZE;
		//
		//			if ( sendbuff[n].pos[2] > BOXSIZE - wid_bnd )
		//				sendbuff[n].pos[2] -= BOXSIZE;
		//		}	
		//




		MPI_Request req[8];
		MPI_Status  sta[8];
		MPI_Status  status[8];

		int nrecv_total = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(&(nsend[m]), 1, MPI_INT, rank_send[m], 2, MPI_COMM_WORLD, &(req[m]) );
		}

		for (m=0; m<8; m++) {
			MPI_Recv( &(nrecv[m]), 1, MPI_INT, rank_recv[m], 2, MPI_COMM_WORLD, &(status[m]));
//			printf("%d:  nrecv[%d] = %d\n", SUDOM_RANK, m, nrecv[m]);
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}

		for (m=0; m<8; m++) {
			nrecv_total += nrecv[m];
		}
		NPART_BND_COM = nrecv_total;

//		for (m=0; m<8; m++) {
//			printf(" [%d] %d  %d: %d \n", SUDOM_RANK, nsend[m], nrecv[m],  rank_send[m]);
//		}
		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(sendbuff+ip, nsend[m]*sizeof(BodyBnd), MPI_BYTE, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			ip += nsend[m];
		}
		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Recv( part_com+ip, nrecv[m]*sizeof(BodyBnd), MPI_BYTE, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			ip += nrecv[m];
//			printf("%d:  nrecv[%d] = %d\n", SUDOM_RANK, m, nrecv[m]);
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}





		myfree(sendbuff, 33);

	}	




}






#else


void download_bnd(cpuType wid_bnd) {
	int n, c, ip,m;
	int nrecv;
	int nsend[NDOMinSUDOM];
	int nsend_total;
	cpuType width = BOXSIZE/NSIDEMESH;
	cpuType x0, x1;
	cpuType x2, x3;

	int x0i, x1i, x2i, x3i;
	MPI_Status  status;	
	MPI_Status  status1;	

	MPI_Request req[NDOMinSUDOM];
	MPI_Status  sta[NDOMinSUDOM];

	MPI_Request req1[NDOMinSUDOM];
	MPI_Status  sta1[NDOMinSUDOM];
	BodyBnd* sendbuff;

			int box_wid_i = (int) ( (BOXSIZE - wid_bnd ) * POS2INT );
			int wid_i = (int) (wid_bnd * POS2INT);
			int box_i = (int) (BOXSIZE * POS2INT);

	if ( isSUDOM ) {	
		sendbuff = (BodyBnd*)mymalloc(sizeof(BodyBnd) * MAXNPART_BND_SUD * 3, 32);
		ip = 0;	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			c =0;
			cpuType disp  = 0.0;

			x0 = dom_grp_mesh_start[m] * width - wid_bnd;
			x1 = dom_grp_mesh_start[m] * width;
			x2 = dom_grp_mesh_end[m] * width;
			x3 = dom_grp_mesh_end[m] * width + wid_bnd;

			x0i = (int) (POS2INT * x0 );
			x1i = (int) (POS2INT * x1 );
			x2i = (int) (POS2INT * x2 );
			x3i = (int) (POS2INT * x3 );


			for (n=0; n< NPART_BND_COM; n++ ){
				if( part_com[n].posi[0] > x0i && part_com[n].posi[0] < x3i ){
					sendbuff[ip].posi[0] = part_com[n].posi[0];
					sendbuff[ip].posi[1] = part_com[n].posi[1];
					sendbuff[ip].posi[2] = part_com[n].posi[2];
					c ++;
					ip ++;
				}
			}




			if (NDOM_HEAD == m) {
				x0i += box_i;
				x1i += box_i;
			}

			if ( m == NDOMinSUDOM - 1 ) {
				x2i -= box_i;
				x3i -= box_i;
			}

			for (n=0; n< NPART_BND; n++ ){
				if (( part_bnd[n].posi[0]>x0i && part_bnd[n].posi[0]<x1i)
				||  ( part_bnd[n].posi[0]>x2i && part_bnd[n].posi[0]<x3i) ) {
					sendbuff[ip].posi[0] = part_bnd[n].posi[0];
					sendbuff[ip].posi[1] = part_bnd[n].posi[1];
					sendbuff[ip].posi[2] = part_bnd[n].posi[2];
					c ++;
					ip ++;
				}
			}
			if (NDOM_HEAD==m)
				for (n=0; n< NPART_BND_COM; n++ ){
					if( part_com[n].posi[0] >  box_wid_i) {
						sendbuff[ip].posi[0] = part_com[n].posi[0]- box_i;
						sendbuff[ip].posi[1] = part_com[n].posi[1];
						sendbuff[ip].posi[2] = part_com[n].posi[2];
						c ++;
						ip ++;
					}
				}

			if (m== NDOMinSUDOM - 1)
				for (n=0; n< NPART_BND_COM; n++ ){
					if ( part_com[n].posi[0] < wid_i ){
						sendbuff[ip].posi[0] = part_com[n].posi[0] + box_i;
						sendbuff[ip].posi[1] = part_com[n].posi[1];
						sendbuff[ip].posi[2] = part_com[n].posi[2];
						c ++;
						ip ++;
					}
				}

			nsend[m] = c;

			MPI_Isend(&(nsend[m]), 1, MPI_INT, dom_grp[m], 1, MPI_COMM_WORLD, &(req[m]));
		}
		nsend_total = ip;
	}
	/////////////

	if ( COM_DOM ) {
		MPI_Recv(&nrecv, 1, MPI_INT, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &status);

		if (nrecv > MAXNPART_BND) {
			printf(" [%d] bd  error download_bnd %d %d\n", PROC_RANK, nrecv, MAXNPART_BND);
			exit(0);
		}
		NPART_BND = nrecv;
	}

	if ( isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req[m]), &(sta[m]));
		}
	}

	//////////////////////////////////////
	if (  isSUDOM ) {
		ip = 0;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
				MPI_Isend(sendbuff+ip, nsend[m]*sizeof(BodyBnd), MPI_BYTE, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			ip += nsend[m];
		}

	}

	if ( COM_DOM ) {
		MPI_Recv(part_b, nrecv*sizeof(BodyBnd), MPI_BYTE, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);

		cpuType x0 = BOX_DOM_L[0] - wid_bnd;	
		cpuType y0 = BOX_DOM_H[0] + wid_bnd;	
		cpuType x1 = BOX_DOM_L[1] - wid_bnd;	
		cpuType y1 = BOX_DOM_H[1] + wid_bnd;	
		cpuType x2 = BOX_DOM_L[2] - wid_bnd;	
		cpuType y2 = BOX_DOM_H[2] + wid_bnd;	


		int x0i = (int) (x0 * POS2INT);
		int x1i = (int) (x1 * POS2INT);
		int x2i = (int) (x2 * POS2INT);
		int y0i = (int) (y0 * POS2INT);
		int y1i = (int) (y1 * POS2INT);
		int y2i = (int) (y2 * POS2INT);
	
		for (n=0; n<nrecv; n++) { 
			if (part_b[n].posi[0] < x0i) part_b[n].posi[0] += box_i;
			if (part_b[n].posi[1] < x1i) part_b[n].posi[1] += box_i;
			if (part_b[n].posi[2] < x2i) part_b[n].posi[2] += box_i;

			if (part_b[n].posi[0] > y0i) part_b[n].posi[0] -= box_i;
			if (part_b[n].posi[1] > y1i) part_b[n].posi[1] -= box_i;
			if (part_b[n].posi[2] > y2i) part_b[n].posi[2] -= box_i;

		}

	}
	if (  isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req1[m]), &(sta1[m]));
		}
	}




	if (  isSUDOM ) {	
		myfree(sendbuff, 32);
	}

}


void global_comm(cpuType wid_bnd) {
	if ( isSUDOM ) {
		int thisidx = SUDOM_RANK;		
		int rank_send[8];
		int rank_recv[8];
		int tr, ts;	
		int x_this, y_this;
		x_this = thisidx / NSIDE1SUDOM;
		y_this = thisidx % NSIDE1SUDOM;
		//printf(" <%d> %d %d\n", thisidx, x_this, y_this);
		int n, m, c;
		int np, mp;
		c=0;
		for (n=-1; n<2; n++) {
			for (m=-1; m<2; m++) {
				if (n==0 && m==0)
					continue;

				np = ( x_this + n + NSIDE0SUDOM ) % NSIDE0SUDOM ;
				mp = ( y_this + m + NSIDE1SUDOM ) % NSIDE1SUDOM ;
				ts = np * NSIDE1SUDOM + mp;

				np = ( x_this - n + NSIDE0SUDOM ) % NSIDE0SUDOM ;
				mp = ( y_this - m + NSIDE1SUDOM ) % NSIDE1SUDOM ;
				tr = np * NSIDE1SUDOM + mp;

				rank_send[c] = proc_rank_sudom[ts];
				rank_recv[c] = proc_rank_sudom[tr];

				//printf(" <%d> %d : %d %d :: %d %d \n",  SUDOM_RANK,c, ts, tr, rank_send[c], rank_recv[c]);

				c++;
			}
		}

		for (n=0; n<8; n++) {
			PROC_RANK_SEND_NGB[n] = rank_send[n];
			PROC_RANK_RECV_NGB[n] = rank_recv[n];

		}

		int nrecv[8];
		int nsend[8];
		int i, j;
		m=0;
		cpuType y0, y1, z0, z1;

		y0 = BOX0L +  wid_bnd;
		y1 = BOX0H -  wid_bnd;
		z0 = BOX1L +  wid_bnd;
		z1 = BOX1H -  wid_bnd;


		int y0i = (int) ( y0 *POS2INT );
		int y1i = (int) ( y1 *POS2INT );
		int z0i = (int) ( z0 *POS2INT );
		int z1i = (int) ( z1 *POS2INT );

		BodyBnd* sendbuff = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND_SUD*3, 33);
		int ip=0;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[1] < y0i && part_bnd[n].posi[2] < z0i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}

		}
		nsend[0] = c;
//		printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[1] < y0i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}
		}
		nsend[1] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[1] < y0i && part_bnd[n].posi[2] > z1i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}
		}
		nsend[2] = c;

		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[2] < z0i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}
		}
		nsend[3] = c;

		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[2] > z1i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}
		}
		nsend[4] = c;

		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[1] > y1i && part_bnd[n].posi[2] < z0i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}

		}
		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		nsend[5] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[1] > y1i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}
		}
		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);
		nsend[6] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd[n].posi[1] > y1i && part_bnd[n].posi[2] > z1i ) {
				sendbuff[ip].posi[0] = part_bnd[n].posi[0];
				sendbuff[ip].posi[1] = part_bnd[n].posi[1];
				sendbuff[ip].posi[2] = part_bnd[n].posi[2];
				c ++;
				ip ++;
			}
		}
		nsend[7] = c;
//		printf(" ip = %d, %d\n", ip, MAXNPART_BND_SUD);

		int nsend_total = ip;

		MPI_Request req[8];
		MPI_Status  sta[8];
		MPI_Status  status[8];

		int nrecv_total = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(&(nsend[m]), 1, MPI_INT, rank_send[m], 2, MPI_COMM_WORLD, &(req[m]) );
		}

		for (m=0; m<8; m++) {
			MPI_Recv( &(nrecv[m]), 1, MPI_INT, rank_recv[m], 2, MPI_COMM_WORLD, &(status[m]));
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}

		for (m=0; m<8; m++) {
			nrecv_total += nrecv[m];
		}
		NPART_BND_COM = nrecv_total;

		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(sendbuff+ip, nsend[m]*sizeof(BodyBnd), MPI_BYTE, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			ip += nsend[m];
		}
		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Recv( part_com+ip, nrecv[m]*sizeof(BodyBnd), MPI_BYTE, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			ip += nrecv[m];
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}





		myfree(sendbuff, 33);

	}	




}


#endif

