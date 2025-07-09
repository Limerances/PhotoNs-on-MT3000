#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "photoNs.h"

typedef struct {
	unsigned long id;
	cpuType pos[3];
	cpuType vel[3];
} BodyBndAnalysis;

MPI_Datatype BodyBndAnlstype;

int MAXNPART_BND_ANLS = 0;
int MAXNPART_BND_ANLS_SUD = 0;
int MAXNPART_COM_ANLS_SUD = 0;
BodyBndAnalysis *part_bnd_anls, *part_com_anls, *part_b_anls;
float ratio = 1.1;

void upload_bnd_anls(cpuType width) ;
void download_bnd_anls(cpuType width);
void global_comm_anls(cpuType width) ;

//void check_and_resize() {
//}

#define MIN_BUFF_LENGTH 5000000

void boundary_analysis(cpuType width_bnd) {

	if (MAXNPART_BND_ANLS == 0) {
		// follow the sim part
		MAXNPART_BND_ANLS = MAXNPART * max_bnd_ratio;


		// qwang 
		if (MAXNPART_BND_ANLS < MIN_BUFF_LENGTH) MAXNPART_BND_ANLS = MIN_BUFF_LENGTH;	
		// qwang


		MAXNPART_BND_ANLS_SUD = MAXNPART_BND_ANLS * NDOMinSUDOM;
		MAXNPART_COM_ANLS_SUD = MAXNPART_BND_ANLS_SUD * 3;

		MPI_Type_contiguous(sizeof(BodyBndAnalysis), MPI_BYTE, &BodyBndAnlstype);
		MPI_Type_commit(&BodyBndAnlstype);

	}

	if ( isSUDOM  ) {
		part_bnd_anls = mymalloc(sizeof(BodyBndAnalysis)* MAXNPART_BND_ANLS_SUD, 210);	
		part_com_anls = mymalloc(sizeof(BodyBndAnalysis)* MAXNPART_COM_ANLS_SUD, 211);	
	}

	if (COM_DOM) {
		part_b_anls = mymalloc(sizeof(BodyBndAnalysis) * MAXNPART_BND_ANLS, 212);
	}

	double t0, t1, t2, t3;

	t0 = dtime();
	upload_bnd_anls( width_bnd );

	t1 = dtime();

	global_comm_anls( width_bnd );

	t2 = dtime();
	download_bnd_anls( width_bnd  );

	t3 = dtime();

	if (isSUDOM ) {
		myfree(part_bnd_anls, 210);
		myfree(part_com_anls, 211);
	}

}

void upload_bnd_anls(cpuType wid_bnd) {
	int n, c;
	int nsend ;
	int nrecv[NDOMinSUDOM];

	MPI_Request req;
	MPI_Status  sta;
	MPI_Status  status[NDOMinSUDOM];

	if ( COM_DOM ) {

#ifndef INTXYZ	
		for (n=0, c=0; n< NPART; n++ ){
			if (
					part[n].pos[0] < BOX_DOM_L[0] + wid_bnd ||
					part[n].pos[1] < BOX_DOM_L[1] + wid_bnd ||
					part[n].pos[2] < BOX_DOM_L[2] + wid_bnd ||
					part[n].pos[0] >=BOX_DOM_H[0] - wid_bnd ||
					part[n].pos[1] >=BOX_DOM_H[1] - wid_bnd ||
					part[n].pos[2] >=BOX_DOM_H[2] - wid_bnd ) 
			{
				part_b_anls[c].pos[0] = part[n].pos[0];
				part_b_anls[c].pos[1] = part[n].pos[1];
				part_b_anls[c].pos[2] = part[n].pos[2];
				part_b_anls[c].vel[0] = part[n].vel[0];
				part_b_anls[c].vel[1] = part[n].vel[1];
				part_b_anls[c].vel[2] = part[n].vel[2];
				part_b_anls[c].id = part[n].id;
				c ++;
				if (c == MAXNPART_BND_ANLS && n < NPART-1) {
					printf(" [%d] bd upload_bnd_anls NPART %d MAXNPART %d send %d MAXNPART_BND_ANLS %d -> %d\n", PROC_RANK, NPART, MAXNPART, c, MAXNPART_BND_ANLS, (int)((float)c * ratio));
					MAXNPART_BND_ANLS = (int)((float)c * ratio);
					part_b_anls = myrealloc(part_b_anls, sizeof(BodyBndAnalysis) * MAXNPART_BND_ANLS, 212);
				}
			}

		}

#else
		for (n=0, c=0; n< NPART; n++ ){
			cpuType ppos[3];
			ppos[0] = part[n].posi[0] * INT2POS;
			ppos[1] = part[n].posi[1] * INT2POS;
			ppos[2] = part[n].posi[2] * INT2POS;

			if (
				ppos[0] < BOX_DOM_L[0] + wid_bnd ||
				ppos[1] < BOX_DOM_L[1] + wid_bnd ||
				ppos[2] < BOX_DOM_L[2] + wid_bnd ||
				ppos[0] >=BOX_DOM_H[0] - wid_bnd ||
				ppos[1] >=BOX_DOM_H[1] - wid_bnd ||
				ppos[2] >=BOX_DOM_H[2] - wid_bnd ) 
			{
				part_b_anls[c].pos[0] =ppos[0];
				part_b_anls[c].pos[1] =ppos[1];
				part_b_anls[c].pos[2] =ppos[2];
				part_b_anls[c].vel[0] = part[n].vel[0];
				part_b_anls[c].vel[1] = part[n].vel[1];
				part_b_anls[c].vel[2] = part[n].vel[2];
				part_b_anls[c].id = part[n].id;
				c ++;
				if (c == MAXNPART_BND_ANLS && n < NPART-1) {
					printf(" [%d] bd upload_bnd_anls NPART %d MAXNPART %d send %d MAXNPART_BND_ANLS %d -> %d\n", PROC_RANK, NPART, MAXNPART, c, MAXNPART_BND_ANLS, (int)((float)c * ratio));
					MAXNPART_BND_ANLS = (int)((float)c * ratio);
					part_b_anls = myrealloc(part_b_anls, sizeof(BodyBndAnalysis) * MAXNPART_BND_ANLS, 212);
				}
			}

		}



#endif

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
		if (nrecv_total > MAXNPART_BND_ANLS_SUD) {
			printf(" [%d] bd  error upload_bnd_anls NPART %d MAXNPART %d recv %d MAXNPART_BND_ANLS_SUD %d\n", PROC_RANK, NPART, MAXNPART, nrecv_total, MAXNPART_BND_ANLS_SUD);
			MAXNPART_BND_ANLS_SUD = (int)((float)nrecv_total * ratio);
			part_bnd_anls = myrealloc(part_bnd_anls, sizeof(BodyBndAnalysis) * MAXNPART_BND_ANLS_SUD, 210);
		}
		NPART_BND = nrecv_total;
	}

	//////////////////////////////////////
	if ( COM_DOM ) {	
		MPI_Isend(part_b_anls, nsend, BodyBndAnlstype, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		int ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(part_bnd_anls+ip, nrecv[n], BodyBndAnlstype, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
			if (ip > MAXNPART_BND_ANLS_SUD) {
				printf(" [%d] bd  error upload_bnd_anls NPART %d MAXNPART %d recv %d MAXNPART_BND_ANLS_SUD %d\n", PROC_RANK, NPART, MAXNPART, ip, MAXNPART_BND_ANLS_SUD);
				exit(0);
			}
		}
	}

	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);		
	}

	/////////////////////////////


}

void download_bnd_anls(cpuType wid_bnd) {
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
	BodyBndAnalysis* sendbuff;

	if ( isSUDOM ) {	
		sendbuff = (BodyBndAnalysis*)mymalloc(sizeof(BodyBndAnalysis) * MAXNPART_COM_ANLS_SUD, 213);
		ip = 0;	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			c =0;
			x0 = dom_grp_mesh_start[m] * width - wid_bnd;
			x1 = dom_grp_mesh_start[m] * width;
			x2 = dom_grp_mesh_end[m] * width;
			x3 = dom_grp_mesh_end[m] * width + wid_bnd;
			cpuType disp  = 0.0;
			for (n=0; n< NPART_BND_COM; n++ ){
				if( part_com_anls[n].pos[0] > x0 && part_com_anls[n].pos[0] < x3 ){
					sendbuff[ip].pos[0] = part_com_anls[n].pos[0];
					sendbuff[ip].pos[1] = part_com_anls[n].pos[1];
					sendbuff[ip].pos[2] = part_com_anls[n].pos[2];
					sendbuff[ip].vel[0] = part_com_anls[n].vel[0];
					sendbuff[ip].vel[1] = part_com_anls[n].vel[1];
					sendbuff[ip].vel[2] = part_com_anls[n].vel[2];
					sendbuff[ip].id = part_com_anls[n].id;
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
				if (( part_bnd_anls[n].pos[0]>x0 && part_bnd_anls[n].pos[0]<x1)
				||  ( part_bnd_anls[n].pos[0]>x2 && part_bnd_anls[n].pos[0]<x3) ) {
					sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
					sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
					sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
					sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
					sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
					sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
					sendbuff[ip].id = part_bnd_anls[n].id;
					c ++;
					ip ++;
				}
			}

			if (NDOM_HEAD==m)
				for (n=0; n< NPART_BND_COM; n++ ){
					if( part_com_anls[n].pos[0] >  BOXSIZE - wid_bnd) {
						sendbuff[ip].pos[0] = part_com_anls[n].pos[0]- BOXSIZE;
						sendbuff[ip].pos[1] = part_com_anls[n].pos[1];
						sendbuff[ip].pos[2] = part_com_anls[n].pos[2];
						sendbuff[ip].vel[0] = part_com_anls[n].vel[0];
						sendbuff[ip].vel[1] = part_com_anls[n].vel[1];
						sendbuff[ip].vel[2] = part_com_anls[n].vel[2];
						sendbuff[ip].id = part_com_anls[n].id;
						c ++;
						ip ++;
					}
				}

			if (m== NDOMinSUDOM - 1)
				for (n=0; n< NPART_BND_COM; n++ ){
					if ( part_com_anls[n].pos[0] < wid_bnd ){
						sendbuff[ip].pos[0] = part_com_anls[n].pos[0] + BOXSIZE;
						sendbuff[ip].pos[1] = part_com_anls[n].pos[1];
						sendbuff[ip].pos[2] = part_com_anls[n].pos[2];
						sendbuff[ip].vel[0] = part_com_anls[n].vel[0];
						sendbuff[ip].vel[1] = part_com_anls[n].vel[1];
						sendbuff[ip].vel[2] = part_com_anls[n].vel[2];
						sendbuff[ip].id = part_com_anls[n].id;
						c ++;
						ip ++;
					}
				}

			nsend[m] = c;

			MPI_Isend(&(nsend[m]), 1, MPI_INT, dom_grp[m], 1, MPI_COMM_WORLD, &(req[m]));
		}
		nsend_total = ip;
		if (nsend_total > MAXNPART_COM_ANLS_SUD) {
			printf(" [%d] bd  error download_bnd_anls NPART %d MAXNPART %d send %d MAXNPART_COM_ANLS_SUD %d\n", PROC_RANK, NPART, MAXNPART, nsend_total, MAXNPART_COM_ANLS_SUD);
			exit(0);
		}
	}
	/////////////

	if ( COM_DOM ) {
		MPI_Recv(&nrecv, 1, MPI_INT, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &status);

		if (nrecv > MAXNPART_BND_ANLS) {
			printf(" [%d] bd download_bnd_anls realloc part_b_anls recv %d MAXNPART_BND_ANLS %d\n", PROC_RANK, nrecv, MAXNPART_BND_ANLS);
			MAXNPART_BND_ANLS = (int)((float)nrecv * ratio);
			part_b_anls = myrealloc(part_b_anls, sizeof(BodyBndAnalysis) * MAXNPART_BND_ANLS, 212);
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
				MPI_Isend(sendbuff+ip, nsend[m], BodyBndAnlstype, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			ip += nsend[m];
		}

	}

	if ( COM_DOM ) {
		MPI_Recv(part_b_anls, nrecv, BodyBndAnlstype, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);

		cpuType x0 = BOX_DOM_L[0] - wid_bnd;	
		cpuType y0 = BOX_DOM_H[0] + wid_bnd;	
		cpuType x1 = BOX_DOM_L[1] - wid_bnd;	
		cpuType y1 = BOX_DOM_H[1] + wid_bnd;	
		cpuType x2 = BOX_DOM_L[2] - wid_bnd;	
		cpuType y2 = BOX_DOM_H[2] + wid_bnd;	
	
		for (n=0; n<nrecv; n++) { 
			if (part_b_anls[n].pos[0] < x0) part_b_anls[n].pos[0] += BOXSIZE;
			if (part_b_anls[n].pos[1] < x1) part_b_anls[n].pos[1] += BOXSIZE;
			if (part_b_anls[n].pos[2] < x2) part_b_anls[n].pos[2] += BOXSIZE;

			if (part_b_anls[n].pos[0] > y0) part_b_anls[n].pos[0] -= BOXSIZE;
			if (part_b_anls[n].pos[1] > y1) part_b_anls[n].pos[1] -= BOXSIZE;
			if (part_b_anls[n].pos[2] > y2) part_b_anls[n].pos[2] -= BOXSIZE;

		}

	}
	if (  isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req1[m]), &(sta1[m]));
		}
	}




	if (  isSUDOM ) {	
		myfree(sendbuff, 213);
	}

}

void global_comm_anls(cpuType wid_bnd) {
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

		BodyBndAnalysis* sendbuff = (BodyBndAnalysis*)mymalloc(sizeof(BodyBndAnalysis)*MAXNPART_COM_ANLS_SUD, 214);
		int ip=0;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[1] < y0 && part_bnd_anls[n].pos[2] < z0 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;

				c ++;
				ip ++;
			}

		}
		nsend[0] = c;
//		printf(" ip = %d, %d\n", ip, MAXNPART_BND_ANLS_SUD);
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[1] < y0 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}
		}
		nsend[1] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[1] < y0 && part_bnd_anls[n].pos[2] > z1 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}
		}
		nsend[2] = c;

		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_ANLS_SUD);
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[2] < z0 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}
		}
		nsend[3] = c;

		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[2] > z1 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}
		}
		nsend[4] = c;

		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[1] > y1 && part_bnd_anls[n].pos[2] < z0 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}

		}
		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_ANLS_SUD);
		nsend[5] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[1] > y1 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}
		}
		//	printf(" ip = %d, %d\n", ip, MAXNPART_BND_ANLS_SUD);
		nsend[6] = c;
		for (n=0, c=0; n< NPART_BND; n++ ){
			if ( part_bnd_anls[n].pos[1] > y1 && part_bnd_anls[n].pos[2] > z1 ) {
				sendbuff[ip].pos[0] = part_bnd_anls[n].pos[0];
				sendbuff[ip].pos[1] = part_bnd_anls[n].pos[1];
				sendbuff[ip].pos[2] = part_bnd_anls[n].pos[2];
				sendbuff[ip].vel[0] = part_bnd_anls[n].vel[0];
				sendbuff[ip].vel[1] = part_bnd_anls[n].vel[1];
				sendbuff[ip].vel[2] = part_bnd_anls[n].vel[2];
				sendbuff[ip].id = part_bnd_anls[n].id;
				c ++;
				ip ++;
			}
		}
		nsend[7] = c;
		printf(" ip = %d, %d\n", ip, MAXNPART_BND_ANLS_SUD);

		int nsend_total = ip;
		if (nsend_total > MAXNPART_COM_ANLS_SUD) {
			printf(" [%d] bd  error global_comm_anls NPART %d MAXNPART %d send %d MAXNPART_COM_ANLS_SUD %d\n", PROC_RANK, NPART, MAXNPART, nsend_total, MAXNPART_COM_ANLS_SUD);
			exit(0);
		}

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
			MPI_Isend(sendbuff+ip, nsend[m], BodyBndAnlstype, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			ip += nsend[m];
		}
		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Recv( part_com_anls+ip, nrecv[m], BodyBndAnlstype, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			ip += nrecv[m];
			if (ip > MAXNPART_COM_ANLS_SUD) {
				printf(" [%d] bd  error global_comm_anls NPART %d MAXNPART %d recv %d MAXNPART_COM_ANLS_SUD %d\n", PROC_RANK, NPART, MAXNPART, ip, MAXNPART_COM_ANLS_SUD);
				exit(0);
			}
//			printf("%d:  nrecv[%d] = %d\n", SUDOM_RANK, m, nrecv[m]);
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}





		myfree(sendbuff, 214);

	}	




}


void particle_combine_analysis() {
	int npart_tot;
	npart_tot = NPART + NPART_BND;


	if (npart_tot > MAXNPART) {
		printf(" [%d] in analysis: need more buff MAXNPART %d: NPART %d + NPART_BND %d\n", PROC_RANK, MAXNPART, NPART, NPART_BND);
		MAXNPART = npart_tot;
		part = (Body*) myrealloc(part, sizeof(Body)* MAXNPART, 9);
	}


	int n, i ;
	for ( n=0; n<NPART; n++ ) {
		part[n].tag=1;
	}
	
	for (i=0, n=NPART; i<NPART_BND; i++, n++) {
		part[n].tag = 0;
#ifndef INTXYZ		
		part[n].pos[0] = part_b_anls[i].pos[0];
		part[n].pos[1] = part_b_anls[i].pos[1];
		part[n].pos[2] = part_b_anls[i].pos[2];
#else

		part[n].posi[0] = (int) POS2INT *  part_b_anls[i].pos[0];
		part[n].posi[1] = (int) POS2INT * part_b_anls[i].pos[1];
		part[n].posi[2] = (int) POS2INT * part_b_anls[i].pos[2];

#endif
		part[n].vel[0] = part_b_anls[i].vel[0];
		part[n].vel[1] = part_b_anls[i].vel[1];
		part[n].vel[2] = part_b_anls[i].vel[2];
		part[n].id = part_b_anls[i].id;
	}

	NPART = npart_tot;
	myfree(part_b_anls, 212);
} 
