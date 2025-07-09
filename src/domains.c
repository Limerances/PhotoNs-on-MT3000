#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "photoNs.h"
#include "measure.h"

#ifndef INTXYZ

void checkpart(int label) {
	return;
	int i, n;
	if (!COM_DOM) return;
	printf("[%d] checkpart %d start.\n", PROC_RANK, label);
	for (i=0; i<NPART; i++) {
	//	if (isnan(acc[0][i])) {printf("[%d] (%d) %d acc0 is nan.\n", PROC_RANK, label, i);exit(0);}
	//	if (isnan(acc[1][i])) {printf("[%d] (%d) %d acc1 is nan.\n", PROC_RANK, label, i);exit(0);}
	//	if (isnan(acc[2][i])) {printf("[%d] (%d) %d acc2 is nan.\n", PROC_RANK, label, i);exit(0);}
	//	if (abs(acc[0][i]) > 1.0e2) {printf("[%d] (%d) %d acc0 is huge.\n", PROC_RANK, label, i);exit(0);}
	//	if (abs(acc[1][i]) > 1.0e2) {printf("[%d] (%d) %d acc1 is huge.\n", PROC_RANK, label, i);exit(0);}
	//	if (abs(acc[2][i]) > 1.0e2) {printf("[%d] (%d) %d acc2 is huge.\n", PROC_RANK, label, i);exit(0);}
	//	if (isnan(acc_pm[0][i])) {printf("[%d] (%d) %d acc_pm0 is nan.\n", PROC_RANK, label, i);exit(0);}
	//	if (isnan(acc_pm[1][i])) {printf("[%d] (%d) %d acc_pm1 is nan.\n", PROC_RANK, label, i);exit(0);}
	//	if (isnan(acc_pm[2][i])) {printf("[%d] (%d) %d acc_pm2 is nan.\n", PROC_RANK, label, i);exit(0);}

		if (isnan(part[i].vel[0])) {printf("[%d] (%d) %d vel0 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].vel[1])) {printf("[%d] (%d) %d vel1 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].vel[2])) {printf("[%d] (%d) %d vel2 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].pos[0])) {printf("[%d] (%d) %d pos0 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].pos[1])) {printf("[%d] (%d) %d pos1 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].pos[2])) {printf("[%d] (%d) %d pos2 is nan.\n", PROC_RANK, label, i);exit(0);}
	}
	printf("[%d] checkpart %d end.\n", PROC_RANK, label);
	return;
}

void checkbound(Body* part, int NPART, int label) {
	int i;
	int flag;
	for (i=0; i<NPART; i++) {
		flag = 0;
		if (part[i].pos[0] < 0.0 || part[i].pos[0] >= BOXSIZE) {printf("[%d] (%d) %d pos0 is out of range vel (%e ).\n", PROC_RANK, i, label, part[i].vel[0]);flag=1;}
		if (part[i].pos[1] < 0.0 || part[i].pos[1] >= BOXSIZE) {printf("[%d] (%d) %d pos1 is out of range vel (%e ).\n", PROC_RANK, i, label, part[i].vel[1]);flag=1;}
		if (part[i].pos[2] < 0.0 || part[i].pos[2] >= BOXSIZE) {printf("[%d] (%d) %d pos2 is out of range vel (%e ).\n", PROC_RANK, i, label, part[i].vel[2]);flag=1;}
		if (flag == 1) {printf("[%d] pos (%d) %d (%e %e %e)\n", PROC_RANK, i, label, part[i].pos[0], part[i].pos[1], part[i].pos[2]);exit(0);}
	}
}

// qwang 2024-10-11 for UMS2-high resolution
void warp(Body* pt, int np) {
	int n;
	for (n=0; n< np; n++ ){
#ifdef CPU_DP
		cpuType eps = 1.0e-15;
#else
		cpuType eps = 1.2e-7;
#endif
		cpuType boxsize = BOXSIZE * (1.0 - eps);

            if (pt[n].pos[0] < 0.0) pt[n].pos[0] += boxsize;
            if (pt[n].pos[1] < 0.0) pt[n].pos[1] += boxsize;
            if (pt[n].pos[2] < 0.0) pt[n].pos[2] += boxsize;
            if (pt[n].pos[0] >=BOXSIZE) pt[n].pos[0] -= boxsize;
            if (pt[n].pos[1] >=BOXSIZE) pt[n].pos[1] -= boxsize;
            if (pt[n].pos[2] >=BOXSIZE) pt[n].pos[2] -= boxsize;
	}
	checkbound(pt, np, 1);
}
/*
void warp(Body* pt, int np) {
	int n;
	for (n=0; n< np; n++ ){
#ifdef CPU_DP
		cpuType eps = 1.0e-15;
#else
		cpuType eps = 1.0e-6;
#endif
		cpuType boxsize = BOXSIZE * (1.0 - eps);
            if (pt[n].pos[0] < 0.0) pt[n].pos[0] += boxsize;
            if (pt[n].pos[1] < 0.0) pt[n].pos[1] += boxsize;
            if (pt[n].pos[2] < 0.0) pt[n].pos[2] += boxsize;
            if (pt[n].pos[0] >=BOXSIZE) pt[n].pos[0] -= boxsize;
            if (pt[n].pos[1] >=BOXSIZE) pt[n].pos[1] -= boxsize;
            if (pt[n].pos[2] >=BOXSIZE) pt[n].pos[2] -= boxsize;
	}
	checkbound(pt, np, 1);
}
*/


void partup() {
	int npart_in;
	int npart_ex;
	int n, c;

	int nrecv[NDOMinSUDOM];
	MPI_Request req;
	MPI_Status  sta;
	MPI_Status  status[NDOMinSUDOM];

	if (COM_DOM) {
		for (n=0, c=0; n<NPART; n++) {
			if (	part[n].pos[0] <  BOX_DOM_L[0] ||
					part[n].pos[0] >= BOX_DOM_H[0] ||
					part[n].pos[1] <  BOX_DOM_L[1] ||
					part[n].pos[1] >= BOX_DOM_H[1] ||
					part[n].pos[2] <  BOX_DOM_L[2] ||
					part[n].pos[2] >= BOX_DOM_H[2] )

			{
				part[n].tag = 0;
				c++;
			}
			else {
				part[n].tag = 1;
			}
		}
		npart_ex = c;
		NPART_IN = NPART - npart_ex;

		// inplace replacement

		Body tpart;
		int p0, p1;

		p0 = 0;
		p1 = NPART-1; 
		while (p0 < p1 ) {
			if ( 0==part[p0].tag ) {
				tpart = part[p0];
				while ( 0==part[p1].tag ) {
					p1--; 
				}
				if (p1 <= p0)
					break;

				part[p0] = part[p1];
				part[p1] = tpart;
			}

			p0++;
		}
		//printf("[%d] p0= %d, p1 = %d , NPART = %d = %d + %d - (%e %e) (%e %e) (%e %e)\n", PROC_RANK, p0, p1, NPART, NPART_IN, npart_ex, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);

		MPI_Isend(&npart_ex, 1, MPI_INT, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);


		//		for (n=0; n<NPART; n++)
		//			part[n].tag = 0;
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
		if (nrecv_total > MAXNPART_DOM_SUD) {
			printf(" error upload_bnd %d %d\n", nrecv_total, MAXNPART_DOM_SUD);
			exit(0);
		}
		NPART_DOM_SUD = nrecv_total;
	}

	if ( COM_DOM ) {	
		MPI_Isend(part+NPART_IN, npart_ex, Bodytype, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
		//	MPI_Isend(part+NPART_IN, npart_ex*sizeof(Body), MPI_BYTE, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		int ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(part_dom_sud+ip, nrecv[n], Bodytype, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			//	MPI_Recv(part_dom_sud+ip, nrecv[n]*sizeof(Body), MPI_BYTE, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}


	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);		
	}


}

void partexch(){
	if ( isSUDOM ) {
		int thisidx = SUDOM_RANK;		
		int rank_send[8];
		int rank_recv[8];
		int tr, ts;	
		int x_this, y_this;
		x_this = thisidx / NSIDE1SUDOM;
		y_this = thisidx % NSIDE1SUDOM;
		//	printf(" <%d> [%d]  %d %d BOX %lf %lf : %lf %lf : %lf %lf\n", thisidx, PROC_RANK, x_this, y_this, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);
		//	printf(" <%d> %d %d\n", thisidx, x_this, y_this);
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

				c++;
			}
		}

		for (n=0; n<8; n++) {
			PROC_RANK_SEND_NGB[n] = rank_send[n];
			PROC_RANK_RECV_NGB[n] = rank_recv[n];

		}
		//	printf(" <%d> send %d %d %d %d %d %d %d %d\n", thisidx, rank_send[0], rank_send[1], rank_send[2], rank_send[3], rank_send[4], rank_send[5], rank_send[6], rank_send[7]);
		//	printf(" <%d> recv %d %d %d %d %d %d %d %d\n", thisidx, rank_recv[0], rank_recv[1], rank_recv[2], rank_recv[3], rank_recv[4], rank_recv[5], rank_recv[6], rank_recv[7]);
		int nrecv[8];
		int nsend[8];
		int i, j;
		m=0;
		cpuType y0, y1, z0, z1;
//		double y0, y1, z0, z1;

		y0 = BOX0L;
		y1 = BOX0H;
		z0 = BOX1L;
		z1 = BOX1H;
		Body* pd = part_dom_sud;

		Body* sendbuff = (Body*)mymalloc(sizeof(Body)*MAXNPART_DOM_SUD*2, 20);
		if (sendbuff == NULL) {
			printf(" ERROR partexch: cannot open sendbuff  \n");
			fflush(stdout);
			exit(0);
		}
		int ip=0;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] < y0 && pd[n].pos[2] < z0 ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[0] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] < y0 && pd[n].pos[2] >= z0 && pd[n].pos[2] < z1 ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[1] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] < y0 && pd[n].pos[2] >= z1  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[2] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] >= y0 &&  pd[n].pos[1] < y1 && pd[n].pos[2] < z0  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[3] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] >= y0 &&  pd[n].pos[1] < y1 && pd[n].pos[2] >= z1  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[4] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] >= y1 && pd[n].pos[2] < z0 ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[5] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] >= y1 && pd[n].pos[2] >= z0 && pd[n].pos[2] < z1 ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[6] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].pos[1] >= y1 && pd[n].pos[2] >= z1  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[7] = c;

		int nsend_total = ip;
		// mengchen
		if (nsend_total > MAXNPART_DOM_SUD * 2) {
			printf(" enlarge MAXNPART_DOM_SUD * 2 ! %d %d\n", nsend_total, MAXNPART_DOM_SUD * 2);
			exit(0);
		}


		MPI_Request req[8];
		MPI_Status  sta[8];
		MPI_Status  status[8];

		int nrecv_total = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(&(nsend[m]), 1, MPI_INT, rank_send[m], 2, MPI_COMM_WORLD, &(req[m]) );
		}
		//////////////////////////////
		for (m=0; m<8; m++) {
			MPI_Recv( &(nrecv[m]), 1, MPI_INT, rank_recv[m], 2, MPI_COMM_WORLD, &(status[m]));
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}

		for (m=0; m<8; m++) {
			nrecv_total += nrecv[m];
		}
		NPART_DOM_COM = nrecv_total;
		//qwang
		if (nrecv_total>MAXNPART_DOM) {
			printf(" enlarge MAXNPART_DOM ! %d %d\n", nrecv_total, MAXNPART_DOM);
			exit(0);
		}

		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(sendbuff+ip, nsend[m], Bodytype, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			//	MPI_Isend(sendbuff+ip, nsend[m]*sizeof(Body), MPI_BYTE, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			ip += nsend[m];
		}
		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Recv( part_dom+ip, nrecv[m], Bodytype, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			//MPI_Recv( part_dom+ip, nrecv[m]*sizeof(Body), MPI_BYTE, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			ip += nrecv[m];
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}

		myfree(sendbuff, 20);


		warp(part_dom, NPART_DOM_COM);
	}
}

void partdown() {
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
	Body* sendbuff;

	if ( isSUDOM ) {	
		sendbuff = (Body*)mymalloc(sizeof(Body) * MAXNPART_DOM_SUD * 2, 21);
		if (sendbuff == NULL) {
			printf(" ERROR partdown : cannot open sendbuff  \n");
			fflush(stdout);
			exit(0);
		}
		Body *pd = part_dom_sud;	

		for (n=0; n< NPART_DOM_SUD; n++ ){
			if (pd[n].pos[0] < 0.0)
				pd[n].pos[0] += BOXSIZE;
			if (pd[n].pos[0] >=BOXSIZE)
				pd[n].pos[0] -= BOXSIZE;
		}

		ip = 0;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			c =0;
			x1 = dom_grp_mesh_start[m] * width;
			x2 = dom_grp_mesh_end[m] * width;
			cpuType disp  = 0.0;
			for (n=0; n< NPART_DOM_COM; n++ ){
				if( part_dom[n].pos[0] >= x1 && part_dom[n].pos[0] < x2 ){
					sendbuff[ip] = part_dom[n];
					c ++;
					ip ++;
				}
			}
			for (n=0; n< NPART_DOM_SUD; n++ ){
				if ( pd[n].pos[0] >= x1 && pd[n].pos[0] < x2
						&& pd[n].pos[1] >= BOX0L && pd[n].pos[1] < BOX0H
						&& pd[n].pos[2] >= BOX1L && pd[n].pos[2] < BOX1H ) {
					sendbuff[ip] = pd[n];
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

		NPART_EXP = nrecv;
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
			MPI_Isend(sendbuff+ip, nsend[m], Bodytype, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			//	MPI_Isend(sendbuff+ip, nsend[m]*sizeof(Body), MPI_BYTE, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			ip += nsend[m];
		}

	}

	if ( COM_DOM ) {
		MPI_Recv(part+NPART_IN, NPART_EXP, Bodytype, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);
		//	MPI_Recv(part+NPART_IN, NPART_EXP*sizeof(Body), MPI_BYTE, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);

	}
	if (  isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req1[m]), &(sta1[m]));
		}
	}

	if (  isSUDOM ) {	
		myfree(sendbuff, 21);
	}


}

void determine_dom_start() {
	int n;
	int r;
	int idx;
	int dom_start[NDOMinSUDOM];
	for (n=0; n<NDOMinSUDOM; n++) {
		dom_start[n] = dom_grp_mesh_start[n];
	}
	cpuType *npart_mesh;
	cpuType *npart_mesh_tot;

	MPI_Request req;
	MPI_Status  sta;
	MPI_Status  status[NDOMinSUDOM];
	MPI_Request  request[NDOMinSUDOM];

	if ( COM_DOM ) {
		int nproc_sud = NDOMinSUDOM - NDOM_HEAD;
		r = dom_grp_rank;
		int len = dom_grp_mesh_size[r];
		npart_mesh = (cpuType*)mymalloc(sizeof(cpuType)*len, 22);
		if (npart_mesh == NULL) {
			printf(" ERROR determine_dom_start: cannot open npart_mesh  \n");
			fflush(stdout);
			exit(0);
		}

		for ( n =0; n<len; n++  ) {
			npart_mesh[n] = 0.0;
		}
		cpuType dx = BOXSIZE/NSIDEMESH;
		int intmp = 1;
		for ( n=0; n<NPART; n++) {			
			idx  = (int)( part[n].pos[0]/dx ) - dom_grp_mesh_start[r] ;
			if (idx >=0 && idx < len )
				npart_mesh[idx] += intmp;
		}
		MPI_Isend(npart_mesh, len, MPI_CPUTYPE, PROC_RANK_SUDOM, 2, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		npart_mesh_tot = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 23);
		if (npart_mesh_tot == NULL) {
			printf(" ERROR determine_dom_start: cannot open npart_mesh_tot  \n");
			fflush(stdout);
			exit(0);
		}
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(npart_mesh_tot+dom_grp_mesh_start[n], dom_grp_mesh_size[n], MPI_CPUTYPE, dom_grp[n], 2, MPI_COMM_WORLD, &(status[n]));
		}
	}

	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);
	}

	if ( isSUDOM ) {
		int n;
		npart_mesh_tot[0] = 0.0;
		for (n=1; n<NSIDEMESH; n++) {
			npart_mesh_tot[n] += npart_mesh_tot[n-1];
		}
		cpuType np_sud = npart_mesh_tot[NSIDEMESH-1];	
		cpuType np_sep = np_sud/ (NDOMinSUDOM - NDOM_HEAD);

		int c=1;
		dom_start[0] = 0;
		for (n=0; n<NSIDEMESH; n++) {
			if (c > NDOMinSUDOM - NDOM_HEAD )
				break;

			if (  npart_mesh_tot[n] >  c*np_sep  ) {
				dom_start[c] = n;
				c++;
			}
		}

		for (n=0; n<NDOM_HEAD; n++) {
			dom_grp_mesh_start[n] = 0;
			dom_grp_mesh_end[n] = 0;
		}

		c =0 ;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			dom_grp_mesh_start[n] = dom_start[c]; 
			c++;
		}

		for (n=NDOM_HEAD; n<NDOMinSUDOM-1; n++) {
			dom_grp_mesh_end[n] = dom_grp_mesh_start[n+1]; 
		}
		dom_grp_mesh_end[NDOMinSUDOM-1] = NSIDEMESH;

		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			dom_grp_mesh_size[n] = dom_grp_mesh_end[n] -  dom_grp_mesh_start[n]; 
		}

		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Isend(dom_grp_mesh_start, NDOMinSUDOM, MPI_INT, dom_grp[n], 3, MPI_COMM_WORLD, &request[n]);
		}

	}

	if (COM_DOM) {
		MPI_Recv(dom_grp_mesh_start, NDOMinSUDOM, MPI_INT, PROC_RANK_SUDOM, 3, MPI_COMM_WORLD, &sta);
		myfree(npart_mesh, 22);
	}


	if (isSUDOM) {


		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Wait(&request[n], &status[n]);
		}

		myfree(npart_mesh_tot, 23);
	}

	for (n=0; n<NDOMinSUDOM-1; n++) {
		dom_grp_mesh_end[n] = dom_grp_mesh_start[n+1];
	}
	dom_grp_mesh_end[NDOMinSUDOM-1] = NSIDEMESH;

	for (n=0; n<NDOMinSUDOM; n++) {
		dom_grp_mesh_size[n] = dom_grp_mesh_end[n] - dom_grp_mesh_start[n];
	}

	cpuType dx = BOXSIZE / NSIDEMESH;

	BOX_DOM_L[0] = dx * dom_grp_mesh_start[dom_grp_rank];
	BOX_DOM_H[0] = dx * dom_grp_mesh_end[dom_grp_rank];

}

void domain_decomp() {

	determine_dom_start();

	if (isSUDOM ) {
		part_dom_sud = (Body*)mymalloc(sizeof(Body)*MAXNPART_DOM_SUD, 24);
		if (part_dom_sud == NULL) {
			printf(" ERROR domain_decomp: cannot open part_dom_sud  \n");
			fflush(stdout);
			exit(0);
		}
		part_dom = mymalloc(sizeof(Body)* MAXNPART_DOM, 25);	
		if (part_dom == NULL) {
			printf(" ERROR domain_decomp: cannot open part_dom  \n");
			fflush(stdout);
			exit(0);
		}
	}



	partup();

	NPART = NPART_IN;

	partexch();
	partdown();

	NPART = NPART_IN + NPART_EXP;

	int n;
	int c = 0;
	for (n=0; n<NPART; n++) {
		part[n].tag = 1;
		if (	part[n].pos[0] <  BOX_DOM_L[0] ||
				part[n].pos[0] >= BOX_DOM_H[0] ||
				part[n].pos[1] <  BOX_DOM_L[1] ||
				part[n].pos[1] >= BOX_DOM_H[1] ||
				part[n].pos[2] <  BOX_DOM_L[2] ||
				part[n].pos[2] >= BOX_DOM_H[2] )

		{

			printf(" [%d] ERROR DOMAIN :%d tag = %d %e %e %e\n", PROC_RANK, n, part[n].tag, part[n].pos[0], part[n].pos[1], part[n].pos[2]);
			printf("      (%lf,%lf) (%lf,%lf) (%lf,%lf)\n", BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);
			exit(0);

		}
		else {
			c++;
		}
	}

//	printf( " dom comp : %d %d\n", c, NPART);
	//realloc

	if ( isSUDOM ) {
		myfree(part_dom_sud, 24);
		myfree(part_dom, 25);
	}
}



#else

void checkpart(int label) {
	return;
	int i, n;
	if (!COM_DOM) return;
	printf("[%d] checkpart %d start.\n", PROC_RANK, label);
	for (i=0; i<NPART; i++) {
		if (isnan(part[i].vel[0])) {printf("[%d] (%d) %d vel0 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].vel[1])) {printf("[%d] (%d) %d vel1 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(part[i].vel[2])) {printf("[%d] (%d) %d vel2 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(INT2POS * part[i].posi[0])) {printf("[%d] (%d) %d pos0 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(INT2POS * part[i].posi[1])) {printf("[%d] (%d) %d pos1 is nan.\n", PROC_RANK, label, i);exit(0);}
		if (isnan(INT2POS * part[i].posi[2])) {printf("[%d] (%d) %d pos2 is nan.\n", PROC_RANK, label, i);exit(0);}
	}
	printf("[%d] checkpart %d end.\n", PROC_RANK, label);
	return;
}

void checkbound(Body* part, int NPART, int label) {
	int i;
	int flag;
	int box_i = (int) (POS2INT * BOXSIZE );
	for (i=0; i<NPART; i++) {
		flag = 0;
		if (part[i].posi[0] < 0 || part[i].posi[0] >= box_i) {printf("[%d] (%d) %d pos0 is out of range vel (%d ).\n", PROC_RANK, i, label, part[i].posi[0]);flag=1;}
		if (part[i].posi[1] < 0 || part[i].posi[1] >= box_i) {printf("[%d] (%d) %d pos1 is out of range vel (%d ).  BOXSIZE = %d\n", PROC_RANK, i, label, part[i].posi[1], (int)(BOXSIZE * POS2INT)  );flag=1;}
		if (part[i].posi[2] < 0 || part[i].posi[2] >= box_i) {printf("[%d] (%d) %d pos2 is out of range vel (%d ).\n", PROC_RANK, i, label, part[i].posi[2]);flag=1;}
		if (flag == 1) {printf("[%d] pos (%d) %d (%d %d %d)\n", PROC_RANK, i, label, part[i].posi[0], part[i].posi[1], part[i].posi[2]);exit(0);}
	}
}

// qwang 2024-10-11 for UMS2-high resolution

void warp(Body* pt, int np) {
	int n;

////printf(" box_i = %d %d %f\n", box_i, (int)(POS2INT * 0.9 *BOXSIZE), BOXSIZE  );
#ifdef CPU_DP
		cpuType eps = 1.0e-15;
#else
		cpuType eps = 1.5e-7;
#endif
	cpuType boxsize = BOXSIZE * (1.0 - eps);
//	boxsize = BOXSIZE  ;

	int box_c = (int) (POS2INT * BOXSIZE );
//	int box_i = (int) (POS2INT * boxsize );
	int box_i = box_c - 3;
	
	for (n=0; n< np; n++ ){
		if ( n == 1092 && PROC_RANK == 2022 ) printf(" posi = %d %d %d(box= %d %d)\n", pt[n].posi[0], pt[n].posi[1], pt[n]. posi[2], box_c, box_i );
            if (pt[n].posi[0] < 0    ) pt[n].posi[0] += box_i;
            if (pt[n].posi[1] < 0    ) pt[n].posi[1] += box_i;
            if (pt[n].posi[2] < 0    ) pt[n].posi[2] += box_i;
            if (pt[n].posi[0] >=box_c) pt[n].posi[0] -= box_i;
            if (pt[n].posi[1] >=box_c) pt[n].posi[1] -= box_i;
            if (pt[n].posi[2] >=box_c) pt[n].posi[2] -= box_i;
	}
	checkbound(pt, np, 1);
}


/*
void warp(Body* pt, int np) {
	int n;
	int box_i = (int) (POS2INT * BOXSIZE );

////printf(" box_i = %d %d %f\n", box_i, (int)(POS2INT * 0.9 *BOXSIZE), BOXSIZE  );
	for (n=0; n< np; n++ ){
#ifdef CPU_DP
		cpuType eps = 1.0e-15;
#else
		cpuType eps = 1.0e-6;
#endif
//	cpuType boxsize = BOXSIZE * (1.0 - eps);
            if (pt[n].posi[0] < 0) pt[n].posi[0] += box_i;
            if (pt[n].posi[1] < 0) pt[n].posi[1] += box_i;
            if (pt[n].posi[2] < 0) pt[n].posi[2] += box_i;
            if (pt[n].posi[0] >=box_i) pt[n].posi[0] -= box_i;
            if (pt[n].posi[1] >=box_i) pt[n].posi[1] -= box_i;
            if (pt[n].posi[2] >=box_i) pt[n].posi[2] -= box_i;
	}
	checkbound(pt, np, 1);
}
*/

void partup() {
	int npart_in;
	int npart_ex;
	int n, c;

	int nrecv[NDOMinSUDOM];
	MPI_Request req;
	MPI_Status  sta;
	MPI_Status  status[NDOMinSUDOM];


	int box_dom_l[3], box_dom_h[3];

	box_dom_l[0] = (int) ( BOX_DOM_L[0]  * POS2INT  );
	box_dom_l[1] = (int) ( BOX_DOM_L[1]  * POS2INT  );
	box_dom_l[2] = (int) ( BOX_DOM_L[2]  * POS2INT  );

	box_dom_h[0] = (int) ( BOX_DOM_H[0]  * POS2INT  );
	box_dom_h[1] = (int) ( BOX_DOM_H[1]  * POS2INT  );
	box_dom_h[2] = (int) ( BOX_DOM_H[2]  * POS2INT  );




	if (COM_DOM) {
		for (n=0, c=0; n<NPART; n++) {
			if (	part[n].posi[0] <  box_dom_l[0] ||
					part[n].posi[0] >=box_dom_h[0] ||
					part[n].posi[1] < box_dom_l[1] ||
					part[n].posi[1] >=box_dom_h[1] ||
					part[n].posi[2] < box_dom_l[2] ||
					part[n].posi[2] >=box_dom_h[2] )

			{
				part[n].tag = 0;
				c++;
			}
			else {
				part[n].tag = 1;
			}
		}
		npart_ex = c;
		NPART_IN = NPART - npart_ex;

		// inplace replacement

		Body tpart;
		int p0, p1;

		p0 = 0;
		p1 = NPART-1; 
		while (p0 < p1 ) {
			if ( 0==part[p0].tag ) {
				tpart = part[p0];
				while ( 0==part[p1].tag ) {
					p1--; 
				}
				if (p1 <= p0)
					break;

				part[p0] = part[p1];
				part[p1] = tpart;
			}

			p0++;
		}
		//printf("[%d] p0= %d, p1 = %d , NPART = %d = %d + %d - (%e %e) (%e %e) (%e %e)\n", PROC_RANK, p0, p1, NPART, NPART_IN, npart_ex, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);

		MPI_Isend(&npart_ex, 1, MPI_INT, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);


		//		for (n=0; n<NPART; n++)
		//			part[n].tag = 0;
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
		if (nrecv_total > MAXNPART_DOM_SUD) {
			printf(" error upload_bnd %d %d\n", nrecv_total, MAXNPART_DOM_SUD);
			exit(0);
		}
		NPART_DOM_SUD = nrecv_total;
	}

	if ( COM_DOM ) {	
		MPI_Isend(part+NPART_IN, npart_ex, Bodytype, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
		//	MPI_Isend(part+NPART_IN, npart_ex*sizeof(Body), MPI_BYTE, PROC_RANK_SUDOM, 1, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		int ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(part_dom_sud+ip, nrecv[n], Bodytype, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			//	MPI_Recv(part_dom_sud+ip, nrecv[n]*sizeof(Body), MPI_BYTE, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}


	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);		
	}


}

void partexch(){
	if ( isSUDOM ) {
		int thisidx = SUDOM_RANK;		
		int rank_send[8];
		int rank_recv[8];
		int tr, ts;	
		int x_this, y_this;
		x_this = thisidx / NSIDE1SUDOM;
		y_this = thisidx % NSIDE1SUDOM;
		//	printf(" <%d> [%d]  %d %d BOX %lf %lf : %lf %lf : %lf %lf\n", thisidx, PROC_RANK, x_this, y_this, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);
		//	printf(" <%d> %d %d\n", thisidx, x_this, y_this);
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

				c++;
			}
		}

		for (n=0; n<8; n++) {
			PROC_RANK_SEND_NGB[n] = rank_send[n];
			PROC_RANK_RECV_NGB[n] = rank_recv[n];

		}
		//	printf(" <%d> send %d %d %d %d %d %d %d %d\n", thisidx, rank_send[0], rank_send[1], rank_send[2], rank_send[3], rank_send[4], rank_send[5], rank_send[6], rank_send[7]);
		//	printf(" <%d> recv %d %d %d %d %d %d %d %d\n", thisidx, rank_recv[0], rank_recv[1], rank_recv[2], rank_recv[3], rank_recv[4], rank_recv[5], rank_recv[6], rank_recv[7]);
		int nrecv[8];
		int nsend[8];
		int i, j;
		m=0;
		cpuType y0, y1, z0, z1;

		y0 = BOX0L;
		y1 = BOX0H;
		z0 = BOX1L;
		z1 = BOX1H;

		int y0i = (int) ( y0 *POS2INT );
		int y1i = (int) ( y1 *POS2INT );
		int z0i = (int) ( z0 *POS2INT );
		int z1i = (int) ( z1 *POS2INT );

		Body* pd = part_dom_sud;

		Body* sendbuff = (Body*)mymalloc(sizeof(Body)*MAXNPART_DOM_SUD*2, 20);
		if (sendbuff == NULL) {
			printf(" ERROR partexch: cannot open sendbuff  \n");
			fflush(stdout);
			exit(0);
		}
		int ip=0;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] < y0i && pd[n].posi[2] < z0i ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[0] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] < y0i && pd[n].posi[2] >= z0i && pd[n].posi[2] < z1i ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[1] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] < y0i && pd[n].posi[2] >= z1i  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[2] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] >= y0i &&  pd[n].posi[1] < y1i && pd[n].posi[2] < z0i  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[3] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] >= y0i &&  pd[n].posi[1] < y1i && pd[n].posi[2] >= z1i  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[4] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] >= y1i && pd[n].posi[2] < z0i ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[5] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] >= y1i && pd[n].posi[2] >= z0i && pd[n].posi[2] < z1i ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[6] = c;
		for (n=0, c=0; n< NPART_DOM_SUD; n++ ){
			if ( pd[n].posi[1] >= y1i && pd[n].posi[2] >= z1i  ) {
				sendbuff[ip] = pd[n];
				c ++;
				ip ++;
			}
		}
		nsend[7] = c;

		int nsend_total = ip;
		// mengchen
		if (nsend_total > MAXNPART_DOM_SUD * 2) {
			printf(" enlarge MAXNPART_DOM_SUD * 2 ! %d %d\n", nsend_total, MAXNPART_DOM_SUD * 2);
			exit(0);
		}


		MPI_Request req[8];
		MPI_Status  sta[8];
		MPI_Status  status[8];

		int nrecv_total = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(&(nsend[m]), 1, MPI_INT, rank_send[m], 2, MPI_COMM_WORLD, &(req[m]) );
		}
		//////////////////////////////
		for (m=0; m<8; m++) {
			MPI_Recv( &(nrecv[m]), 1, MPI_INT, rank_recv[m], 2, MPI_COMM_WORLD, &(status[m]));
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}

		for (m=0; m<8; m++) {
			nrecv_total += nrecv[m];
		}
		NPART_DOM_COM = nrecv_total;
		//qwang
		if (nrecv_total>MAXNPART_DOM) {
			printf(" enlarge MAXNPART_DOM ! %d %d\n", nrecv_total, MAXNPART_DOM);
			exit(0);
		}

		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Isend(sendbuff+ip, nsend[m], Bodytype, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			//	MPI_Isend(sendbuff+ip, nsend[m]*sizeof(Body), MPI_BYTE, rank_send[m], 3, MPI_COMM_WORLD, &(req[m]) );
			ip += nsend[m];
		}
		ip = 0;
		for (m=0; m<8; m++) {
			MPI_Recv( part_dom+ip, nrecv[m], Bodytype, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			//MPI_Recv( part_dom+ip, nrecv[m]*sizeof(Body), MPI_BYTE, rank_recv[m], 3, MPI_COMM_WORLD, &(status[m]));
			ip += nrecv[m];
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&req[m], &sta[m]);		
		}

		myfree(sendbuff, 20);


		warp(part_dom, NPART_DOM_COM);
	}
}

void partdown() {
	int n, c, ip,m;
	int nrecv;
	int nsend[NDOMinSUDOM];
	int nsend_total;
	cpuType width = BOXSIZE/NSIDEMESH;
//	double width = BOXSIZE/NSIDEMESH;
	cpuType x0, x1;
	cpuType x2, x3;
//	double x0, x1;
//	double x2, x3;
	MPI_Status  status;	
	MPI_Status  status1;	

	MPI_Request req[NDOMinSUDOM];
	MPI_Status  sta[NDOMinSUDOM];

	MPI_Request req1[NDOMinSUDOM];
	MPI_Status  sta1[NDOMinSUDOM];
	Body* sendbuff;

	int box_i = (int) (POS2INT * BOXSIZE );

	if ( isSUDOM ) {	
		sendbuff = (Body*)mymalloc(sizeof(Body) * MAXNPART_DOM_SUD * 2, 21);
		if (sendbuff == NULL) {
			printf(" ERROR partdown : cannot open sendbuff  \n");
			fflush(stdout);
			exit(0);
		}
		Body *pd = part_dom_sud;	

		for (n=0; n< NPART_DOM_SUD; n++ ){
			if (pd[n].posi[0] < 0)
				pd[n].posi[0] += box_i;
			if (pd[n].posi[0] >=box_i)
				pd[n].posi[0] -= box_i;
		}

		ip = 0;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			c =0;
			x1 = dom_grp_mesh_start[m] * width;
			x2 = dom_grp_mesh_end[m] * width;

			int x1i  = (int) (POS2INT * x1 );
			int x2i  = (int) (POS2INT * x2 );
			int y1i  = (int) (POS2INT * BOX0L );
			int y2i  = (int) (POS2INT * BOX0H );
			int z1i  = (int) (POS2INT * BOX1L );
			int z2i  = (int) (POS2INT * BOX1H );

			cpuType disp  = 0.0;
			for (n=0; n< NPART_DOM_COM; n++ ){
				if( part_dom[n].posi[0] >= x1i && part_dom[n].posi[0] < x2i ){
					sendbuff[ip] = part_dom[n];
					c ++;
					ip ++;
				}
			}
			for (n=0; n< NPART_DOM_SUD; n++ ){
				if ( pd[n].posi[0] >= x1i && pd[n].posi[0] < x2i
						&& pd[n].posi[1] >= y1i && pd[n].posi[1] < y2i
						&& pd[n].posi[2] >= z1i && pd[n].posi[2] < z2i ) {
					sendbuff[ip] = pd[n];
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

		NPART_EXP = nrecv;
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
			MPI_Isend(sendbuff+ip, nsend[m], Bodytype, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			//	MPI_Isend(sendbuff+ip, nsend[m]*sizeof(Body), MPI_BYTE, dom_grp[m], 4, MPI_COMM_WORLD, &(req1[m]));
			ip += nsend[m];
		}

	}

	if ( COM_DOM ) {
		MPI_Recv(part+NPART_IN, NPART_EXP, Bodytype, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);
		//	MPI_Recv(part+NPART_IN, NPART_EXP*sizeof(Body), MPI_BYTE, PROC_RANK_SUDOM, 4, MPI_COMM_WORLD, &status1);

	}
	if (  isSUDOM ) {	
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait(&(req1[m]), &(sta1[m]));
		}
	}

	if (  isSUDOM ) {	
		myfree(sendbuff, 21);
	}


}

void determine_dom_start() {
	int n;
	int r;
	int idx;
	int dom_start[NDOMinSUDOM];
	for (n=0; n<NDOMinSUDOM; n++) {
		dom_start[n] = dom_grp_mesh_start[n];
	}
	cpuType *npart_mesh;
	cpuType *npart_mesh_tot;

	MPI_Request req;
	MPI_Status  sta;
	MPI_Status  status[NDOMinSUDOM];
	MPI_Request  request[NDOMinSUDOM];

	if ( COM_DOM ) {
		int nproc_sud = NDOMinSUDOM - NDOM_HEAD;
		r = dom_grp_rank;
		int len = dom_grp_mesh_size[r];
		npart_mesh = (cpuType*)mymalloc(sizeof(cpuType)*len, 22);
		if (npart_mesh == NULL) {
			printf(" ERROR determine_dom_start: cannot open npart_mesh  \n");
			fflush(stdout);
			exit(0);
		}

		for ( n =0; n<len; n++  ) {
			npart_mesh[n] = 0.0;
		}
		cpuType dx = BOXSIZE/NSIDEMESH;
	//	double dx = BOXSIZE/NSIDEMESH;
		int intmp = 1;
		for ( n=0; n<NPART; n++) {		
	//		double pnx =(double) ( part[n].posi[0] * INT2POS ) / dx;	
			cpuType pnx =(float) ( part[n].posi[0] * INT2POS ) / dx;	
			idx  = (int)( pnx ) - dom_grp_mesh_start[r] ;
			if (idx >=0 && idx < len )
				npart_mesh[idx] += intmp;
		}
		MPI_Isend(npart_mesh, len, MPI_CPUTYPE, PROC_RANK_SUDOM, 2, MPI_COMM_WORLD, &req);
	}

	if ( isSUDOM ) {
		npart_mesh_tot = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 23);
		if (npart_mesh_tot == NULL) {
			printf(" ERROR determine_dom_start: cannot open npart_mesh_tot  \n");
			fflush(stdout);
			exit(0);
		}
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Recv(npart_mesh_tot+dom_grp_mesh_start[n], dom_grp_mesh_size[n], MPI_CPUTYPE, dom_grp[n], 2, MPI_COMM_WORLD, &(status[n]));
		}
	}

	if ( COM_DOM ) {	
		MPI_Wait(&req, &sta);
	}

	if ( isSUDOM ) {
		int n;
		npart_mesh_tot[0] = 0.0;
		for (n=1; n<NSIDEMESH; n++) {
			npart_mesh_tot[n] += npart_mesh_tot[n-1];
		}
		cpuType np_sud = npart_mesh_tot[NSIDEMESH-1];	
		cpuType np_sep = np_sud/ (NDOMinSUDOM - NDOM_HEAD);

		int c=1;
		dom_start[0] = 0;
		for (n=0; n<NSIDEMESH; n++) {
			if (c > NDOMinSUDOM - NDOM_HEAD )
				break;

			if (  npart_mesh_tot[n] >  c*np_sep  ) {
				dom_start[c] = n;
				c++;
			}
		}

		for (n=0; n<NDOM_HEAD; n++) {
			dom_grp_mesh_start[n] = 0;
			dom_grp_mesh_end[n] = 0;
		}

		c =0 ;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			dom_grp_mesh_start[n] = dom_start[c]; 
			c++;
		}

		for (n=NDOM_HEAD; n<NDOMinSUDOM-1; n++) {
			dom_grp_mesh_end[n] = dom_grp_mesh_start[n+1]; 
		}
		dom_grp_mesh_end[NDOMinSUDOM-1] = NSIDEMESH;

		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			dom_grp_mesh_size[n] = dom_grp_mesh_end[n] -  dom_grp_mesh_start[n]; 
		}

		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Isend(dom_grp_mesh_start, NDOMinSUDOM, MPI_INT, dom_grp[n], 3, MPI_COMM_WORLD, &request[n]);
		}

	}

	if (COM_DOM) {
		MPI_Recv(dom_grp_mesh_start, NDOMinSUDOM, MPI_INT, PROC_RANK_SUDOM, 3, MPI_COMM_WORLD, &sta);
		myfree(npart_mesh, 22);
	}


	if (isSUDOM) {


		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Wait(&request[n], &status[n]);
		}

		myfree(npart_mesh_tot, 23);
	}

	for (n=0; n<NDOMinSUDOM-1; n++) {
		dom_grp_mesh_end[n] = dom_grp_mesh_start[n+1];
	}
	dom_grp_mesh_end[NDOMinSUDOM-1] = NSIDEMESH;

	for (n=0; n<NDOMinSUDOM; n++) {
		dom_grp_mesh_size[n] = dom_grp_mesh_end[n] - dom_grp_mesh_start[n];
	}

	cpuType dx = BOXSIZE / NSIDEMESH;
//	double dx = BOXSIZE / NSIDEMESH;

	BOX_DOM_L[0] = dx * dom_grp_mesh_start[dom_grp_rank];
	BOX_DOM_H[0] = dx * dom_grp_mesh_end[dom_grp_rank];

}

void domain_decomp() {

	determine_dom_start();

	if (isSUDOM ) {
		part_dom_sud = (Body*)mymalloc(sizeof(Body)*MAXNPART_DOM_SUD, 24);
		if (part_dom_sud == NULL) {
			printf(" ERROR domain_decomp: cannot open part_dom_sud  \n");
			fflush(stdout);
			exit(0);
		}
		part_dom = mymalloc(sizeof(Body)* MAXNPART_DOM, 25);	
		if (part_dom == NULL) {
			printf(" ERROR domain_decomp: cannot open part_dom  \n");
			fflush(stdout);
			exit(0);
		}
	}



	partup();

	NPART = NPART_IN;

	partexch();
	partdown();

	NPART = NPART_IN + NPART_EXP;

	int box_dom_l[3], box_dom_h[3];

	box_dom_l[0] = (int) ( BOX_DOM_L[0]  * POS2INT  );
	box_dom_l[1] = (int) ( BOX_DOM_L[1]  * POS2INT  );
	box_dom_l[2] = (int) ( BOX_DOM_L[2]  * POS2INT  );

	box_dom_h[0] = (int) ( BOX_DOM_H[0]  * POS2INT  );
	box_dom_h[1] = (int) ( BOX_DOM_H[1]  * POS2INT  );
	box_dom_h[2] = (int) ( BOX_DOM_H[2]  * POS2INT  );

	domain_center_i[0] = BOX_DOM_L[0]/2 + BOX_DOM_H[0]/2  ;
	domain_center_i[1] = BOX_DOM_L[1]/2 + BOX_DOM_H[1]/2  ;
	domain_center_i[2] = BOX_DOM_L[2]/2 + BOX_DOM_H[2]/2  ;


	int n;
	int c = 0;
	for (n=0; n<NPART; n++) {
		part[n].tag = 1;
		if (	part[n].posi[0] <   box_dom_l[0] ||
			part[n].posi[0] >=  box_dom_h[0] ||
			part[n].posi[1] <   box_dom_l[1] ||
			part[n].posi[1] >=  box_dom_h[1] ||
			part[n].posi[2] <   box_dom_l[2] ||
			part[n].posi[2] >=  box_dom_h[2] )

		{

			printf(" [%d] ERROR DOMAIN :%d tag = %d %d %d %d\n", PROC_RANK, n, part[n].tag, part[n].posi[0], part[n].posi[1], part[n].posi[2]);
			printf("      (%lf,%lf) (%lf,%lf) (%lf,%lf)\n", BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);
			printf("      (%d,%d) (%d,%d) (%d,%d)\n", box_dom_l[0], box_dom_h[0], box_dom_l[1], box_dom_h[1], box_dom_l[2], box_dom_h[2]);
			fflush(stdout);
//			exit(0);

		}
		else {
			c++;
		}

//		if (PROC_RANK == 1 && n< 5) {
//			printf(" %d parti: %f %f %f\n",n,part[n].posi[0]*INT2POS, part[n].posi[1]*INT2POS, part[n].posi[2]*INT2POS);
//			printf(" %d part : %f %f %f\n",n, part[n].pos[0], part[n].pos[1], part[n].pos[2]);
//		}

	}

	if ( isSUDOM ) {
		myfree(part_dom_sud, 24);
		myfree(part_dom, 25);
	}
}


#endif
