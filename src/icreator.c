/* use Zel'dovich approximation to generate the IC, modified from 2LPTic */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include <limits.h>
#include "photoNs.h"
#include "typedef.h"
static double Norm;
static double Dplus;
long int seedi;


int MODE_IC;


cpuType *rp1, *rp2;

static double UnitLength_in_cm  = 3.085678e21;  

void setup_ic_ps();
void clean_ic_ps();

void pot_grad();
void part_grad(int direct_ic);
double PowerSpec(double k) ;  /* Eisenstein & Hu */


void ic_ZA(cpuType ai, int NSIDEMESH_IC) {
	InitialTime = ai;

	int vproc[2];

	vproc[0] = NSIDE0PROC;
	vproc[1] = NSIDE1PROC;

#ifdef TEST_FLOAT
	float dispf = 2000000.0;
#endif


	NPART_TOTAL = (unsigned long long)(NSIDEMESH_IC) * NSIDEMESH_IC * NSIDEMESH_IC;

	int nside[3] = { NSIDEMESH_IC, NSIDEMESH_IC, NSIDEMESH_IC };


	setup_ic_ps();

	//	printf(" cal + %e\n", PowerSpec(0.02 * 0.001));

	int x_n = NSIDEMESH_IC;
	int y_n = NSIDEMESH_IC/NSIDE0PROC;
	int z_n = NSIDEMESH_IC/NSIDE1PROC;
	int meshsize_proc_pad = ( x_n + NDOMinSUDOM *2)  * y_n * z_n;
	int meshsize_proc =  x_n * y_n * z_n;
	cpuType param[3] = { 1.2*BOXSIZE/NSIDEMESH_IC, BOXSIZE,0};

	cpuType pmass=(OmegaM0*3*0.01)/(8.0*M_PI*GravConst)*pow(BOXSIZE/NSIDEMESH_IC,3.0);
	MASSPART = pmass;
	int n, m, d, npart = meshsize_proc;
	//	Body *pic = (Body*)mymalloc(sizeof(Body)*npart, 6);

	cpuType hubble_ai =  0.1*sqrt(OmegaM0/(ai*ai*ai)+OmegaX0);
	cpuType F_Omega_ai = pow(OmegaM0/(OmegaM0+ai*ai*ai*OmegaX0), 0.6);
	cpuType vel_prefac = ai*hubble_ai*F_Omega_ai;
	vel_prefac /= sqrt(ai);
	if ( 0 == PROC_RANK) {
		printf("\n\n Initial Condition generating...\n");
		printf(" vel_prefac = %lf, hubble_a = %lf, F_Omega = %lf\n", vel_prefac, hubble_ai, F_Omega_ai);
	}
	// to co-vel
	vel_prefac *= sqrt(ai) * ai;


	cpuType dsep = BOXSIZE / NSIDEMESH_IC;


	int iy = PROC_RANK / NSIDE0PROC;
	int iz = PROC_RANK % NSIDE0PROC;
//	int iy = PROC_RANK / NSIDE1PROC;
//	int iz = PROC_RANK % NSIDE1PROC;
	cpuType maxdisp = 0.0;
	cpuType meandisp = 0.0;

	int i, j, k;
	int nrecv[NDOMinSUDOM];	
	int nsend[NDOMinSUDOM];	
	int ndisp[NDOMinSUDOM];	

//	int ns_x = NSIDEMESH_IC / NSIDE0PROC;
//	int ns_y = NSIDEMESH_IC / NSIDE1PROC;

	for (n=0; n<NDOMinSUDOM; n++) {
		nsend[n] = dom_grp_mesh_size[n] * y_n * z_n;
		//	nsend[n] = dom_grp_mesh_size[n] * ns_x * ns_y;
		//		printf("ic_ZA : [%d] (%d) <%d> dom_grp_mesh_size = %d, (%d, %d)\n", PROC_RANK, dom_grp[n], dom_grp_rank,  dom_grp_mesh_size[n], dom_grp_mesh_start[n], dom_grp_mesh_size[n] );
	}

	MPI_Request req[NDOMinSUDOM];
	MPI_Status  sta[NDOMinSUDOM];
	MPI_Status  status[NDOMinSUDOM];

	for (n=0; n<NDOMinSUDOM; n++) {
		MPI_Isend(&nsend[n], 1, MPI_INT, dom_grp[n], 1, MPI_COMM_WORLD, &req[n]);
	}

	for (n=0; n<NDOMinSUDOM; n++) {
		MPI_Recv(&(nrecv[n]), 1, MPI_INT, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
	}
	for (n=0; n<NDOMinSUDOM; n++) {
		MPI_Wait(&req[n], &sta[n]);		
	}


	int npart_proc = 0;
	for (n=0; n<NDOMinSUDOM; n++) {
		npart_proc += nrecv[n];	
	}

	NPART = npart_proc; 

	int buff_len = NPART;
	if (meshsize_proc > buff_len) buff_len = meshsize_proc;

	cpuType *mic = (cpuType*)mymalloc(sizeof(cpuType) * buff_len, 5);
	cpuType *pic = (cpuType*)mymalloc(sizeof(cpuType) * buff_len, 6);


	if (isSUDOM) {
		NPART = (int) ( NPART_TOTAL /( PROC_SIZE - NDOM_HEAD * NSIDE0SUDOM * NSIDE1SUDOM ) );
	}
	if (COM_DOM)
		NPART = (int) ( npart_proc );  

	MAXNPART = (int) (NPART * max_npart_ratio ) ; 

	MAXNPART_BND = (int)  MAXNPART * max_bnd_ratio;
	MAXNPART_BND_SUD = MAXNPART_BND * NDOMinSUDOM;

	MAXNPART_DOM = MAXNPART_BND_SUD;
	MAXNPART_DOM_SUD = MAXNPART_BND_SUD;

	int x_start = dom_grp_mesh_start[dom_grp_rank];
	int x_size = dom_grp_mesh_size[dom_grp_rank];

	int dsx = x_start;

	int ns_X = x_size;
	int ns_Y = NSIDEMESH_IC / NSIDE0SUDOM;
	int ns_Z = NSIDEMESH_IC / NSIDE1SUDOM;
	int dsy = (SUDOM_RANK / NSIDE1SUDOM) * ns_Y;
	int dsz = (SUDOM_RANK % NSIDE1SUDOM) * ns_Z;

//	int iy = PROC_RANK / NSIDE0PROC;
//	int iz = PROC_RANK % NSIDE0PROC;

	unsigned long pid_x_start = (unsigned long) x_start; 
	unsigned long pid_y_start = (unsigned long) dsy ;
	unsigned long pid_z_start = (unsigned long) dsz ;


	if ( COM_DOM  ) {
		part = (Body*)mymalloc(sizeof(Body) * MAXNPART, 9);

		//		unsigned long nstart = PROC_RANK * NSIDEMESH * y_n * z_n; //mengchen add

		unsigned long pid_start =  ns_Y * ns_Z;
		pid_start *= SUDOM_RANK * NSIDEMESH ;
		pid_start += x_start * ns_Y * ns_Z ;

		cpuType dp = BOXSIZE / NSIDEMESH;
		for (n=0, i=0; i<x_size; i++) {
			for (j=0; j<ns_Y; j++) {
				for (k=0; k<ns_Z; k++) {
//					cpuType px = (dsx + i + 0.5) * dp;
//					cpuType py = (dsy + j + 0.5) * dp;
//					cpuType pz = (dsz + k + 0.5) * dp;

					double px = (dsx + i + 0.5) * dp;
					double py = (dsy + j + 0.5) * dp;
					double pz = (dsz + k + 0.5) * dp;

#ifndef INTXYZ
					part[n].pos[0] = px;
					part[n].pos[1] = py;
					part[n].pos[2] = pz;


#ifdef TEST_FLOAT
					part[n].pos[0] += dispf;
					part[n].pos[1] += dispf;
					part[n].pos[2] += dispf;
#endif

#else
					part[n].posi[0] = (int)((double) px*POS2INT+0.5);
					part[n].posi[1] = (int)((double) py*POS2INT+0.5);
					part[n].posi[2] = (int)((double) pz*POS2INT+0.5);

				//	double ddx[3] ;
				//	ddx[0] = part[n].pos[0] - (double)(part[n].posi[0] *INT2POS);
				//	ddx[1] = part[n].pos[1] - (double)(part[n].posi[1] *INT2POS);
				//	ddx[2] = part[n].pos[2] - (double)(part[n].posi[2] *INT2POS);
				//	ddx[0] =(double)(part[n].posi[0] *INT2POS);
				//	ddx[1] =(double)(part[n].posi[1] *INT2POS);
				//	ddx[2] =(double)(part[n].posi[2] *INT2POS);
#endif					
				//	if (fabs(ddx[0]/px-1.0) > 1e-4) printf(" div4: %e %e %d\n", px, (double)(part[n].posi[0] * INT2POS), part[n].posi[0] );
				//	if (fabs(ddx[1]/py-1.0) > 1e-5) printf(" div5: %e %e %d\n", py, (double)(part[n].posi[1] * INT2POS), part[n].posi[1] );
				//	if (fabs(ddx[2]/pz-1.0) > 1e-6) printf(" div6: %e %e %d\n", pz, (double)(part[n].posi[2] * INT2POS), part[n].posi[2] );
					part[n].tag = 1;
					part[n].act = 0;
					part[n].id = pid_start + n; //mengchen add // first version ok 
					part[n].id = ( ( pid_x_start + i ) * NSIDEMESH_IC + pid_y_start + j ) * NSIDEMESH_IC + pid_z_start + k;
					n++;

				}
			}
		}

		NPART = n;

//printf("[%d] id start = %d, npart = %d ,disp = %d\n", PROC_RANK, pid_start, NPART, pid_start +NPART); // first version ok
		if (part == NULL ){
			printf(" error part malloc %d\n", PROC_RANK);
			exit(0);
		}
	}

	rp1 = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH*(NSIDEMESH+1), 7);
	rp2 = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH*(NSIDEMESH+1), 8);

	seedi = SEED_IC ;
	for (n=0; n<NSIDEMESH * (NSIDEMESH + 1); n++ ) {	
		cpuType x1 = ran2(&seedi);
		cpuType x2 = ran2(&seedi);

		while (x1 == 0.0) {
			x1 = ran2(&seedi);
		}
		rp1[n] = x1;
		rp2[n] = x2;
	}
	int disp_int;
	double disp_max, disp_mean;
	double velo_max, velo_mean;
	disp_max = disp_mean = 0.0;
	velo_max = velo_mean = 0.0;
	for (d=0; d<3; d++) {		
		seedi =  SEED_IC - 7*PROC_RANK ;
	
		param[2] = d+1;
		ic_grad_fft(mic, nside, param);

		for (n=0; n<npart; n++) {
			meandisp += fabs(mic[n]);

			if ( isnan( mic[n] ) ) {printf("[%d] d %d nan ", PROC_RANK, d);exit(0);}
			if (mic[n] > maxdisp)
				maxdisp = mic[n] ;
		}
		int idx, c=0;
		for (n=0; n<NDOMinSUDOM; n++) {
			for (i=dom_grp_mesh_start[n]; i<dom_grp_mesh_end[n]; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						idx  = (i * y_n + j ) * z_n + k ; 
						pic[c] = mic[idx];
						c++;
					}
				}
			}
		}

		cpuType *psend = pic;	
		cpuType *precv = mic;
		for (n=0; n<NDOMinSUDOM; n++) {
			MPI_Isend(psend, nsend[n], MPI_FLOAT, dom_grp[n], 1, MPI_COMM_WORLD, &req[n]);
			psend += nsend[n];
		}

		for (n=0; n<NDOMinSUDOM; n++) {
			MPI_Recv(precv, nrecv[n], MPI_FLOAT, dom_grp[n], 1, MPI_COMM_WORLD, &(status[n]));
			precv += nrecv[n];
		}
		for (n=0; n<NDOMinSUDOM; n++) {
			MPI_Wait(&req[n], &sta[n]);		
		}

		if ( COM_DOM  ) {
			c=0;
			int dx, dy, dz;
			dx = dy = dz = 0;

			int dproc = NSIDE1PROC/NSIDE1SUDOM;
			for (n=0; n<NDOMinSUDOM; n++) {
				dy = y_n * (n / dproc );
				dz = z_n * (n % dproc );
				for (i=0; i<ns_X; i++) {
					for (j=dy; j<dy+y_n; j++) {
						for (k=dz; k<dz+z_n; k++) {
							idx  =  (i * ns_Y + j ) * ns_Z + k; 
							pic[idx] = mic[c] ;
							c++;
						}
					}
				}
			}
			for (n=0; n<NPART; n++) {

#ifndef INTXYZ 
				part[n].pos[d] += pic[n];
#ifdef TEST_FLOAT
				part[n].pos[d] -= dispf;
#endif

#else
				part[n].posi[d] += (int) (pic[n]*POS2INT + 0.5);
//				if ( pic[n] ) printf (" mic = %lf posi = %d\n", pic[n], part[n].posi[1] );
#endif
				part[n].vel[d]  = pic[n] * vel_prefac;


				if (fabs(pic[n]) > disp_max) {disp_max = fabs(pic[n]); velo_max = fabs(vel_prefac * pic[n]);disp_int =(int)( pic[n] * POS2INT); }
				disp_mean += (double) (fabs(pic[n]));
				velo_mean += (double) (fabs(pic[n]*vel_prefac));
			}

		//	disp_mean = sqrt(disp_mean);
		//	velo_mean = sqrt(velo_mean);

			disp_mean /= NPART;
			velo_mean /= NPART;
//			printf(" [%d]  d=%d: disp = %lf (%lf), vel = %lf (%lf), %d\n", PROC_RANK, d ,disp_mean, disp_max, velo_mean, velo_max, disp_int);

		}
	}

	if ( COM_DOM  ) {
		part_b = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND, 10);
	}

	clean_ic_ps() ;

	myfree(pic, 6);
	myfree(rp1, 7);
	myfree(rp2, 8);
	myfree(mic, 5);
}

////////////////////
double tk_eh(double k)          /* from Martin White */
{
	double q, theta, ommh2, a, s, gamma, L0, C0;
	double tmp;
	double omegam, ombh2, hubble;

	/* other input parameters */
	hubble = Hubble0;

	omegam = OmegaM0;
	ombh2 = OmegaB0 * Hubble0 * Hubble0;
	if ( 0 == OmegaB0)
		ombh2 = 0.04 * Hubble0 * Hubble0;

	k *= (3.085678e24 / UnitLength_in_cm);        /* convert to h/Mpc */

	theta = 2.728 / 2.7;
	ommh2 = omegam * hubble * hubble;
	s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
	a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
		+ 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
	gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
	gamma *= omegam * hubble;
	q = k * theta * theta / gamma;
	L0 = log(2. * exp(1.) + 1.8 * q);
	C0 = 14.2 + 731. / (1. + 62.5 * q);
	tmp = L0 / (L0 + C0 * q * q);
	return (tmp);
}


double PowerSpec_EH(double k)   /* Eisenstein & Hu */
{
	return Norm * k * pow(tk_eh(k) , 2);
}

static int NPowerTable;

static struct pow_table
{
	double logk, logD;
}
*PowerTable;



void read_power_table(char FileWithInputSpectrum[])
{
	FILE *fd;
	char buf[500];
	double k, p;
	if (0 == PROC_RANK) {
		printf(" read power spectrum from '%s'\n", FileWithInputSpectrum);
		fflush(stdout);
	}

	sprintf(buf, FileWithInputSpectrum);
	if(!(fd = fopen(buf, "r")))
	{
		printf("can't read input spectrum in file '%s' on task %d\n", buf, PROC_RANK);
		exit(0);  
	}
	NPowerTable = 0;
	do {
		if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
			NPowerTable++;
		else
			break;
	}  while(1);
	fclose(fd);
	if(PROC_RANK == 0)
	{
		printf("found %d pairs of values in input spectrum table\n", NPowerTable);
		fflush(stdout);
	}
	PowerTable = mymalloc(NPowerTable * sizeof(struct pow_table), 11);
	//  sprintf(buf, FileWithInputSpectrum);
	if(!(fd = fopen(buf, "r")))
	{
		printf("can't read input spectrum in file '%s' on task %d\n", buf, PROC_RANK);
		exit(0); 
	}
	NPowerTable = 0;
	do
	{
		if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
		{
			PowerTable[NPowerTable].logk = log10(k);
			PowerTable[NPowerTable].logD = log10(p);
			NPowerTable++;
		}
		else
			break;
	}
	while(1);
	fclose(fd);
	//  qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);
}



double PowerSpec_Tabulated(double k)
{
	double logk, logD, P, kold, u, dlogk, Delta2;
	int binlow, binhigh, binmid;

	kold = k;

	//  k *= 1000.0 ;  // MPC 2 KPC
	// k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);     /* convert to h/Mpc */
	double unit = (3.085678e24 / UnitLength_in_cm);        /* convert to h/Mpc */
	k *= unit;
	logk = log10(k);

	if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
		return 0;

	binlow = 0;
	binhigh = NPowerTable - 1;

	while(binhigh - binlow > 1)
	{
		binmid = (binhigh + binlow) / 2;
		if(logk < PowerTable[binmid].logk)
			binhigh = binmid;
		else
			binlow = binmid;
	}

	dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

	if(dlogk == 0) {
		printf(" error dlog == 0 \n");
		//    FatalError(777);
		exit(0);
	}


	u = (logk - PowerTable[binlow].logk) / dlogk;

	logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;

	Delta2 = pow(10.0, logD);
	P = Norm * Delta2;



	return P;
}





double PowerSpec(double k)
{
	double power, alpha, Tf;

	if ( 1 == MODE_IC){
		power =  PowerSpec_EH(k);
		if (isnan(power))
			printf("power %lf\n", power);

		power *= pow(k, PrimordialIndex - 1.0);
		if (isnan(power))
			printf("k %lf Pri %lf pow %lf check %lf\n", k, PrimordialIndex - 1.0, pow(k, PrimordialIndex - 1.0), pow(0, -0.03));
	}
	if ( 2 == MODE_IC)
		power = PowerSpec_Tabulated(k);


	return power;
}

double r_tophat;

double sigma2_int(double k)
{
	double kr, kr3, kr2, w, x;

	kr = r_tophat * k;
	kr2 = kr * kr;
	kr3 = kr2 * kr;

	if(kr < 1e-8)
		return 0;

	w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
	x = 4 * M_PI * k * k * w * w * PowerSpec(k);

	return x;
}


double qromb(double (*func)(double), double a, double b);

double TopHatSigma2(double R)
{
	r_tophat = R;

	return qromb(sigma2_int, 0, 500.0 * 1 / R);   /* note: 500/R is here chosen as 
							 integration boundary (infinity) */
}


//////// 2LPTic ///////
double growth_int(double a) {
	return pow(a / (OmegaM0 + (1 - OmegaM0 - OmegaX0) * a + OmegaX0 * a * a * a), 1.5);
}

double growth(double a) {
	double hubble_a;

	hubble_a = sqrt(OmegaM0 / (a * a * a) + (1 - OmegaM0 - OmegaX0) / (a * a) + OmegaX0);

	return hubble_a * qromb(growth_int, 0, a);
}


void clean_ic_ps() {
	if ( 2 == MODE_IC ) {
		myfree(PowerTable, 11);
	}

}

void setup_ic_ps()
{

	if ( 2 == MODE_IC ) {
	//	read_power_table("../demo/pk18_hyper_z127.dat");
		read_power_table("../demo/planck_2018_linear_powspec_z_127.txt");
//		read_power_table("../demo/planck13_z149.dat");
//		read_power_table("../demo/plik18_z149.dat");
		Norm = 1.0e9/(8*M_PI*M_PI*M_PI);
		Dplus = 1.0;
	} 
	if (1 == MODE_IC ) {
		//		res = TopHatSigma2(R8);
		double R8 = 8 * (3.085678e24 / UnitLength_in_cm);    /* 8 Mpc/h */
		Norm = Sigma8 * Sigma8 / TopHatSigma2(R8);
		Dplus = growth(1.0)/growth(InitialTime); // flat lcdm
	}	
}


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran_ic(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	long mj,mk;
	int i,ii,k;
	static int iff=0;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
//////// NR subrutine //////
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=(double *)mymalloc((size_t) ((n+1)*sizeof(double)), 12);
	d=(double *)mymalloc((size_t) ((n+1)*sizeof(double)), 13);

	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {printf("Error in routine polint");exit(111);}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	myfree(d, 13);
	myfree(c, 12);
}

#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define EPS 1.0e-8
#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double), double a, double b, int n);
	//        void nrerror(char error_text[]);
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	printf("Too many steps in routine qromb\n");
	exit(222);
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
//////// NR subrutine //////


void ps_ic_k_(double *kn, double *dre, double *dim) {

	double fac = pow(2*M_PI/BOXSIZE,1.5);
	double k2 = (*kn);
	double phase, k = sqrt(k2);
	double p_of_k = PowerSpec(k);

	//printf(" Dplus = %lf, %lf %lf\n", Dplus, p_of_k, k2);
	double rn;

	double x1 = ran2(&seedi);
	double x2 = ran2(&seedi);

	if (x1 == 0.0) {
		x1 = ran2(&seedi);
	}

	p_of_k *= -log(x1);


	double amp = fac*sqrt(p_of_k)/Dplus ;

	double y1 = ( amp * sin(2*M_PI*x2) );
	double y2 = ( amp * cos(2*M_PI*x2) );
	*dre = y1;
	*dim = y2;

}

void ps_ic_lmn_(int *l, int *m, int *n, double *dre, double *dim) {
	double nor = 2*M_PI / BOXSIZE;
	double fac = pow(2*M_PI/BOXSIZE, 1.5);
	double kx = (*l) * nor;
	double ky = (*m) * nor;
	double kz = (*n) * nor;
	double k2 = kx*kx + ky *ky  + kz *kz;
	double phase, kp = sqrt(k2);
	double p_of_k = PowerSpec(kp);
	double pk;
	int nhalf = NSIDEMESH / 2;

	double rn;
	double x1 = ran2(&seedi);
	double x2 = ran2(&seedi);

	while (x1 == 0.0) {
		x1 = ran2(&seedi);
	}
///printf(" (%d) l,m,n = %d, %d, %d\n", PROC_RANK, *l, *m, *n);
	pk = p_of_k * ( -log(x1) );
#ifdef FIXPK
	pk = p_of_k;
#endif
	double amp = fac*sqrt(pk)/Dplus ;
	double y1 = ( amp * sin(2*M_PI*x2) );
	double y2 = ( amp * cos(2*M_PI*x2) );


	int si, sj, sk;
	si = *l;
	sj = *m;
	sk = *n;
	double sign = 1.0;

	if ( si == 0 ) {
		if ( sj == 0 ) {
			if ( sk < 0 ) {
				sk = -sk; 	
			
				sign = - 1.0;
			}
		}

		if ( sj < 0 ) {
			sj = - sj;
			sk = - sk;
	
			sign = -1.0;
		}

		int idx =  sj * NSIDEMESH  + sk + nhalf;

		x1 = rp1[idx];
		x2 = rp2[idx];

		pk = p_of_k * ( -log(x1) );

#ifdef FIXPK
		pk = p_of_k;
#endif
		amp = fac * sqrt(pk) / Dplus ;

		y1 = sign * amp * sin(2*M_PI*x2) ;
		y2 = amp * cos(2*M_PI*x2);
	}

	*dre = y1;
	*dim = y2;

}

