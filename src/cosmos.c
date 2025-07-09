#include "photoNs.h"
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/////////////////////////////////////////////////////////////////////////////////////////////////

cpuType a_flat_lcdm_t(cpuType time)
{
	cpuType t_star = 3.0*sqrt(OmegaX0)/20.0;
	cpuType kernel = sinh(t_star * time);
	cpuType a = pow(kernel*kernel*OmegaM0/OmegaX0, 0.33333333f);
	return a;
}

cpuType t_flat_lcdm_a(cpuType a)
{
	cpuType t_star = 3.0*sqrt(OmegaX0)/2.0;
//	cpuType t_star = 3.0*sqrt(OmegaX0)/20.0;
	cpuType a3 = a*a*a;
	cpuType f = OmegaX0/OmegaM0;
	cpuType time = log( sqrt(f*a3) + sqrt(1.0+f*a3) )/t_star;
	return 9.769946*time/Hubble0;
}

cpuType kick_loga(cpuType loga_i, cpuType loga_f) {

	int n;
	int Nblock = 1024;
	cpuType kick_time = 0.0;
	cpuType dloga = (loga_f - loga_i)/Nblock;
	cpuType a_f = exp(loga_f);
	cpuType a_i = exp(loga_i);
	cpuType z1 = 1.0/(a_i);
	cpuType h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time = dloga*z1/h;
	for (n=1; n<Nblock; n++) {
		z1 = 1.0/(exp(loga_i+dloga*n));
		h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
		kick_time += 2.0*(1+n%2)*dloga*z1/h;
	}
	z1 = 1.0/(a_f);
	h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time += dloga*z1/h;
	kick_time /= (3.0);
	return kick_time;
}

cpuType drift_loga(cpuType loga_i, cpuType loga_f)
{
	int n;
	int Nblock = 1024;
	cpuType kick_time = 0.0;
	cpuType dloga = (loga_f - loga_i)/Nblock;
	cpuType a_f = exp(loga_f);
	cpuType a_i = exp(loga_i);
	cpuType z1 = 1.0/(a_i);
	cpuType h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time = dloga*z1*z1/h;
	for (n=1; n<Nblock; n++) {
		z1 = 1.0/(exp(loga_i+dloga*n));
		h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
		kick_time += 2.0*(1+n%2)*dloga*z1*z1/h;
	}
	z1 = 1.0/(a_f);
	h = 0.1*sqrt(OmegaM0*z1*z1*z1 + OmegaX0);
	kick_time += dloga*z1*z1/h;
	kick_time /= (3.0);
	return kick_time;
}

////////////////// NR ran3 ///////////////

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

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

/////////////////// NR ran2 ///////////////
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
#define RAN2NTAB 32
#define NDIV (1+IMM1/RAN2NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum) {
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[RAN2NTAB];
	float temp;

	if (*idum <=0 ) {
		if (-(*idum)<1) *idum = 1;
		else *idum = -(*idum);

		idum2=(*idum);
		for (j=RAN2NTAB+7;j>=0; j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j<RAN2NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1) - k*IR1;
	if (*idum < 0) *idum += IM1;
	k = idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	
	if (idum2<0) idum2+= IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] =*idum;
	
	if (iy < 1) iy += IMM1;
	if (( temp = AM * iy ) > RNMX) return RNMX;
	else return temp;	

}







