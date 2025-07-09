#include "photoNs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bitset>
#include <fstream>
#include <iostream>

using namespace std;

#define NBITS 20
#define ERR 1.0e-6

#define MAXINT (1<<NBITS)

#define bLEN (256 * NBITS * 3)
#define SZ (sizeof(unsigned long) * 8)

//#define W2F
#define CHECK

template<typename T>
struct CmprMeta {
	T init;
	uint8_t nbits;
};

template<typename T>
void bit_fwrite(T* i_part[], int np, CmprMeta<T>* meta, ofstream& fout, int dim) {
	int p, d, i;
	bitset<bLEN> buf;
	i = 0;
	for (d=0; d<dim; ++d) {
		for (p=0;p<np;++p) {
			i_part[d][p] -= meta[d].init;
			for (int j=0; j<meta[d].nbits; ++j) {
				buf.set(i++, i_part[d][p]&((T)1));
				i_part[d][p] >>= 1;
			}
		}
	}
	bitset<SZ> bitvec;
	unsigned long ul;
	fout.write((char*)&np, sizeof(int));
	for (d=0; d<dim; ++d) {
		fout.write((char*)&meta[d], sizeof(CmprMeta<T>));
	}
	for (int j=0; j<(i+SZ-1)/SZ; ++j) {
		for (int k =0; k<SZ; ++k) {
			bitvec[k] = buf[j*SZ + k];
		}
		ul = bitvec.to_ulong();
		fout.write((char*)&ul, sizeof(unsigned long));	
	}
	return;
}

template<typename T>
void bit_fread(T* i_part[3], int& np, CmprMeta<T>* meta, ifstream& fin, int dim) {
	int p, d, i, blen;
	bitset<bLEN> buf;
	fin.read((char*)&np, sizeof(int));
	blen = 0;
	for (d=0; d<dim; ++d) {
		fin.read((char*)&meta[d], sizeof(CmprMeta<T>));
		blen += meta[d].nbits * np;
	}
	unsigned long ul;
	i = 0;
	for (int j=0; j<(blen+SZ-1)/SZ; ++j) {
		fin.read((char*)&ul, sizeof(unsigned long));	
		bitset<SZ> bitvec(ul);
		for (int k =0; k<SZ; ++k) {
			buf[i++] = bitvec[k];
		}
	}
	i = blen - 1;
	for (d=dim-1; d>=0; --d) {
		for (p=np-1;p>=0;--p) {
			i_part[d][p] = 0;
			for (int j=0; j<meta[d].nbits; ++j) {
				i_part[d][p] <<= 1;
				i_part[d][p] |= buf[i--];
			}
			i_part[d][p] += meta[d].init;
		}
	}
	return;
}

extern "C" void build_localtree(int, cpuType*, cpuType*);
extern "C" void fmm_deconstruct();
extern "C" int cmprs_fread(Body* part, Pack* leaf, char* fn);

void cmprs_fwrite_pos(Body* part, int NPART, ofstream& fout) {
	int n, p, d, i, cnt;
	int ip, np;
	float* f_part[3];
	unsigned int* i_part[3]; 
	unsigned int minv[3], maxv[3];
	CmprMeta<unsigned int> meta[3];	
	for (d=0; d<3; ++d) {
		f_part[d] = (float*)mymalloc(sizeof(float) * MAXLEAF, 301+d);
		if (f_part[d] == NULL) {
			printf(" ERROR cmprs: cannot open f_part[%d]  \n", d);
			fflush(stdout);
			exit(0);
		}
		i_part[d] = (unsigned int*)mymalloc(sizeof(unsigned int) * MAXLEAF, 304+d);
		if (i_part[d] == NULL) {
			printf(" ERROR cmprs: cannot open i_part[%d]  \n", d);
			fflush(stdout);
			exit(0);
		}
	}
	// Positions
	cnt = 0;
	for (n=first_leaf; n<last_leaf; ++n) {
		for (d=0; d<3; ++d) {
			minv[d] = MAXINT;
			maxv[d] = 0;
		}
		ip = leaf[n].ipart;
		np = leaf[n].npart;
		for (d=0; d<3; ++d) {
			for (p=ip;p<ip + np;++p) {

				if (part[p].pos[d] < 0.0 || part[p].pos[d] >= BOXSIZE) {
					if (PROC_RANK == 1)
						printf("[%d] check1 pos [%d] part %lf\n", PROC_RANK, d, part[p].pos[d]);
					return;
				}
				f_part[d][p-ip] = part[p].pos[d] * ((float)MAXINT / BOXSIZE);

				i_part[d][p-ip] = (unsigned int)f_part[d][p-ip];
				minv[d] = minv[d]<i_part[d][p-ip]?minv[d]:i_part[d][p-ip];
				maxv[d] = maxv[d]>i_part[d][p-ip]?maxv[d]:i_part[d][p-ip];
			}
		        meta[d].init = minv[d];	
			int lg = log2(maxv[d]-minv[d]);
			meta[d].nbits = lg + 1;
			cnt += meta[d].nbits;
		}
		bit_fwrite(i_part, np, meta, fout, 3);
	}
	for (d=0; d<3; ++d) {
		myfree(f_part[d], 301+d);
		myfree(i_part[d], 304+d);
	}
	if (PROC_RANK == 1)
		printf("[%d] output compress pos done bits %lf\n", PROC_RANK, (float)cnt / 3.0 / (last_leaf - first_leaf));
}

void cmprs_fwrite_vel(Body* part, int NPART, ofstream& fout) {
	int n, p, d, i, cnt;
	int ip, np;
	float* f_part[3];
	unsigned int* i_part[3]; 
	unsigned int minv[3], maxv[3];
	CmprMeta<unsigned int> meta[3];	
	for (d=0; d<3; ++d) {
		f_part[d] = (float*)mymalloc(sizeof(float) * MAXLEAF, 301+d);
		if (f_part[d] == NULL) {
			printf(" ERROR cmprs: cannot open f_part[%d]  \n", d);
			fflush(stdout);
			exit(0);
		}
		i_part[d] = (unsigned int*)mymalloc(sizeof(unsigned int) * MAXLEAF, 304+d);
		if (i_part[d] == NULL) {
			printf(" ERROR cmprs: cannot open i_part[%d]  \n", d);
			fflush(stdout);
			exit(0);
		}
	}
	// Velocities
	float vmax = -1.0e6, vmin = 1.0e6;
	for (n=first_leaf; n<last_leaf; ++n) {
		ip = leaf[n].ipart;
		np = leaf[n].npart;

		for (p=ip;p<ip + np;++p) {
			for (d=0; d<3; ++d) {
				if (part[p].vel[d] > vmax)
					vmax = part[p].vel[d];
				if (part[p].vel[d] < vmin)
					vmin = part[p].vel[d];
			}
		}
	}
	fout.write((char*)&vmin, sizeof(float));
	fout.write((char*)&vmax, sizeof(float));
	cnt = 0;
	for (n=first_leaf; n<last_leaf; ++n) {
		for (d=0; d<3; ++d) {
			minv[d] = MAXINT;
			maxv[d] = 0;
		}
		ip = leaf[n].ipart;
		np = leaf[n].npart;

		for (d=0; d<3; ++d) {
			for (p=ip;p<ip + np;++p) {
				f_part[d][p-ip] = (part[p].vel[d] - vmin) * ((float)MAXINT / (vmax - vmin));

				i_part[d][p-ip] = (unsigned int)f_part[d][p-ip];
				minv[d] = minv[d]<i_part[d][p-ip]?minv[d]:i_part[d][p-ip];
				maxv[d] = maxv[d]>i_part[d][p-ip]?maxv[d]:i_part[d][p-ip];
			}
		        meta[d].init = minv[d];	
			int lg = log2(maxv[d]-minv[d]);
			meta[d].nbits = lg + 1;
			cnt += meta[d].nbits;
		}
		bit_fwrite(i_part, np, meta, fout, 3);
	}
	for (d=0; d<3; ++d) {
		myfree(f_part[d], 301+d);
		myfree(i_part[d], 304+d);
	}
	if (PROC_RANK == 1)
		printf("[%d] output compress vel done bits %lf\n", PROC_RANK, (float)cnt / 3.0 / (last_leaf - first_leaf));
}

void cmprs_fwrite_id(Body* part, int NPART, ofstream& fout) {
	int n, p, d, i, cnt;
	int ip, np;
	// IDs
	unsigned long* i64_part[1];
	i64_part[0] = (unsigned long*)mymalloc(sizeof(unsigned long) * MAXLEAF, 301);
	unsigned long* id_part = i64_part[0];
	unsigned long idmin, idmax;
	CmprMeta<unsigned long> meta64;
	if (id_part == NULL) {
		printf(" ERROR cmprs: cannot open id_part.\n");
		fflush(stdout);
		exit(0);
	}

	cnt = 0;
	for (n=first_leaf; n<last_leaf; ++n) {
		ip = leaf[n].ipart;
		np = leaf[n].npart;
		idmin = part[ip].id;
		idmax = part[ip].id;
		for (p=ip;p<ip + np;++p) {
			id_part[p-ip] = part[p].id;
			idmin = idmin<id_part[p-ip]?idmin:id_part[p-ip];
			idmax = idmax>id_part[p-ip]?idmax:id_part[p-ip];
		}
		meta64.init = idmin;	
		int lg = log2(idmax-idmin);
		meta64.nbits = lg + 1;
		cnt += meta64.nbits;
		bit_fwrite(i64_part, np, &meta64, fout, 1);
	}
	myfree(id_part, 301);

	if (PROC_RANK == 1)
		printf("[%d] output compress ids done bits %lf\n", PROC_RANK, (float)cnt / (last_leaf - first_leaf));
}

extern "C" void compress_fwrite(Body* part, int NPART, char* fn) {
	int direct = 0;
	cpuType bdl[3];
	cpuType bdr[3];
	build_localtree(direct, bdl, bdr);

        ofstream fout(fn, ios::binary);

	int nleafs = last_leaf - first_leaf;
	fout.write((char*)&nleafs, sizeof(int));

	cmprs_fwrite_pos(part, NPART, fout);
	cmprs_fwrite_vel(part, NPART, fout);
	cmprs_fwrite_id(part, NPART, fout);

	fout.close();

#ifdef CHECK
	int npart = cmprs_fread(part, leaf, fn);
	if (npart != NPART) {
		if (PROC_RANK == 1)
			printf("[%d] check npart %d != NPART %d\n", PROC_RANK, npart, NPART);
		return;

	}
#endif
	fmm_deconstruct();
}

extern "C" int cmprs_fread(Body* part, Pack* leaf, char* fn) {
	int npart, nleafs;
	int cnt;
	int p, d, np, i, j;
	float* f_part[3];
	unsigned int* i_part[3]; 
	CmprMeta<unsigned int> meta[3];	
	for (d=0; d<3; ++d) {
		f_part[d] = (float*)mymalloc(sizeof(float) * MAXLEAF, 307+d);
		if (f_part[d] == NULL) {
			printf(" ERROR cmprs: cannot open f_part[%d]  \n", d);
			fflush(stdout);
			exit(0);
		}
		i_part[d] = (unsigned int*)mymalloc(sizeof(unsigned int) * MAXLEAF, 310+d);
		if (i_part[d] == NULL) {
			printf(" ERROR cmprs: cannot open i_part[%d]  \n", d);
			fflush(stdout);
			exit(0);
		}
	}
        ifstream fin(fn, ios::binary);
	if (!fin.is_open()) {
		printf("cannot open file %s\n", fn);
		exit(0);
	}
#ifdef W2F
	char fn_d[64];
	sprintf(fn_d, "%s.d", fn);
        ofstream fout(fn_d, ios::binary);
	if (!fout.is_open()) {
		printf("cannot open file %s\n", fn_d);
		exit(0);
	}
#endif
	fin.read((char*)&nleafs, sizeof(int));
	npart = 0;
	j = 0;
	cnt = 0;
	for (i=0; i<nleafs; ++i) {
		bit_fread(i_part, np, meta, fin, 3);
		for (p=0;p<np;p++) {
			for (d=0; d<3; d++) {
				f_part[d][p] = (float)i_part[d][p];
				f_part[d][p] /= (float)MAXINT;
#ifdef CHECK
				if (leaf && fabs(part[leaf[i+first_leaf].ipart+p].pos[d] / BOXSIZE - f_part[d][p]) > ERR) {
					cnt++;
				}
#endif
				f_part[d][p] *= BOXSIZE;
				if (!leaf) {
					part[j].pos[d] = f_part[d][p];
				}
#ifdef W2F
				fout.write((char*)&f_part[d][p], sizeof(float));
#endif
			}
			j++;
		}
		npart += np;
	}
	fin.close();
#ifdef W2F
	fout.close();
#endif
#ifdef CHECK
	if (leaf && PROC_RANK == 1)
		printf("[%d] input compress done error %d/%d\n", PROC_RANK, cnt, npart*3);
#endif
	for (d=0; d<3; ++d) {
		myfree(f_part[d], 307+d);
		myfree(i_part[d], 310+d);
	}
	return npart;
}
