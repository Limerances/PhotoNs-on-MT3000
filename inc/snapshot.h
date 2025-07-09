#ifndef SNAPSHOT_H
#define SNAPSHOT_H


void input_snapshot(char fname[], int timestamp, int rank ) ;
void output_snapshot(cpuType redshift,  int timestamp, int rank ); 
void output_snapshot_new(cpuType redshift, int timestamp, int rank ) ;

void input_snapshot_compress(char fname[], int timestamp, int rank ) ;
void output_snapshot_compress(cpuType redshift,  int timestamp, int rank ); 

void read_GadgetHeader(const char filename[], int *npart_infile);
void read_Particle_Gadget2(char filename[], int n_start, int n_count);
void read_Particle(char filename[], int n_start, int n_count);
void read_Particle_text(char filename[], int n_start, int n_count);
void read_Particle_Gadget2_mfile(char filename[], int n_start, int n_end, int ip);
void read_Particle_Gadget2_domain(char filename[],int *part_start);

void npart_infile(const char filename[], int ith, int np[]);
void write_snapshot(int type, int timestamp, int rank ); 
//void write_Particle_Gadget2(char filename[], int n_start, int n_count);
void write_Particle_Gadget2(char filename[], int n_start, int n_count, cpuType Redshift_Time);

#endif
