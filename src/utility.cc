#include "utility.h"

#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

double dtime() // return double style timestamp using gettimeofday().
{
	double tseconds;
	struct timeval mytime;

	gettimeofday(&mytime, NULL);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);

	return (tseconds);

} /* dtime() */

void make_dir(const char* path, const char* name, int r_snap){
	struct stat st = {0};
	char name_dir[128];
	sprintf(name_dir, "%s/%s_%d", path, name, r_snap);
	if (stat(name_dir, &st) == -1 ) {
//		printf(" mkdir `%s'\n", name_dir);
		mkdir(name_dir, 0700);
	}
}
