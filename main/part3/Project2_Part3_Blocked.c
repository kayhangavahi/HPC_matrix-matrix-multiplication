#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define BILLION 1000000000L

double MAXf(double *c, int n){
	
	double max=c[0];
	int i;
	for (i=1;i<n*n; i++)
	{
		if (c[i]>max){
			max=c[i];
		}
	}
	
	return max;
}

double MINf(double *c, int n){
	
	int i;
	double min=c[0];
	for (i=1;i<n*n; i++)
	{
		if (c[i]<min){
			min=c[i];
		}
	}
	
	
	return min;
}



//ijk algorithm no blocking
// This part is only for checking the correctness of the blocking version algorithms
void dgemm_1 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (i=0; i<n; i++)  
		for (j=0; j<n; j++) { 
			register double r = c[i*n+j] ; 
			for (k=0; k<n; k++) 
				r += a[i*n+k] * b[k*n+j]; 
			c[i*n+j] = r;


	}
	
				
	
 }
 
 //ijk algorithm blocked
void dgemm_1B (double *a, double *b, double *c, int n,int B) {
	
	
	int i,j,k,i1,j1,k1;
	for(i=0; i<n; i+=B)
		for(j=0; j<n; j+=B)
			for(k=0; k<n; k+=B)
				for(i1=i; i1<i+B; i1++)
					for(j1=j; j1<j+B; j1++){
						register double r = c[i1*n+j1];
						for(k1=k; k1<k+B; k1++)
							r+=a[i1*n+k1]*b[k1*n+j1];
						c[i1*n+j1]=r;
					}

}

 //jik algorithm blocked
void dgemm_2B (double *a, double *b, double *c, int n,int B) {
	
	
	int i,j,k,i1,j1,k1;
	for(j=0; j<n; j+=B)
		for(i=0; i<n; i+=B)
			for(k=0; k<n; k+=B)
				for(j1=j; j1<j+B; j1++)
					for(i1=i; i1<i+B; i1++){
						register double r = c[i1*n+j1];
						for(k1=k; k1<k+B; k1++)
							r+=a[i1*n+k1]*b[k1*n+j1];
						c[i1*n+j1]=r;
					}

}

 //kij algorithm blocked
void dgemm_3B (double *a, double *b, double *c, int n,int B) {
	
	
	int i,j,k,i1,j1,k1;
	for(k=0; k<n; k+=B)
		for(i=0; i<n; i+=B)
			for(j=0; j<n; j+=B)
				for(k1=k; k1<k+B; k1++)
					for(i1=i; i1<i+B; i1++){
						register double r = a[i1*n+k1];
						for(j1=j; j1<j+B; j1++)
							c[i1*n+j1]+=r*b[k1*n+j1];
						
					}

}

 //ikj algorithm blocked
void dgemm_4B (double *a, double *b, double *c, int n,int B) {
	
	
	int i,j,k,i1,j1,k1;
	for(i=0; i<n; i+=B)
		for(k=0; k<n; k+=B)
			for(j=0; j<n; j+=B)
				for(i1=i; i1<i+B; i1++)
					for(k1=k; k1<k+B; k1++){
						register double r = a[i1*n+k1];
						for(j1=j; j1<j+B; j1++)
							c[i1*n+j1]+=r*b[k1*n+j1];
						
					}

}

 //jki algorithm blocked
void dgemm_5B (double *a, double *b, double *c, int n,int B) {
	
	
	int i,j,k,i1,j1,k1;
	for(j=0; j<n; j+=B)
		for(k=0; k<n; k+=B)
			for(i=0; i<n; i+=B)
				for(j1=j; j1<j+B; j1++)
					for(k1=k; k1<k+B; k1++){
						register double r = b[k1*n+j1];
						for(i1=i; i1<i+B; i1++)
							c[i1*n+j1]+=r*a[i1*n+k1];
						
					}

}

 //kji algorithm blocked
void dgemm_6B (double *a, double *b, double *c, int n,int B) {
	
	
	int i,j,k,i1,j1,k1;
	for(k=0; k<n; k+=B)
		for(j=0; j<n; j+=B)
			for(i=0; i<n; i+=B)
				for(k1=k; k1<k+B; k1++)
					for(j1=j; j1<j+B; j1++){
						register double r = b[k1*n+j1];
						for(i1=i; i1<i+B; i1++)
							c[i1*n+j1]+=r*a[i1*n+k1];
						
					}

}
				







void main()
{

	int n,B;

	n=2048;
	
	double *a, *b, *c1, *c1B, *c2B, *c3B, *c4B, *c5B, *c6B;
	double t1,t1B,t2B,t3B,t4B,t5B,t6B;
	struct timespec start, end;
	
	printf ("matrix size: %d\n", n);
	printf("-------------------------------|code is running|---------------------------------------\n");
	
	a = (double *) malloc(n*n*sizeof(double));
	b = (double *) malloc(n*n*sizeof(double));
	c1 = (double *) malloc(n*n*sizeof(double));	
	c1B = (double *) malloc(n*n*sizeof(double));	
	c2B = (double *) malloc(n*n*sizeof(double));	
	c3B = (double *) malloc(n*n*sizeof(double));	
	c4B = (double *) malloc(n*n*sizeof(double));	
	c5B = (double *) malloc(n*n*sizeof(double));	
	c6B = (double *) malloc(n*n*sizeof(double));	


	double *diff1;
	diff1 = (double *) malloc(n*n*sizeof(double));
	double *diff2;
	diff2 = (double *) malloc(n*n*sizeof(double));
	double *diff3;
	diff3 = (double *) malloc(n*n*sizeof(double));
	double *diff4;
	diff4 = (double *) malloc(n*n*sizeof(double));
	double *diff5;
	diff5 = (double *) malloc(n*n*sizeof(double));
	double *diff6;
	diff6 = (double *) malloc(n*n*sizeof(double));	
	
	double random1;
	double random2;	
	
	int i;
	for(i=0;i<(n*n);i++)
	{	
		random1=(double)rand()/RAND_MAX;
		
		random2=(double)rand()/RAND_MAX;
		a[i]=random1; 
		b[i]=random2;	
		c1[i]=0;

	}
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_1(a, b, c1, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t1 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t1);
	printf (" and %f seconds\n", (double) t1/1000000000.0);
	printf("------------------------------|ijk algorithm done no block|---------------------------------------\n");
	
	for(B=8; B<=2048; B=B*2){
		
	int i;
	for(i=0;i<(n*n);i++)
	{	
		
		c1B[i]=0;
		c2B[i]=0;	
		c3B[i]=0;	
		c4B[i]=0;	
		c5B[i]=0;
		c6B[i]=0;

	}
	
	printf("  \n");
	printf(" B= %d\n", B);
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_1B(a, b, c1B, n, B);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t1B = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t1B);
	printf (" and %f seconds\n", (double) t1B/1000000000.0);
	printf("------------------------------|ijk algorithm done|---------------------------------------\n");
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_2B(a, b, c2B, n, B);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t2B = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t2B);
	printf (" and %f seconds\n", (double) t2B/1000000000.0);
	printf("------------------------------|jik algorithm done|---------------------------------------\n");
	
    clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_3B(a, b, c3B, n, B);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t3B = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t3B);
	printf (" and %f seconds\n", (double) t3B/1000000000.0);
	printf("------------------------------|kij algorithm done|---------------------------------------\n");
	
	
    clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_4B(a, b, c4B, n, B);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t4B = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t4B);
	printf (" and %f seconds\n", (double) t4B/1000000000.0);
	printf("------------------------------|ikj algorithm done|---------------------------------------\n");
	
    clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_5B(a, b, c5B, n, B);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t5B = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t5B);
	printf (" and %f seconds\n", (double) t5B/1000000000.0);
	printf("------------------------------|jki algorithm done|---------------------------------------\n");	

    clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_6B(a, b, c6B, n, B);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t6B = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t6B);
	printf (" and %f seconds\n", (double) t6B/1000000000.0);
	printf("------------------------------|kji algorithm done|---------------------------------------\n");


	

	
	
	int u;
	for (u=0;u<n*n;u++){
		
		diff1[u]=c1[u]-c1B[u];
		diff2[u]=c1B[u]-c2B[u];
		diff3[u]=c2B[u]-c3B[u];
		diff4[u]=c3B[u]-c4B[u];
		diff5[u]=c4B[u]-c5B[u];
		diff6[u]=c5B[u]-c6B[u];

		
	}
	
	printf ("Max Difference dgemm_1 & dgemm_1B: %f\n", MAXf(diff1,n));
	printf ("Max Difference dgemm_1B & dgemm_2B: %f\n", MAXf(diff2,n));
	printf ("Max Difference dgemm_2B & dgemm_3B: %f\n", MAXf(diff3,n));
	printf ("Max Difference dgemm_3B & dgemm_4B: %f\n", MAXf(diff4,n));
	printf ("Max Difference dgemm_4B & dgemm_5B: %f\n", MAXf(diff5,n));
	printf ("Max Difference dgemm_5B & dgemm_6B: %f\n", MAXf(diff6,n));

}
	
  	
}





