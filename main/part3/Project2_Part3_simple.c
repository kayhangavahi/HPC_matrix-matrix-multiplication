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



//ijk algorithm
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
 
 //jik algorithm
void dgemm_2 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (j=0; j<n; j++)  {
		for (i=0; i<n; i++) { 
			register double r = c[i*n+j] ; 
			for (k=0; k<n; k++) 
				r += a[i*n+k] * b[k*n+j]; 
			c[i*n+j] = r;


	}
}
				
	
 }
  //kij algorithm
void dgemm_3 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (k=0; k<n; k++)  
		for (i=0; i<n; i++) { 
			register double r = a[i*n+k] ; 
			for (j=0; j<n; j++) 
				c[i*n+j] += r * b[k*n+j]; 
			


	}
	
				
	
 }
   //ikj algorithm
void dgemm_4 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (i=0; i<n; i++)  
		for (k=0; k<n; k++) { 
			register double r = a[i*n+k] ; 
			for (j=0; j<n; j++) 
				c[i*n+j] += r * b[k*n+j]; 
			


	}
	
				
	
 }
 
    //jki algorithm
void dgemm_5 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (j=0; j<n; j++)  
		for (k=0; k<n; k++) { 
			register double r = b[k*n+j] ; 
			for (i=0; i<n; i++) 
				c[i*n+j] += r * a[i*n+k]; 
			


	}
	
				
	
 }
 
     //kji algorithm
void dgemm_6 (double *a, double *b, double *c, int n) {
	
	
	int i,j,k;
	for (k=0; k<n; k++)  
		for (j=0; j<n; j++) { 
			register double r = b[k*n+j] ; 
			for (i=0; i<n; i++) 
				c[i*n+j] += r * a[i*n+k]; 
			


	}
	
				
	
 }
 




void main()
{

	int n;

	n=2048;
	
	double *a, *b, *c1,*c2,*c3,*c4,*c5,*c6;
	double t1,t2,t3,t4,t5,t6;
	struct timespec start, end;
	
	
	printf ("\n");
	printf ("Matrix Size is: %d\n", n);
	printf("-------------------------------|code is running|---------------------------------------\n");
	
	a = (double *) malloc(n*n*sizeof(double));
	b = (double *) malloc(n*n*sizeof(double));
	c1 = (double *) malloc(n*n*sizeof(double));	
	c2 = (double *) malloc(n*n*sizeof(double));	
	c3 = (double *) malloc(n*n*sizeof(double));	
	c4 = (double *) malloc(n*n*sizeof(double));	
	c5 = (double *) malloc(n*n*sizeof(double));	
	c6 = (double *) malloc(n*n*sizeof(double));		
	
	double random1;
	double random2;	
	
	int i;
	for(i=0;i<(n*n);i++)
	{	
		random1=(double)rand()/RAND_MAX;
		
		random2=(double)rand()/RAND_MAX;
		a[i]=random1; 
		b[i]=random2;		

	}
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_1(a, b, c1, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t1 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t1);
	printf (" and %f seconds\n", (double) t1/1000000000.0);
	printf("------------------------------|ijk algorithm done|---------------------------------------\n");
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_2(a, b, c2, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t2 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;	
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t2);
	printf (" and %f seconds\n", (double) t2/1000000000.0);
	printf("------------------------------|jik algorithm done|---------------------------------------\n");
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_3(a, b, c3, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t3 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t3);
	printf (" and %f seconds\n", (double) t3/1000000000.0);
	printf("------------------------------|kij algorithm done|---------------------------------------\n");	
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_4(a, b, c4, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t4 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t4);
	printf (" and %f seconds\n", (double) t4/1000000000.0);
	printf("------------------------------|ikj algorithm done|---------------------------------------\n");		
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_5(a, b, c5, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t5 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t5);
	printf (" and %f seconds\n", (double) t5/1000000000.0);
	printf("------------------------------|jki algorithm done|---------------------------------------\n");		
	
	clock_gettime(CLOCK_MONOTONIC, &start);
	dgemm_6(a, b, c6, n);
	clock_gettime(CLOCK_MONOTONIC, &end);
	t6 = (end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec;
	printf ("runtime: %llu nanoseconds", (long long unsigned int) t6);
	printf (" and %f seconds\n", (double) t6/1000000000.0);
	printf("------------------------------|kji algorithm done|---------------------------------------\n");	
	
	
	
	
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

	
	
	int u;
	for (u=0;u<n*n;u++){
		
		diff1[u]=c1[u]-c2[u];
		diff2[u]=c2[u]-c3[u];
		diff3[u]=c3[u]-c4[u];	
		diff4[u]=c4[u]-c5[u];
		diff5[u]=c5[u]-c6[u];
		
	}
	printf ("Max Difference dgemm_1 and dgemm_2: %f\n", MAXf(diff1,n));
	printf ("Max Difference dgemm_2 and dgemm_3: %f\n", MAXf(diff2,n));
	printf ("Max Difference dgemm_3 and dgemm_4: %f\n", MAXf(diff3,n));	
	printf ("Max Difference dgemm_4 and  dgemm_5: %f\n", MAXf(diff4,n));	
	printf ("Max Difference dgemm_5 and dgemm_6: %f\n", MAXf(diff5,n));	
		
	
	
  	
}





