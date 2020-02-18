#include "ising.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>

#define SWAP(x, y)  { typeof(x) SWAP = x; x = y; y = SWAP; } 

int val(double value,int example,int *condition) 
{ 
	int e = (int ) (value*(1/(fabs(value)-DBL_EPSILON)));
	int r = e+example;
    int result =e+r%2;
    *condition+=fabs(result);
	return e+r%2;
} 

//! Ising model evolution
/*!
  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]
  NOTE: Both matrices G and w are stored in row-major format.
*/
void ising(int *G,double *w,int k, int n){
    int *g =(int *) malloc( n*n*sizeof(int ) );
    int condition =0;
    struct timeval startwtime, endwtime;
double t = 0;
    for(int iter = 0;iter<k;iter++){
    gettimeofday (&startwtime, NULL); 
        for (int z = 0;z<n*n;z++){
            int i = z/n;
            int j = z%n;
            double temp = 0.0;
            for(int l = 0;l<25;l++){
                if (l==12)l++;
                int r = l/5;
                int c = l%5;
                int y = ((r-2) + i + n) % n;
                int x = ((c-2) + j + n) % n;
                temp += w[r*5+c]*G[y*n+x];
            }
            if(temp>DBL_EPSILON){
                g[i*n+j]=1;
                condition =1;
            }
            else if (temp<-DBL_EPSILON){
                g[i*n+j]=-1;
                condition = 1;
            }
            else{
                g[i*n+j]=G[i*n+j];
            }
            //g[i*n+j]=val(temp,G[i*n+j],&condition);
        }
    gettimeofday (&endwtime, NULL); 
	t += (double)((endwtime.tv_usec - startwtime.tv_usec)
				/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
        if(condition ==0){
            printf("No changes made terminating...\n");
            break;
        }
        SWAP(g,G);
        condition=0;
    }
        printf("Elapsed: %f seconds\n",t);
    if(k%2==1){
        memcpy(g,G,n*n*sizeof(int));
        SWAP(g,G);
    }
    free(g); 
}
