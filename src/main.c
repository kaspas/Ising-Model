#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ising.h"
#include <sys/time.h>

int randomNumber(int min,int max){
    int range = max - min + 1;

    // Largest value that when multiplied by "range"
    // is less than or equal to RAND_MAX
    int chunkSize = (2 + 1) / range; 
    int endOfLastChunk = chunkSize * range;

    int r = rand();
    while(r >= endOfLastChunk){
        r = rand();
    }
    return min + r / chunkSize;
}

int *generateSpins(int n){
    int *G=(int *)calloc(n*n,sizeof(int));
    for(int i =0;i<n*n;i++){
        int r = (int)rand()%2;
        G[i]=r?1:-1;
    }
    return G;
}

int main(int argc, char *argv[]){
    int serial;

    if(argc!=3 && argc !=4){
        printf("Error Executing... \nUsage: %s N k\nWhere:\nN   Number of elements on NxN array\nk   Number of iterations\n",argv[0]);
	printf("For serial execution you can give an extra argument what ever\n");
        exit(1);
    }
    if(argc == 3){
        serial =-1;
    }
    else{
        serial =0;
    }
    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    struct timeval startwtime, endwtime;
    int *G = generateSpins(n);
    
    double w[] ={0.004, 0.016, 0.026, 0.016, 0.004,
    		     0.016, 0.071, 0.117, 0.071, 0.016,
			     0.026, 0.117,  0 , 0.117, 0.026,
			     0.016, 0.071, 0.117, 0.071, 0.016,
			     0.004, 0.016, 0.026, 0.016, 0.004};
    gettimeofday (&startwtime, NULL); 
    ising(G,w,k,n);
    gettimeofday (&endwtime, NULL); 
    double t = (double)((endwtime.tv_usec - startwtime.tv_usec)
				/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    if (serial ==0)
        printf("Elapsed: %f seconds\n",t);

    return 0;
}
    
