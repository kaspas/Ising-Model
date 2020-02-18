#include "ising.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <sys/time.h>

#define SWAP(x, y)  { typeof(x) SWAP = x; x = y; y = SWAP; } 

__global__ void ising_cuda(int *g_new,int *G,float *w,int n,int *C){

    int id = (blockIdx.y * blockDim.y + threadIdx.y)*n+(blockIdx.x * blockDim.x + threadIdx.x);
    
    if((id)<(n*n)){
        int i = id/n;
        int j = id%n;
        float temp=0.0f;
        temp+= (w[0]  *  G[ ((i-2+n)%n)*n+(j-2+n)%n]);
        temp+= (w[1]  *  G[ ((i-2+n)%n)*n+(j-1+n)%n]);
        temp+= (w[2]  *  G[ ((i-2+n)%n)*n+(j+n)%n]);
        temp+= (w[3]  *  G[ ((i-2+n)%n)*n+(j+1+n)%n]);
        temp+= (w[4]  *  G[ ((i-2+n)%n)*n+(j+2+n)%n]);
        temp+= (w[5]  *  G[ ((i-1+n)%n)*n+(j-2+n)%n]);
        temp+= (w[6]  *  G[ ((i-1+n)%n)*n+(j-1+n)%n]);
        temp+= (w[7]  *  G[ ((i-1+n)%n)*n+(j+n)%n]);
        temp+= (w[8]  *  G[ ((i-1+n)%n)*n+(j+1+n)%n]);
        temp+= (w[9]  *  G[ ((i-1+n)%n)*n+(j+2+n)%n]);
        temp+= (w[10] *  G[ ((i+n)%n)*n+(j-2+n)%n]);
        temp+= (w[11] *  G[ ((i+n)%n)*n+(j-1+n)%n]);
        temp+= (w[13] *  G[ ((i+n)%n)*n+(j+1+n)%n]);
        temp+= (w[14] *  G[ ((i+n)%n)*n+(j+2+n)%n]);
        temp+= (w[15] *  G[ ((i+1+n)%n)*n+(j-2+n)%n]);
        temp+= (w[16] *  G[ ((i+1+n)%n)*n+(j-1+n)%n]);
        temp+= (w[17] *  G[ ((i+1+n)%n)*n+(j+n)%n]);
        temp+= (w[18] *  G[ ((i+1+n)%n)*n+(j+1+n)%n]);
        temp+= (w[19] *  G[ ((i+1+n)%n)*n+(j+2+n)%n]);
        temp+= (w[20] *  G[ ((i+2+n)%n)*n+(j-2+n)%n]);
        temp+= (w[21] *  G[ ((i+2+n)%n)*n+(j-1+n)%n]);
        temp+= (w[22] *  G[ ((i+2+n)%n)*n+(j+n)%n]);
        temp+= (w[23] *  G[ ((i+2+n)%n)*n+(j+1+n)%n]);
        temp+= (w[24] *  G[ ((i+2+n)%n)*n+(j+2+n)%n]);

        if(temp>FLT_EPSILON){
            g_new[i*n+j]=1;
            *C=1;
        }
        else if ( temp<-FLT_EPSILON){
            *C=1;
            g_new[i*n+j]=-1;
        }
        else{
            g_new[i*n+j]=G[i*n+j];
        }
      
    }
    
}


//! Ising model evolution
/*!
  \param G      Spins on the square la,ttice             [n-by-n]
  \param s_w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]
  NOTE: Both matrices G and w are stored in row-major format.
*/

void ising(int *G,double *w,int k, int n){
    //int *g =(int *) malloc( n*n*sizeof(int ) );
    int *dev_g,*dev_G;
    float *dev_w;
    int *cond,*dev_C;
	float W[25];
    struct timeval startwtime, endwtime;
for(int i=0;i<25;i++) W[i]=(float)w[i];
    int blockSize;      // The launch configurator returned block size 
   // int minGridSize;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch 
    int gridSize;       // The actual grid size needed, based on input size

    cudaMalloc(&dev_g,n*n*sizeof(int));
    cudaMalloc(&dev_G,n*n*sizeof(int));
    cudaMalloc(&dev_C,n*n*sizeof(int));
    cudaMalloc(&dev_w,25*sizeof(float));
    cudaMemcpy(dev_G,G,n*n*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_w,W,25*sizeof(float),cudaMemcpyHostToDevice);
    
    int condition =0;
    //cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, ising_cuda, 0, n*n); 
    cudaMalloc(&cond,sizeof(int));
    cudaMemcpy(cond,&condition,sizeof(int),cudaMemcpyHostToDevice);
    blockSize=8*8;
    gridSize = (n*n + blockSize - 1) / blockSize;
    //printf("GridSize : %d\n, BlockSize : %d\n",gridSize,blockSize);
    double t=0;
    for(int iter = 0;iter<k;iter++){
        gettimeofday (&startwtime, NULL); 
        ising_cuda<<<gridSize,blockSize>>>(dev_g,dev_G,dev_w,n,cond);
        cudaDeviceSynchronize();
        gettimeofday (&endwtime, NULL); 
        t += (double)((endwtime.tv_usec - startwtime.tv_usec)
				/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
        cudaMemcpy(&condition,cond,sizeof(int),cudaMemcpyDeviceToHost);
        
        SWAP(dev_G,dev_g);
        if(condition!=0){
            condition=0;
            cudaMemcpy(cond,&condition,sizeof(int),cudaMemcpyHostToDevice);
        }
    }
    printf("Elapsed: %f seconds\n",t);

    /*if(k%2==1){
        cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDeviceToHost);
        
    }
    else{
        cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDeviceToHost);
    }*/
    cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDeviceToHost);

    cudaFree(dev_g);
    cudaFree(dev_G);
    cudaFree(dev_w);
    //free(g); 
}
