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
    int stride =  blockDim.x * gridDim.x;
    
    if((id)<(n*n)){
        for(int ind =id;ind<n*n;ind+=stride){
            int i = ind/n;
            int j = ind%n;
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
            if(temp>0.0001){
                g_new[i*n+j]=1;
                *C=1;
            }
            else if ( temp<-0.0001){
                *C=1;
                g_new[i*n+j]=-1;
            }
            else{
                g_new[i*n+j]=G[i*n+j];
            }
          //  int e = (int)(temp*(1/(fabs(temp)-DBL_EPSILON)));
            //cond=cond%2+e*e;
            //g_new[i*n+j] = e +(e+ G[i*n+j] )%2;
        }
        
        //g_new[i*n+j] =(int)(__dmul_rn(temp,__drcp_ru(fabs(temp)-1e-12)))+( (int)(__dmul_rn(temp,__drcp_ru(fabs(temp)-1e-12))) + G[i*n+j])%2;
    }
    
}


//! Ising model evolution
/*!
  \param G      Spins on the square la,ttice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]
  NOTE: Both matrices G and w are stored in row-major format.
*/

void ising(int *G,double *w,int k, int n){
    //int *g =(int *) malloc( n*n*sizeof(int ) );
    int *dev_g,*dev_G;
    float *dev_w;
    int *cond;
    struct timeval startwtime, endwtime;

    int blockSize;      // The launch configurator returned block size 
    int minGridSize;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch 
    int gridSize;       // The actual grid size needed, based on input size
	float W[25];
for(int i =0;i<25;i++) W[i] = (float)w[i];
    cudaMalloc(&dev_g,n*n*sizeof(int));
    cudaMalloc(&dev_G,n*n*sizeof(int));
    cudaMalloc(&dev_w,25*sizeof(float));
    cudaMemcpy(dev_G,G,n*n*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_w,W,25*sizeof(float),cudaMemcpyHostToDevice);
    
    int condition =0;
    //cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, ising_cuda, 0, n*n); 
    
   // cudaOccupancyMaxActiveBlocksPerMultiprocessor(&gridSize,ising_cuda,blockSize,25*sizeof(double)+sizeof(int));
    //blockSize =256;
    cudaMalloc(&cond,sizeof(int));
    cudaMemcpy(cond,&condition,sizeof(int),cudaMemcpyHostToDevice);
    //gridSize = (n*n/32 + blockSize - 1) / blockSize;
    //printf("GridSize : %d\n, BlockSize : %d\n",minGridSize,blockSize);
    double t=0;
    blockSize=16*16;
    gridSize= (n*n/4+blockSize-1)/(blockSize);
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

    cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDeviceToHost);

    cudaFree(dev_g);
    cudaFree(dev_G);
    cudaFree(dev_w);
}
