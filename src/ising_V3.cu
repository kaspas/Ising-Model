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
#define CUDA_SAFE_CALL(call)                                          \
do {                                                                  \
    cudaError_t err = call;                                           \
    if (cudaSuccess != err) {                                         \
        fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
                 __FILE__, __LINE__, cudaGetErrorString(err) );       \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
} while (0)
// Macro to catch CUDA errors in kernel launches
#define CHECK_LAUNCH_ERROR()                                          \
do {                                                                  \
    /* Check synchronous errors, i.e. pre-launch */                   \
    cudaError_t err = cudaGetLastError();                             \
    if (cudaSuccess != err) {                                         \
        fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
                 __FILE__, __LINE__, cudaGetErrorString(err) );       \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
    /* Check asynchronous errors, i.e. kernel failed (ULF) */         \
    err = cudaDeviceSynchronize();                                    \
    if (cudaSuccess != err) {                                         \
        fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
                 __FILE__, __LINE__, cudaGetErrorString( err) );      \
        exit(EXIT_FAILURE);                                           \
    }                                                                 \
} while (0)

/*


    int dest = threadIdx.x + (threadIdx.y * blockDim.x);
    int destY = dest / (blockDim.x+4);
    int destX = dest % (blockDim.x+4);

    int srcY = destY + (blockIdx.y * blockDim.y) - 2;
    int srcX = destX + (blockIdx.x * blockDim.x) - 2;
    int src = (srcX+n)%n + (((srcY+n)%n) * n);

   s_G[destY*(blockDim.x+4)+destX]=G[src];
        dest = threadIdx.x + (threadIdx.y * blockDim.x) + blockDim.x * blockDim.y;
        destY = dest / (blockDim.x+4);
        destX = dest % (blockDim.x+4);

        srcY = destY + (blockIdx.y * blockDim.y) - 2;
        srcX = destX + (blockIdx.x * blockDim.x) - 2;
        src = (srcX+n)%n + (((srcY+n)%n) * n);
        if (destY<(blockDim.y+4) )  s_G[destY*(blockDim.x+4)+destX]=G[src];
*/
 
        /*s_G[threadIdx.x+threadIdx.y*(blockDim.x+4)]=G[(x-2+n)%n+((y-2+n)%n)*n];
        //check right bottom and right and bottom
        //right side with left threads
        if(threadIdx.x<4)
            s_G[blockDim.x+threadIdx.x+threadIdx.y*(blockDim.x+4)]=G[(blockDim.x+x-2+n)%n+((y-2+n)%n)*n];
        //bottom side with up threads
        if(threadIdx.y<4)
            s_G[threadIdx.x+(threadIdx.y+blockDim.y)*(blockDim.x+4)]=G[(x-2+n)%n+((y+blockDim.y-2+n)%n)*n];
        //trinagle bottom right side that left with upper left sid of threads
        if(threadIdx.x<4 && threadIdx.y<4)
            s_G[blockDim.x+threadIdx.x+(threadIdx.y+blockDim.y)*(blockDim.x+4)]=G[(blockDim.x+x-2+n)%n+((y+blockDim.y-2+n)%n)*n];
        */

__constant__ float s_w[25];

__global__ void ising_cuda(int  *g_new,const int  *G,const int n,int *C){

    int dest;    
    int destY;
    int destX;
    int x;
    int y;
    int srcY;
    int srcX;
    int src;
    int i=0,j=0;
    float temp;

    extern __shared__ int8_t s_G[];
    
    for (i =0; i<(n+blockDim.x-1)/(blockDim.x);i+=gridDim.x){
        for(j=0;j<(n+blockDim.y-1)/(blockDim.y);j+=gridDim.y){
            __syncthreads();
	    for(int l=0;l<=(blockDim.x+4)*(blockDim.y+4)/(blockDim.x*blockDim.y);l++){
                dest = threadIdx.x + (threadIdx.y * blockDim.x)  + l*blockDim.x * (blockDim.y);
                destY = dest / (blockDim.x+4);
                destX = dest % (blockDim.x+4);

                srcY = destY + (blockIdx.y * blockDim.y)+j*blockDim.y - 2;
                srcX = destX + (blockIdx.x * blockDim.x)+i*blockDim.x - 2;
                src = (srcX+n)%n + (((srcY+n)%n) * n);
                if (destY<(blockDim.y+4) )  
                    s_G[destY*(blockDim.x+4)+destX]=(int8_t)G[src];
            }
            x = threadIdx.x + (blockIdx.x +i)*blockDim.x;
            y = threadIdx.y + (blockIdx.y +j)*blockDim.y;
          
            __syncthreads();
            if(x<n && y<n){       
                temp = 0.0f;

                temp+= (s_w[0]  * s_G[(threadIdx.y   ) * (blockDim.x+4)  +  (threadIdx.x  )] );
                temp+= (s_w[1]  * s_G[(threadIdx.y   ) * (blockDim.x+4)  +  (threadIdx.x+1)] );
                temp+= (s_w[2]  * s_G[(threadIdx.y   ) * (blockDim.x+4)  +  (threadIdx.x+2)] );
                temp+= (s_w[3]  * s_G[(threadIdx.y   ) * (blockDim.x+4)  +  (threadIdx.x+3)] );
                temp+= (s_w[4]  * s_G[(threadIdx.y   ) * (blockDim.x+4)  +  (threadIdx.x+4)] );
                temp+= (s_w[5]  * s_G[(threadIdx.y+1 ) * (blockDim.x+4)  +  (threadIdx.x  )] );
                temp+= (s_w[6]  * s_G[(threadIdx.y+1 ) * (blockDim.x+4)  +  (threadIdx.x+1)] );
                temp+= (s_w[7]  * s_G[(threadIdx.y+1 ) * (blockDim.x+4)  +  (threadIdx.x+2)] );
                temp+= (s_w[8]  * s_G[(threadIdx.y+1 ) * (blockDim.x+4)  +  (threadIdx.x+3)] );
                temp+= (s_w[9]  * s_G[(threadIdx.y+1 ) * (blockDim.x+4)  +  (threadIdx.x+4)] );
                temp+= (s_w[10] * s_G[(threadIdx.y+2 ) * (blockDim.x+4)  +  (threadIdx.x  )] );
                temp+= (s_w[11] * s_G[(threadIdx.y+2 ) * (blockDim.x+4)  +  (threadIdx.x+1)] );
                temp+= (s_w[13] * s_G[(threadIdx.y+2 ) * (blockDim.x+4)  +  (threadIdx.x+3)] );
                temp+= (s_w[14] * s_G[(threadIdx.y+2 ) * (blockDim.x+4)  +  (threadIdx.x+4)] );
                temp+= (s_w[15] * s_G[(threadIdx.y+3 ) * (blockDim.x+4)  +  (threadIdx.x  )] );
                temp+= (s_w[16] * s_G[(threadIdx.y+3 ) * (blockDim.x+4)  +  (threadIdx.x+1)] );
                temp+= (s_w[17] * s_G[(threadIdx.y+3 ) * (blockDim.x+4)  +  (threadIdx.x+2)] );
                temp+= (s_w[18] * s_G[(threadIdx.y+3 ) * (blockDim.x+4)  +  (threadIdx.x+3)] );
                temp+= (s_w[19] * s_G[(threadIdx.y+3 ) * (blockDim.x+4)  +  (threadIdx.x+4)] );
                temp+= (s_w[20] * s_G[(threadIdx.y+4 ) * (blockDim.x+4)  +  (threadIdx.x  )] );
                temp+= (s_w[21] * s_G[(threadIdx.y+4 ) * (blockDim.x+4)  +  (threadIdx.x+1)] );
                temp+= (s_w[22] * s_G[(threadIdx.y+4 ) * (blockDim.x+4)  +  (threadIdx.x+2)] );
                temp+= (s_w[23] * s_G[(threadIdx.y+4 ) * (blockDim.x+4)  +  (threadIdx.x+3)] );
                temp+= (s_w[24] * s_G[(threadIdx.y+4 ) * (blockDim.x+4)  +  (threadIdx.x+4)] );
                        
                if(temp>FLT_EPSILON){
                    g_new[x+y*n]=1;
                    *C = 1;
                }
                else if(temp<-FLT_EPSILON){
                    g_new[x+y*n]=-1;
                    *C=1;
                }
                else{
                    g_new[x+y*n]=s_G[(threadIdx.y+2)*(blockDim.x+4)+(threadIdx.x+2)];
                }   
                
            }
        }
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
    int *dev_g,*dev_G;
    int *cond;
    struct timeval startwtime, endwtime;
 	float W[25];
for(int i =0;i<25;i++) W[i] = (float)w[i];
    //int blockSize;      // The launch configurator returned block size 
    //int gridSize;       // The actual grid size needed, based on input size
    //int minGridSize;
    CUDA_SAFE_CALL( cudaMalloc(&dev_g,n*n*sizeof(int)));
    CUDA_SAFE_CALL (cudaMalloc(&dev_G,n*n*sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpy(dev_G,G,n*n*sizeof(int),cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpyToSymbol    (  s_w,  W,   sizeof(float)*25  ));
    int condition =0;

    CUDA_SAFE_CALL(cudaMalloc(&cond,sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpy(cond,&condition,sizeof(int),cudaMemcpyHostToDevice));
    int bx =32;
    int by =16;
    //blockSize=bx*by;
    int gx = 16;//(n/4+bx-1)/bx; //for 1024 size best solution
    int gy = 8;//(n/4+by-1)/by;

    int memsize = (bx+4)*(by+4)*sizeof(int8_t);
    //cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, ising_cuda, 0, n*n); 
  //  printf("Grid Size : %d , Block Size : %d\n",minGridSize,blockSize);

    dim3 block(bx,by);
    dim3 grid(gx,gy);
    double t=0;
    for(int iter = 0;iter<k;iter++){
        gettimeofday (&startwtime, NULL); 
        ising_cuda<<<grid,block,memsize>>>(dev_g,dev_G,n,cond);
        CUDA_SAFE_CALL(cudaDeviceSynchronize());
        gettimeofday (&endwtime, NULL); 
        t += (double)((endwtime.tv_usec - startwtime.tv_usec)
				/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
        CUDA_SAFE_CALL(cudaMemcpy(&condition,cond,sizeof(unsigned int),cudaMemcpyDeviceToHost));

        SWAP(dev_G,dev_g);
        if(condition!=0){
            condition=0;
            CUDA_SAFE_CALL(cudaMemcpy(cond,&condition,sizeof(unsigned int),cudaMemcpyHostToDevice));
        }
    }
            printf("Elapsed: %f seconds\n",t);

    CUDA_SAFE_CALL(cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaFree(cond));
    CUDA_SAFE_CALL(cudaFree(dev_g));
    CUDA_SAFE_CALL(cudaFree(dev_G));
}
