#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ising.h"


int test(double *w,int k);

int main(int argc, char *argv[]){
   
    double w[] ={0.004, 0.016, 0.026, 0.016, 0.004,
    		     0.016, 0.071, 0.117, 0.071, 0.016,
			     0.026, 0.117,  0 , 0.117, 0.026,
			     0.016, 0.071, 0.117, 0.071, 0.016,
			     0.004, 0.016, 0.026, 0.016, 0.004};
    int result = test(w,1);
    if(result==1) printf("Test succeed for k=1\n");
    if(result==-1) printf("Test failed for k=1\n");

    result = test(w,4);

    if(result==1) printf("Test succeed for k=4\n");
    if(result==-1) printf("Test failed for k=4\n");

    result = test(w,11);

    if(result==1) printf("Test succeed for k=11\n");
    if(result==-1) printf("Test failed for k=11\n");
}

int test(double w[],int k){
    int n=517;
    FILE *fp;
    fp = fopen("test/conf-init.bin","rb");
    int *G = (int *)malloc(n*n*sizeof(int));
    fread(G,sizeof(int), n*n,fp);
    fclose(fp);
    ising(G,w,k,n);

    char file[100];
    sprintf(file,"test/conf-%d.bin",k);
    fp = fopen(file,"rb");
    int *ans=(int *)malloc(n*n*sizeof(int));
    fread(ans,sizeof(int),n*n,fp);
    fclose(fp);
    
    int c=0;
    for(int i =0;i<(n*n);i++) {
        if(G[i]!=ans[i]){
            printf("Failed on %d,%d value %d, real value %d\n",i/n,i%n,G[i],ans[i]);
            return -1;
        }
        else{
            c++;
        }
    }

    free(G);
    free(ans);
    return 1;
}