#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ising.h"
#include <sys/time.h>


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

    if(argc !=4){
        printf("Error Executing... \nUsage: %s N k e\nWhere:\nN   Number of elements on NxN array\nk   Number of iterations\nWhere e = 0 for serial 1 for parallel execution\n",argv[0]);
	printf("serial execution generates the conf binaries for parallel executions.\n");        
exit(1);
    }

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    int e = atoi(argv[3]);
    if(e!=0 && e!=1)printf("Wrong third parameter check again\n");
    char file[100];
     FILE *fp;
    int *G;
    if(e==0){
    G = generateSpins(n);

    sprintf(file,"conf-init-%d",n);
    fp= fopen(file,"wb");
    fwrite(G,sizeof(int),n*n,fp);
    fclose(fp);
}
if(e==1){

    sprintf(file,"conf-init-%d",n);
    fp= fopen(file,"rb");
    G=(int *) malloc(n*n*sizeof(int));
    fread(G,sizeof(int),n*n,fp);
    fclose(fp);
}
    double w[] ={0.004, 0.016, 0.026, 0.016, 0.004,
    		     0.016, 0.071, 0.117, 0.071, 0.016,
			     0.026, 0.117,  0 , 0.117, 0.026,
			     0.016, 0.071, 0.117, 0.071, 0.016,
			     0.004, 0.016, 0.026, 0.016, 0.004};
    ising(G,w,k,n);
    sprintf(file,"conf-final-%d",n);
    if(e==0){
    fp=fopen(file,"wb");
    fwrite(G,sizeof(int),n*n,fp);
}
    if(e==1){
    fp=fopen(file,"rb");
    int *ans = (int*)malloc(n*n*sizeof(int));
    fread(ans,sizeof(int),n*n,fp);
    
    int c=0;
    int r =0;
    for(int i =0;i<(n*n);i++) {
        if(G[i]!=ans[i]){
            printf("Failed on %d,%d value %d, real value %d\n",i/n,i%n,G[i],ans[i]);
            r=-1;
	    break;
        }
        else{
            c++;
        }
    }
    free(ans);
if(r ==0)printf("Test succeed for size n = %d and k = %d\n",n,k);

}
    free(G);
    fclose(fp);
 
    return 0;
}
    
