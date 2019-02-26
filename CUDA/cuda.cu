#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#define Inf 9999
#define inf 9999
#define INF 9999

float *dist;
float *kernel1;
float *kernel2;
float *kernel3;
struct timeval startwtime, endwtime;
double seq_time;

__global__ void floydWarshellKernel2 (float* dist, int k, int n);
__global__ void floydWarshellKernel1(float *dist, int k, int n);
__global__ void floydWarshellKernel3(float *dist, int k, int n);
void floydWarshellSerial (float* graph, float* result, int n);
void printSolution(float* dist, int n);
void Check(float* array1, float* array2, int n);

int main(int argc, char*argv[]){


  int n;

/*Check Arguments*/
  if (argc!=3){
    printf("Error, two arguments are needed. arg1 =  full path of the input txt file"
    " which contains the matrix, arg2 = n, where n X n is matrix dimension \n");
    exit(1);
  }

/*Open file*/
  FILE *inputMatrix;
  inputMatrix=fopen(argv[1], "r+");

/*Check if success*/
  if (inputMatrix==NULL){
    printf("Error opening file. Check file permissions\n");
    exit(1);
  }

  n= 1<<atoi(argv[2]);

  float *graph;
  graph=(float*)malloc(n*n*sizeof(float));

  //printf("Initial Distance-Matrix between vertices is:\n");

  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      fscanf(inputMatrix, "%f", &graph[j+i*n]);
    }
  }
  printf("\n");
  fclose(inputMatrix);

  kernel1=(float*)malloc(n*n*sizeof(float));
  kernel2=(float*)malloc(n*n*sizeof(float));
  kernel3=(float*)malloc(n*n*sizeof(float));

  //////////////////Serial algorithm/////////////////////

  gettimeofday (&startwtime, NULL);
  float *result;
  result=(float*)malloc(n*n*sizeof(float));

  floydWarshellSerial(graph, result, n);

  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);

  printf("Serial time: %f\n", seq_time);

  ///////////////Cuda kernel 1 algorithm///////////////////


  int blocksize= 4;
  dim3 dimBlock( blocksize, blocksize );
  dim3 dimGrid( n/dimBlock.x, n/dimBlock.y );

  gettimeofday (&startwtime, NULL);
  cudaMalloc((void**)&dist, n*n*sizeof(float));
  cudaMemcpy(dist, graph, n*n*sizeof(float), cudaMemcpyHostToDevice);
  for (int k=0; k<n; k++) floydWarshellKernel1<<<dimGrid, dimBlock>>>(dist, k, n);
  cudaMemcpy(kernel1, dist, n*n*sizeof(float), cudaMemcpyDeviceToHost);
  gettimeofday (&endwtime, NULL);
  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);

  Check(kernel1, result, n);
  printf("Cuda time, kernel 1: %f \n", seq_time);

  ///////////////Cuda kernel 2 algorithm///////////////////
  blocksize= 4;
  dim3 dimBlock2( blocksize, blocksize );
  dim3 dimGrid2(n, (n+blocksize-1)/blocksize);

  gettimeofday (&startwtime, NULL);
  cudaMemcpy(dist, graph, n*n*sizeof(float), cudaMemcpyHostToDevice);
  for (int k=0; k<n; k++) {
    floydWarshellKernel2<<<dimGrid2, dimBlock2>>>(dist, k, n);
  }
  cudaMemcpy(kernel2, dist, n*n*sizeof(float), cudaMemcpyDeviceToHost);
  gettimeofday (&endwtime, NULL);
  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);

  Check(kernel2, result, n);
  printf("Cuda time, kernel 2: %f\n", seq_time);

  ///////////////Cuda kernel 3 algorithm///////////////////


  blocksize= 128;
  int gridsize=n/blocksize;

  gettimeofday (&startwtime, NULL);
  cudaMalloc((void**)&dist, n*n*sizeof(float));
  cudaMemcpy(dist, graph, n*n*sizeof(float), cudaMemcpyHostToDevice);
  for (int k=0; k<n; k++) {
    floydWarshellKernel3<<<gridsize, blocksize>>>(dist, k, n);
  }
  cudaMemcpy(kernel3, dist, n*n*sizeof(float), cudaMemcpyDeviceToHost);
  gettimeofday (&endwtime, NULL);
  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);

  Check(kernel3, result, n);
  printf("Cuda time, kernel 3: %f\n", seq_time);

}

__global__ void floydWarshellKernel2(float *dist, int k, int n)
{

  int j=blockIdx.y*blockDim.y + threadIdx.y;
	if(j>=n) return;
	int idx=n*blockIdx.x+j;

	__shared__ float best;

	if(threadIdx.y==0) best=dist[n*blockIdx.x+k];
	__syncthreads();

	if(dist[k*n+j]+best<dist[idx]){
		dist[idx]=dist[k*n+j]+best;
	}
}

__global__ void floydWarshellKernel1(float *dist, int k, int n)
{


    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i*n + j;

    if (i<n && j<n){
      if (dist[k+i*n] + dist[j+k*n] < dist[index]){
        dist[index] = dist[k+i*n]+dist[j+k*n];
      }
    }
    __syncthreads();

}

void floydWarshellSerial(float *graph, float *result, int n)
{
    for (int i = 0; i<n; i++){
      for (int j = 0; j<n; j++){
        result[j+i*n] = graph[j+i*n];
      }
    }

    for (int k=0; k<n; k++)
    {
        // Pick all vertices as source one by one
        for (int i=0; i<n; i++)
        {
            // Pick all vertices as destination for the
            // above picked source
            for (int j=0; j<n; j++)
            {

                // If vertex k is on the shortest path from
                // i to j, then update the value of Distance-Matrix[i][j]
                if (result[k+i*n] + result[j+k*n] < result[j+i*n])
                    result[j+i*n] = result[k+i*n]+result[j+k*n];
            }
        }
    }
}

__global__ void floydWarshellKernel3(float *dist, int k, int n)
{

    // int i=blockIdx.x * blockDim.x + threadIdx.x;
    // int index=0;
    // for (int j=0; j<n; j++){
    //   if(i<n){
    //     index=i*n+j;
    //     if (dist[k+i*n] + dist[j+k*n] < dist[index]){
    //       dist[index] = dist[k+i*n]+dist[j+k*n];
    //     }
    //   }
    // }
    int j= blockDim.x*blockIdx.x+threadIdx.x;
    if (j>=n) return;

    __shared__ float best;

    for (int i=0; i<n; i++){
      if (threadIdx.x==0) best=dist[n*i+k];
      __syncthreads();

      if(best+dist[k*n+j]<dist[n*i+j]){
        dist[n*i+j]=best+dist[k*n+j];
      }
    }

}

void printSolution(float* dist, int n){
  printf ("Following matrix shows the shortest distances"
            " between every pair of vertices \n");
  printf("\n");
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
                printf ("%f ", dist[j+i*n]);
        }
        printf("\n");
    }
}

void Check(float* array1, float* array2, int n){
  for (int i=0; i<n*n; i++){
    if (array1[i]!=array2[i]){
      printf("Incorrect Solution\n");
      exit(1);
    }
  }
  printf("Correct Solution\n");
}
