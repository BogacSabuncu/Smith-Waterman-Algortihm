#include <cstdio>
#include <cstdlib>
#include <math.h>
#include "SW_GPU.h"

#define GAP_PEN 3

#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"CUDA_SAFE_CALL: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

void InitNucArray (const char* RefSeq, SeqData OutSeq) {
  	int i = 1;
  	*(OutSeq.Seq + 0) = 0;
  	char CurNuc = 'B'; // Input checked by sequence input function
  	while (CurNuc != '\0')
  	{
    	CurNuc = *(RefSeq + i-1);
    	*(OutSeq.Seq + i) = (int8_t) CurNuc;
		//printf("%d ", *(OutSeq.Seq + i));
    	i++;
  	}
	printf("\n");
}

void InitSubGPMat (SeqData RefSeq, SeqData ReadSeq, int* DataMat) {
  uint32_t ReadLength = ReadSeq.SeqLength;
  uint32_t RefLength = RefSeq.SeqLength;
  int32_t i,j;
  for (i = 1; i < ReadLength; i++) {
    for (j = 1; j < RefLength; j++)
	{
      *(DataMat + i*RefLength + j) = (*(RefSeq.Seq + j) == *(ReadSeq.Seq + i)) ? MATCH_HIT : MATCH_MISS;
    }

  }
}

__global__ void kernal_scoring(SeqData *d_read, SeqData *d_ref, int * d_sub_matrix, int *d_score_matrix, int step, int ReadLength, int RefLength)
{
	int i = threadIdx.x +1; //+ blockIdx.x * blockDim.x;	
	int j = step - i;



	int16_t Score1 = 0;
	int16_t Score2 = 0;
 	int16_t Score3 = 0;	
	
	//printf("%d %d \n", i, j);	

        Score1 = *(d_score_matrix + ((i-1)*RefLength + (j-1))) + *(d_sub_matrix + (i*RefLength + j));
        Score2 = *(d_score_matrix + ((i-1)*RefLength + j)) - GAP_PEN;
        Score3 = *(d_score_matrix + (i*RefLength + (j-1))) - GAP_PEN;
	
	//printf(" %d %d %d \n ", Score1, Score2, Score3 );
	
	int16_t MaxScore = 0;
	MaxScore = (Score1 > MaxScore) ? Score1 : MaxScore;
        MaxScore = (Score2 > MaxScore) ? Score2 : MaxScore;
        MaxScore = (Score3 > MaxScore) ? Score3 : MaxScore;
	
	//printf("%d \n ", MaxScore);

	*(d_score_matrix + (i*RefLength + j)) = MaxScore;
}

int iDivUp(int a, int b){

return (a%b != 0) ? (a/b + 1) : (a/b);
}

int main(){
	//// GPU Timing variables///////
	cudaEvent_t gpu_start, gpu_stop, cpu_start, cpu_stop;
	float elapsed_gpu, elapsed_cpu;
	

	///// initializing ///////////

	uint32_t RefLength = 40;  // Includes additional row needed for SW
	uint32_t ReadLength = 16; // Includes additional row needed for SW
	SeqData ReadSeq, RefSeq;
	ReadSeq.SeqLength = ReadLength;
	RefSeq.SeqLength = RefLength;

	  /* Default Sequences for testing purposes */
	  const char ref_seq_def[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A','C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T','C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
	  const char read_seq_def[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};
	
	  /* Initialize input sequences as int8_t arrays */
	  RefSeq.Seq = (uint8_t*) malloc(RefSeq.SeqLength);	// Change to dynamic malloc depending on input
	  ReadSeq.Seq = (uint8_t*) malloc(ReadSeq.SeqLength);  // Change to dynamic malloc depending on input
	  InitNucArray (ref_seq_def, RefSeq);
	  InitNucArray (read_seq_def, ReadSeq);
	  
	  //substution matrix
 	  int* DataMat = (int*) calloc((RefLength)*(ReadLength),sizeof(int));
  	  InitSubGPMat (RefSeq, ReadSeq, DataMat);
	  
	  //scoring matrix
	  
	  int* score_matrix = (int*) calloc(((RefSeq.SeqLength)+1)*((ReadSeq.SeqLength)+1),sizeof(uint32_t));
	  
	  ////// cuda initializing ///////

	 	// Create the cuda events
       	cudaEventCreate(&gpu_start);
	cudaEventCreate(&gpu_stop);
	// Record event  on the default stream
	cudaEventRecord(gpu_start, 0);
	  
	  CUDA_SAFE_CALL(cudaSetDevice(0));
	  
	  SeqData * d_read, * d_ref; 
	  int* d_score_matrix, * d_sub_matrix;
	  
	  //setting the arrays
	  CUDA_SAFE_CALL(cudaMalloc((void **)&d_read, ((ReadSeq.SeqLength)*sizeof(SeqData))));
      CUDA_SAFE_CALL(cudaMalloc((void **)&d_ref, ((RefSeq.SeqLength)*sizeof(SeqData))));
	  CUDA_SAFE_CALL(cudaMalloc((void **)&d_score_matrix, (((RefSeq.SeqLength)+1)*((ReadSeq.SeqLength)+1)*sizeof(uint32_t))));
	  CUDA_SAFE_CALL(cudaMalloc((void **)&d_sub_matrix,((RefSeq.SeqLength)*(ReadSeq.SeqLength))*sizeof(int)));
	  
	  //transfering the arrays
	  CUDA_SAFE_CALL(cudaMemcpy(d_read, ReadSeq.Seq, ((ReadSeq.SeqLength)*sizeof(SeqData)), cudaMemcpyHostToDevice));
	  CUDA_SAFE_CALL(cudaMemcpy(d_ref, RefSeq.Seq, ((RefSeq.SeqLength)*sizeof(SeqData)), cudaMemcpyHostToDevice));
	  CUDA_SAFE_CALL(cudaMemcpy(d_sub_matrix, DataMat, ((RefSeq.SeqLength)*(ReadSeq.SeqLength))*sizeof(int), cudaMemcpyHostToDevice));
	  CUDA_SAFE_CALL(cudaMemcpy(d_score_matrix, score_matrix, (((RefSeq.SeqLength)+1)*((ReadSeq.SeqLength)+1)*sizeof(uint32_t)), cudaMemcpyHostToDevice ));	  

	
	  //kernal function call

	int step;
	int num_iters = RefSeq.SeqLength+ReadSeq.SeqLength-1;
	int max_thread = RefSeq.SeqLength-ReadSeq.SeqLength+1;

	int k = 0;
	//for(k = 0; k < 1600; k++){
	for(step = 1; step< num_iters; step++)
	{
		dim3 dimGrid(iDivUp(step,max_thread));

		kernal_scoring<<< dimGrid, 1>>>(d_read, d_ref, d_sub_matrix, d_score_matrix, step, ReadSeq.SeqLength, RefSeq.SeqLength);
	}
	//}
	//CUDA return
	 CUDA_SAFE_CALL(cudaPeekAtLastError());
	  CUDA_SAFE_CALL(cudaMemcpy(ReadSeq.Seq, d_read, ((ReadSeq.SeqLength)*sizeof(uint32_t)), cudaMemcpyDeviceToHost));
	  CUDA_SAFE_CALL(cudaMemcpy(RefSeq.Seq, d_ref, ((RefSeq.SeqLength)*sizeof(uint32_t)), cudaMemcpyDeviceToHost));
	  CUDA_SAFE_CALL(cudaMemcpy(score_matrix, d_score_matrix, (((RefSeq.SeqLength)+1)*((ReadSeq.SeqLength)+1)*sizeof(uint32_t)), cudaMemcpyDeviceToHost ));	 

	

	// Stop and destroy the timer
	cudaEventRecord(gpu_stop,0);
	cudaEventSynchronize(gpu_stop);
	cudaEventElapsedTime(&elapsed_gpu, gpu_start, gpu_stop);
	printf("\nGPU time: %f (msec)\n", elapsed_gpu);	

	int i,j;
//	  for (i = 0; i < RefLength; i++) {printf(" %3c ",*(RefSeq.Seq + i));}
 /* //printf("\n");
  for (i = 1; i < ReadLength; i++) {
    //printf("%c ",*(ReadSeq.Seq + i));
    for (j = 1; j < RefLength; j++) {
      printf(" %3d ", *(score_matrix + (i*RefLength + j)));
    }
	printf("\n");
	}
*/
	return 0;  

	

}
