/* gcc -o EX Example_MT.c SW.c -msse2 -msse4.1 -pthread -fopenmp -lrt -mavx 
 * -mavx2 for AVX2, needs newer version of GCC than PHO307 version
 */
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include "SW.h"

/* Sequence Lengths */
//#define REF_SEQ_LENGTH 10000001			// 10M
//#define REF_SEQ_LENGTH 1000001			// 1M
#define READ_SEQ_LENGTH 1001
#define REF_SEQ_LENGTH 100001				// 100K
//#define READ_SEQ_LENGTH 33

/* Sequence Filenames */
//char RefFilename[] = "10M.fa";
char RefFilename[] = "100k.fa";
//char RefFilename[] = "1M.fa";
char ReadFilename[] = "1KQuery.fa";
//char ReadFilename[] = "query1.fa";


/* Timing Infastructure Declerations */
#define GIG 1000000000
#define CPG 2.90           // Cycles per GHz
#define OPTIONS 4
#define ITERS 20
#define NUM_THREADS 1
#define NUMSEQ NUM_THREADS
#define NUMEL REF_SEQ_LENGTH * READ_SEQ_LENGTH

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

/* Global PThread barrier */

int main (int argc, char* const argv[]) {
  int i,j,k,ii,jj = 0;
  int Option = 0;
  
  /* Timing Analysis Variables */
  struct timespec diff(struct timespec start, struct timespec end);
  struct timespec time1, time2;
  struct timespec time_stamp[OPTIONS][ITERS];
  int clock_gettime(clockid_t clk_id, struct timespec *tp);

  /* PThread Declerations */
  pthread_t id[NUM_THREADS];

  /* Initialize Input Sequences */
  SeqData* RefSeq_1 = InitArrayExt(RefFilename, REF_SEQ_LENGTH);
  SeqData* ReadSeq_1 = InitArrayExt(ReadFilename, READ_SEQ_LENGTH);
  int ReadLength_1 = ReadSeq_1->SeqLength;
  int RefLength_1 = RefSeq_1->SeqLength;

  /* Declare 8 bit HMat */
  int HMatSize = (RefLength_1 * ReadLength_1) * sizeof(uint8_t);		
  uint8_t* HMat8 = (uint8_t*) _mm_malloc(HMatSize, VEC_ALIGN);
  assert(HMat8 != NULL);

  /* Declare 16 bit HMat */
  HMatSize = (RefLength_1 * ReadLength_1) * sizeof(uint16_t);			
  uint16_t* HMat16 = (uint16_t*) _mm_malloc(HMatSize, VEC_ALIGN);
  assert(HMat16 != NULL);
  
  /* Declare Serial HEF Mat */
  int HEF_MatSize = (RefLength_1 * ReadLength_1) * sizeof(HEF_Mat);	       
  HEF_Mat* HMatSerial = (HEF_Mat*) _mm_malloc(HEF_MatSize, VEC_ALIGN);
  assert(HMatSerial != NULL);


  /* Initialize Query Sequence */
  Vec_t* VQuery_B = InitQuerySegments_B(*ReadSeq_1);					// NOTE: VQuery should be iterpreted as unsigned
  Vec_t* VQuery_HW = InitQuerySegments_HW(*ReadSeq_1);


  /* Initialize PT_BParam so pthread parameters can be passed
   * NOTE: If overflow occurs InitHMat_B_MT will call InitHMat_HW_MT. After which InitHMat_HW_MT exits.
   */
  int ParamSize = NUM_THREADS * sizeof(PT_HMat);
  PT_HMat* PT_BParam = (PT_HMat*) _mm_malloc(ParamSize, VEC_ALIGN);
  PT_HMat* PT_HWParam = (PT_HMat*) _mm_malloc(ParamSize, VEC_ALIGN);

  for (i = 0; i < NUM_THREADS; i++) {
    (PT_HWParam + i)->ReadSeq = ReadSeq_1;
    (PT_HWParam + i)->RefSeq  = RefSeq_1;
    (PT_HWParam + i)->QueryData = VQuery_B;
    (PT_HWParam + i)->HMat = (Vec_t*) HMat16;
    (PT_HWParam + i)->HWParam = NULL ;				// HW function will exit if overflow occurs
  }

  for (i = 0; i < NUM_THREADS; i++) {
    (PT_BParam + i)->ReadSeq = ReadSeq_1;
    (PT_BParam + i)->RefSeq  = RefSeq_1;
    (PT_BParam + i)->QueryData = VQuery_B;
    (PT_BParam + i)->HMat = (Vec_t*) HMat8;
    (PT_BParam + i)->HWParam = (void*) (PT_HWParam +  i);
  }
  
  /* Error Return Variables for SSW and SW */
  int Err = 0;
  int PTRet = 0;
  int OFCount = 0;

  /* SW parallized with PThreads  */
  Option = 0;
  for (i = 0; i < ITERS; i++) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (j = 0; j < NUM_THREADS; j++) {
      (PT_BParam + j)->TID = j;					// HW ID is the same, assigned in InitHMat_B_MT
      assert((PT_BParam + j)->ReadSeq != NULL);
      assert((PT_BParam + j)->RefSeq != NULL);
      PTRet = pthread_create(&id[j], NULL, InitHMat_B_MT_128, (void*) (PT_BParam + j));
      if (PTRet) {printf("Error creating thread:%d\n",j); exit(19);}
    }
    for (j = 0; j < NUM_THREADS; j++) {
      if (pthread_join(id[j], NULL)) {
	printf("\n Error on join:%d \n", j); exit(19);
      }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[Option][i] = diff(time1,time2);
  }

  /* SW parallized w/ OMP   */
  Option++;
  omp_set_num_threads(NUM_THREADS);
  for (i = 0; i < ITERS; i++) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    #pragma omp parallel for private(Err)
    for (j = 0; j < NUMSEQ; j++) {
      Err  = 0;						  // Err = 1 if max value returned InitHMat >=255 - BIAS, e.g. need more precision
      Err = InitHMat_B_128(VQuery_B, HMat8, *ReadSeq_1, *RefSeq_1);

      if (Err) {
	printf("Overflow count: %d\n", ++OFCount);
	Err = InitHMat_HW_128(VQuery_HW, HMat16, *ReadSeq_1, *RefSeq_1);
      }
    }
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
      time_stamp[Option][i] = diff(time1,time2);
  }


  Option++;
  for (i = 0; i < ITERS; i++) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    Err  = 0;						  // Err = 1 if max value returned InitHMat >=255 - BIAS, e.g. need more precision
    Err = InitHMat_B_128(VQuery_B, HMat8, *ReadSeq_1, *RefSeq_1);
    if (Err) {
      printf("Overflow count: %d\n", ++OFCount);
      Err = InitHMat_HW_128(VQuery_HW, HMat16, *ReadSeq_1, *RefSeq_1);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[Option][i] = diff(time1,time2);
  }

  Option++;
  for (i = 0; i < ITERS; i++) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (j = 0; j < NUMSEQ; j++) {
      InitHMat_Serial(HMatSerial, *ReadSeq_1, *RefSeq_1);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[Option][i] = diff(time1,time2);
  }

  printf("Iter, NumElem, NumSeq, Striped_128_PThreads, Striped_128_OMP, Serial");
  for (i = 0; i < ITERS; i++) {
    printf("\n%d, %ld, %d, ",i, NUMEL, NUMSEQ);
    for (j = 0; j < OPTIONS; j++) {
      if (j != 0) printf(", ");
      printf("%ld", (long int)((double)(CPG)*(double)
		 (GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));
    }
  }
  printf("\n");

  /* Free all memory */
  _mm_free(VQuery_B);
  _mm_free(VQuery_HW);
  _mm_free(HMatSerial);
  _mm_free(HMat8);
  _mm_free(HMat16);
  _mm_free(PT_BParam);
  _mm_free(PT_HWParam);
  _mm_free(RefSeq_1->Seq);
  _mm_free(ReadSeq_1->Seq);
  _mm_free(RefSeq_1);
  _mm_free(ReadSeq_1);
  return 0;
}
