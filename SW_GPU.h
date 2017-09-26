#ifndef SW_H
#define SW_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#define MATCH_HIT 3
#define MATCH_MISS -3
#define THREAD_COUNT 16


#define SEQ_GAP 64000 // Value if sequence alignment goes up a row
#define SEQ_HORZ 64001 // Value if sequence alignment move horizontally
/* Data Structure Declerations */
typedef struct {
  int8_t SubScore;
  int8_t GP; // Gap Penalty, may need larger integer.
} ScoreElem, *SE_ptr;

typedef struct {
  int32_t i;
  int32_t j;
  int8_t Match; // 1 if nucleotides corresponding to i and j match
} MaxElem;

typedef struct {
  uint8_t* Seq;
  uint32_t SeqLength;

} SeqData;

/* Number of bytes in SSE vector */
#define VBYTES 32

/* Number of elements in SSE vector */
#define V_MAT_SIZE  VBYTES/sizeof(ScoreMat)
#define V_GP_SIZE VBYTES/sizeof(int8_t)
#define V_SUB_SIZE VBYTES/sizeof(int8_t)


/*
 * InitNucArray initializes an array of integers corresponding to the inputted nucleotide sequence
 * 65 = A
 * 84 = T
 * 71 = G
 * 67 = C
 */
void InitNucArray (const char* RefSeq, SeqData OutSeq);

/* Length refers to the length of RefSeq e.g. number of columns
 * SubMat[i*length + j] = 1 if(RefSeq[i] == ReadSequence[j])
 * SubMat[i*length + j] = 1 if(RefSeq[i] == ReadSequence[j])
 *
 * Gap Penalty scheme is affine scheme
 * k is defined as the length of current gap
 * GP[i*length + j] = GAP_EXT(k - 1) + GAP_START
 */
void InitSubGPMat (SeqData RefSeq, SeqData ReadSeq, ScoreElem* SubMat);

/*
 *
 *
 */


MaxElem* InitScoreMat (ScoreElem* SubMat, SeqData RefSeq, SeqData ReadSeq, int16_t* ScoreMat);
int16_t MaxGPScore (int16_t Score1, int16_t Score2, int16_t Score3);


/* 
 * MaxTBScore => Max Traceback Score: Finds maximum value of three adjacent elements
 * Returns 1 if Elem(i,j-1) is greatest
 * Returns 2 if Elem(i-1,j-1) is greatest
 * Returns 3 if Elem(i-1,j) is greatest
 * If tie between 2 and 1 or 3 or both 2 will win.
 * Tie between 1 and 3 will have 1 win.
 */
int8_t MaxTBScore (int16_t Score1, int16_t Score2, int16_t Score3);

/*
 *
 *
 *
 */
void Traceback (MaxElem* AlignSeq, int16_t* ScoreMat, SeqData RefSeq, SeqData ReadSeq, MaxElem* MaxIndex);

void PrintAlignment (MaxElem* AlignSeq, SeqData RefSeq, SeqData ReadSeq);
#endif
