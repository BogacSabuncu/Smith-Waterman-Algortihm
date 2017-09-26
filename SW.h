#ifndef SW_H
#define SW_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>

/* Parallel Processing Declertaions */

/* Number of bytes in SSE vector(128 and 256 bit) */
#define VEC_SIZE 0	 // VEC_SIZE = 1 (256 bit), VEC_SIZE = 0 (128 bit)	
#define VBYTES0 16	 // Number of bytes in 128 bit vector, VEC_SIZE = 0
#define VBYTES1 32	 // Number of bytes in 256 bit vector, VEC_SIZE = 1

/* Set byte and vector size to 128 bit SSE vector */
#if VEC_SIZE == 0
#define VBYTES VBYTES0
typedef __m128i Vec_t;

/* Set byte and vector size to 256 bit SSE vector */
#elif VEC_SIZE == 1
#define VBYTES VBYTES1
typedef __m256i Vec_t;

#endif 

//#define VEC_ALIGN	VBYTES
#define VEC_ALIGN	64				// Align to 64 bytes because of cache block size
#define VSIZE_SEQ	(VBYTES/sizeof(uint8_t))
#define VSIZE_SEQHW	(VBYTES/sizeof(uint16_t))

/* SW Constants*/
#define MATCH_HIT 2
#define MATCH_MISS -2
#define BIAS 2	       // Bias = abs(lowest square in matching matrix), see InitQuerySequence for use
#define GAP_EXT 1      // Gap extension penalty
#define GAP_INIT 3     // Gap initialization penalty, or gap start penalty

/* Constants for types of alignments */
#define SP 0	        // If SP == 1, query segments for all amino acids will be generated (NUM_AA), else query for all nucleotides will be generated(NUM_NT)
#define NUM_AA 20	// Number of amino acids
#define NUM_NT 4	// Number of nucleotides      

/* Scoring matrix struct for serial implementation */
typedef struct {
  uint32_t E;    
  uint32_t F;    
  uint32_t H; 
  uint32_t pad;
} HEF_Mat;

/* SeqData stores data relavent to the respective sequence */
typedef struct {
  uint8_t* Seq;
  uint32_t SeqLength;
  uint32_t PadLength; 
} SeqData;

/* PT_HMat refers to the struct used to pass parameters to PThread versions of InitHMat_B() and InitHMat_HW */
typedef struct {
  SeqData* ReadSeq;
  SeqData* RefSeq;
  Vec_t* QueryData;
  void* HWParam;
  int TID;
  Vec_t* HMat;
  uint8_t Pad[20];		// Pad to 64 byte boundary
} PT_HMat;


/* Vec_t* InitQuerySegments_B(SeqData ReadSeq)
 * Intialiize Query Segments with the elements of H Matrix being 1 byte
 * VSIZE_SEQ represents the number of elements processed per vector for this case
 */
Vec_t* InitQuerySegments_B(SeqData ReadSeq);

/* Vec_t* InitQuerySegments_HW(SeqData ReadSeq)
 * Intialiize Query Segments with the size of the H Matrix = 2 BYTES
 * VSIZE_SEQHW represents the number of elements processed per vector for this case
 */
Vec_t* InitQuerySegments_HW(SeqData ReadSeq);

/* int InitHMat_B(Vec_t* VQuery, uint8_t* HMat, SeqData ReadSeq, SeqData RefSeq)
 * Initializes the scoring matrix using uint8_t integers
 */
int InitHMat_B_128(Vec_t* VQuery, uint8_t* HMat, SeqData ReadSeq, SeqData RefSeq);

/* int InitHMat_HW(Vec_t* VQuery, uint16_t* HMat, SeqData ReadSeq, SeqData RefSeq)
 * Initializes the scoring matrix using uint16_t integers
 * Used when 8 bit precision is unable to provide accurate results
 * When 8 bits are inadequate 1 is returned
 */
int InitHMat_HW_128(Vec_t* VQuery, uint16_t* HMat, SeqData ReadSeq, SeqData RefSeq);

/* void* InitHMat_B_MT(void* Parameters)
 * Pthread implementation of SW
 * Only change is OF is returned as a pointer to OF
 *
 */
void* InitHMat_B_MT_128(void* Parameters);

/* void* InitHMat_HW_MT(void* Parameters)
 * Pthread implementation of SW
 * Only change is OF is returned as a pointer to OF
 *
 */
void* InitHMat_HW_MT_128(void* Parameters);

/* int InitHMat_B(Vec_t* VQuery, uint8_t* HMat, SeqData ReadSeq, SeqData RefSeq)
 * Initializes the scoring matrix using uint8_t integers
 * 256 bit SSE vectors used
 */
int InitHMat_B_256(Vec_t* VQuery, uint8_t* HMat, SeqData ReadSeq, SeqData RefSeq);

/* int InitHMat_HW(Vec_t* VQuery, uint16_t* HMat, SeqData ReadSeq, SeqData RefSeq)
 * Initializes the scoring matrix using uint16_t integers
 * Used when 8 bit precision is unable to provide accurate results
 * When 8 bits are inadequate 1 is returned
 * 256 bit vectors are used
 */
int InitHMat_HW_256(Vec_t* VQuery, uint16_t* HMat, SeqData ReadSeq, SeqData RefSeq);

/* void* InitHMat_B_MT(void* Parameters)
 * Pthread implementation of SW
 * 256 bit SSE vectors used
 *
 */
void* InitHMat_B_MT_256(void* Parameters);

/* void* InitHMat_HW_MT_256(void* Parameters)
 * Pthread implementation of SW
 * 256 SSE Vectors used
 *
 */
void* InitHMat_HW_MT_256(void* Parameters);


/* SeqData* InitArrayExt(const char* Filename, int SeqSize)
 * Inputs text file of nucleotides or protein into integer values correseponding to respective base or protein.
 * Function is parameterized in SW.h for simplicity
 * Data is padded to the an integer multiple of VSIZE_SEQ(Size of vector registers) to ensure data fits is alligned with vectors
 *
 */
SeqData* InitArrayExt (const char* Filename, int SeqSize);

/* uint8_t* InitHMat_SerialB(uint8_t* HMat, SeqData ReadSeq, SeqData RefSeq)
 * Serial implementation of SW. 
 * Same scoring scheme used as vectorised versions
 */
void InitHMat_Serial(HEF_Mat* HMat, SeqData ReadSeq, SeqData RefSeq);

void InitWMat_Vec(SeqData ReadSeq, SeqData RefSeq, uint8_t* WMat);

inline int Max_2(int A, int B);

inline int Max_3(int A, int B, int C);

inline void StoreHEF(HEF_Mat* Loc, int EMax, int HMax, int FMax);

#endif

