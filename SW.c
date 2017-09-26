/*
 * gcc -o EX Example.c SW.c -msse2 -lrt -msse4.1 -mavx -omp -pthreads
 * "immintrin.h" 
 * 
 */
#include "SW.h"
#include <assert.h>

/* Serial SW Functions */
void InitHMat_Serial(HEF_Mat* HMat, SeqData ReadSeq, SeqData RefSeq) {
  int i,j,k,ii,jj;
  int ReadLength = ReadSeq.SeqLength;
  int ReadPadLength = ReadSeq.PadLength;
  int RefLength = RefSeq.SeqLength;
  int RefPadLength = RefSeq.PadLength;
  int GapExt = GAP_EXT;
  int GapInit = GAP_INIT;
  
  /* Score calculation values */
  int EMax = 0;
  int HMax = 0;
  int FMax = 0;
  int Mismatch = 0;

    /* See SW.h for score calcuation */
  for (i = 1; i < ReadLength-ReadPadLength; i++) {
    for (j = 1; j < RefLength-RefPadLength; j++) {
      Mismatch = (*(ReadSeq.Seq + i) == *(RefSeq.Seq + j)) ? MATCH_HIT : MATCH_MISS;
      EMax = Max_2((HMat + i*RefLength + (j-1))->E - GapExt, (HMat + i*RefLength + (j-1))->H - GapInit);
      FMax = Max_2((HMat + (i-1)*RefLength + j)->F - GapExt, (HMat + (i-1)*RefLength + j)->H - GapInit);
      HMax = Max_3(EMax, FMax, (HMat + (i-1)*RefLength + (j-1))->H - Mismatch);
      StoreHEF((HMat + i*RefLength + j), EMax, FMax, HMax);
    }
  }
}

inline int Max_2(int A, int B) {
  return (A >= B) ? A : B;
}

inline int Max_3(int A, int B, int C) {
  int Temp = 0;
  Temp = (A >= B) ? A : B;
  Temp = (Temp >= C) ? Temp : C;
  Temp = (Temp >= 0) ? Temp : 0;
  return Temp;
}

inline void StoreHEF(HEF_Mat* Loc, int EMax, int HMax, int FMax) {
  Loc->E = EMax;
  Loc->H = HMax;
  Loc->F = FMax;
  return;
}

/* SW Array Initialization */
SeqData* InitArrayExt (const char* Filename, int SeqSize) {
  int i = 0;
  char InChar;
  int  InCount = 0;          // Number of characters read
  int  PadCount = 0;
  int  PadLength = VSIZE_SEQ - SeqSize % VSIZE_SEQ;
  uint8_t Temp = 4;
  
  /* Initialize Output Sequence */
  //SeqData* OutSeq =  (SeqData*) _mm_malloc(sizeof(SeqData), VSIZE_SEQ);
  SeqData* OutSeq =  (SeqData*) _mm_malloc(sizeof(SeqData), VEC_ALIGN);
  OutSeq->PadLength = 0;
  OutSeq->SeqLength = 0;
  OutSeq->Seq = (uint8_t*) _mm_malloc((SeqSize + PadLength) * sizeof(uint8_t), VEC_ALIGN);
  
  /* Open Data File */
  FILE* fp;
  fp = fopen(Filename, "r");
  if (fp == NULL) {
    printf("Error opening file \n");
    return OutSeq;
  }
  InChar = getc(fp);
  while (i < (SeqSize) || InChar != -1) {
    switch (InChar) {
     case 'A': {Temp = 0; break;}	// Change
     case 'R': {Temp = 2; break;}
     case 'N': {Temp = 3; break;}
     case 'D': {Temp = 4; break;}
     case 'C': {Temp = 1; break;}	// Change
     case 'Q': {Temp = 6; break;}
     case 'E': {Temp = 7; break;}
     case 'G': {Temp = 2; break;}	// Change
     case 'H': {Temp = 9; break;}
     case 'I': {Temp = 10; break;}
     case 'L': {Temp = 11; break;}
     case 'K': {Temp = 12; break;}
     case 'M': {Temp = 13; break;}
     case 'F': {Temp = 14; break;}
     case 'P': {Temp = 15; break;}
     case 'S': {Temp = 16; break;}
     case 'T': {Temp = 3; break;}	// change
     case 'W': {Temp = 18; break;}
     case 'Y': {Temp = 19; break;}
     case 'V': {Temp = 20; break;}
     default:  {Temp = 4; break;}
    }
    if (Temp != 4) {
      *(OutSeq->Seq + i) = (uint8_t) Temp; // A = 1, B = 2 etc.
      i++;
    }
    InChar = getc(fp);
  }
  
  /* Pad OutSeq if not an integer multiple of VSIZE_SEQ(Vector Size) */
  while (i % VSIZE_SEQ != 0) {
    *(OutSeq->Seq+ i) = (int8_t) 4;  // Set to NOP value, e.g. 4 for nucleotide, 20 for protein
    PadCount++;
    i++;
  }
  OutSeq->PadLength = PadLength;
  OutSeq->SeqLength = SeqSize + PadLength;
  return OutSeq;
}


/* SSW Query Initializations */
Vec_t* InitQuerySegments_B(SeqData ReadSeq) {
  int i,j,k,h,p;
  int QueryNum = (SP) ? NUM_AA+1 : NUM_NT+1;			// Number of query sequences being generated, +1 to account for non-reconizable proteins & i|j = 0
  int ReadLength = ReadSeq.SeqLength;				// -1 to account for the additonal pad bit needed in initWMat()

  int PadLength = ReadSeq.PadLength;
  int SegLength = (ReadLength + (VSIZE_SEQ-1)) / VSIZE_SEQ;	// Size of query sequences being generated in accordance with the Farrar striped Smith-Waterman
  
  int QuerySize = QueryNum * SegLength * VSIZE_SEQ * sizeof(uint8_t);			// Number of bytes needed in query(all protein/nucleotides)
  int BQuerySize = SegLength * VSIZE_SEQ * sizeof(uint8_t);				// Size of one query in bytes
  //Vec_t* VQuery = (Vec_t*) _mm_malloc(QuerySize, VEC_ALIGN);
  Vec_t* VQuery = (Vec_t*) malloc(QuerySize);
  
  uint8_t* QPtr = (uint8_t*) VQuery;				 // Unsigned b/c bias is added
  
  int8_t ScoreMat[QueryNum*QueryNum];

  /* Generate scoring matrix for all possible combinations e.g A := A, A != T etc.
   * First row corresponds to the matching of the first protein/nucleotide denoted as 0 in InitNucArray()
   *	 A    T   C   G   NOP
   * A   MH   MM  MH  MM  MM
   * T   MM   MH  MM  MM  MM
   * C   MM   MM  MH  MM  MM
   * G   MM   MM  MM  MH  MM
   * NOP MM   MM  MM  MM  MM
   * 
   * MH = Match Hit
   * MM = Match Miss
   */
  for (i = 0; i < QueryNum; i++) {
    for (j = 0; j < QueryNum; j++) {
      *(ScoreMat + i*QueryNum + j) = (i == j) ? MATCH_HIT : MATCH_MISS;
      if (i == QueryNum-1 || j == QueryNum-1) {*(ScoreMat + i*QueryNum + j) = MATCH_MISS;}
    }
  }


  /* Generate Query Profile 
   * Query profile generates is n (Num. Proteins/Nucleotides + 1(Non-OP characters)) rows with length := ReadLength
   * For each for element(i,j) i denotes the item being compared to the readsequence. Each item's number(given in InitNucArray) corresponds 
   * with the row in the query where those comparisons are located.
   * For example, A is = 0 so all elements in the first row are equal to MATCH_HIT + BIAS if they are equivalent.
   * ex.
   *  A  T  C  G  T  A  A
   * MH MM MM MM MM MH MH
  */

  for (i = 0; i < QueryNum; i++) {
    h = 0;
    for (j = 0; j < SegLength; j++) {
      p = j;
      for (k = 0; k < VSIZE_SEQ; k++) {
	*(QPtr + i*BQuerySize + h) = (p >= ReadLength-PadLength) ? BIAS : *(ScoreMat + i*QueryNum + *(ReadSeq.Seq + p)) + BIAS;
	h++;
	p += SegLength; 
      }
    }
  }
  return VQuery;
}

Vec_t* InitQuerySegments_HW(SeqData ReadSeq) {
  int i,j,k,h,p;
  int QueryNum = (SP) ? NUM_AA+1 : NUM_NT+1;			// Number of query sequences being generated, +1 to account for non-reconizable proteins & i|j = 0
  int ReadLength = ReadSeq.SeqLength;				// -1 to account for the additonal pad bit needed in initWMat()
  int PadLength = ReadSeq.PadLength;
  int SegLength = (ReadLength + (VSIZE_SEQHW - 1));		// Size of query sequences being generated in accordance with the Farrar striped Smith-Waterman
  SegLength /= VSIZE_SEQHW;

  int QuerySize = QueryNum * SegLength * VSIZE_SEQHW * sizeof(uint16_t);	// Number of bytes needed in query(all protein/nucleotides)
  int BQuerySize = SegLength * VSIZE_SEQHW * sizeof(uint16_t);			// Size of one query in bytes
  Vec_t* VQuery = (Vec_t*) _mm_malloc(QuerySize, VEC_ALIGN);
  
  int16_t* QPtr = (int16_t*) VQuery;
  int8_t ScoreMat[QueryNum*QueryNum];

  /* Generate scoring matrix for all possible combinations e.g A := A, A != T etc.
   * First row corresponds to the matching of the first protein/nucleotide denoted as 0 in InitNucArray()
   *	 A    T   C   G   NOP
   * A   MH   MM  MH  MM  MM
   * T   MM   MH  MM  MM  MM
   * C   MM   MM  MH  MM  MM
   * G   MM   MM  MM  MH  MM
   * NOP MM   MM  MM  MM  MM
   * 
   * MH = Match Hit
   * MM = Match Miss
   */
  for (i = 0; i < QueryNum; i++) {
    for (j = 0; j < QueryNum; j++) {
      *(ScoreMat + i*QueryNum + j) = (i == j) ? MATCH_HIT : MATCH_MISS;
      if (i == QueryNum-1 || j == QueryNum-1) {*(ScoreMat + i*QueryNum + j) = MATCH_MISS;}
    }
  }

  /* Generate Query Profile 
   * Query profile generates is n (Num. Proteins/Nucleotides + 1(Non-OP characters)) rows with length := ReadLength
   * For each for element(i,j) i denotes the item being compared to the readsequence. Each item's number(given in InitNucArray) corresponds 
   * with the row in the query where those comparisons are located.
   * For example, A is = 0 so all elements in the first row are equal to MATCH_HIT + BIAS if they are equivalent.
   * ex.
   *  A  T  C  G  T  A  A
   * MH MM MM MM MM MH MH
  */
  for (i = 0; i < QueryNum; i++) {
    h = 0;
    for (j = 0; j < SegLength; j++) {
      p = j;
      for (k = 0; k < VSIZE_SEQHW; k++) {
	*(QPtr + i*BQuerySize + h) = (p >= ReadLength-PadLength) ? BIAS : *(ScoreMat + i*QueryNum + *(ReadSeq.Seq + p)) + BIAS;
	h++;
	p += SegLength; 
      }
    }
  }
  return VQuery;
}

/* SSW vectorised and MT vectorised implementation with 128 bit SSE vectors used */
int InitHMat_B_128(Vec_t* VQuery, uint8_t* HMat, SeqData ReadSeq, SeqData RefSeq) {
  int i,j,k,ii,jj;
  int ReadLength = ReadSeq.SeqLength;
  int ReadPadLength = ReadSeq.PadLength;
  int RefLength = RefSeq.SeqLength;
  int RefPadLength = RefSeq.PadLength;
  int SegLen = (ReadLength + (VSIZE_SEQ-1))/ VSIZE_SEQ; 
  __m128i* HMatStore = (__m128i*) HMat;

  uint8_t HMax = 0;						// Maximum H Value, if HMax >= 255 - BIAS will return 1(need higher precision)

  /* Algorithmic constants */
  uint8_t GInit = GAP_INIT;			// Gap initialization penalty
  uint8_t GExt  = GAP_EXT;			// Gap extension penalty
  uint8_t Bias  = BIAS;				// Bias, see SW.h
  uint16_t VFMask = 0;				// Stores results of max(F(i,j)) == zero for vector operations, (0xFFFF) if all zeros, else non-zero values. Used for Lazy-F evaluation
  uint8_t Overflow = 0xFF - Bias;

  /* Constant Vector Declerations */
  Vec_t VZero = _mm_set1_epi8((char) 0);
  Vec_t VGapI = _mm_set1_epi8((char) GInit);
  Vec_t VGapE = _mm_set1_epi8((char) GExt);
  Vec_t VBias = _mm_set1_epi8((char) Bias);
  Vec_t VOverflow = _mm_set1_epi8((char) Overflow);

  /* Dynamic Vector Declerations */
  int VecSize = SegLen * VBYTES;							// Size of vectors for one iteration in bytes
  Vec_t* VHStore = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for previous iteration
  Vec_t* VHLoad  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for current iteration
  Vec_t* VHSwap  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Temp vector for swapping H vectors
  Vec_t* VE      = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains E values for current iteration

   /* Loop Vector Declerations */
  Vec_t VHTemp;			// VH stores H scores for current iteration
  Vec_t VETemp;			// Holds values from previous E calculation
  Vec_t VFTemp;			// Holds values needed to calculate F for current iteration
  Vec_t VProfData;			// VProf data corresponds to the vector of the query line needed for the current iteration. e.g. first vector of A
  Vec_t VColMax;		// Maximum value for current reference protein/nucleotide
  Vec_t VTemp;		// Temp vector for inner loop calculations
  Vec_t* VProf;		// VProf corresponds to the VQuery row needed for the current iteration, e.g. for A, T, C etc.
  
  /* Integer used to determine if overflow has occured */
  uint32_t OF = 0;
  
  
  for (i = 0; i < (RefLength-RefPadLength); i++) {
    /* Set F to zero. If errors occur they will be fixed in lazy F loop. */
    VFTemp = VZero;

    /* Set VColMax to zero, max only valid for one reference protein/nucleotide */
    VColMax = VZero;
    
    /* Adjust H Value to be used in next segment */
    VHTemp = *(VHStore + SegLen - 1);	// Store last H vector of prev. col. as 1st inner loop iteration is dependant on this vector.
    VHTemp = _mm_slli_si128(VHTemp,1);	// Shift VH left as the first element in current iteration is dependant on [i-1,j-1]. This always references i = 0(H always 0)
    
    /* Swap H Vectors */
    *VHSwap = *VHLoad;
    *VHLoad = *VHStore;
    *VHStore = *VHSwap;

    /* Store VQuery data corresponding to current iteration */
    VProf = VQuery + *(RefSeq.Seq + i) * SegLen;
    
    for (j = 0; j < SegLen; j++) {

      /* Load Query Data in H */
      VProfData = _mm_load_si128(VProf + j);		// Vector corresponding to the current iteration's profile data
      VHTemp = _mm_adds_epu8(VHTemp, VProfData);	// Load query data into current iteration's VH vector
      VHTemp = _mm_subs_epu8(VHTemp, VBias);		// Remove bias from query b/c

      /* Calculate H for current column. Note no comparison w/ VZero needed as saturation operations are used */
      VETemp  = _mm_load_si128(VE + j);			// E(i,j-1) for all relavent comparisons, e.g. previous col. value
      VHTemp  = _mm_max_epu8(VHTemp, VETemp);	       	// Max(VH,VE)
      VHTemp  = _mm_max_epu8(VHTemp, VFTemp);		// Max(VH,VF)
      VColMax = _mm_max_epu8(VHTemp, VColMax);		// Store column maxes for later traversal

      /* Store VH values for next i iteration */
      _mm_store_si128((VHStore + j), VHTemp);


      /* Calcualte VE & VH for current iteration and store. Values to be used on next iteration */
      VHTemp  = _mm_subs_epu8(VHTemp, VGapI);		// VH[i-1,j-1] - GapInit
      VETemp  = _mm_subs_epu8(VETemp, VGapE);		// VE[i-i,j] - GapExt
      VETemp  = _mm_max_epu8 (VETemp, VHTemp);		// Max(VH,VE)

      /* Store E values for next iteration */
      _mm_store_si128((VE + j), VETemp);		// Store E for next iteration

      /* Calculate VF for next iteration. F is based upon previous rows F so if we subtract GapExt then F can propogate*/
      VFTemp = _mm_subs_epu8(VFTemp, VGapE);		// VF[i,j-1] - GapExt
      VFTemp = _mm_max_epu8(VFTemp, VHTemp);		// Max(VF,VH)

      /* Load next VH */
      VHTemp = _mm_load_si128(VHLoad + j);

    }
    /* Lazy F Evaluation */
    j = 0;

    /* Load first H vector for current column */
    VHTemp = _mm_load_si128(VHStore + j);	
    
    /* Last F vector is used to calculate first H vector. The first element is depedent with F[i-1,j]. Since i-1 for 1st element always 0, F is always zero  */
    VFTemp = _mm_slli_si128(VFTemp,1);
    
    /* VH - GapInit */
    VHTemp = _mm_subs_epu8(VHTemp, VGapI);		// VH[j] - GapInit
    VTemp = _mm_subs_epu8(VFTemp, VHTemp);		// Subtract (VF-VH), VH larger VHTemp = 0 b/c sub_saturation 
    VTemp = _mm_cmpeq_epi8(VTemp, VZero);		// If all zeros evaluation was correct, and F values will propogate
    VFMask = _mm_movemask_epi8(VTemp);

    /* Evaluate F is VF is non-zero for any value and store results in VHStore */
    while (VFMask != 0xFFFF) {
      VHTemp = _mm_max_epu8(VFTemp,VHTemp);		// New VH = Max(VH,VF)
      VColMax = _mm_max_epu8(VHTemp, VColMax);		// Store new column max
      _mm_store_si128((VHStore + j), VHTemp);		// Store newly calculated values
      VFTemp = _mm_subs_epu8(VFTemp, VGapE);		// Calculate VF for next vector, VF will propogate vertically

      /* If segment is processed shift VF values to next segment */
      if (++j >= SegLen) {
	VFTemp = _mm_slli_si128(VFTemp,1);
	j = 0;
      }
      
      /* Test VF again */
      VHTemp = _mm_load_si128(VHStore + j);
      VTemp = _mm_subs_epu8(VHTemp, VGapI);
      VTemp = _mm_subs_epu8(VFTemp,VTemp);
      VTemp = _mm_cmpeq_epi8(VTemp, VZero); 
      VFMask = _mm_movemask_epi8(VTemp);
    }

    /* Check if maximum precision has been reached */
    VTemp = VOverflow;					// VOverflow = (255-Bias)
    VTemp = _mm_cmpeq_epi8(VTemp, VColMax);
    OF = _mm_movemask_epi8(VTemp);

    /* Exit if overflow detected */
    if (OF >= 1) {
      _mm_free((void*) VHStore);
      _mm_free((void*) VHLoad);
      _mm_free((void*) VHSwap);
      _mm_free((void*) VE);
      return 1;
    }	

    /* Store HMat with non-temporal hit */
    for (ii = 0; ii < SegLen; ii++) 
      _mm_stream_si128(HMatStore + ii , *(VHStore + ii));
  }   
  _mm_free((void*) VHStore);
  _mm_free((void*) VHLoad);
  _mm_free((void*) VHSwap);
  _mm_free((void*) VE);
  return 0;
}



int InitHMat_HW_128(Vec_t* VQuery, uint16_t* HMat, SeqData ReadSeq, SeqData RefSeq) {
  int i,j,k,ii,jj;
  int ReadLength = ReadSeq.SeqLength;
  int ReadPadLength = ReadSeq.PadLength;
  int RefLength = RefSeq.SeqLength;
  int RefPadLength = RefSeq.PadLength;
  __m128i* HMatStore = (__m128i*) HMat;

  /* DONT CHANGE
   * SegLen is calculated weird because if calculated in one statement SegLen would return as 1.
   */
  int SegLen = (ReadLength + (VSIZE_SEQHW-1));
  SegLen /= VSIZE_SEQHW; 

  uint8_t HMax = 0;						// Maximum H Value, if HMax >= 255 - BIAS will return 1(need higher precision)

  /* Algorithmic constants */
  uint16_t GInit = GAP_INIT;	// Gap initialization penalty
  uint16_t GExt  = GAP_EXT;	// Gap extension penalty
  uint16_t Bias  = BIAS;	// Bias, see SW.h
  uint16_t Overflow = 0xFFFF - Bias;
  uint16_t VFMask = 0;		// Stores results of max(F(i,j)) == zero for vector operations, (0xFFFF) if all zeros, else non-zero values. Used for Lazy-F evaluation

  /* Constant Vector Declerations */
  Vec_t VZero = _mm_set1_epi16(0);
  Vec_t VGapI = _mm_set1_epi16(GInit);
  Vec_t VGapE = _mm_set1_epi16(GExt);
  Vec_t VBias = _mm_set1_epi16(Bias);
  Vec_t VOverflow = _mm_set1_epi16(Overflow);

  /* Dynamic Vector Declerations */
  int VecSize = SegLen * VBYTES;							// Size of vectors for one iteration in bytes
  Vec_t* VHStore = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for previous iteration
  Vec_t* VHLoad  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for current iteration
  Vec_t* VHSwap  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Temp vector for swapping H vectors
  Vec_t* VE      = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains E values for current iteration

   /* Loop Vector Declerations */
  Vec_t VHTemp;			// VH stores H scores for current iteration
  Vec_t VETemp;			// Holds values from previous E calculation
  Vec_t VFTemp;			// Holds values needed to calculate F for current iteration
  Vec_t VProfData;			// VProf data corresponds to the vector of the query line needed for the current iteration. e.g. first vector of A
  Vec_t VColMax;		// Maximum value for current reference protein/nucleotide
  Vec_t VTemp;		// Temp vector for inner loop calculations
  Vec_t* VProf;		// VProf corresponds to the VQuery row needed for the current iteration, e.g. for A, T, C etc.
  
  /* Integer used to determine if overflow has occured */
  uint32_t OF = 0;
  
  
  for (i = 0; i < (RefLength-RefPadLength); i++) {
    /* Set F to zero. If errors occur they will be fixed in lazy F loop. VColMax only valid for one reference protein or nucleotide */
    VFTemp = VZero; VETemp = VZero; VColMax = VZero;

    /* Adjust H Value to be used in next segment */
    VHTemp = *(VHStore + SegLen - 1);	// Store last H vector of prev. col. as 1st inner loop iteration is dependant on this vector.
    VHTemp = _mm_slli_epi16(VHTemp,1);	// Shift VH left as the first element in current iteration is dependant on [i-1,j-1]. This always references i = 0(H always 0)
    
    /* Swap H Vectors */
    *VHSwap = *VHLoad;
    *VHLoad = *VHStore;
    *VHStore = *VHSwap;

    /* Store VQuery data corresponding to current iteration */
    VProf = VQuery + *(RefSeq.Seq + i) * SegLen;
    
    for (j = 0; j < SegLen; j++) {

      /* Load Query Data in H */
      VProfData = _mm_load_si128(VProf + j);		// Vector corresponding to the current iteration's profile data
      VHTemp = _mm_adds_epu16(VHTemp, VProfData);	// Load query data into current iteration's VH vector
      VHTemp = _mm_subs_epu16(VHTemp, VBias);		// Remove bias from query b/c

      /* Calculate H for current column. Note no comparison w/ VZero needed as saturation operations are used */
      VETemp  = _mm_load_si128(VE + j);			// E(i,j-1) for all relavent comparisons, e.g. previous col. value
      VHTemp  = _mm_max_epu16(VHTemp, VETemp);	       	// Max(VH,VE)
      VHTemp  = _mm_max_epu16(VHTemp, VFTemp);		// Max(VH,VF)
      VColMax = _mm_max_epu16(VHTemp, VColMax);		// Store column maxes for later traversal

      /* Store VH values for next i iteration */
      _mm_store_si128((VHStore + j), VHTemp);

      /* Calcualte VE & VH for current iteration and store. Values to be used on next iteration */
      VHTemp  = _mm_subs_epu16(VHTemp, VGapI);		// VH[i-1,j-1] - GapInit
      VETemp  = _mm_subs_epu16(VETemp, VGapE);		// VE[i-i,j] - GapExt
      VETemp  = _mm_max_epu16 (VETemp, VHTemp);		// Max(VH,VE)
      _mm_store_si128((VE + j), VETemp);		// Store E for next iteration

      /* Calculate VF for next iteration. F is based upon previous rows F so if we subtract GapExt then F can propogate*/
      VFTemp = _mm_subs_epu16(VFTemp, VGapE);		// VF[i,j-1] - GapExt
      VFTemp = _mm_max_epu16(VFTemp, VHTemp);		// Max(VF,VH)

      /* Load next VH */
      VHTemp = _mm_load_si128(VHLoad + j);

    }
    /* Lazy F Evaluation */
    j = 0;

    /* Load first H vector for current column */
    VHTemp = _mm_load_si128(VHStore + j);	
    
    /* Last F vector is used to calculate first H vector. The first element is depedent with F[i-1,j]. Since i-1 for 1st element always 0, F is always zero  */
    VFTemp = _mm_slli_si128(VFTemp,1);
    
    /* VH - GapInit */
    VHTemp = _mm_subs_epu16(VHTemp, VGapI);		// VH[j] - GapInit
    VTemp = _mm_subs_epu16(VFTemp, VHTemp);		// Subtract (VF-VH), VH larger VHTemp = 0 b/c sub_saturation 
    VTemp = _mm_cmpeq_epi16(VTemp, VZero);		// If all zeros evaluation was correct, and F values will propogate
    VFMask = _mm_movemask_epi8(VTemp);

    /* Evaluate F is VF is non-zero for any value and store results in VHStore */
    while (VFMask != 0xFFFF) {
      VHTemp = _mm_max_epu16(VFTemp,VHTemp);		// New VH = Max(VH,VF)
      VColMax = _mm_max_epu16(VHTemp, VColMax);		// Store new column max
      _mm_store_si128((VHStore + j), VHTemp);		// Store newly calculated values
      VFTemp = _mm_subs_epu16(VFTemp, VGapE);		// Calculate VF for next vector, VF will propogate vertically

      /* If segment is processed shift VF values to next segment */
      if (++j >= SegLen) {
	VFTemp = _mm_slli_epi16(VFTemp,1);
	j = 0;
      }
      
      /* If segment is not finished revaluate for next vector */
      
      /* Test VF again */
      VHTemp = _mm_load_si128(VHStore + j);
      VTemp = _mm_subs_epu16(VHTemp, VGapI);
      VTemp = _mm_subs_epu16(VFTemp,VTemp);
      VTemp = _mm_cmpeq_epi16(VTemp, VZero); 
      VFMask = _mm_movemask_epi8(VTemp);
    }

    /* Check if maximum precision has been reached */
    VTemp = VOverflow;					// VOverflow = (2^16 - Bias)
    VTemp = _mm_cmpeq_epi16(VTemp, VColMax);
    OF = _mm_movemask_epi8(VTemp);
    assert(OF == 0);

    /* Exit if overflow is detected */
    if (OF >= 1) {
      _mm_free((void*) VHStore);
      _mm_free((void*) VHLoad);
      _mm_free((void*) VHSwap);
      _mm_free((void*) VE);
      return 1;
    }				
    
    /* Store HMat with non-temporal hit */
    for (ii = 0; ii < SegLen; ii++) 
      _mm_stream_si128((HMatStore + ii), *(VHStore + ii));
  }
  _mm_free((void*) VHStore);
  _mm_free((void*) VHLoad);
  _mm_free((void*) VHSwap);
  _mm_free((void*) VE);
  return 0;
}


void* InitHMat_B_MT_128(void* H_PT_Struct) {
  /* Recast Struct Pointer */
  PT_HMat* DataIn = (PT_HMat*) H_PT_Struct;
  
  /* Check Data Integrity */
  assert(DataIn->ReadSeq != NULL);
  assert(DataIn->RefSeq != NULL);
  assert(DataIn->QueryData != NULL);
  assert(DataIn->HMat != NULL);
  assert(DataIn->HWParam != NULL);

  /* Unpack Input Data */
  SeqData ReadSeq = *(DataIn->ReadSeq);
  SeqData RefSeq  = *(DataIn->RefSeq);
  Vec_t* VQuery = DataIn->QueryData;
  int TID = DataIn->TID;
  PT_HMat* HWData = (PT_HMat*) DataIn->HWParam;			// Parameters for HW function. Used if overflow occurs
  __m128i* HMatStore = DataIn->HMat;
  HWData->TID = TID;						// Set thread of ID of 16 bit data parameter struct

  /* Loop/Function Overhead Variables */
  int i,j,k,ii,jj;
  int ReadLength = ReadSeq.SeqLength;
  int ReadPadLength = ReadSeq.PadLength;
  int RefLength = RefSeq.SeqLength;
  int RefPadLength = RefSeq.PadLength;
  int SegLen = (ReadLength + (VSIZE_SEQ-1))/ VSIZE_SEQ; 

  uint8_t HMax = 0;						// Maximum H Value, if HMax >= 255 - BIAS will return 1(need higher precision)

  /* Algorithmic constants */
  uint8_t GInit = GAP_INIT;	// Gap initialization penalty
  uint8_t GExt  = GAP_EXT;	// Gap extension penalty
  uint8_t Bias  = BIAS;		// Bias, see SW.h
  uint16_t VFMask = 0;		// Stores results of max(F(i,j)) == zero for vector operations, (0xFFFF) if all zeros, else non-zero values. Used for Lazy-F evaluation

  /* Constant Vector Declerations */
  Vec_t VZero = _mm_set1_epi8((char) 0);
  Vec_t VGapI = _mm_set1_epi8((char) GInit);
  Vec_t VGapE = _mm_set1_epi8((char) GExt);
  Vec_t VBias = _mm_set1_epi8((char) Bias);
  Vec_t VOverflow = _mm_set1_epi8((char) 255-Bias);
  

  /* Dynamic Vector Declerations */
  int VecSize = SegLen * VBYTES;							// Size of vectors for one iteration in bytes
  Vec_t* VHStore = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for previous iteration
  Vec_t* VHLoad  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for current iteration
  Vec_t* VHSwap  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Temp vector for swapping H vectors
  Vec_t* VE      = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains E values for current iteration
  
   /* Loop Vector Declerations */
  Vec_t VHTemp;			// VH stores H scores for current iteration
  Vec_t VETemp;			// Holds values from previous E calculation
  Vec_t VFTemp;			// Holds values needed to calculate F for current iteration
  Vec_t VProfData;			// VProf data corresponds to the vector of the query line needed for the current iteration. e.g. first vector of A
  Vec_t VColMax;		// Maximum value for current reference protein/nucleotide
  Vec_t VTemp;		// Temp vector for inner loop calculations
  Vec_t* VProf;		// VProf corresponds to the VQuery row needed for the current iteration, e.g. for A, T, C etc.
  
  /* Integer used to determine if overflow has occured */
  uint32_t OF = 0;
  
  
  for (i = 0; i < (RefLength-RefPadLength); i++) {
    /* Set F to zero. If errors occur they will be fixed in lazy F loop. */
    VFTemp = VZero;

    /* Set VColMax to zero, max only valid for one reference protein/nucleotide */
    VColMax = VZero;
    
    /* Adjust H Value to be used in next segment */
    VHTemp = *(VHStore + SegLen - 1);	// Store last H vector of prev. col. as 1st inner loop iteration is dependant on this vector.
    VHTemp = _mm_slli_si128(VHTemp,1);	// Shift VH left as the first element in current iteration is dependant on [i-1,j-1]. This always references i = 0(H always 0)
    
    /* Swap H Vectors */
    *VHSwap = *VHLoad;
    *VHLoad = *VHStore;
    *VHStore = *VHSwap;

    /* Store VQuery data corresponding to current iteration */
    VProf = VQuery + *(RefSeq.Seq + i) * SegLen;
    
    for (j = 0; j < SegLen; j++) {

      /* Load Query Data in H */
      VProfData = _mm_load_si128(VProf + j);		// Vector corresponding to the current iteration's profile data
      VHTemp = _mm_adds_epu8(VHTemp, VProfData);	// Load query data into current iteration's VH vector
      VHTemp = _mm_subs_epu8(VHTemp, VBias);		// Remove bias from query b/c

      /* Calculate H for current column. Note no comparison w/ VZero needed as saturation operations are used */
      VETemp  = _mm_load_si128(VE + j);			// E(i,j-1) for all relavent comparisons, e.g. previous col. value
      VHTemp  = _mm_max_epu8(VHTemp, VETemp);	       	// Max(VH,VE)
      VHTemp  = _mm_max_epu8(VHTemp, VFTemp);		// Max(VH,VF)
      VColMax = _mm_max_epu8(VHTemp, VColMax);		// Store column maxes for later traversal

      /* Store VH values for next i iteration */
      _mm_store_si128((VHStore + j), VHTemp);


      /* Calcualte VE & VH for current iteration and store. Values to be used on next iteration */
      VHTemp  = _mm_subs_epu8(VHTemp, VGapI);		// VH[i-1,j-1] - GapInit
      VETemp  = _mm_subs_epu8(VETemp, VGapE);		// VE[i-i,j] - GapExt
      VETemp  = _mm_max_epu8 (VETemp, VHTemp);		// Max(VH,VE)
      _mm_store_si128((VE + j), VETemp);		// Store E for next iteration

      /* Calculate VF for next iteration. F is based upon previous rows F so if we subtract GapExt then F can propogate*/
      VFTemp = _mm_subs_epu8(VFTemp, VGapE);		// VF[i,j-1] - GapExt
      VFTemp = _mm_max_epu8(VFTemp, VHTemp);		// Max(VF,VH)

      /* Load next VH */
      VHTemp = _mm_load_si128(VHLoad + j);

    }
    /* Lazy F Evaluation */
    j = 0;

    /* Load first H vector for current column */
    VHTemp = _mm_load_si128(VHStore + j);	
    
    /* Last F vector is used to calculate first H vector. The first element is depedent with F[i-1,j]. Since i-1 for 1st element always 0, F is always zero  */
    VFTemp = _mm_slli_si128(VFTemp,1);
    
    /* VH - GapInit */
    VHTemp = _mm_subs_epu8(VHTemp, VGapI);		// VH[j] - GapInit
    VTemp = _mm_subs_epu8(VFTemp, VHTemp);		// Subtract (VF-VH), VH larger VHTemp = 0 b/c sub_saturation 
    VTemp = _mm_cmpeq_epi8(VTemp, VZero);		// If all zeros evaluation was correct, and F values will propogate
    VFMask = _mm_movemask_epi8(VTemp);

    /* Evaluate F is VF is non-zero for any value and store results in VHStore */
    while (VFMask != 0xFFFF) {
      VHTemp = _mm_max_epu8(VFTemp,VHTemp);		// New VH = Max(VH,VF)
      VColMax = _mm_max_epu8(VHTemp, VColMax);		// Store new column max
      _mm_store_si128((VHStore + j), VHTemp);		// Store newly calculated values
      VFTemp = _mm_subs_epu8(VFTemp, VGapE);		// Calculate VF for next vector, VF will propogate vertically

      /* If segment is processed shift VF values to next segment */
      if (++j >= SegLen) {
	VFTemp = _mm_slli_si128(VFTemp,1);
	j = 0;
      }
      
      /* If segment is not finished revaluate for next vector */
      
      /* Test VF again */
      VHTemp = _mm_load_si128(VHStore + j);
      VTemp = _mm_subs_epu8(VHTemp, VGapI);
      VTemp = _mm_subs_epu8(VFTemp,VTemp);
      VTemp = _mm_cmpeq_epi8(VTemp, VZero); 
      VFMask = _mm_movemask_epi8(VTemp);
    }

    /* Check if maximum precision has been reached */
    VTemp = VOverflow;					// VOverflow = (255-Bias)
    VTemp = _mm_cmpeq_epi8(VTemp, VColMax);
    OF = _mm_movemask_epi8(VTemp);

    /* Do higher precision calculation if overflow detected */
    assert(OF == 0);
    if (OF >= 1) {
      _mm_free((void*) VHStore);
      _mm_free((void*) VHLoad);
      _mm_free((void*) VHSwap);
      _mm_free((void*) VE);
      InitHMat_HW_MT_128((void*) HWData);
    }     

    /* Store HMat with non-temporal hit */
    for (ii = 0; ii < SegLen; ii++) 
      _mm_stream_si128((HMatStore + ii), *(VHStore + ii));
  }
  _mm_free((void*) VHStore);
  _mm_free((void*) VHLoad);
  _mm_free((void*) VHSwap);
  _mm_free((void*) VE);
  pthread_exit(NULL);
}

void* InitHMat_HW_MT_128(void* H_PT_Struct) {
  /* Recast Parameter Input Struct */
  static int OFCount = 0;
  printf("MT Overflow count: %d\n", ++OFCount);
  PT_HMat* DataIn = (PT_HMat*) H_PT_Struct;

  /* Test Data Integrity */
  assert(DataIn->ReadSeq != NULL);
  assert(DataIn->RefSeq != NULL);
  assert(DataIn->QueryData != NULL);
  assert(DataIn->HMat != NULL);
  assert(DataIn->HWParam == NULL);

  /* Unpack Input Data */
  SeqData ReadSeq = *(DataIn->ReadSeq);
  SeqData RefSeq  = *(DataIn->RefSeq);
  Vec_t* VQuery = DataIn->QueryData;
  __m128i* HMatStore = DataIn->HMat;
  int TID = DataIn->TID;

  /* Loop/Function Overhead Variables */
  int i,j,k,ii,jj;
  int ReadLength = ReadSeq.SeqLength;
  int ReadPadLength = ReadSeq.PadLength;
  int RefLength = RefSeq.SeqLength;
  int RefPadLength = RefSeq.PadLength;

  /* DONT CHANGE
   * SegLen is calculated weird because if calculated in one statement SegLen would return as 1.
   */
  int SegLen = (ReadLength + (VSIZE_SEQHW-1));
  SegLen /= VSIZE_SEQHW; 

  uint8_t HMax = 0;						// Maximum H Value, if HMax >= 255 - BIAS will return 1(need higher precision)

  /* Algorithmic constants */
  uint16_t GInit = GAP_INIT;	// Gap initialization penalty
  uint16_t GExt  = GAP_EXT;	// Gap extension penalty
  uint16_t Bias  = BIAS;	// Bias, see SW.h
  uint16_t Overflow = 0xFFFF - Bias;
  uint16_t VFMask = 0;		// Stores results of max(F(i,j)) == zero for vector operations, (0xFFFF) if all zeros, else non-zero values. Used for Lazy-F evaluation

  /* Constant Vector Declerations */
  Vec_t VZero = _mm_set1_epi16(0);
  Vec_t VGapI = _mm_set1_epi16(GInit);
  Vec_t VGapE = _mm_set1_epi16(GExt);
  Vec_t VBias = _mm_set1_epi16(Bias);
  Vec_t VOverflow = _mm_set1_epi16(Overflow);
  

  /* Dynamic Vector Declerations */
  int VecSize = SegLen * VBYTES;							// Size of vectors for one iteration in bytes
  Vec_t* VHStore = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for previous iteration
  Vec_t* VHLoad  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains H values for current iteration
  Vec_t* VHSwap  = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Temp vector for swapping H vectors
  Vec_t* VE      = (Vec_t*) _mm_malloc(VecSize, VSIZE_SEQ);				// Contains E values for current iteration
   
   /* Loop Vector Declerations */
  Vec_t VHTemp;			// VH stores H scores for current iteration
  Vec_t VETemp;			// Holds values from previous E calculation
  Vec_t VFTemp;			// Holds values needed to calculate F for current iteration
  Vec_t VProfData;			// VProf data corresponds to the vector of the query line needed for the current iteration. e.g. first vector of A
  Vec_t VColMax;		// Maximum value for current reference protein/nucleotide
  Vec_t VTemp;		// Temp vector for inner loop calculations
  Vec_t* VProf;		// VProf corresponds to the VQuery row needed for the current iteration, e.g. for A, T, C etc.
  
  /* Integer used to determine if overflow has occured */
  uint32_t OF = 0;

  for (i = 0; i < (RefLength-RefPadLength); i++) {
    /* Set F to zero. If errors occur they will be fixed in lazy F loop. VColMax only valid for one reference protein or nucleotide */
    VFTemp = VZero; VETemp = VZero; VColMax = VZero;

    /* Adjust H Value to be used in next segment */
    VHTemp = *(VHStore + SegLen - 1);	// Store last H vector of prev. col. as 1st inner loop iteration is dependant on this vector.
    VHTemp = _mm_slli_epi16(VHTemp,1);	// Shift VH left as the first element in current iteration is dependant on [i-1,j-1]. This always references i = 0(H always 0)
    
    /* Swap H Vectors */
    *VHSwap = *VHLoad;
    *VHLoad = *VHStore;
    *VHStore = *VHSwap;

    /* Store VQuery data corresponding to current iteration */
    VProf = VQuery + *(RefSeq.Seq + i) * SegLen;
    
    for (j = 0; j < SegLen; j++) {

      /* Load Query Data in H */
      VProfData = _mm_load_si128(VProf + j);		// Vector corresponding to the current iteration's profile data
      VHTemp = _mm_adds_epu16(VHTemp, VProfData);	// Load query data into current iteration's VH vector
      VHTemp = _mm_subs_epu16(VHTemp, VBias);		// Remove bias from query b/c

      /* Calculate H for current column. Note no comparison w/ VZero needed as saturation operations are used */
      VETemp  = _mm_load_si128(VE + j);			// E(i,j-1) for all relavent comparisons, e.g. previous col. value
      VHTemp  = _mm_max_epu16(VHTemp, VETemp);	       	// Max(VH,VE)
      VHTemp  = _mm_max_epu16(VHTemp, VFTemp);		// Max(VH,VF)
      VColMax = _mm_max_epu16(VHTemp, VColMax);		// Store column maxes for later traversal

      /* Store VH values for next i iteration */
      _mm_store_si128((VHStore + j), VHTemp);

      /* Calcualte VE & VH for current iteration and store. Values to be used on next iteration */
      VHTemp  = _mm_subs_epu16(VHTemp, VGapI);		// VH[i-1,j-1] - GapInit
      VETemp  = _mm_subs_epu16(VETemp, VGapE);		// VE[i-i,j] - GapExt
      VETemp  = _mm_max_epu16 (VETemp, VHTemp);		// Max(VH,VE)
      _mm_store_si128((VE + j), VETemp);		// Store E for next iteration

      /* Calculate VF for next iteration. F is based upon previous rows F so if we subtract GapExt then F can propogate*/
      VFTemp = _mm_subs_epu16(VFTemp, VGapE);		// VF[i,j-1] - GapExt
      VFTemp = _mm_max_epu16(VFTemp, VHTemp);		// Max(VF,VH)

      /* Load next VH */
      VHTemp = _mm_load_si128(VHLoad + j);

    }
    /* Lazy F Evaluation */
    j = 0;

    /* Load first H vector for current column */
    VHTemp = _mm_load_si128(VHStore + j);	
    
    /* Last F vector is used to calculate first H vector. The first element is depedent with F[i-1,j]. Since i-1 for 1st element always 0, F is always zero  */
    VFTemp = _mm_slli_si128(VFTemp,1);
    
    /* VH - GapInit */
    VHTemp = _mm_subs_epu16(VHTemp, VGapI);		// VH[j] - GapInit
    VTemp = _mm_subs_epu16(VFTemp, VHTemp);		// Subtract (VF-VH), VH larger VHTemp = 0 b/c sub_saturation 
    VTemp = _mm_cmpeq_epi16(VTemp, VZero);		// If all zeros evaluation was correct, and F values will propogate
    VFMask = _mm_movemask_epi8(VTemp);

    /* Evaluate F is VF is non-zero for any value and store results in VHStore */
    while (VFMask != 0xFFFF) {
      VHTemp = _mm_max_epu16(VFTemp,VHTemp);		// New VH = Max(VH,VF)
      VColMax = _mm_max_epu16(VHTemp, VColMax);		// Store new column max
      _mm_store_si128((VHStore + j), VHTemp);		// Store newly calculated values
      VFTemp = _mm_subs_epu16(VFTemp, VGapE);		// Calculate VF for next vector, VF will propogate vertically

      /* If segment is processed shift VF values to next segment */
      if (++j >= SegLen) {
	VFTemp = _mm_slli_epi16(VFTemp,1);
	j = 0;
      }
      
      /* If segment is not finished revaluate for next vector */
      
      /* Test VF again */
      VHTemp = _mm_load_si128(VHStore + j);
      VTemp = _mm_subs_epu16(VHTemp, VGapI);
      VTemp = _mm_subs_epu16(VFTemp,VTemp);
      VTemp = _mm_cmpeq_epi16(VTemp, VZero); 
      VFMask = _mm_movemask_epi8(VTemp);
    }

    /* Check if maximum precision has been reached */
    VTemp = VOverflow;					// VOverflow = (2^16 - Bias)
    VTemp = _mm_cmpeq_epi16(VTemp, VColMax);
    OF = _mm_movemask_epi8(VTemp);
    
    /* If overflow detected exit */
    assert(OF == 0);
    if (OF >= 1) {printf("\n 16 bit overflow detected \n"); exit(-1);}			       
    
    /* Store HMat with non-temporal hit */
    for (ii = 0; ii < SegLen; ii++) 
      _mm_stream_si128((HMatStore + ii), *(VHStore + ii));
  }
  
  _mm_free((void*) VHStore);
  _mm_free((void*) VHLoad);
  _mm_free((void*) VHSwap);
  _mm_free((void*) VE);
  pthread_exit(NULL);
}

void InitWMat_Vec(SeqData ReadSeq, SeqData RefSeq, uint8_t* WMat){
  uint32_t i,j;
  uint32_t ReadLength = ReadSeq.SeqLength;
  uint32_t RefLength = RefSeq.SeqLength;
  uint32_t ReadPadLength = ReadSeq.PadLength;
  uint32_t RefPadLength = RefSeq.PadLength;
  uint32_t RefVecNum = RefLength / VSIZE_SEQ;

  /* Values to broadcasted for matches and mismatches */
  int8_t MatchHit = MATCH_HIT;
  int8_t MatchMiss = MATCH_MISS;
  
  /* Calculation/Storage Vectors */
  __m128i RefVec0;
  __m128i ReadVec0;
  __m128i ResultVec0;
  __m128i TempVec0;

  /* Vector storing protein matches */
  __m128i EQMaskVec0;
  __m128i EQMaskVec1;
  
  /* Vectors for match and mismatch penalities */
  __m128i MatchHitVector; 
  __m128i MatchMissVector;

  MatchHitVector = _mm_set1_epi8((char) MatchHit);
  MatchMissVector = _mm_set1_epi8((char) MatchMiss);

  for (i = 0; i < RefVecNum; i++) {
    RefVec0 = _mm_loadu_si128((__m128i*) (RefSeq.Seq + i*VSIZE_SEQ + 1)); // +1 to avoid 1st column zero pad
    for (j = 1; j < (ReadLength - ReadPadLength); j++) {
	/* Broadcast single protein*/ 
       ReadVec0 = _mm_set1_epi8((char) *(ReadSeq.Seq + j));

       /* Record protein matches */ 
       EQMaskVec0 = _mm_cmpeq_epi8(ReadVec0, RefVec0);

       /* Set protein matches to MATCH_HIT */
       ResultVec0 = _mm_and_si128(EQMaskVec0, MatchHitVector);

       /* Set protein mismatches to MATCH_MISS */
       TempVec0 = _mm_andnot_si128(EQMaskVec0, MatchMissVector);

       /* Combine protein matches and mismatches in one vector */
       ResultVec0 = _mm_add_epi8(ResultVec0, TempVec0);

       /* Store hits and misses in WMat */
       _mm_storeu_si128 ((__m128i*) (WMat + j*RefLength + i*VSIZE_SEQ + 1), ResultVec0);

    }
  }
}

