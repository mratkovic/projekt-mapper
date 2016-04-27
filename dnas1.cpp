#include <bits/stdc++.h>
using namespace std;

// If you submit, comment the following line
#define LOCAL_ONLY

namespace hawker {
uint32_t tot_len = 0;
uint32_t kmer_cnt = 0;

uint32_t pos_tot_len = 0;
uint32_t pos_kmer_cnt = 0;

uint32_t neg_tot_len = 0;
uint32_t neg_kmer_cnt = 0;

const bool LOAD = true;
const bool SAVE = false;
const int BREAK_CNT = INT32_MAX;

const int KMER = 52;
const int MAX_KMER = 144;
const int LO_CNT = 1;
const int HI_CNT = 1;
const int THREADS = 1;
const double KEEP_F = 1.8;
const int KEEP_NUM = 8;
const int MAX_EDIT = 35;

const double TRESH = -1;

const bool ALIGN = false;
const bool FIND_STARTS = false;

#ifndef EDLIB_H
#define EDLIB_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 The MIT License (MIT)

 Copyright (c) 2014 Martin Sosic

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */
// Status codes
#define EDLIB_STATUS_OK 0
#define EDLIB_STATUS_ERROR 1

// Alignment modes
#define EDLIB_MODE_HW  0
#define EDLIB_MODE_NW  1
#define EDLIB_MODE_SHW 2
#define EDLIB_MODE_OV  3

// Cigar formats
#define EDLIB_CIGAR_EXTENDED 0
#define EDLIB_CIGAR_STANDARD 1

/**
 * Calculates Levenshtein distance of query and target
 * using Myers's fast bit-vector algorithm and Ukkonen's algorithm.
 * In Levenshtein distance mismatch and indel have cost of 1, while match has cost of 0.
 * Query and target are represented as arrays of numbers, where each number is
 * index of corresponding letter in alphabet. So for example if alphabet is ['A','C','T','G']
 * and query string is "AACG" and target string is "GATTCGG" then our input query should be
 * [0,0,1,3] and input target should be [3,0,2,2,1,3,3] (and alphabetLength would be 4).
 * @param [in] query  Array of alphabet indices.
 * @param [in] queryLength
 * @param [in] target  Array of alphabet indices.
 * @param [in] targetLength
 * @param [in] alphabetLength
 * @param [in] k  Non-negative number, constraint for Ukkonen.
 *     Only best score <= k will be searched for.
 *     If k is smaller then calculation is faster.
 *     If you are interested in score only if it is <= K, set k to K.
 *     If k is negative then k will be auto-adjusted (increased) until score is found.
 * @param [in] mode  Mode that determines alignment algorithm.
 *     EDLIB_MODE_NW: global (Needleman-Wunsch)
 *     EDLIB_MODE_HW: semi-global. Gaps before and after query are not penalized.
 *     EDLIB_MODE_SHW: semi-global. Gap after query is not penalized.
 *     EDLIB_MODE_OV: semi-global. Gaps before and after query and target are not penalized.
 * @param [in] findStartLocations  If true, start locations are returned.
 *                                 May somewhat slow down the calculation.
 *                                 If findAlignment is true, start locations will also be found.
 * @param [in] findAlignment  If true and if score != -1, reconstruction of alignment will be performed
 *                            and alignment will be returned.
 *                            Notice: Finding aligment will increase execution time
 *                                    and could take large amount of memory.
 * @param [out] bestScore  Best score (smallest edit distance) or -1 if there is no score <= k.
 * @param [out] endLocations  Array of zero-based positions in target where
 *     query ends (position of last character) with the best score.
 *     If gap after query is penalized, gap counts as part of query (NW), otherwise not.
 *     If there is no score <= k, endLocations is set to NULL.
 *     Otherwise, array is returned and it is on you to free it with free().
 * @param [out] startLocations  Array of zero-based positions in target where
 *     query starts, they correspond to endLocations.
 *     If gap before query is penalized, gap counts as part of query (NW), otherwise not.
 *     If there is no score <= k, startLocations is set to NULL.
 *     Otherwise, array is returned and it is on you to free it with free().
 * @param [out] numLocations  Number of positions returned.
 * @param [out] alignment  Alignment is found for first position returned.
 *                         Will contain alignment if findAlignment is true and score != -1.
 *                         Otherwise it will be set NULL.
 *                         Alignment is sequence of numbers: 0, 1, 2, 3.
 *                         0 stands for match.
 *                         1 stands for insertion to target.
 *                         2 stands for insertion to query.
 *                         3 stands for mismatch.
 *                         Alignment aligns query to target from begining of query till end of query.
 *                         Alignment ends at @param positions[0] in target.
 *                         If gaps are not penalized, they are not in alignment.
 *                         Needed memory is allocated and given pointer is set to it.
 *                         Important: Do not forget to free memory allocated for alignment!
 *                                    Use free().
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
int edlibCalcEditDistance(const unsigned char* query, int queryLength,
                          const unsigned char* target, int targetLength,
                          int alphabetLength, int k, int mode,
                          bool findStartLocations,
                          bool findAlignment, int* bestScore,
                          int** endLocations, int** startLocations,
                          int* numLocations, unsigned char** alignment,
                          int* alignmentLength);

/**
 * Builds cigar string from given alignment sequence.
 * @param [in] alignment  Alignment sequence.
 *     0 stands for match.
 *     1 stands for insertion to target.
 *     2 stands for insertion to query.
 *     3 stands for mismatch.
 * @param [in] alignmentLength
 * @param [in] cigarFormat
 *     If EDLIB_CIGAR_EXTENDED, extended cigar is returned.
 *     If EDLIB_CIGAR_STANDARD, standard cigar is returned (contains only I, D and M).
 * @param [out] cigar  Will contain cigar string.
 *     I stands for insertion.
 *     D stands for deletion.
 *     X stands for mismatch. (used only in extended format)
 *     = stands for match. (used only in extended format)
 *     M stands for (mis)match. (used only in standard format)
 *     String is null terminated.
 *     Needed memory is allocated and given pointer is set to it.
 *     Do not forget to free it later using free()!
 * @return Status code.
 */
int edlibAlignmentToCigar(unsigned char* alignment, int alignmentLength,
                          int cigarFormat, char** cigar);

#ifdef __cplusplus
}
#endif

#endif // EDLIB_H

#include <stdint.h>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cassert>

using namespace std;

typedef uint64_t Word;
static const int WORD_SIZE = sizeof(Word) * 8;  // Size of Word in bits
static const Word WORD_1 = (Word) 1;
static const Word HIGH_BIT_MASK = WORD_1 << (WORD_SIZE - 1);  // 100..00

// Data needed to find alignment.
struct AlignmentData {
    Word* Ps;
    Word* Ms;
    int* scores;
    int* firstBlocks;
    int* lastBlocks;

    AlignmentData(int maxNumBlocks, int targetLength) {
        Ps = new Word[maxNumBlocks * targetLength];
        Ms = new Word[maxNumBlocks * targetLength];
        scores = new int[maxNumBlocks * targetLength];
        firstBlocks = new int[targetLength];
        lastBlocks = new int[targetLength];
    }

    ~AlignmentData() {
        delete[] Ps;
        delete[] Ms;
        delete[] scores;
        delete[] firstBlocks;
        delete[] lastBlocks;
    }
};

struct Block {
    Word P;  // Pvin
    Word M;  // Mvin
    int score;  // score of last cell in block;
};

static int myersCalcEditDistanceSemiGlobal(Block* blocks, Word* Peq, int W,
                                           int maxNumBlocks,
                                           const unsigned char* query,
                                           int queryLength,
                                           const unsigned char* target,
                                           int targetLength, int alphabetLength,
                                           int k, int mode, int* bestScore,
                                           int** positions, int* numPositions);

static int myersCalcEditDistanceNW(Block* blocks, Word* Peq, int W,
                                   int maxNumBlocks, const unsigned char* query,
                                   int queryLength, const unsigned char* target,
                                   int targetLength, int alphabetLength, int k,
                                   int* bestScore, int* position,
                                   bool findAlignment,
                                   AlignmentData** alignData);

static void obtainAlignment(int maxNumBlocks, int queryLength, int targetLength,
                            int W, int bestScore, int position,
                            AlignmentData* alignData, unsigned char** alignment,
                            int* alignmentLength);

static inline int ceilDiv(int x, int y);

static inline unsigned char* createReverseCopy(const unsigned char* seq,
                                               int length);

static inline Word* buildPeq(int alphabetLength, const unsigned char* query,
                             int queryLength);

/**
 * Entry function.
 */
int edlibCalcEditDistance(const unsigned char* query, int queryLength,
                          const unsigned char* target, int targetLength,
                          int alphabetLength, int k, int mode,
                          bool findStartLocations,
                          bool findAlignment, int* bestScore,
                          int** endLocations, int** startLocations,
                          int* numLocations, unsigned char** alignment,
                          int* alignmentLength) {

    *alignment = NULL;
    /*--------------------- INITIALIZATION ------------------*/
    int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);  // bmax in Myers
    int W = maxNumBlocks * WORD_SIZE - queryLength;  // number of redundant cells in last level blocks

    Block* blocks = new Block[maxNumBlocks];
    Word* Peq = buildPeq(alphabetLength, query, queryLength);
    /*-------------------------------------------------------*/

    /*------------------ MAIN CALCULATION -------------------*/
    // TODO: Store alignment data only after k is determined? That could make things faster.
    *bestScore = -1;
    *endLocations = *startLocations = NULL;
    *numLocations = 0;
    int positionNW;  // Used only when mode is NW.
    AlignmentData* alignData = NULL;
    bool dynamicK = false;
    if(k < 0) {  // If valid k is not given, auto-adjust k until solution is found.
        dynamicK = true;
        k = WORD_SIZE;  // Gives better results then smaller k.
    }

    do {
        if(alignData) delete alignData;
        if(mode == EDLIB_MODE_HW || mode == EDLIB_MODE_SHW) {
            myersCalcEditDistanceSemiGlobal(blocks, Peq, W, maxNumBlocks, query,
                                            queryLength, target, targetLength,
                                            alphabetLength, k, mode, bestScore,
                                            endLocations, numLocations);
        } else {  // mode == EDLIB_MODE_NW
            myersCalcEditDistanceNW(blocks, Peq, W, maxNumBlocks, query,
                                    queryLength, target, targetLength,
                                    alphabetLength, k, bestScore, &positionNW,
                                    findAlignment, &alignData);
        }
        k *= 2;
    } while (dynamicK && *bestScore == -1);

    if(*bestScore >= 0) {  // If there is solution.
        // If NW mode, set end location explicitly.
        if(mode == EDLIB_MODE_NW) {
            *endLocations = (int *) malloc(sizeof(int) * 1);
            (*endLocations)[0] = targetLength - 1;
            *numLocations = 1;
        }

        // Find starting locations.
        if(findStartLocations || findAlignment) {
            *startLocations = (int*) malloc((*numLocations) * sizeof(int));
            if(mode == EDLIB_MODE_HW) {  // If HW, I need to calculate start locations.
                const unsigned char* rTarget = createReverseCopy(target,
                                                                 targetLength);
                const unsigned char* rQuery = createReverseCopy(query,
                                                                queryLength);
                Word* rPeq = buildPeq(alphabetLength, rQuery, queryLength);  // Peq for reversed query
                for (int i = 0; i < *numLocations; i++) {
                    int endLocation = (*endLocations)[i];
                    int bestScoreSHW, numPositionsSHW;
                    int* positionsSHW;
                    myersCalcEditDistanceSemiGlobal(
                            blocks, rPeq, W, maxNumBlocks, rQuery, queryLength,
                            rTarget + targetLength - endLocation - 1,
                            endLocation + 1, alphabetLength, *bestScore,
                            EDLIB_MODE_SHW,
                            &bestScoreSHW, &positionsSHW, &numPositionsSHW);
                    // Taking last location as start ensures that alignment will not start with insertions
                    // if it can start with mismatches instead.
                    if(numPositionsSHW != 0) {
                        (*startLocations)[i] = endLocation
                                - positionsSHW[numPositionsSHW - 1];
                        delete[] positionsSHW;
                    }
                }
                delete[] rTarget;
                delete[] rQuery;
                delete[] rPeq;
            } else {  // If mode is SHW or NW
                for (int i = 0; i < *numLocations; i++) {
                    (*startLocations)[i] = 0;
                }
            }
        }

        // Find alignment -> all comes down to finding alignment for NW.
        // Currently we return alignment only for first pair of locations.
        if(findAlignment) {
            int alnStartLocation = (*startLocations)[0];
            int alnEndLocation = (*endLocations)[0];
            if(mode != EDLIB_MODE_NW) {  // Calculate align data.
                int score_, endLocation_;  // Used only to call function.
                myersCalcEditDistanceNW(blocks, Peq, W, maxNumBlocks, query,
                                        queryLength, target + alnStartLocation,
                                        alnEndLocation - alnStartLocation + 1,
                                        alphabetLength, *bestScore, &score_,
                                        &endLocation_, true, &alignData);
                //assert(score_ == *bestScore);
                //assert(endLocation_ == alnEndLocation - alnStartLocation);
            }
            obtainAlignment(maxNumBlocks, queryLength,
                            alnEndLocation - alnStartLocation + 1, W,
                            *bestScore, alnEndLocation - alnStartLocation,
                            alignData, alignment, alignmentLength);
        }
    }
    /*-------------------------------------------------------*/

    //--- Free memory ---//
    delete[] blocks;
    delete[] Peq;
    if(alignData) delete alignData;
    //-------------------//

    return EDLIB_STATUS_OK;
}

int edlibAlignmentToCigar(unsigned char* alignment, int alignmentLength,
                          int cigarFormat, char** cigar_) {
    *cigar_ = NULL;
    if(cigarFormat != EDLIB_CIGAR_EXTENDED
            && cigarFormat != EDLIB_CIGAR_STANDARD) {
        return EDLIB_STATUS_ERROR;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = { '=', 'I', 'D', 'X' };
    if(cigarFormat == EDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    vector<char>* cigar = new vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if(i == alignmentLength
                || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if(i < alignmentLength) {
                // Check if alignment has valid values.
                if(alignment[i] < 0 || alignment[i] > 3) {
                    delete cigar;
                    return EDLIB_STATUS_ERROR;
                }
                numOfSameMoves = 0;
            }
        }
        if(i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0);  // Null character termination.
    *cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    memcpy(*cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return EDLIB_STATUS_OK;
}

/**
 * Build Peq table for given query and alphabet.
 * Peq is table of dimensions alphabetLength+1 x maxNumBlocks.
 * Bit i of Peq[s * maxNumBlocks + b] is 1 if i-th symbol from block b of query equals symbol s, otherwise it is 0.
 * NOTICE: free returned array with delete[]!
 */
static inline Word* buildPeq(int alphabetLength, const unsigned char* query,
                             int queryLength) {
    int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
    int W = maxNumBlocks * WORD_SIZE - queryLength;  // number of redundant cells in last level blocks
    // table of dimensions alphabetLength+1 x maxNumBlocks. Last symbol is wildcard.
    Word* Peq = new Word[(alphabetLength + 1) * maxNumBlocks];

    // Build Peq (1 is match, 0 is mismatch). NOTE: last column is wildcard(symbol that matches anything) with just 1s
    for (int symbol = 0; symbol <= alphabetLength; symbol++) {
        for (int b = 0; b < maxNumBlocks; b++) {
            if(symbol < alphabetLength) {
                Peq[symbol * maxNumBlocks + b] = 0;
                for (int r = (b + 1) * WORD_SIZE - 1; r >= b * WORD_SIZE; r--) {
                    Peq[symbol * maxNumBlocks + b] <<= 1;
                    // NOTE: We pretend like query is padded at the end with W wildcard symbols
                    if(r >= queryLength || query[r] == symbol) Peq[symbol
                            * maxNumBlocks + b] += 1;
                }
            } else {  // Last symbol is wildcard, so it is all 1s
                Peq[symbol * maxNumBlocks + b] = (Word) -1;
            }
        }
    }

    return Peq;
}

/**
 * Returns new sequence that is reverse of given sequence.
 */
static inline unsigned char* createReverseCopy(const unsigned char* seq,
                                               int length) {
    unsigned char* rSeq = new unsigned char[length];
    for (int i = 0; i < length; i++) {
        rSeq[i] = seq[length - i - 1];
    }
    return rSeq;
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word &PvOut, Word &MvOut) {
    // hin can be 1, -1 or 0.
    // 1  -> 00...01
    // 0  -> 00...00
    // -1 -> 11...11 (2-complement)

    Word hinIsNeg = (Word) (hin >> 2) & WORD_1;  // 00...001 if hin is -1, 00...000 if 0 or 1

    Word Xv = Eq | Mv;
    // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
    Eq |= hinIsNeg;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
    hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
    // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
    hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);

    Ph <<= 1;
    Mh <<= 1;

    // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
    Mh |= hinIsNeg;
    // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
    Ph |= (Word) ((hin + 1) >> 1);

    PvOut = Mh | ~(Xv | Ph);
    MvOut = Ph & Xv;

    return hout;
}

/**
 * Does ceiling division x / y.
 * Note: x and y must be non-negative and x + y must not overflow.
 */
static inline int ceilDiv(int x, int y) {
    return x % y ? x / y + 1 : x / y;
}

static inline int min(int x, int y) {
    return x < y ? x : y;
}

static inline int max(int x, int y) {
    return x > y ? x : y;
}

/**
 * @param [in] block
 * @return Values of cells in block, starting with bottom cell in block.
 */
static inline vector<int> getBlockCellValues(Block* const block) {
    vector<int> scores(WORD_SIZE);
    int score = block->score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; i++) {
        scores[i] = score;
        if(block->P & mask) score--;
        if(block->M & mask) score++;
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
    return scores;
}

/**
 * @param [in] block
 * @param [in] k
 * @return True if all cells in block have value larger than k, otherwise false.
 */
static inline bool allBlockCellsLarger(Block* const block, const int k) {
    vector<int> scores = getBlockCellValues(block);
    for (int i = 0; i < WORD_SIZE; i++) {
        if(scores[i] <= k) return false;
    }
    return true;
}

/**
 * @param [in] mode  EDLIB_MODE_HW or EDLIB_MODE_SHW or EDLIB_MODE_OV
 */
static int myersCalcEditDistanceSemiGlobal(Block* const blocks, Word* const Peq,
                                           const int W, const int maxNumBlocks,
                                           const unsigned char* const query,
                                           const int queryLength,
                                           const unsigned char* const target,
                                           const int targetLength,
                                           const int alphabetLength, int k,
                                           const int mode, int* bestScore_,
                                           int** positions_,
                                           int* numPositions_) {
    *positions_ = NULL;
    *numPositions_ = 0;

    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of last block in Ukkonen band.
    int firstBlock = 0;
    int lastBlock = min(ceilDiv(k + 1, WORD_SIZE), maxNumBlocks) - 1;  // y in Myers
    Block *bl;  // Current block

    // For HW, solution will never be larger then queryLength.
    if(mode == EDLIB_MODE_HW) {
        k = min(queryLength, k);
    }

    // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
    // This gives speed up of about 2 times for small k.
    const int STRONG_REDUCE_NUM = 2048;

    // Initialize P, M and score
    bl = blocks;
    for (int b = 0; b <= lastBlock; b++) {
        bl->score = (b + 1) * WORD_SIZE;
        bl->P = (Word) -1;  // All 1s
        bl->M = (Word) 0;
        bl++;
    }

    int bestScore = -1;
    vector<int> positions;  // TODO: Maybe put this on heap?
    const int startHout = mode == EDLIB_MODE_HW ? 0 : 1;  // If 0 then gap before query is not penalized;
    const unsigned char* targetChar = target;
    for (int c = 0; c < targetLength; c++) {  // for each column
        const Word* Peq_c = Peq + (*targetChar) * maxNumBlocks;

        //----------------------- Calculate column -------------------------//
        int hout = startHout;
        bl = blocks + firstBlock;
        Peq_c += firstBlock;
        for (int b = firstBlock; b <= lastBlock; b++) {
            hout = calculateBlock(bl->P, bl->M, *Peq_c, hout, bl->P, bl->M);
            bl->score += hout;
            bl++;
            Peq_c++;
        }
        bl--;
        Peq_c--;
        //------------------------------------------------------------------//

        //---------- Adjust number of blocks according to Ukkonen ----------//
        if((lastBlock < maxNumBlocks - 1) && (bl->score - hout <= k)  // bl is pointing to last block
                && ((*(Peq_c + 1) & WORD_1) || hout < 0)) {  // Peq_c is pointing to last block
            // If score of left block is not too big, calculate one more block
            lastBlock++;
            bl++;
            Peq_c++;
            bl->P = (Word) -1;  // All 1s
            bl->M = (Word) 0;
            bl->score = (bl - 1)->score - hout + WORD_SIZE
                    + calculateBlock(bl->P, bl->M, *Peq_c, hout, bl->P, bl->M);
        } else {
            while (lastBlock >= firstBlock && bl->score >= k + WORD_SIZE) {
                lastBlock--;
                bl--;
                Peq_c--;
            }
        }

        // Every some columns, do some expensive but also more efficient block reducing -> this is important!
        if(c % STRONG_REDUCE_NUM == 0) {
            while (lastBlock >= firstBlock && allBlockCellsLarger(bl, k)) {
                lastBlock--;
                bl--;
                Peq_c--;
            }
        }

        if(mode != EDLIB_MODE_HW) {
            while (firstBlock <= lastBlock
                    && blocks[firstBlock].score >= k + WORD_SIZE) {
                firstBlock++;
            }
            if(c % STRONG_REDUCE_NUM == 0) {  // Do strong reduction every some blocks
                while (firstBlock <= lastBlock
                        && allBlockCellsLarger(blocks + firstBlock, k)) {
                    firstBlock++;
                }
            }
        }

        // For HW, even if all cells are > k, there still may be solution in next
        // column because starting conditions at upper boundary are 0.
        // That means that first block is always candidate for solution,
        // and we can never end calculation before last column.
        if(mode == EDLIB_MODE_HW) {
            lastBlock = max(0, lastBlock);
        }

        // If band stops to exist finish
        if(lastBlock < firstBlock) {
            *bestScore_ = bestScore;
            if(bestScore != -1) {
                *positions_ = (int *) malloc(sizeof(int) * positions.size());
                *numPositions_ = positions.size();
                copy(positions.begin(), positions.end(), *positions_);
            }
            return EDLIB_STATUS_OK;
        }
        //------------------------------------------------------------------//

        //------------------------- Update best score ----------------------//
        if(lastBlock == maxNumBlocks - 1) {
            int colScore = bl->score;
            if(colScore <= k) {  // Scores > k dont have correct values (so we cannot use them), but are certainly > k.
                // NOTE: Score that I find in column c is actually score from column c-W
                if(bestScore == -1 || colScore <= bestScore) {
                    if(colScore != bestScore) {
                        positions.clear();
                        bestScore = colScore;
                        // Change k so we will look only for equal or better
                        // scores then the best found so far.
                        k = bestScore;
                    }
                    positions.push_back(c - W);
                }
            }
        }
        //------------------------------------------------------------------//

        targetChar++;
    }

    // Obtain results for last W columns from last column.
    if(lastBlock == maxNumBlocks - 1) {
        vector<int> blockScores = getBlockCellValues(bl);
        for (int i = 0; i < W; i++) {
            int colScore = blockScores[i + 1];
            if(colScore <= k && (bestScore == -1 || colScore <= bestScore)) {
                if(colScore != bestScore) {
                    positions.clear();
                    k = bestScore = colScore;
                }
                positions.push_back(targetLength - W + i);
            }
        }
    }

    *bestScore_ = bestScore;
    if(bestScore != -1) {
        *positions_ = (int *) malloc(sizeof(int) * positions.size());
        *numPositions_ = positions.size();
        copy(positions.begin(), positions.end(), *positions_);
    }

    return EDLIB_STATUS_OK;
}

/**
 * @param alignData  Data generated during calculation, that is needed for reconstruction of alignment.
 *                   I it is allocated with new, so free it with delete.
 *                   Data is generated only if findAlignment is true.
 */
static int myersCalcEditDistanceNW(Block* blocks, Word* Peq, int W,
                                   int maxNumBlocks, const unsigned char* query,
                                   int queryLength, const unsigned char* target,
                                   int targetLength, int alphabetLength, int k,
                                   int* bestScore_, int* position_,
                                   bool findAlignment,
                                   AlignmentData** alignData) {

    // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
    const int STRONG_REDUCE_NUM = 2048;  // TODO: Choose this number dinamically (based on query and target lengths?), so it does not affect speed of computation

    if(k < abs(targetLength - queryLength)) {
        *bestScore_ = *position_ = -1;
        return EDLIB_STATUS_OK;
    }

    k = min(k, max(queryLength, targetLength));  // Upper bound for k

    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of last block in Ukkonen band.
    int firstBlock = 0;
    // This is optimal now, by my formula.
    int lastBlock = min(
            maxNumBlocks,
            ceilDiv(min(k, (k + queryLength - targetLength) / 2) + 1,
                    WORD_SIZE)) - 1;
    Block* bl;  // Current block

    // Initialize P, M and score
    bl = blocks;
    for (int b = 0; b <= lastBlock; b++) {
        bl->score = (b + 1) * WORD_SIZE;
        bl->P = (Word) -1;  // All 1s
        bl->M = (Word) 0;
        bl++;
    }

    // If we want to find alignment, we have to store needed data.
    if(findAlignment)
        *alignData = new AlignmentData(maxNumBlocks, targetLength);
    else *alignData = NULL;

    const unsigned char* targetChar = target;
    for (int c = 0; c < targetLength; c++) {  // for each column
        Word* Peq_c = Peq + *targetChar * maxNumBlocks;

        //----------------------- Calculate column -------------------------//
        int hout = 1;
        bl = blocks + firstBlock;
        for (int b = firstBlock; b <= lastBlock; b++) {
            hout = calculateBlock(bl->P, bl->M, Peq_c[b], hout, bl->P, bl->M);
            bl->score += hout;
            bl++;
        }
        bl--;
        //------------------------------------------------------------------//
        // bl now points to last block

        // Update k. I do it only on end of column because it would slow calculation too much otherwise.
        // NOTICE: I add W when in last block because it is actually result from W cells to the left and W cells up.
        k = min(k,
                bl->score
                        + max(targetLength - c - 1,
                              queryLength - ((1 + lastBlock) * WORD_SIZE - 1)
                                      - 1)
                        + (lastBlock == maxNumBlocks - 1 ? W : 0));

        //---------- Adjust number of blocks according to Ukkonen ----------//
        //--- Adjust last block ---//
        // If block is not beneath band, calculate next block. Only next because others are certainly beneath band.
        if(lastBlock + 1 < maxNumBlocks && !(  //score[lastBlock] >= k + WORD_SIZE ||  // NOTICE: this condition could be satisfied if above block also!
                ((lastBlock + 1) * WORD_SIZE - 1
                        > k - bl->score + 2 * WORD_SIZE - 2 - targetLength + c
                                + queryLength))) {
            lastBlock++;
            bl++;
            bl->P = (Word) -1;  // All 1s
            bl->M = (Word) 0;
            int newHout = calculateBlock(bl->P, bl->M, Peq_c[lastBlock], hout,
                                         bl->P, bl->M);
            bl->score = (bl - 1)->score - hout + WORD_SIZE + newHout;
            hout = newHout;
        }

        // While block is out of band, move one block up. - This is optimal now, by my formula.
        // NOTICE: I added + W, and now it works! This has to be added because query is padded with W cells.
        while (lastBlock >= firstBlock
                && (bl->score >= k + WORD_SIZE
                        || ((lastBlock + 1) * WORD_SIZE - 1
                                > k - bl->score + 2 * WORD_SIZE - 2
                                        - targetLength + c + queryLength + W))) {
            lastBlock--;
            bl--;
        }
        //-------------------------//

        //--- Adjust first block ---//
        // While outside of band, advance block
        while (firstBlock <= lastBlock
                && (blocks[firstBlock].score >= k + WORD_SIZE
                        || ((firstBlock + 1) * WORD_SIZE - 1
                                < blocks[firstBlock].score - k - targetLength
                                        + queryLength + c))) {
            firstBlock++;
        }
        //--------------------------/

        // TODO: consider if this part is useful, it does not seem to help much
        if(c % STRONG_REDUCE_NUM == 0) {  // Every some columns do more expensive but more efficient reduction
            while (lastBlock >= firstBlock) {
                // If all cells outside of band, remove block
                vector<int> scores = getBlockCellValues(bl);
                int r = (lastBlock + 1) * WORD_SIZE - 1;
                bool reduce = true;
                for (int i = 0; i < WORD_SIZE; i++) {
                    // TODO: Does not work if do not put +1! Why???
                    if(scores[i] <= k
                            && r
                                    <= k - scores[i] - targetLength + c
                                            + queryLength + W + 1) {
                        reduce = false;
                        break;
                    }
                    r--;
                }
                if(!reduce) break;
                lastBlock--;
                bl--;
            }

            while (firstBlock <= lastBlock) {
                // If all cells outside of band, remove block
                vector<int> scores = getBlockCellValues(blocks + firstBlock);
                int r = (firstBlock + 1) * WORD_SIZE - 1;
                bool reduce = true;
                for (int i = 0; i < WORD_SIZE; i++) {
                    if(scores[i] <= k
                            && r
                                    >= scores[i] - k - targetLength + c
                                            + queryLength) {
                        reduce = false;
                        break;
                    }
                    r--;
                }
                if(!reduce) break;
                firstBlock++;
            }
        }

        // If band stops to exist finish
        if(lastBlock < firstBlock) {
            *bestScore_ = *position_ = -1;
            return EDLIB_STATUS_OK;
        }
        //------------------------------------------------------------------//

        //---- Save column so it can be used for reconstruction ----//
        if(findAlignment && c < targetLength) {
            bl = blocks + firstBlock;
            for (int b = firstBlock; b <= lastBlock; b++) {
                (*alignData)->Ps[maxNumBlocks * c + b] = bl->P;
                (*alignData)->Ms[maxNumBlocks * c + b] = bl->M;
                (*alignData)->scores[maxNumBlocks * c + b] = bl->score;
                (*alignData)->firstBlocks[c] = firstBlock;
                (*alignData)->lastBlocks[c] = lastBlock;
                bl++;
            }
        }
        //----------------------------------------------------------//

        targetChar++;
    }

    if(lastBlock == maxNumBlocks - 1) {  // If last block of last column was calculated
    // Obtain best score from block -> it is complicated because query is padded with W cells
        int bestScore = getBlockCellValues(blocks + lastBlock)[W];
        if(bestScore <= k) {
            *bestScore_ = bestScore;
            *position_ = targetLength - 1;
            return EDLIB_STATUS_OK;
        }
    }

    *bestScore_ = *position_ = -1;
    return EDLIB_STATUS_OK;
}

/**
 * Finds one possible alignment that gives optimal score.
 * @param [in] maxNumBlocks
 * @param [in] queryLength  Normal length, without W.
 * @param [in] targetLength  Normal length, without W.
 * @param [in] W  Padding.
 * @param [in] bestScore  Best score.
 * @param [in] position  Position in target where best score was found.
 * @param [in] alignData  Data obtained during finding best score that is useful for finding alignment.
 * @param [out] alignment  Alignment.
 * @param [out] alignmentLength  Length of alignment.
 */
static void obtainAlignment(int maxNumBlocks, int queryLength, int targetLength,
                            int W, int bestScore, int position,
                            AlignmentData* alignData, unsigned char** alignment,
                            int* alignmentLength) {
    if(alignData == NULL) return;
    *alignment = (unsigned char*) malloc(
            (queryLength + targetLength - 1) * sizeof(unsigned char));
    *alignmentLength = 0;
    // TODO(martinsos): is this position needed? If this function is used for NW only,
    // then I can assume value of position is equal to targetLength - 1 and remove position argument.
    int c = position;  // index of column
    int b = maxNumBlocks - 1;  // index of block in column
    int currScore = bestScore;  // Score of current cell
    int lScore = -1;  // Score of left cell
    int uScore = -1;  // Score of upper cell
    int ulScore = -1;  // Score of upper left cell
    Word currP = alignData->Ps[c * maxNumBlocks + b];  // P of current block
    Word currM = alignData->Ms[c * maxNumBlocks + b];  // M of current block
    // True if block to left exists and is in band
    bool thereIsLeftBlock = c > 0 && b >= alignData->firstBlocks[c - 1]
            && b <= alignData->lastBlocks[c - 1];
    Word lP, lM;
    if(thereIsLeftBlock) {
        lP = alignData->Ps[(c - 1) * maxNumBlocks + b];  // P of block to the left
        lM = alignData->Ms[(c - 1) * maxNumBlocks + b];  // M of block to the left
    }
    currP <<= W;
    currM <<= W;
    int blockPos = WORD_SIZE - W - 1;  // 0 based index of current cell in blockPos
    if(c == 0) {
        thereIsLeftBlock = true;
        lScore = b * WORD_SIZE + blockPos + 1;
        ulScore = lScore - 1;
    }
    while (true) {
        // TODO: improvement: calculate only those cells that are needed,
        //       for example if I calculate upper cell and can move up,
        //       there is no need to calculate left and upper left cell
        //---------- Calculate scores ---------//
        if(lScore == -1 && thereIsLeftBlock) {
            lScore = alignData->scores[(c - 1) * maxNumBlocks + b];  // score of block to the left
            for (int i = 0; i < WORD_SIZE - blockPos - 1; i++) {
                if(lP & HIGH_BIT_MASK) lScore--;
                if(lM & HIGH_BIT_MASK) lScore++;
                lP <<= 1;
                lM <<= 1;
            }
        }
        if(ulScore == -1) {
            if(lScore != -1) {
                ulScore = lScore;
                if(lP & HIGH_BIT_MASK) ulScore--;
                if(lM & HIGH_BIT_MASK) ulScore++;
            } else if(c > 0 && b - 1 >= alignData->firstBlocks[c - 1]
                    && b - 1 <= alignData->lastBlocks[c - 1]) {
                // This is the case when upper left cell is last cell in block,
                // and block to left is not in band so lScore is -1.
                ulScore = alignData->scores[(c - 1) * maxNumBlocks + b - 1];
            }
        }
        if(uScore == -1) {
            uScore = currScore;
            if(currP & HIGH_BIT_MASK) uScore--;
            if(currM & HIGH_BIT_MASK) uScore++;
            currP <<= 1;
            currM <<= 1;
        }
        //-------------------------------------//

        // TODO: should I check if there is upper block?

        //-------------- Move --------------//
        // Move up - insertion to target - deletion from query
        if(uScore != -1 && uScore + 1 == currScore) {
            currScore = uScore;
            lScore = ulScore;
            uScore = ulScore = -1;
            if(blockPos == 0) {  // If entering new (upper) block
                if(b == 0) {  // If there are no cells above (only boundary cells)
                    (*alignment)[(*alignmentLength)++] = 1;  // Move up
                    for (int i = 0; i < c + 1; i++)  // Move left until end
                        (*alignment)[(*alignmentLength)++] = 2;
                    break;
                } else {
                    blockPos = WORD_SIZE - 1;
                    b--;
                    currP = alignData->Ps[c * maxNumBlocks + b];
                    currM = alignData->Ms[c * maxNumBlocks + b];
                    if(c > 0 && b >= alignData->firstBlocks[c - 1]
                            && b <= alignData->lastBlocks[c - 1]) {
                        thereIsLeftBlock = true;
                        lP = alignData->Ps[(c - 1) * maxNumBlocks + b];  // TODO: improve this, too many operations
                        lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
                    } else {
                        thereIsLeftBlock = false;
                    }
                }
            } else {
                blockPos--;
                lP <<= 1;
                lM <<= 1;
            }
            // Mark move
            (*alignment)[(*alignmentLength)++] = 1;  // TODO: enumeration?
        }
        // Move left - deletion from target - insertion to query
        else if(lScore != -1 && lScore + 1 == currScore) {
            currScore = lScore;
            uScore = ulScore;
            lScore = ulScore = -1;
            c--;
            if(c == -1) {  // If there are no cells to the left (only boundary cells)
                (*alignment)[(*alignmentLength)++] = 2;  // Move left
                int numUp = b * WORD_SIZE + blockPos + 1;
                for (int i = 0; i < numUp; i++)  // Move up until end
                    (*alignment)[(*alignmentLength)++] = 1;
                break;
            }
            currP = lP;
            currM = lM;
            if(c > 0 && b >= alignData->firstBlocks[c - 1]
                    && b <= alignData->lastBlocks[c - 1]) {
                thereIsLeftBlock = true;
                lP = alignData->Ps[(c - 1) * maxNumBlocks + b];
                lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
            } else {
                if(c == 0) {  // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = true;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = false;
                }
            }
            // Mark move
            (*alignment)[(*alignmentLength)++] = 2;
        }
        // Move up left - (mis)match
        else if(ulScore != -1) {
            unsigned char moveCode = ulScore == currScore ? 0 : 3;  // 0 for match, 3 for mismatch.
            currScore = ulScore;
            uScore = lScore = ulScore = -1;
            c--;
            if(c == -1) {  // If there are no cells to the left (only boundary cells)
                (*alignment)[(*alignmentLength)++] = moveCode;  // Move left
                int numUp = b * WORD_SIZE + blockPos;
                for (int i = 0; i < numUp; i++)  // Move up until end
                    (*alignment)[(*alignmentLength)++] = 1;
                break;
            }
            if(blockPos == 0) {  // If entering upper left block
                if(b == 0) {  // If there are no more cells above (only boundary cells)
                    (*alignment)[(*alignmentLength)++] = moveCode;  // Move up left
                    for (int i = 0; i < c + 1; i++)  // Move left until end
                        (*alignment)[(*alignmentLength)++] = 2;
                    break;
                }
                blockPos = WORD_SIZE - 1;
                b--;
                currP = alignData->Ps[c * maxNumBlocks + b];
                currM = alignData->Ms[c * maxNumBlocks + b];
            } else {  // If entering left block
                blockPos--;
                currP = lP;
                currM = lM;
                currP <<= 1;
                currM <<= 1;
            }
            // Set new left block
            if(c > 0 && b >= alignData->firstBlocks[c - 1]
                    && b <= alignData->lastBlocks[c - 1]) {
                thereIsLeftBlock = true;
                lP = alignData->Ps[(c - 1) * maxNumBlocks + b];
                lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
            } else {
                if(c == 0) {  // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = true;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = false;
                }
            }
            // Mark move
            (*alignment)[(*alignmentLength)++] = moveCode;
        } else {
            // Reached end - finished!
            break;
        }
        //----------------------------------//
    }

    *alignment = (unsigned char*) realloc(
            *alignment, (*alignmentLength) * sizeof(unsigned char));
    reverse(*alignment, *alignment + (*alignmentLength));
    return;
}

/*
 * bioutil.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 *
 *  Header contains useful functions for working with DNA sequences.
 */

#ifndef BIOUTIL_H_
#define BIOUTIL_H_

#include <cstring>
#include <cctype>
#include <cassert>
#include <cstdio>
#include <cstdint>
#include <climits>
#include <map>

/**
 * Function tests if base passed as argument is valid.
 * Valid bases are A, C, G, T and N.
 *
 * @param base base that is tested
 * @return true if base is valid, false otherwise
 */
static inline bool isValidBaseACGT(char base) {
    base = toupper(base);
    return base == 'A' || base == 'C' || base == 'G' || base == 'T'
            || base == 'N';

}

/**
 * Function transforms base to small int number.
 *  '0' - 'A'
 *  '1' - 'C'
 *  '2' - 'G'
 *  '3' - 'T'
 *  '4' - 'N'
 *
 * @param base base
 * @return number assigned to base
 */
static inline uint8_t baseToInt(char base) {
    base = toupper(base);
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
            return 4;
        default:
            return UCHAR_MAX;
    }

}
/**
 * Function transforms base from small int number to ACTG base.
 *  '0' - 'A'
 *  '1' - 'C'
 *  '2' - 'G'
 *  '3' - 'T'
 *  '4' - 'N'
 *
 * @param base base
 * @return number assigned to base
 */
inline char intToBase(const uint8_t num) {
    switch (num) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 4:
            return 'N';
        default:
            return 0;
    }
}
/**
 * Function that calculates complement of given base.
 * Valid bases are A, C, G, T and N.
 *
 * @param base base whose complement is being calculated
 * @return complement of given base
 */
inline char getACGTComplement(char base) {
    base = toupper(base);
    assert(isValidBaseACGT(base));

    switch (base) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
    }
}

/**
 * Function that calculates complement of base as small int number.
 * Valid bases are A, C, G, T and N given as small int numbers 0, 1, 2, 3
 *
 * @param base base whose complement is being calculated given as small int number
 * @return complement of given base also as small int
 */
inline uint8_t getACGTComplementAsSmallInt(uint8_t base) {

    switch (base) {
        case 0:
            // A to T
            return 3;
        case 1:
            // C to G
            return 2;
        case 2:
            // G to C
            return 1;
        case 3:
            // T to A
            return 0;
        default:
            return 4;
    }

}

#endif /* UTIL_H_ */
/*
 * utility_functions.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 *
 */
#ifndef SRC_UTIL_UTILITY_H_
#define SRC_UTIL_UTILITY_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <dirent.h>
#include <unistd.h>

inline bool validateOutputFile(char*filePath) {
    fprintf(stderr, "W %s\n", filePath);
    FILE* f = fopen(filePath, "w");
    if(f != NULL) {
        fclose(f);
        return true;
    }
    return false;
}

inline bool validateInputFile(char*filePath) {
    fprintf(stderr, "R %s\n", filePath);
    FILE* f = fopen(filePath, "r");

    if(f != NULL) {
        fclose(f);
        return true;
    }
    return false;
}

inline bool isFolder(char*folderPath) {
    DIR* dir = opendir(folderPath);
    if(dir != NULL) {
        closedir(dir);
        return true;
    }
    return false;
}

template<typename T, typename Pred = std::less<T> >
struct ptr_compare : Pred {
    ptr_compare(Pred const & p = Pred())
            : Pred(p) {
    }

    bool operator()(T const * p1, T const * p2) const {
        return Pred::operator()(*p1, *p2);
    }
};

#endif
/*
 * util_structures.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_UTIL_STRUCTURES_H_
#define SRC_UTIL_STRUCTURES_H_

struct event_t {
    uint32_t first;
    uint32_t second;bool isStart;
    uint32_t index;

    event_t(uint32_t first, uint32_t second, bool isStart, uint32_t index)
            : first(first),
              second(second),
              isStart(isStart),
              index(index) {
    }

    bool operator<(const event_t &b) const {
        if(first != b.first) return first < b.first;

        if(second != b.second) return second < b.second;

        return isStart < b.isStart;
    }

};

struct eventK_t {
    uint32_t first;
    uint32_t second;
    uint32_t k;bool isStart;
    uint32_t index;

    eventK_t(uint32_t first, uint32_t second, uint32_t k, bool isStart,
             uint32_t index)
            : first(first),
              second(second),
              k(k),
              isStart(isStart),
              index(index) {
    }

};

template<typename T>
struct triplet_t {
    T first;
    T second;
    T third;

    triplet_t(T first, T second, T third)
            : first(first),
              second(second),
              third(third) {
    }

};

inline bool operator<(const triplet_t<uint32_t> &a,
                      const triplet_t<uint32_t> &b) {
    if(a.first != b.first) {
        return a.first < b.first;
    }

    return a.second < b.second;

}

#endif /* SRC_UTIL_STRUCTURES_H_ */
/*
 * Fenwick.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 *
 *  Template class that implements Fenwick tree data structure also known as logaritmic structure.
 */

#ifndef SRC_FENWICK_H_
#define SRC_FENWICK_H_

#include <vector>
#include <utility>

template<class T>
class Fenwick {

 public:
    /**
     * Constructor.
     * @param n number of elements or positions
     */
    Fenwick(int n) {
        elements_ = std::vector<T>(n + 1, T());
    }
    ~Fenwick() {
    }

    /**
     * Function that calculates value of lowest significant bit in given value
     * @param a value
     * @return lowest significant bit of a
     */
    inline int lobitOf(const int& a) {
        return a & -a;
    }

    /**
     * Method that inserts value val at position pos and updates data structure with new maximum element.
     * @param pos position on wich the value in inserted
     * @param val value that is inserted in data structure
     */
    void updateMax(int pos, const T& val) {
        ++pos;
        for (; pos < elements_.size(); pos += lobitOf(pos)) {
            elements_[pos] = std::max(elements_[pos], val);
        }
    }
    /**
     * Method that returns maximum value stored in data structure in position range of [0, pos].
     * @param pos upper bound of range
     * @return maximum element in given range
     */
    T getMax(int pos) {
        ++pos;
        T ret = T();
        for (; pos > 0; pos -= lobitOf(pos)) {
            ret = std::max(ret, elements_[pos]);
        }
        return ret;
    }

 private:
    std::vector<T> elements_;

};

#endif /* SRC_FENWICK_H_ */
/*
 * suffix_array.h
 *
 *  Created on: Dec 27, 2014
 *      Author: marko
 */

#ifndef DTRA_SUFFIX_ARRAY_DATABASE
#define DTRA_SUFFIX_ARRAY_DATABASE

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <vector>
#include <stdint.h>

typedef int32_t saint_t;
typedef int32_t saidx_t;
typedef uint8_t sauchar_t;

class SuffixArray {

 public:
    SuffixArray(const char* text, size_t textLen);
    SuffixArray(FILE *in, const char* text, size_t textLen);
    ~SuffixArray();

    /**
     * @param pattern to search for
     * @param length length of pattern
     * @param numOfSolutions number of found patterns in the text
     * @returns the pointer to the index of the first occurrence
     */
    const int* search(const char* pattern, int length, int* numOfSolutions);

    const int* iterativeSearch(const char* pattern, int length, int startLen,
                               int* numOfSolutions, const int solUpperLimit,
                               const int solLowerLimit, int* finalLen);
    const int* iterativeSearchDev(const char* pattern, int length, int startLen,
                                  int* numOfSolutions, const int solUpperLimit,
                                  const int solLowerLimit, int* finalLen);

    void saveSuffixArray(FILE *out);

    uint32_t size();
    const char* text();

 private:
    const char* text_;
    uint32_t textLen_;
    std::vector<saidx_t> suffix_array_;
};

#endif

/*
 * position.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_POSITION_H_
#define SRC_UTIL_POSITION_H_

#include <stdint.h>
#include <cstdlib>

class Position {

 public:
    Position();
    Position(uint32_t score, uint32_t secondaryScore, uint32_t start,
             uint32_t end, bool complemented, const char* cigar,
             uint32_t cigarLen);

    virtual ~Position();bool isComplement();
    uint32_t score() const;
    uint32_t secondaryScore() const;
    uint32_t end();
    uint32_t start() const;

    void end(uint32_t end);
    void start(uint32_t start);
    void cigar(const char* cigar, uint32_t len);
    void setComplement(bool complemented);
    void score(uint32_t score);
    void secondaryScore(uint32_t secondary_score);

    const char* cigar();bool operator <(const Position &other) const {
        if(score_ != other.score()) {
            return score_ < other.score();
        }

        if(abs(start_ - other.start()) > 20) {
            // iste pozicije ako su +- 20
            return start_ - other.start();

        }

        return 0;

    }

 private:
    uint32_t score_;
    uint32_t secondary_score_;
    uint32_t start_;
    uint32_t end_;

    bool complement_;
    char* cigar_;

};

#endif /* SRC_UTIL_POSITION_H_ */
/*
 * sequence.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <vector>

#include <cstdlib>
#include <cstdio>
#include <cstdint>

class Sequence {
 public:
    Sequence();
    ~Sequence();

    void clear();
    //void readSequencesFromFASTA(FILE* fastaIn);
    //void readSingleSequenceFromFASTA(FILE* fastaIn);
    void allBasesToSmallInt();
    void allBasesToLetters();

    const char* data();
    uint32_t dataLen();
    uint32_t numOfSequences();
    uint32_t seqLen(uint32_t index);

    uint32_t positionInSeq(uint32_t positionGlobal);
    uint32_t sequenceIndex(uint32_t positionGlobal);

    bool basesInt();

    char* data_;
    uint32_t dataLen_;
    uint32_t numOfSequences_;
    std::vector<uint32_t> seqEndIndex_;
    std::vector<int> info_;bool basesAsInt_;

};
#endif /* SEQUENCE_H_ */
/*
 * read.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SRC_CORE_READ_H_
#define SRC_CORE_READ_H_

#include <zlib.h>
#include <set>
#include <string>

#define MAX_KEEP 10
#define KEEP_FACTOR 1.50f

using namespace std;

class Read {
 public:
    Read(float keepRatio = KEEP_FACTOR, uint32_t maxPositions = MAX_KEEP);
    ~Read();

    void setData(string data, string id);
    void clear();
    //bool readNextFromFASTQ(kseq_t *seq);
    void allBasesToSmallInt();
    void allBasesToLetters();

    Read* getReverseComplement();
    void printReadSAM(FILE* outFile, Sequence* seq);

    string getReadMiniSAM(Sequence* seq);
    string generateMiniSAM(Position* pos, double score, Sequence* seq);
    vector<string> getReadMiniSAMs(Sequence* seq);

    const char* data();
    const char* id();
    uint32_t dataLen();

    void addPosition(uint32_t score, uint32_t start, uint32_t end,
    bool isComplement = false,
                     const char* cigar = NULL, uint32_t cigarLen = 0,
                     uint32_t secondaryScore = 0);
    void addPosition(Position* position);

    std::set<Position*, ptr_compare<Position> >& positions();
    Position* bestPosition(uint32_t index);
    uint32_t positionsSize();
    void keepRatio(float keepRatio);
    void maxPositions(uint32_t maxPositions);bool basesInt();

    char* id_;
    char* data_;
    uint32_t dataLen_;

    char* optional_identifier_;
    char* quality_;

    bool basesAsInt_;
    std::set<Position*, ptr_compare<Position> > positions_;

    float keepRatio_;
    uint32_t maxPositions_;

    uint32_t kmerCnt[2] = { 0, 0 };
    uint32_t kmerTotalLen[2] = { 0, 0 };

};

#endif /* SRC_CORE_READ_H_ */
/*
 * lcskpp.h
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */
#ifndef SRC_METRICS_ALGORITHM_LCSKPP_H_
#define SRC_METRICS_ALGORITHM_LCSKPP_H_

#include <vector>
#include <stdint.h>

class LCSkpp {
 public:

    /**

     * @param k parameter k of LCSk++
     * @result  vector for keeping reconstruction of LCSk++. It contains pairs of integers that
     *          correspond with positions in strings.
     * @elements vector filled with pairs (i, j) that represent starting positions of common substrings
     *            of length k. i is starting position in first string, j starting position in second string
     * @return  value of LCSk++ between two strings
     */
    static uint32_t calcLCSkpp(
            uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >& result,
            std::vector<std::pair<uint32_t, uint32_t> >& elements);

    static uint32_t estimateBeginingPosFromLCSkpp(
            std::vector<std::pair<uint32_t, uint32_t> >& reconstruction);

    static uint32_t calcLCSkppSlow(
            uint32_t k,
            std::vector<std::pair<uint32_t, uint32_t>> &reconstruction,
            std::vector<std::pair<uint32_t, uint32_t>> &matches);

    static void reconstructLCSkpp(
            std::vector<std::pair<uint32_t, uint32_t>> &elements, uint32_t k,
            std::vector<int> &prevIndex, int lastIndex, int lcskLen,
            std::vector<std::pair<uint32_t, uint32_t>> &reconstruction);
};

#endif /* SRC_METRICS_ALGORITHM_LCSKPP_H_ */
/*
 * lcskppV2.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_METRICS_ALGORITHM_LCSKPPV2_H_
#define SRC_METRICS_ALGORITHM_LCSKPPV2_H_

#include <vector>
#include <utility>
#include <cstdint>

class LCSkppV2 {
 public:

    static uint32_t calcLCSkpp(
            std::vector<std::pair<uint32_t, uint32_t>> &result,
            std::vector<triplet_t<uint32_t>> &elements);

    static uint32_t estimateBeginingPosFromLCSkpp(
            std::vector<std::pair<uint32_t, uint32_t> >& reconstruction);

 private:
    static void reconstructLCSpp(
            std::vector<triplet_t<uint32_t>> &elements,
            std::vector<int> &prevIndex, int lastIndex, int lcskLen,
            std::vector<std::pair<uint32_t, uint32_t>> &reconstruction);

};
#endif /* SRC_METRICS_ALGORITHM_LCSKPPV2_H_ */
/*
 * solver.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_SOLVER_H_
#define SRC_CORE_SOLVER_H_

class Solver {
 public:
    Solver() {
    }
    virtual ~Solver() {
    }

    virtual void printInfo() = 0;
    virtual void findReadPosition(Read* read) = 0;
    virtual void findReadPosition(Read* read, bool orientation) = 0;
};

#endif /* SRC_CORE_SOLVER_H_ */
/*
 * lcsk_solver.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_LCSK_SOLVER_H_
#define SRC_CORE_LCSK_SOLVER_H_

#include <cstdint>
#include <vector>

#define KMER_K 15
#define WINDOW_SIZE 1.08

#define SSW_MATCH 5
#define SSW_MISMATCH 4
#define SSW_GAP_OPEN 8
#define SSW_GAP_EXTEND 8

class LCSkSolver : public Solver {
 public:
    LCSkSolver(Sequence* seq);
    virtual ~LCSkSolver();

    virtual void findReadPosition(Read* read);
    virtual void findReadPosition(Read* read, bool orientation);
    void readSuffixArrayFromFile(const char* saInPath);
    virtual void printInfo();

    Sequence* seq();
    SuffixArray* sa();

    uint32_t kmerK_;
    double windowSize_;

    virtual void runLCSkpp(int startIndex, int endIndex,
                           std::vector<std::pair<uint32_t, uint32_t>> &pos,
                           Read* read);
    virtual void fillPositions(Read* read);
    virtual void getKmerPositions(
            Read* read, std::vector<std::pair<uint32_t, uint32_t>> &positions,
            int kmerStart);

    Sequence* seq_;
    SuffixArray* sa_;
};

#endif /* SRC_CORE_LCSK_SOLVER_H_ */
/*
 * incremental_lcsk_solver.h
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#ifndef SRC_CORE_INCREMENTAL_LCSK_SOLVER_H_
#define SRC_CORE_INCREMENTAL_LCSK_SOLVER_H_

#define MAX_MATCH_NUM 45
#define MIN_MATCH_NUM 20

#include <cstdint>
#include <vector>

class IncrementalLCSkSolver : public LCSkSolver {

 public:
    IncrementalLCSkSolver(Sequence* seq);
    virtual ~IncrementalLCSkSolver();
    virtual void printInfo();

    uint32_t maxMatchNum_;
    uint32_t minMatchNum_;

 private:
    virtual void runLCSkpp(int startIndex, int endIndex,
                           std::vector<triplet_t<uint32_t>> &pos, Read* read);
    virtual void fillPositions(Read* read);
    virtual uint32_t getKmerPositions(
            Read* read, std::vector<triplet_t<uint32_t>> &positions,
            int kmerStart, uint32_t inittialLen);

};

#endif /* SRC_CORE_INCREMENTAL_LCSK_SOLVER_H_ */
/*
 * mapper.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_MAPPER_H_
#define SRC_UTIL_MAPPER_H_

#include <stdint.h>
#include <omp.h>

#define MAX_TMP_NAME_LEN 50

class Mapper {

 public:

    Mapper(Sequence* seq, Solver* solver, uint32_t threadNum =
                   omp_get_num_procs(),
           float minKeepScoreRatio = KEEP_FACTOR, uint32_t maxPositionsPerRead =
           MAX_KEEP);
    void mapAllReads(char* readsInPath, char* solutionOutPath);

 private:
    // TODO REMOVE OR REFACTOR
    static void runLIS(int startIndex, int endIndex,
                       std::vector<std::pair<uint32_t, uint32_t> > &pos,
                       Read* read);

    void fillSAMHeader(FILE* out);

    void copyFromTmpFileAndDelete(char* tmpFileName, FILE* src, FILE *dest);
    void createTmpFiles(FILE* tempFiles[],
                        char tmpFileNames[][MAX_TMP_NAME_LEN],
                        uint32_t numOfFiles);
    void mergeTmpFiles(char fileNames[][MAX_TMP_NAME_LEN], FILE* tmpFiles[],
                       FILE* solutionFile, int numberOfFiles);

    Sequence* seq_;
    Solver* solver_;
    uint32_t threadNum_;

    float minKeepScoreRatio_;
    uint32_t maxPositionsPerRead_;
};
#endif

/*
 * mapping.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <cstring>

Position::Position() {
    score_ = secondary_score_ = start_ = end_ = 0;
    complement_ = false;
    cigar_ = 0;
}

Position::Position(uint32_t score, uint32_t secondaryScore, uint32_t start,
                   uint32_t end, bool complemented, const char* cigarStr,
                   uint32_t cigarLen)
        : score_(score),
          secondary_score_(secondaryScore),
          start_(start),
          end_(end),
          complement_(complemented),
          cigar_(0) {

    cigar(cigarStr, cigarLen);
}

bool Position::isComplement() {
    return complement_;
}
void Position::setComplement(bool complemented) {
    complement_ = complemented;
}

uint32_t Position::score() const {
    return score_;
}
uint32_t Position::secondaryScore() const {
    return secondary_score_;
}

uint32_t Position::end() {
    return end_;
}

uint32_t Position::start() const {
    return start_;
}

void Position::end(uint32_t end) {
    end_ = end;
}
void Position::start(uint32_t start) {
    start_ = start;
}

void Position::score(uint32_t score) {
    score_ = score;
}
void Position::secondaryScore(uint32_t sscore) {
    secondary_score_ = sscore;
}

const char* Position::cigar() {
    return cigar_;
}
void Position::cigar(const char* cigar, uint32_t len) {
    if(cigar == NULL) {
        return;
    }

    cigar_ = new char[len + 1];
    memcpy(cigar_, cigar, len);
    cigar_[len] = 0;
}

Position::~Position() {
    if(cigar_) {
        delete[] cigar_;
        cigar_ = 0;
    }
}
/*
 * sequence.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>
#include <vector>

//KSEQ_INIT(int, read)

Sequence::Sequence() {
    data_ = NULL;
    dataLen_ = 0;
    basesAsInt_ = false;
    numOfSequences_ = 0;
}
Sequence::~Sequence() {
    clear();
}

void Sequence::clear() {
    if(data_) {
        delete[] data_;
        data_ = 0;
    }

    info_.clear();
    dataLen_ = 0;
}

const char* Sequence::data() {
    return data_;
}

uint32_t Sequence::dataLen() {
    return dataLen_;
}

uint32_t Sequence::numOfSequences() {
    return numOfSequences_;
}

inline bool Sequence::basesInt() {
    return basesAsInt_;
}

void Sequence::allBasesToSmallInt() {
    if(basesAsInt_) {
        return;
    }

    for (uint32_t i = 0; i < dataLen_; ++i) {
        data_[i] = baseToInt(data_[i]);
    }
    basesAsInt_ = true;
}
void Sequence::allBasesToLetters() {
    if(!basesAsInt_) {
        return;
    }

    for (uint32_t i = 0; i < dataLen_; ++i) {
        data_[i] = intToBase(data_[i]);
    }
    basesAsInt_ = false;
}

uint32_t Sequence::sequenceIndex(uint32_t positionGlobal) {
    uint32_t lo, hi, mid;
    lo = 0;
    hi = numOfSequences_ - 1;

    while (lo < hi) {
        mid = lo + (hi - lo) / 2;
        uint32_t lastPos = seqEndIndex_[mid];

        if(positionGlobal < lastPos) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    return lo;
}

uint32_t Sequence::seqLen(uint32_t index) {
    if(index == 0) {
        return seqEndIndex_[0];
    } else {
        return seqEndIndex_[index] - seqEndIndex_[index - 1];
    }
}

uint32_t Sequence::positionInSeq(uint32_t positionGlobal) {
    uint32_t seqIndex = sequenceIndex(positionGlobal);

    if(seqIndex == 0) {
        return positionGlobal;
    } else {
        return positionGlobal - seqEndIndex_[seqIndex - 1];
    }
    return 0;
}

/*
 * read.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <sstream>

#define LINE_SIZE 30000

Read::Read(float keepRatio, uint32_t maxPositions)
        : keepRatio_(keepRatio),
          maxPositions_(maxPositions) {
    id_ = data_ = optional_identifier_ = quality_ = 0;
    dataLen_ = 0;
    basesAsInt_ = false;

}
Read::~Read() {
    clear();
}

void Read::setData(string data, string id) {
    dataLen_ = data.size();

    data_ = (char*) malloc(dataLen_ * sizeof(char));
    id_ = (char*) malloc((id.size() + 1) * sizeof(char));
    id_[id.size()] = 0;

    memcpy(this->data_, data.c_str(), dataLen_);
    memcpy(this->id_, id.c_str(), id.size());
}

void Read::addPosition(uint32_t score, uint32_t start, uint32_t end,
bool isComplement,
                       const char* cigar, uint32_t cigarLen,
                       uint32_t secondaryScore) {

    Position* p = new Position(score, secondaryScore, start, end, isComplement,
                               cigar, cigarLen);
    addPosition(p);

}

void Read::addPosition(Position* p) {
    positions_.insert(p);
    while (positions_.size() > 1
            && ((((double) (*positions_.rbegin())->score()
                    / (*positions_.begin())->score()) > keepRatio_)
                    || positions_.size() > maxPositions_)) {

        delete *positions_.begin();
        positions_.erase(positions_.begin());
    }
}

Position* Read::bestPosition(uint32_t index) {
    std::multiset<Position*, ptr_compare<Position> >::reverse_iterator it =
            positions_.rbegin();
    uint32_t cntr = 0;

    for (; it != positions_.rend(); ++it, ++cntr) {
        if(cntr == index) {
            return *it;
        }
    }
    return NULL;
}
std::set<Position*, ptr_compare<Position> >& Read::positions() {
    return positions_;
}
inline uint32_t Read::positionsSize() {
    return positions_.size();
}

inline const char* Read::data() {
    return data_;
}
inline const char* Read::id() {
    return id_;
}

inline void Read::keepRatio(float keepRatio) {
    keepRatio_ = keepRatio;
}
inline void Read::maxPositions(uint32_t maxPositions) {
    maxPositions_ = maxPositions;
}

inline uint32_t Read::dataLen() {
    return dataLen_;
}

inline bool Read::basesInt() {
    return basesAsInt_;
}

void Read::clear() {
    if(data_) {
        delete[] data_;
        data_ = 0;
    }
    if(id_) {
        delete[] id_;
        id_ = 0;
    }
    if(optional_identifier_) {
        delete[] optional_identifier_;
        optional_identifier_ = 0;
    }
    if(quality_) {
        delete[] quality_;
        quality_ = 0;
    }

    uint32_t cntr = 0;
    for (auto it = positions_.begin(); it != positions_.end(); ++it, ++cntr) {
        delete *it;
    }
    dataLen_ = 0;

}

void Read::allBasesToSmallInt() {
    if(basesAsInt_) {
        return;
    }

    for (uint32_t i = 0; i < dataLen_; ++i) {
        data_[i] = baseToInt(data_[i]);
    }

    basesAsInt_ = true;
}
void Read::allBasesToLetters() {
    if(!basesAsInt_) {
        return;
    }
    for (uint32_t i = 0; i < dataLen_; ++i) {
        data_[i] = intToBase(data_[i]);
    }
    basesAsInt_ = false;
}

Read* Read::getReverseComplement() {
    Read *rev = new Read();
    rev->dataLen_ = dataLen_;
    rev->data_ = new char[dataLen_ + 1];

    for (uint32_t i = 0; i < dataLen_; ++i) {
        if(basesAsInt_) {
            rev->data_[i] = getACGTComplementAsSmallInt(
                    data_[dataLen_ - 1 - i]);
        } else {
            rev->data_[i] = getACGTComplement(data_[dataLen_ - 1 - i]);
        }
    }
    rev->data_[dataLen_] = 0;
    return rev;

}
void Read::printReadSAM(FILE* outFile, Sequence* seq) {
    Position* best = bestPosition(0);

    if(best != NULL) {
        allBasesToLetters();

        uint32_t seqIndex = seq->sequenceIndex(best->start());
        uint32_t start = seq->positionInSeq(best->start());

    }

}

string Read::getReadMiniSAM(Sequence* seq) {
    Position* best = bestPosition(0);
    char line[LINE_MAX];

    if(best != NULL) {
        return generateMiniSAM(best, 0.99, seq);
    }
    sprintf(line, "%s,%d,%d,%d,%c,%.2f", id_, seq->info_[0], 0, 0, '+', 0.00);
    return string(line);

}

string Read::generateMiniSAM(Position* pos, double score, Sequence* seq) {
    char line[LINE_MAX];
    uint32_t seqIndex = seq->sequenceIndex(pos->start());
    uint32_t start = seq->positionInSeq(pos->start());
    sprintf(line, "%s,%d,%d,%d,%c,", id_, seq->info_[seqIndex], start + 1,
            start + this->dataLen() + 1, pos->isComplement() ? '-' : '+');
    //cerr << "CIG" << string(pos->cigar()) << " " << pos->score() << endl;
    return string(line);

}

vector<string> Read::getReadMiniSAMs(Sequence* seq) {
    vector<string> sams;
    if(this->positionsSize() == 0) {
        //cerr << "Not mapped" << endl;
        char line[LINE_MAX];
        sprintf(line, "%s,%d,%d,%d,%c,", id_, seq->info_[0], 0 + 1,
                0 + this->dataLen() + 1, true ? '-' : '+');
        sams.push_back(string(line));
        return sams;

    }
    int hits = this->positionsSize();

    for (int i = 0; i < 1; ++i) {
        Position* best = bestPosition(i);
        if(best != NULL) {
            sams.push_back(generateMiniSAM(best, 0.0 / hits, seq));
        }

    }
    return sams;

}

/*
 * suffix_array.cpp
 *
 *  Created on: Dec 27, 2014
 *      Author: marko
 */

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SINGLE FILE LIBDIVSUFSORT WRAPPER
/*
 * divsufsort.c for libdivsufsort
 * Copyright (c) 2003-2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef int32_t saint_t;
typedef int32_t saidx_t;
typedef uint8_t sauchar_t;

using namespace std;

saidx_t divbwt(const sauchar_t *T, sauchar_t *U, saidx_t *A, saidx_t n);
saint_t divsufsort(const sauchar_t *T, saidx_t *SA, saidx_t n);

///////////////////////////////////////////////////////////////////////
// header file
#ifndef _DIVSUFSORT_PRIVATE_H
#define _DIVSUFSORT_PRIVATE_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if HAVE_CONFIG_H
# include "config.h"
#endif
#include <assert.h>
#include <stdio.h>
#if HAVE_STRING_H
# include <string.h>
#endif
#if HAVE_STDLIB_H
# include <stdlib.h>
#endif
#if HAVE_MEMORY_H
# include <memory.h>
#endif
#if HAVE_STDDEF_H
# include <stddef.h>
#endif
#if HAVE_STRINGS_H
# include <strings.h>
#endif
#if HAVE_INTTYPES_H
# include <inttypes.h>
#else
# if HAVE_STDINT_H
#  include <stdint.h>
# endif
#endif
#if defined(BUILD_DIVSUFSORT64)
# ifndef SAIDX_T
#  define SAIDX_T
#  define saidx_t saidx64_t
# endif /* SAIDX_T */
# ifndef PRIdSAIDX_T
#  define PRIdSAIDX_T PRIdSAIDX64_T
# endif /* PRIdSAIDX_T */
# define divsufsort divsufsort64
# define divbwt divbwt64
# define divsufsort_version divsufsort64_version
# define bw_transform bw_transform64
# define inverse_bw_transform inverse_bw_transform64
# define sufcheck sufcheck64
# define sa_search sa_search64
# define sa_simplesearch sa_simplesearch64
# define sssort sssort64
# define trsort trsort64
#endif

/*- Constants -*/
#if !defined(UINT8_MAX)
# define UINT8_MAX (255)
#endif /* UINT8_MAX */
#if defined(ALPHABET_SIZE) && (ALPHABET_SIZE < 1)
# undef ALPHABET_SIZE
#endif
#if !defined(ALPHABET_SIZE)
# define ALPHABET_SIZE (UINT8_MAX + 1)
#endif
/* for divsufsort.c */
#define BUCKET_A_SIZE (ALPHABET_SIZE)
#define BUCKET_B_SIZE (ALPHABET_SIZE * ALPHABET_SIZE)
/* for sssort.c */
#if defined(SS_INSERTIONSORT_THRESHOLD)
# if SS_INSERTIONSORT_THRESHOLD < 1
#  undef SS_INSERTIONSORT_THRESHOLD
#  define SS_INSERTIONSORT_THRESHOLD (1)
# endif
#else
# define SS_INSERTIONSORT_THRESHOLD (8)
#endif
#if defined(SS_BLOCKSIZE)
# if SS_BLOCKSIZE < 0
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (0)
# elif 32768 <= SS_BLOCKSIZE
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (32767)
# endif
#else
# define SS_BLOCKSIZE (1024)
#endif
/* minstacksize = log(SS_BLOCKSIZE) / log(3) * 2 */
#if SS_BLOCKSIZE == 0
# if defined(BUILD_DIVSUFSORT64)
#  define SS_MISORT_STACKSIZE (96)
# else
#  define SS_MISORT_STACKSIZE (64)
# endif
#elif SS_BLOCKSIZE <= 4096
# define SS_MISORT_STACKSIZE (16)
#else
# define SS_MISORT_STACKSIZE (24)
#endif
#if defined(BUILD_DIVSUFSORT64)
# define SS_SMERGE_STACKSIZE (64)
#else
# define SS_SMERGE_STACKSIZE (32)
#endif
/* for trsort.c */
#define TR_INSERTIONSORT_THRESHOLD (8)
#if defined(BUILD_DIVSUFSORT64)
# define TR_STACKSIZE (96)
#else
# define TR_STACKSIZE (64)
#endif

/*- Macros -*/
#ifndef SWAP
# define SWAP(_a, _b) do { t = (_a); (_a) = (_b); (_b) = t; } while(0)
#endif /* SWAP */
#ifndef MIN
# define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif /* MIN */
#ifndef MAX
# define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif /* MAX */
#define STACK_PUSH(_a, _b, _c, _d)\
  do {\
    assert(ssize < STACK_SIZE);\
    stack[ssize].a = (_a), stack[ssize].b = (_b),\
    stack[ssize].c = (_c), stack[ssize++].d = (_d);\
  } while(0)
#define STACK_PUSH5(_a, _b, _c, _d, _e)\
  do {\
    assert(ssize < STACK_SIZE);\
    stack[ssize].a = (_a), stack[ssize].b = (_b),\
    stack[ssize].c = (_c), stack[ssize].d = (_d), stack[ssize++].e = (_e);\
  } while(0)
#define STACK_POP(_a, _b, _c, _d)\
  do {\
    assert(0 <= ssize);\
    if(ssize == 0) { return; }\
    (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
    (_c) = stack[ssize].c, (_d) = stack[ssize].d;\
  } while(0)
#define STACK_POP5(_a, _b, _c, _d, _e)\
  do {\
    assert(0 <= ssize);\
    if(ssize == 0) { return; }\
    (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
    (_c) = stack[ssize].c, (_d) = stack[ssize].d, (_e) = stack[ssize].e;\
  } while(0)
/* for divsufsort.c */
#define BUCKET_A(_c0) bucket_A[(_c0)]
#if ALPHABET_SIZE == 256
#define BUCKET_B(_c0, _c1) (bucket_B[((_c1) << 8) | (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[((_c0) << 8) | (_c1)])
#else
#define BUCKET_B(_c0, _c1) (bucket_B[(_c1) * ALPHABET_SIZE + (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[(_c0) * ALPHABET_SIZE + (_c1)])
#endif

/*- Private Prototypes -*/
/* sssort.c */
void
sssort(const sauchar_t *Td, const saidx_t *PA, saidx_t *first, saidx_t *last,
       saidx_t *buf, saidx_t bufsize, saidx_t depth, saidx_t n,
       saint_t lastsuffix);
/* trsort.c */
void
trsort(saidx_t *ISA, saidx_t *SA, saidx_t n, saidx_t depth);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _DIVSUFSORT_PRIVATE_H */

/////////////////////////////////////////////////////////////////////////////////////////
// UTILS

/*- Private Function -*/

/* Binary search for inverse bwt. */
static saidx_t binarysearch_lower(const saidx_t *A, saidx_t size,
                                  saidx_t value) {
    saidx_t half, i;
    for (i = 0, half = size >> 1; 0 < size; size = half, half >>= 1) {
        if(A[i + half] < value) {
            i += half + 1;
            half -= (size & 1) ^ 1;
        }
    }
    return i;
}

/*- Functions -*/

/* Burrows-Wheeler transform. */
saint_t bw_transform(const sauchar_t *T, sauchar_t *U, saidx_t *SA, saidx_t n,
                     saidx_t *idx) {
    saidx_t *A, i, j, p, t;
    saint_t c;

    /* Check arguments. */
    if((T == NULL) || (U == NULL) || (n < 0) || (idx == NULL)) {
        return -1;
    }
    if(n <= 1) {
        if(n == 1) {
            U[0] = T[0];
        }
        *idx = n;
        return 0;
    }

    if((A = SA) == NULL) {
        i = divbwt(T, U, NULL, n);
        if(0 <= i) {
            *idx = i;
            i = 0;
        }
        return (saint_t) i;
    }

    /* BW transform. */
    if(T == U) {
        t = n;
        for (i = 0, j = 0; i < n; ++i) {
            p = t - 1;
            t = A[i];
            if(0 <= p) {
                c = T[j];
                U[j] = (j <= p) ? T[p] : (sauchar_t) A[p];
                A[j] = c;
                j++;
            } else {
                *idx = i;
            }
        }
        p = t - 1;
        if(0 <= p) {
            c = T[j];
            U[j] = (j <= p) ? T[p] : (sauchar_t) A[p];
            A[j] = c;
        } else {
            *idx = i;
        }
    } else {
        U[0] = T[n - 1];
        for (i = 0; A[i] != 0; ++i) {
            U[i + 1] = T[A[i] - 1];
        }
        *idx = i + 1;
        for (++i; i < n; ++i) {
            U[i] = T[A[i] - 1];
        }
    }

    if(SA == NULL) {
        /* Deallocate memory. */
        free(A);
    }

    return 0;
}

/* Inverse Burrows-Wheeler transform. */
saint_t inverse_bw_transform(const sauchar_t *T, sauchar_t *U, saidx_t *A,
                             saidx_t n, saidx_t idx) {
    saidx_t C[ALPHABET_SIZE];
    sauchar_t D[ALPHABET_SIZE];
    saidx_t *B;
    saidx_t i, p;
    saint_t c, d;

    /* Check arguments. */
    if((T == NULL) || (U == NULL) || (n < 0) || (idx < 0) || (n < idx)
            || ((0 < n) && (idx == 0))) {
        return -1;
    }
    if(n <= 1) {
        return 0;
    }

    if((B = A) == NULL) {
        /* Allocate n*sizeof(saidx_t) bytes of memory. */
        if((B = (saidx_t *) malloc((size_t) n * sizeof(saidx_t))) == NULL) {
            return -2;
        }
    }

    /* Inverse BW transform. */
    for (c = 0; c < ALPHABET_SIZE; ++c) {
        C[c] = 0;
    }
    for (i = 0; i < n; ++i) {
        ++C[T[i]];
    }
    for (c = 0, d = 0, i = 0; c < ALPHABET_SIZE; ++c) {
        p = C[c];
        if(0 < p) {
            C[c] = i;
            D[d++] = (sauchar_t) c;
            i += p;
        }
    }
    for (i = 0; i < idx; ++i) {
        B[C[T[i]]++] = i;
    }
    for (; i < n; ++i) {
        B[C[T[i]]++] = i + 1;
    }
    for (c = 0; c < d; ++c) {
        C[c] = C[D[c]];
    }
    for (i = 0, p = idx; i < n; ++i) {
        U[i] = D[binarysearch_lower(C, d, p)];
        p = B[p - 1];
    }

    if(A == NULL) {
        /* Deallocate memory. */
        free(B);
    }

    return 0;
}

/* Checks the suffix array SA of the string T. */
saint_t sufcheck(const sauchar_t *T, const saidx_t *SA, saidx_t n,
                 saint_t verbose) {
    saidx_t C[ALPHABET_SIZE];
    saidx_t i, p, q, t;
    saint_t c;

    if(verbose) {
        fprintf(stderr, "sufcheck: ");
    }

    /* Check arguments. */
    if((T == NULL) || (SA == NULL) || (n < 0)) {
        if(verbose) {
            fprintf(stderr, "Invalid arguments.\n");
        }
        return -1;
    }
    if(n == 0) {
        if(verbose) {
            fprintf(stderr, "Done.\n");
        }
        return 0;
    }

    /* check range: [0..n-1] */
    for (i = 0; i < n; ++i) {
        if((SA[i] < 0) || (n <= SA[i])) {
            return -2;
        }
    }

    /* check first characters. */
    for (i = 1; i < n; ++i) {
        if(T[SA[i - 1]] > T[SA[i]]) {
            return -3;
        }
    }

    /* check suffixes. */
    for (i = 0; i < ALPHABET_SIZE; ++i) {
        C[i] = 0;
    }
    for (i = 0; i < n; ++i) {
        ++C[T[i]];
    }
    for (i = 0, p = 0; i < ALPHABET_SIZE; ++i) {
        t = C[i];
        C[i] = p;
        p += t;
    }

    q = C[T[n - 1]];
    C[T[n - 1]] += 1;
    for (i = 0; i < n; ++i) {
        p = SA[i];
        if(0 < p) {
            c = T[--p];
            t = C[c];
        } else {
            c = T[p = n - 1];
            t = q;
        }
        if((t < 0) || (p != SA[t])) {
            return -4;
        }
        if(t != q) {
            ++C[c];
            if((n <= C[c]) || (T[SA[C[c]]] != c)) {
                C[c] = -1;
            }
        }
    }

    if(1 <= verbose) {
        fprintf(stderr, "Done.\n");
    }
    return 0;
}

static
int _compare(const sauchar_t *T, saidx_t Tsize, const sauchar_t *P,
             saidx_t Psize, saidx_t suf, saidx_t *match) {
    saidx_t i, j;
    saint_t r;
    for (i = suf + *match, j = *match, r = 0;
            (i < Tsize) && (j < Psize) && ((r = T[i] - P[j]) == 0); ++i, ++j) {
    }
    *match = j;
    return (r == 0) ? -(j != Psize) : r;
}

/* Search for the pattern P in the string T. */
saidx_t sa_search(const sauchar_t *T, saidx_t Tsize, const sauchar_t *P,
                  saidx_t Psize, const saidx_t *SA, saidx_t SAsize,
                  saidx_t *idx) {
    saidx_t size, lsize, rsize, half;
    saidx_t match, lmatch, rmatch;
    saidx_t llmatch, lrmatch, rlmatch, rrmatch;
    saidx_t i, j, k;
    saint_t r;

    if(idx != NULL) {
        *idx = -1;
    }
    if((T == NULL) || (P == NULL) || (SA == NULL) || (Tsize < 0) || (Psize < 0)
            || (SAsize < 0)) {
        return -1;
    }
    if((Tsize == 0) || (SAsize == 0)) {
        return 0;
    }
    if(Psize == 0) {
        if(idx != NULL) {
            *idx = 0;
        }
        return SAsize;
    }

    for (i = j = k = 0, lmatch = rmatch = 0, size = SAsize, half = size >> 1;
            0 < size; size = half, half >>= 1) {
        match = MIN(lmatch, rmatch);
        r = _compare(T, Tsize, P, Psize, SA[i + half], &match);
        if(r < 0) {
            i += half + 1;
            half -= (size & 1) ^ 1;
            lmatch = match;
        } else if(r > 0) {
            rmatch = match;
        } else {
            lsize = half, j = i, rsize = size - half - 1, k = i + half + 1;

            /* left part */
            for (llmatch = lmatch, lrmatch = match, half = lsize >> 1;
                    0 < lsize; lsize = half, half >>= 1) {
                lmatch = MIN(llmatch, lrmatch);
                r = _compare(T, Tsize, P, Psize, SA[j + half], &lmatch);
                if(r < 0) {
                    j += half + 1;
                    half -= (lsize & 1) ^ 1;
                    llmatch = lmatch;
                } else {
                    lrmatch = lmatch;
                }
            }

            /* right part */
            for (rlmatch = match, rrmatch = rmatch, half = rsize >> 1;
                    0 < rsize; rsize = half, half >>= 1) {
                rmatch = MIN(rlmatch, rrmatch);
                r = _compare(T, Tsize, P, Psize, SA[k + half], &rmatch);
                if(r <= 0) {
                    k += half + 1;
                    half -= (rsize & 1) ^ 1;
                    rlmatch = rmatch;
                } else {
                    rrmatch = rmatch;
                }
            }

            break;
        }
    }

    if(idx != NULL) {
        *idx = (0 < (k - j)) ? j : i;
    }
    return k - j;
}

/* Search for the character c in the string T. */
saidx_t sa_simplesearch(const sauchar_t *T, saidx_t Tsize, const saidx_t *SA,
                        saidx_t SAsize, saint_t c, saidx_t *idx) {
    saidx_t size, lsize, rsize, half;
    saidx_t i, j, k, p;
    saint_t r;

    if(idx != NULL) {
        *idx = -1;
    }
    if((T == NULL) || (SA == NULL) || (Tsize < 0) || (SAsize < 0)) {
        return -1;
    }
    if((Tsize == 0) || (SAsize == 0)) {
        return 0;
    }

    for (i = j = k = 0, size = SAsize, half = size >> 1; 0 < size;
            size = half, half >>= 1) {
        p = SA[i + half];
        r = (p < Tsize) ? T[p] - c : -1;
        if(r < 0) {
            i += half + 1;
            half -= (size & 1) ^ 1;
        } else if(r == 0) {
            lsize = half, j = i, rsize = size - half - 1, k = i + half + 1;

            /* left part */
            for (half = lsize >> 1; 0 < lsize; lsize = half, half >>= 1) {
                p = SA[j + half];
                r = (p < Tsize) ? T[p] - c : -1;
                if(r < 0) {
                    j += half + 1;
                    half -= (lsize & 1) ^ 1;
                }
            }

            /* right part */
            for (half = rsize >> 1; 0 < rsize; rsize = half, half >>= 1) {
                p = SA[k + half];
                r = (p < Tsize) ? T[p] - c : -1;
                if(r <= 0) {
                    k += half + 1;
                    half -= (rsize & 1) ^ 1;
                }
            }

            break;
        }
    }

    if(idx != NULL) {
        *idx = (0 < (k - j)) ? j : i;
    }
    return k - j;
}

//////////////////////////////////////////////////////////////////////////////////////////
// TRSORT

/*- Private Functions -*/
static const saint_t lg_table[256] = { -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3,
        3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7 };

static inline saint_t tr_ilg(saidx_t n) {
#if defined(BUILD_DIVSUFSORT64)
    return (n >> 32) ?
    ((n >> 48) ?
            ((n >> 56) ?
                    56 + lg_table[(n >> 56) & 0xff] :
                    48 + lg_table[(n >> 48) & 0xff]) :
            ((n >> 40) ?
                    40 + lg_table[(n >> 40) & 0xff] :
                    32 + lg_table[(n >> 32) & 0xff])) :
    ((n & 0xffff0000) ?
            ((n & 0xff000000) ?
                    24 + lg_table[(n >> 24) & 0xff] :
                    16 + lg_table[(n >> 16) & 0xff]) :
            ((n & 0x0000ff00) ?
                    8 + lg_table[(n >> 8) & 0xff] :
                    0 + lg_table[(n >> 0) & 0xff]));
#else
    return (n & 0xffff0000) ?
            ((n & 0xff000000) ?
                    24 + lg_table[(n >> 24) & 0xff] :
                    16 + lg_table[(n >> 16) & 0xff]) :
            ((n & 0x0000ff00) ?
                    8 + lg_table[(n >> 8) & 0xff] :
                    0 + lg_table[(n >> 0) & 0xff]);
#endif
}

/*---------------------------------------------------------------------------*/

/* Simple insertionsort for small size groups. */
static
void tr_insertionsort(const saidx_t *ISAd, saidx_t *first, saidx_t *last) {
    saidx_t *a, *b;
    saidx_t t, r;

    for (a = first + 1; a < last; ++a) {
        for (t = *a, b = a - 1; 0 > (r = ISAd[t] - ISAd[*b]);) {
            do {
                *(b + 1) = *b;
            } while ((first <= --b) && (*b < 0));
            if(b < first) {
                break;
            }
        }
        if(r == 0) {
            *b = ~*b;
        }
        *(b + 1) = t;
    }
}

/*---------------------------------------------------------------------------*/

static inline
void tr_fixdown(const saidx_t *ISAd, saidx_t *SA, saidx_t i, saidx_t size) {
    saidx_t j, k;
    saidx_t v;
    saidx_t c, d, e;

    for (v = SA[i], c = ISAd[v]; (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
        d = ISAd[SA[k = j++]];
        if(d < (e = ISAd[SA[j]])) {
            k = j;
            d = e;
        }
        if(d <= c) {
            break;
        }
    }
    SA[i] = v;
}

/* Simple top-down heapsort. */
static
void tr_heapsort(const saidx_t *ISAd, saidx_t *SA, saidx_t size) {
    saidx_t i, m;
    saidx_t t;

    m = size;
    if((size % 2) == 0) {
        m--;
        if(ISAd[SA[m / 2]] < ISAd[SA[m]]) {
            SWAP(SA[m], SA[m / 2]);
        }
    }

    for (i = m / 2 - 1; 0 <= i; --i) {
        tr_fixdown(ISAd, SA, i, m);
    }
    if((size % 2) == 0) {
        SWAP(SA[0], SA[m]);
        tr_fixdown(ISAd, SA, 0, m);
    }
    for (i = m - 1; 0 < i; --i) {
        t = SA[0], SA[0] = SA[i];
        tr_fixdown(ISAd, SA, 0, i);
        SA[i] = t;
    }
}

/*---------------------------------------------------------------------------*/

/* Returns the median of three elements. */
static inline saidx_t *
tr_median3(const saidx_t *ISAd, saidx_t *v1, saidx_t *v2, saidx_t *v3) {
    saidx_t *t;
    if(ISAd[*v1] > ISAd[*v2]) {
        SWAP(v1, v2);
    }
    if(ISAd[*v2] > ISAd[*v3]) {
        if(ISAd[*v1] > ISAd[*v3]) {
            return v1;
        } else {
            return v3;
        }
    }
    return v2;
}

/* Returns the median of five elements. */
static inline saidx_t *
tr_median5(const saidx_t *ISAd, saidx_t *v1, saidx_t *v2, saidx_t *v3,
           saidx_t *v4, saidx_t *v5) {
    saidx_t *t;
    if(ISAd[*v2] > ISAd[*v3]) {
        SWAP(v2, v3);
    }
    if(ISAd[*v4] > ISAd[*v5]) {
        SWAP(v4, v5);
    }
    if(ISAd[*v2] > ISAd[*v4]) {
        SWAP(v2, v4);SWAP(v3, v5);
    }
    if(ISAd[*v1] > ISAd[*v3]) {
        SWAP(v1, v3);
    }
    if(ISAd[*v1] > ISAd[*v4]) {
        SWAP(v1, v4);SWAP(v3, v5);
    }
    if(ISAd[*v3] > ISAd[*v4]) {
        return v4;
    }
    return v3;
}

/* Returns the pivot element. */
static inline saidx_t *
tr_pivot(const saidx_t *ISAd, saidx_t *first, saidx_t *last) {
    saidx_t *middle;
    saidx_t t;

    t = last - first;
    middle = first + t / 2;

    if(t <= 512) {
        if(t <= 32) {
            return tr_median3(ISAd, first, middle, last - 1);
        } else {
            t >>= 2;
            return tr_median5(ISAd, first, first + t, middle, last - 1 - t,
                              last - 1);
        }
    }
    t >>= 3;
    first = tr_median3(ISAd, first, first + t, first + (t << 1));
    middle = tr_median3(ISAd, middle - t, middle, middle + t);
    last = tr_median3(ISAd, last - 1 - (t << 1), last - 1 - t, last - 1);
    return tr_median3(ISAd, first, middle, last);
}

/*---------------------------------------------------------------------------*/

typedef struct _trbudget_t trbudget_t;
struct _trbudget_t {
    saidx_t chance;
    saidx_t remain;
    saidx_t incval;
    saidx_t count;
};

static inline
void trbudget_init(trbudget_t *budget, saidx_t chance, saidx_t incval) {
    budget->chance = chance;
    budget->remain = budget->incval = incval;
}

static inline saint_t trbudget_check(trbudget_t *budget, saidx_t size) {
    if(size <= budget->remain) {
        budget->remain -= size;
        return 1;
    }
    if(budget->chance == 0) {
        budget->count += size;
        return 0;
    }
    budget->remain += budget->incval - size;
    budget->chance -= 1;
    return 1;
}

/*---------------------------------------------------------------------------*/

static inline
void tr_partition(const saidx_t *ISAd, saidx_t *first, saidx_t *middle,
                  saidx_t *last, saidx_t **pa, saidx_t **pb, saidx_t v) {
    saidx_t *a, *b, *c, *d, *e, *f;
    saidx_t t, s;
    saidx_t x = 0;

    for (b = middle - 1; (++b < last) && ((x = ISAd[*b]) == v);) {
    }
    if(((a = b) < last) && (x < v)) {
        for (; (++b < last) && ((x = ISAd[*b]) <= v);) {
            if(x == v) {
                SWAP(*b, *a);
                ++a;
            }
        }
    }
    for (c = last; (b < --c) && ((x = ISAd[*c]) == v);) {
    }
    if((b < (d = c)) && (x > v)) {
        for (; (b < --c) && ((x = ISAd[*c]) >= v);) {
            if(x == v) {
                SWAP(*c, *d);
                --d;
            }
        }
    }
    for (; b < c;) {
        SWAP(*b, *c);
        for (; (++b < c) && ((x = ISAd[*b]) <= v);) {
            if(x == v) {
                SWAP(*b, *a);
                ++a;
            }
        }
        for (; (b < --c) && ((x = ISAd[*c]) >= v);) {
            if(x == v) {
                SWAP(*c, *d);
                --d;
            }
        }
    }

    if(a <= d) {
        c = b - 1;
        if((s = a - first) > (t = b - a)) {
            s = t;
        }
        for (e = first, f = b - s; 0 < s; --s, ++e, ++f) {
            SWAP(*e, *f);
        }
        if((s = d - c) > (t = last - d - 1)) {
            s = t;
        }
        for (e = b, f = last - s; 0 < s; --s, ++e, ++f) {
            SWAP(*e, *f);
        }
        first += (b - a), last -= (d - c);
    }
    *pa = first, *pb = last;
}

static
void tr_copy(saidx_t *ISA, const saidx_t *SA, saidx_t *first, saidx_t *a,
             saidx_t *b, saidx_t *last, saidx_t depth) {
    /* sort suffixes of middle partition
     by using sorted order of suffixes of left and right partition. */
    saidx_t *c, *d, *e;
    saidx_t s, v;

    v = b - SA - 1;
    for (c = first, d = a - 1; c <= d; ++c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *++d = s;
            ISA[s] = d - SA;
        }
    }
    for (c = last - 1, e = d + 1, d = b; e < d; --c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *--d = s;
            ISA[s] = d - SA;
        }
    }
}

static
void tr_partialcopy(saidx_t *ISA, const saidx_t *SA, saidx_t *first, saidx_t *a,
                    saidx_t *b, saidx_t *last, saidx_t depth) {
    saidx_t *c, *d, *e;
    saidx_t s, v;
    saidx_t rank, lastrank, newrank = -1;

    v = b - SA - 1;
    lastrank = -1;
    for (c = first, d = a - 1; c <= d; ++c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *++d = s;
            rank = ISA[s + depth];
            if(lastrank != rank) {
                lastrank = rank;
                newrank = d - SA;
            }
            ISA[s] = newrank;
        }
    }

    lastrank = -1;
    for (e = d; first <= e; --e) {
        rank = ISA[*e];
        if(lastrank != rank) {
            lastrank = rank;
            newrank = e - SA;
        }
        if(newrank != rank) {
            ISA[*e] = newrank;
        }
    }

    lastrank = -1;
    for (c = last - 1, e = d + 1, d = b; e < d; --c) {
        if((0 <= (s = *c - depth)) && (ISA[s] == v)) {
            *--d = s;
            rank = ISA[s + depth];
            if(lastrank != rank) {
                lastrank = rank;
                newrank = d - SA;
            }
            ISA[s] = newrank;
        }
    }
}

static
void tr_introsort(saidx_t *ISA, const saidx_t *ISAd, saidx_t *SA,
                  saidx_t *first, saidx_t *last, trbudget_t *budget) {
#define STACK_SIZE TR_STACKSIZE
    struct {
        const saidx_t *a;
        saidx_t *b, *c;
        saint_t d, e;
    } stack[STACK_SIZE];
    saidx_t *a, *b, *c;
    saidx_t t;
    saidx_t v, x = 0;
    saidx_t incr = ISAd - ISA;
    saint_t limit, next;
    saint_t ssize, trlink = -1;

    for (ssize = 0, limit = tr_ilg(last - first);;) {

        if(limit < 0) {
            if(limit == -1) {
                /* tandem repeat partition */
                tr_partition(ISAd - incr, first, first, last, &a, &b,
                             last - SA - 1);

                /* update ranks */
                if(a < last) {
                    for (c = first, v = a - SA - 1; c < a; ++c) {
                        ISA[*c] = v;
                    }
                }
                if(b < last) {
                    for (c = a, v = b - SA - 1; c < b; ++c) {
                        ISA[*c] = v;
                    }
                }

                /* push */
                if(1 < (b - a)) {
                    STACK_PUSH5(NULL, a, b, 0, 0);
                    STACK_PUSH5(ISAd - incr, first, last, -2, trlink);
                    trlink = ssize - 2;
                }
                if((a - first) <= (last - b)) {
                    if(1 < (a - first)) {
                        STACK_PUSH5(ISAd, b, last, tr_ilg(last - b), trlink);
                        last = a, limit = tr_ilg(a - first);
                    } else if(1 < (last - b)) {
                        first = b, limit = tr_ilg(last - b);
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                } else {
                    if(1 < (last - b)) {
                        STACK_PUSH5(ISAd, first, a, tr_ilg(a - first), trlink);
                        first = b, limit = tr_ilg(last - b);
                    } else if(1 < (a - first)) {
                        last = a, limit = tr_ilg(a - first);
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                }
            } else if(limit == -2) {
                /* tandem repeat copy */
                a = stack[--ssize].b, b = stack[ssize].c;
                if(stack[ssize].d == 0) {
                    tr_copy(ISA, SA, first, a, b, last, ISAd - ISA);
                } else {
                    if(0 <= trlink) {
                        stack[trlink].d = -1;
                    }
                    tr_partialcopy(ISA, SA, first, a, b, last, ISAd - ISA);
                }
                STACK_POP5(ISAd, first, last, limit, trlink);
            } else {
                /* sorted partition */
                if(0 <= *first) {
                    a = first;
                    do {
                        ISA[*a] = a - SA;
                    } while ((++a < last) && (0 <= *a));
                    first = a;
                }
                if(first < last) {
                    a = first;
                    do {
                        *a = ~*a;
                    } while (*++a < 0);
                    next = (ISA[*a] != ISAd[*a]) ? tr_ilg(a - first + 1) : -1;
                    if(++a < last) {
                        for (b = first, v = a - SA - 1; b < a; ++b) {
                            ISA[*b] = v;
                        }
                    }

                    /* push */
                    if(trbudget_check(budget, a - first)) {
                        if((a - first) <= (last - a)) {
                            STACK_PUSH5(ISAd, a, last, -3, trlink);
                            ISAd += incr, last = a, limit = next;
                        } else {
                            if(1 < (last - a)) {
                                STACK_PUSH5(ISAd + incr, first, a, next,
                                            trlink);
                                first = a, limit = -3;
                            } else {
                                ISAd += incr, last = a, limit = next;
                            }
                        }
                    } else {
                        if(0 <= trlink) {
                            stack[trlink].d = -1;
                        }
                        if(1 < (last - a)) {
                            first = a, limit = -3;
                        } else {
                            STACK_POP5(ISAd, first, last, limit, trlink);
                        }
                    }
                } else {
                    STACK_POP5(ISAd, first, last, limit, trlink);
                }
            }
            continue;
        }

        if((last - first) <= TR_INSERTIONSORT_THRESHOLD) {
            tr_insertionsort(ISAd, first, last);
            limit = -3;
            continue;
        }

        if(limit-- == 0) {
            tr_heapsort(ISAd, first, last - first);
            for (a = last - 1; first < a; a = b) {
                for (x = ISAd[*a], b = a - 1; (first <= b) && (ISAd[*b] == x);
                        --b) {
                    *b = ~*b;
                }
            }
            limit = -3;
            continue;
        }

        /* choose pivot */
        a = tr_pivot(ISAd, first, last);
        SWAP(*first, *a);
        v = ISAd[*first];

        /* partition */
        tr_partition(ISAd, first, first + 1, last, &a, &b, v);
        if((last - first) != (b - a)) {
            next = (ISA[*a] != v) ? tr_ilg(b - a) : -1;

            /* update ranks */
            for (c = first, v = a - SA - 1; c < a; ++c) {
                ISA[*c] = v;
            }
            if(b < last) {
                for (c = a, v = b - SA - 1; c < b; ++c) {
                    ISA[*c] = v;
                }
            }

            /* push */
            if((1 < (b - a)) && (trbudget_check(budget, b - a))) {
                if((a - first) <= (last - b)) {
                    if((last - b) <= (b - a)) {
                        if(1 < (a - first)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            STACK_PUSH5(ISAd, b, last, limit, trlink);
                            last = a;
                        } else if(1 < (last - b)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            first = b;
                        } else {
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else if((a - first) <= (b - a)) {
                        if(1 < (a - first)) {
                            STACK_PUSH5(ISAd, b, last, limit, trlink);
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            last = a;
                        } else {
                            STACK_PUSH5(ISAd, b, last, limit, trlink);
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else {
                        STACK_PUSH5(ISAd, b, last, limit, trlink);
                        STACK_PUSH5(ISAd, first, a, limit, trlink);
                        ISAd += incr, first = a, last = b, limit = next;
                    }
                } else {
                    if((a - first) <= (b - a)) {
                        if(1 < (last - b)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            STACK_PUSH5(ISAd, first, a, limit, trlink);
                            first = b;
                        } else if(1 < (a - first)) {
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            last = a;
                        } else {
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else if((last - b) <= (b - a)) {
                        if(1 < (last - b)) {
                            STACK_PUSH5(ISAd, first, a, limit, trlink);
                            STACK_PUSH5(ISAd + incr, a, b, next, trlink);
                            first = b;
                        } else {
                            STACK_PUSH5(ISAd, first, a, limit, trlink);
                            ISAd += incr, first = a, last = b, limit = next;
                        }
                    } else {
                        STACK_PUSH5(ISAd, first, a, limit, trlink);
                        STACK_PUSH5(ISAd, b, last, limit, trlink);
                        ISAd += incr, first = a, last = b, limit = next;
                    }
                }
            } else {
                if((1 < (b - a)) && (0 <= trlink)) {
                    stack[trlink].d = -1;
                }
                if((a - first) <= (last - b)) {
                    if(1 < (a - first)) {
                        STACK_PUSH5(ISAd, b, last, limit, trlink);
                        last = a;
                    } else if(1 < (last - b)) {
                        first = b;
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                } else {
                    if(1 < (last - b)) {
                        STACK_PUSH5(ISAd, first, a, limit, trlink);
                        first = b;
                    } else if(1 < (a - first)) {
                        last = a;
                    } else {
                        STACK_POP5(ISAd, first, last, limit, trlink);
                    }
                }
            }
        } else {
            if(trbudget_check(budget, last - first)) {
                limit = tr_ilg(last - first), ISAd += incr;
            } else {
                if(0 <= trlink) {
                    stack[trlink].d = -1;
                }
                STACK_POP5(ISAd, first, last, limit, trlink);
            }
        }
    }
#undef STACK_SIZE
}

/*---------------------------------------------------------------------------*/

/*- Function -*/

/* Tandem repeat sort */
void trsort(saidx_t *ISA, saidx_t *SA, saidx_t n, saidx_t depth) {
    saidx_t *ISAd;
    saidx_t *first, *last;
    trbudget_t budget;
    saidx_t t, skip, unsorted;

    trbudget_init(&budget, tr_ilg(n) * 2 / 3, n);
    /*  trbudget_init(&budget, tr_ilg(n) * 3 / 4, n); */
    for (ISAd = ISA + depth; -n < *SA; ISAd += ISAd - ISA) {
        first = SA;
        skip = 0;
        unsorted = 0;
        do {
            if((t = *first) < 0) {
                first -= t;
                skip += t;
            } else {
                if(skip != 0) {
                    *(first + skip) = skip;
                    skip = 0;
                }
                last = SA + ISA[t] + 1;
                if(1 < (last - first)) {
                    budget.count = 0;
                    tr_introsort(ISA, ISAd, SA, first, last, &budget);
                    if(budget.count != 0) {
                        unsorted += budget.count;
                    } else {
                        skip = first - last;
                    }
                } else if((last - first) == 1) {
                    skip = -1;
                }
                first = last;
            }
        } while (first < (SA + n));
        if(skip != 0) {
            *(first + skip) = skip;
        }
        if(unsorted == 0) {
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// sssort

/*- Private Functions -*/

#if (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE)

static inline saint_t ss_ilg(saidx_t n) {
#if SS_BLOCKSIZE == 0
# if defined(BUILD_DIVSUFSORT64)
    return (n >> 32) ?
    ((n >> 48) ?
            ((n >> 56) ?
                    56 + lg_table[(n >> 56) & 0xff] :
                    48 + lg_table[(n >> 48) & 0xff]) :
            ((n >> 40) ?
                    40 + lg_table[(n >> 40) & 0xff] :
                    32 + lg_table[(n >> 32) & 0xff])) :
    ((n & 0xffff0000) ?
            ((n & 0xff000000) ?
                    24 + lg_table[(n >> 24) & 0xff] :
                    16 + lg_table[(n >> 16) & 0xff]) :
            ((n & 0x0000ff00) ?
                    8 + lg_table[(n >> 8) & 0xff] :
                    0 + lg_table[(n >> 0) & 0xff]));
# else
    return (n & 0xffff0000) ?
    ((n & 0xff000000) ?
            24 + lg_table[(n >> 24) & 0xff] :
            16 + lg_table[(n >> 16) & 0xff]) :
    ((n & 0x0000ff00) ?
            8 + lg_table[(n >> 8) & 0xff] :
            0 + lg_table[(n >> 0) & 0xff]);
# endif
#elif SS_BLOCKSIZE < 256
    return lg_table[n];
#else
    return (n & 0xff00) ?
            8 + lg_table[(n >> 8) & 0xff] : 0 + lg_table[(n >> 0) & 0xff];
#endif
}

#endif /* (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE) */

#if SS_BLOCKSIZE != 0

static const saint_t sqq_table[256] = { 0, 16, 22, 27, 32, 35, 39, 42, 45, 48,
        50, 53, 55, 57, 59, 61, 64, 65, 67, 69, 71, 73, 75, 76, 78, 80, 81, 83,
        84, 86, 87, 89, 90, 91, 93, 94, 96, 97, 98, 99, 101, 102, 103, 104, 106,
        107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
        122, 123, 124, 125, 126, 128, 128, 129, 130, 131, 132, 133, 134, 135,
        136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 145, 146, 147, 148,
        149, 150, 150, 151, 152, 153, 154, 155, 155, 156, 157, 158, 159, 160,
        160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168, 169, 170, 170,
        171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178, 179, 180, 181,
        181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188, 189, 189, 190,
        191, 192, 192, 193, 193, 194, 195, 195, 196, 197, 197, 198, 199, 199,
        200, 201, 201, 202, 203, 203, 204, 204, 205, 206, 206, 207, 208, 208,
        209, 209, 210, 211, 211, 212, 212, 213, 214, 214, 215, 215, 216, 217,
        217, 218, 218, 219, 219, 220, 221, 221, 222, 222, 223, 224, 224, 225,
        225, 226, 226, 227, 227, 228, 229, 229, 230, 230, 231, 231, 232, 232,
        233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238, 239, 240, 240,
        241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246, 246, 247, 247,
        248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253, 253, 254, 254,
        255 };

static inline saidx_t ss_isqrt(saidx_t x) {

    saidx_t y, e;

    if(x >= (SS_BLOCKSIZE * SS_BLOCKSIZE)) {
        return SS_BLOCKSIZE;
    }
    e = (x & 0xffff0000) ?
            ((x & 0xff000000) ?
                    24 + lg_table[(x >> 24) & 0xff] :
                    16 + lg_table[(x >> 16) & 0xff]) :
            ((x & 0x0000ff00) ?
                    8 + lg_table[(x >> 8) & 0xff] :
                    0 + lg_table[(x >> 0) & 0xff]);

    if(e >= 16) {
        y = sqq_table[x >> ((e - 6) - (e & 1))] << ((e >> 1) - 7);
        if(e >= 24) {
            y = (y + 1 + x / y) >> 1;
        }
        y = (y + 1 + x / y) >> 1;
    } else if(e >= 8) {
        y = (sqq_table[x >> ((e - 6) - (e & 1))] >> (7 - (e >> 1))) + 1;
    } else {
        return sqq_table[x] >> 4;
    }

    return (x < (y * y)) ? y - 1 : y;
}

#endif /* SS_BLOCKSIZE != 0 */

/*---------------------------------------------------------------------------*/

/* Compares two suffixes. */
static inline saint_t ss_compare(const sauchar_t *T, const saidx_t *p1,
                                 const saidx_t *p2, saidx_t depth) {
    const sauchar_t *U1, *U2, *U1n, *U2n;

    for (U1 = T + depth + *p1, U2 = T + depth + *p2, U1n = T + *(p1 + 1) + 2, U2n =
            T + *(p2 + 1) + 2; (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
            ++U1, ++U2) {
    }

    return U1 < U1n ? (U2 < U2n ? *U1 - *U2 : 1) : (U2 < U2n ? -1 : 0);
}

/*---------------------------------------------------------------------------*/

#if (SS_BLOCKSIZE != 1) && (SS_INSERTIONSORT_THRESHOLD != 1)

/* Insertionsort for small size groups */
static
void ss_insertionsort(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
                      saidx_t *last, saidx_t depth) {
    saidx_t *i, *j;
    saidx_t t;
    saint_t r;

    for (i = last - 2; first <= i; --i) {
        for (t = *i, j = i + 1; 0 < (r = ss_compare(T, PA + t, PA + *j, depth));
                ) {
            do {
                *(j - 1) = *j;
            } while ((++j < last) && (*j < 0));
            if(last <= j) {
                break;
            }
        }
        if(r == 0) {
            *j = ~*j;
        }
        *(j - 1) = t;
    }
}

#endif /* (SS_BLOCKSIZE != 1) && (SS_INSERTIONSORT_THRESHOLD != 1) */

/*---------------------------------------------------------------------------*/

#if (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE)

static inline
void ss_fixdown(const sauchar_t *Td, const saidx_t *PA, saidx_t *SA, saidx_t i,
                saidx_t size) {
    saidx_t j, k;
    saidx_t v;
    saint_t c, d, e;

    for (v = SA[i], c = Td[PA[v]]; (j = 2 * i + 1) < size;
            SA[i] = SA[k], i = k) {
        d = Td[PA[SA[k = j++]]];
        if(d < (e = Td[PA[SA[j]]])) {
            k = j;
            d = e;
        }
        if(d <= c) {
            break;
        }
    }
    SA[i] = v;
}

/* Simple top-down heapsort. */
static
void ss_heapsort(const sauchar_t *Td, const saidx_t *PA, saidx_t *SA,
                 saidx_t size) {
    saidx_t i, m;
    saidx_t t;

    m = size;
    if((size % 2) == 0) {
        m--;
        if(Td[PA[SA[m / 2]]] < Td[PA[SA[m]]]) {
            SWAP(SA[m], SA[m / 2]);
        }
    }

    for (i = m / 2 - 1; 0 <= i; --i) {
        ss_fixdown(Td, PA, SA, i, m);
    }
    if((size % 2) == 0) {
        SWAP(SA[0], SA[m]);
        ss_fixdown(Td, PA, SA, 0, m);
    }
    for (i = m - 1; 0 < i; --i) {
        t = SA[0], SA[0] = SA[i];
        ss_fixdown(Td, PA, SA, 0, i);
        SA[i] = t;
    }
}

/*---------------------------------------------------------------------------*/

/* Returns the median of three elements. */
static inline saidx_t *
ss_median3(const sauchar_t *Td, const saidx_t *PA, saidx_t *v1, saidx_t *v2,
           saidx_t *v3) {
    saidx_t *t;
    if(Td[PA[*v1]] > Td[PA[*v2]]) {
        SWAP(v1, v2);
    }
    if(Td[PA[*v2]] > Td[PA[*v3]]) {
        if(Td[PA[*v1]] > Td[PA[*v3]]) {
            return v1;
        } else {
            return v3;
        }
    }
    return v2;
}

/* Returns the median of five elements. */
static inline saidx_t *
ss_median5(const sauchar_t *Td, const saidx_t *PA, saidx_t *v1, saidx_t *v2,
           saidx_t *v3, saidx_t *v4, saidx_t *v5) {
    saidx_t *t;
    if(Td[PA[*v2]] > Td[PA[*v3]]) {
        SWAP(v2, v3);
    }
    if(Td[PA[*v4]] > Td[PA[*v5]]) {
        SWAP(v4, v5);
    }
    if(Td[PA[*v2]] > Td[PA[*v4]]) {
        SWAP(v2, v4);SWAP(v3, v5);
    }
    if(Td[PA[*v1]] > Td[PA[*v3]]) {
        SWAP(v1, v3);
    }
    if(Td[PA[*v1]] > Td[PA[*v4]]) {
        SWAP(v1, v4);SWAP(v3, v5);
    }
    if(Td[PA[*v3]] > Td[PA[*v4]]) {
        return v4;
    }
    return v3;
}

/* Returns the pivot element. */
static inline saidx_t *
ss_pivot(const sauchar_t *Td, const saidx_t *PA, saidx_t *first,
         saidx_t *last) {
    saidx_t *middle;
    saidx_t t;

    t = last - first;
    middle = first + t / 2;

    if(t <= 512) {
        if(t <= 32) {
            return ss_median3(Td, PA, first, middle, last - 1);
        } else {
            t >>= 2;
            return ss_median5(Td, PA, first, first + t, middle, last - 1 - t,
                              last - 1);
        }
    }
    t >>= 3;
    first = ss_median3(Td, PA, first, first + t, first + (t << 1));
    middle = ss_median3(Td, PA, middle - t, middle, middle + t);
    last = ss_median3(Td, PA, last - 1 - (t << 1), last - 1 - t, last - 1);
    return ss_median3(Td, PA, first, middle, last);
}

/*---------------------------------------------------------------------------*/

/* Binary partition for substrings. */
static inline saidx_t *
ss_partition(const saidx_t *PA, saidx_t *first, saidx_t *last, saidx_t depth) {
    saidx_t *a, *b;
    saidx_t t;
    for (a = first - 1, b = last;;) {
        for (; (++a < b) && ((PA[*a] + depth) >= (PA[*a + 1] + 1));) {
            *a = ~*a;
        }
        for (; (a < --b) && ((PA[*b] + depth) < (PA[*b + 1] + 1));) {
        }
        if(b <= a) {
            break;
        }
        t = ~*b;
        *b = *a;
        *a = t;
    }
    if(first < a) {
        *first = ~*first;
    }
    return a;
}

/* Multikey introsort for medium size groups. */
static
void ss_mintrosort(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
                   saidx_t *last, saidx_t depth) {
#define STACK_SIZE SS_MISORT_STACKSIZE
    struct {
        saidx_t *a, *b, c;
        saint_t d;
    } stack[STACK_SIZE];
    const sauchar_t *Td;
    saidx_t *a, *b, *c, *d, *e, *f;
    saidx_t s, t;
    saint_t ssize;
    saint_t limit;
    saint_t v, x = 0;

    for (ssize = 0, limit = ss_ilg(last - first);;) {

        if((last - first) <= SS_INSERTIONSORT_THRESHOLD) {
#if 1 < SS_INSERTIONSORT_THRESHOLD
            if(1 < (last - first)) {
                ss_insertionsort(T, PA, first, last, depth);
            }
#endif
            STACK_POP(first, last, depth, limit);
            continue;
        }

        Td = T + depth;
        if(limit-- == 0) {
            ss_heapsort(Td, PA, first, last - first);
        }
        if(limit < 0) {
            for (a = first + 1, v = Td[PA[*first]]; a < last; ++a) {
                if((x = Td[PA[*a]]) != v) {
                    if(1 < (a - first)) {
                        break;
                    }
                    v = x;
                    first = a;
                }
            }
            if(Td[PA[*first] - 1] < v) {
                first = ss_partition(PA, first, a, depth);
            }
            if((a - first) <= (last - a)) {
                if(1 < (a - first)) {
                    STACK_PUSH(a, last, depth, -1);
                    last = a, depth += 1, limit = ss_ilg(a - first);
                } else {
                    first = a, limit = -1;
                }
            } else {
                if(1 < (last - a)) {
                    STACK_PUSH(first, a, depth + 1, ss_ilg(a - first));
                    first = a, limit = -1;
                } else {
                    last = a, depth += 1, limit = ss_ilg(a - first);
                }
            }
            continue;
        }

        /* choose pivot */
        a = ss_pivot(Td, PA, first, last);
        v = Td[PA[*a]];
        SWAP(*first, *a);

        /* partition */
        for (b = first; (++b < last) && ((x = Td[PA[*b]]) == v);) {
        }
        if(((a = b) < last) && (x < v)) {
            for (; (++b < last) && ((x = Td[PA[*b]]) <= v);) {
                if(x == v) {
                    SWAP(*b, *a);
                    ++a;
                }
            }
        }
        for (c = last; (b < --c) && ((x = Td[PA[*c]]) == v);) {
        }
        if((b < (d = c)) && (x > v)) {
            for (; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
                if(x == v) {
                    SWAP(*c, *d);
                    --d;
                }
            }
        }
        for (; b < c;) {
            SWAP(*b, *c);
            for (; (++b < c) && ((x = Td[PA[*b]]) <= v);) {
                if(x == v) {
                    SWAP(*b, *a);
                    ++a;
                }
            }
            for (; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
                if(x == v) {
                    SWAP(*c, *d);
                    --d;
                }
            }
        }

        if(a <= d) {
            c = b - 1;

            if((s = a - first) > (t = b - a)) {
                s = t;
            }
            for (e = first, f = b - s; 0 < s; --s, ++e, ++f) {
                SWAP(*e, *f);
            }
            if((s = d - c) > (t = last - d - 1)) {
                s = t;
            }
            for (e = b, f = last - s; 0 < s; --s, ++e, ++f) {
                SWAP(*e, *f);
            }

            a = first + (b - a), c = last - (d - c);
            b = (v <= Td[PA[*a] - 1]) ? a : ss_partition(PA, a, c, depth);

            if((a - first) <= (last - c)) {
                if((last - c) <= (c - b)) {
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    STACK_PUSH(c, last, depth, limit);
                    last = a;
                } else if((a - first) <= (c - b)) {
                    STACK_PUSH(c, last, depth, limit);
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    last = a;
                } else {
                    STACK_PUSH(c, last, depth, limit);
                    STACK_PUSH(first, a, depth, limit);
                    first = b, last = c, depth += 1, limit = ss_ilg(c - b);
                }
            } else {
                if((a - first) <= (c - b)) {
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    STACK_PUSH(first, a, depth, limit);
                    first = c;
                } else if((last - c) <= (c - b)) {
                    STACK_PUSH(first, a, depth, limit);
                    STACK_PUSH(b, c, depth + 1, ss_ilg(c - b));
                    first = c;
                } else {
                    STACK_PUSH(first, a, depth, limit);
                    STACK_PUSH(c, last, depth, limit);
                    first = b, last = c, depth += 1, limit = ss_ilg(c - b);
                }
            }
        } else {
            limit += 1;
            if(Td[PA[*first] - 1] < v) {
                first = ss_partition(PA, first, last, depth);
                limit = ss_ilg(last - first);
            }
            depth += 1;
        }
    }
#undef STACK_SIZE
}

#endif /* (SS_BLOCKSIZE == 0) || (SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE) */

/*---------------------------------------------------------------------------*/

#if SS_BLOCKSIZE != 0

static inline
void ss_blockswap(saidx_t *a, saidx_t *b, saidx_t n) {
    saidx_t t;
    for (; 0 < n; --n, ++a, ++b) {
        t = *a, *a = *b, *b = t;
    }
}

static inline
void ss_rotate(saidx_t *first, saidx_t *middle, saidx_t *last) {
    saidx_t *a, *b, t;
    saidx_t l, r;
    l = middle - first, r = last - middle;
    for (; (0 < l) && (0 < r);) {
        if(l == r) {
            ss_blockswap(first, middle, l);
            break;
        }
        if(l < r) {
            a = last - 1, b = middle - 1;
            t = *a;
            do {
                *a-- = *b, *b-- = *a;
                if(b < first) {
                    *a = t;
                    last = a;
                    if((r -= l + 1) <= l) {
                        break;
                    }
                    a -= 1, b = middle - 1;
                    t = *a;
                }
            } while (1);
        } else {
            a = first, b = middle;
            t = *a;
            do {
                *a++ = *b, *b++ = *a;
                if(last <= b) {
                    *a = t;
                    first = a + 1;
                    if((l -= r + 1) <= r) {
                        break;
                    }
                    a += 1, b = middle;
                    t = *a;
                }
            } while (1);
        }
    }
}

/*---------------------------------------------------------------------------*/

static
void ss_inplacemerge(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
                     saidx_t *middle, saidx_t *last, saidx_t depth) {
    const saidx_t *p;
    saidx_t *a, *b;
    saidx_t len, half;
    saint_t q, r;
    saint_t x;

    for (;;) {
        if(*(last - 1) < 0) {
            x = 1;
            p = PA + ~*(last - 1);
        } else {
            x = 0;
            p = PA + *(last - 1);
        }
        for (a = first, len = middle - first, half = len >> 1, r = -1; 0 < len;
                len = half, half >>= 1) {
            b = a + half;
            q = ss_compare(T, PA + ((0 <= *b) ? *b : ~*b), p, depth);
            if(q < 0) {
                a = b + 1;
                half -= (len & 1) ^ 1;
            } else {
                r = q;
            }
        }
        if(a < middle) {
            if(r == 0) {
                *a = ~*a;
            }
            ss_rotate(a, middle, last);
            last -= middle - a;
            middle = a;
            if(first == middle) {
                break;
            }
        }
        --last;
        if(x != 0) {
            while (*--last < 0) {
            }
        }
        if(middle == last) {
            break;
        }
    }
}

/*---------------------------------------------------------------------------*/

/* Merge-forward with internal buffer. */
static
void ss_mergeforward(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
                     saidx_t *middle, saidx_t *last, saidx_t *buf,
                     saidx_t depth) {
    saidx_t *a, *b, *c, *bufend;
    saidx_t t;
    saint_t r;

    bufend = buf + (middle - first) - 1;
    ss_blockswap(buf, first, middle - first);

    for (t = *(a = first), b = buf, c = middle;;) {
        r = ss_compare(T, PA + *b, PA + *c, depth);
        if(r < 0) {
            do {
                *a++ = *b;
                if(bufend <= b) {
                    *bufend = t;
                    return;
                }
                *b++ = *a;
            } while (*b < 0);
        } else if(r > 0) {
            do {
                *a++ = *c, *c++ = *a;
                if(last <= c) {
                    while (b < bufend) {
                        *a++ = *b, *b++ = *a;
                    }
                    *a = *b, *b = t;
                    return;
                }
            } while (*c < 0);
        } else {
            *c = ~*c;
            do {
                *a++ = *b;
                if(bufend <= b) {
                    *bufend = t;
                    return;
                }
                *b++ = *a;
            } while (*b < 0);

            do {
                *a++ = *c, *c++ = *a;
                if(last <= c) {
                    while (b < bufend) {
                        *a++ = *b, *b++ = *a;
                    }
                    *a = *b, *b = t;
                    return;
                }
            } while (*c < 0);
        }
    }
}

/* Merge-backward with internal buffer. */
static
void ss_mergebackward(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
                      saidx_t *middle, saidx_t *last, saidx_t *buf,
                      saidx_t depth) {
    const saidx_t *p1, *p2;
    saidx_t *a, *b, *c, *bufend;
    saidx_t t;
    saint_t r;
    saint_t x;

    bufend = buf + (last - middle) - 1;
    ss_blockswap(buf, middle, last - middle);

    x = 0;
    if(*bufend < 0) {
        p1 = PA + ~*bufend;
        x |= 1;
    } else {
        p1 = PA + *bufend;
    }
    if(*(middle - 1) < 0) {
        p2 = PA + ~*(middle - 1);
        x |= 2;
    } else {
        p2 = PA + *(middle - 1);
    }
    for (t = *(a = last - 1), b = bufend, c = middle - 1;;) {
        r = ss_compare(T, p1, p2, depth);
        if(0 < r) {
            if(x & 1) {
                do {
                    *a-- = *b, *b-- = *a;
                } while (*b < 0);
                x ^= 1;
            }
            *a-- = *b;
            if(b <= buf) {
                *buf = t;
                break;
            }
            *b-- = *a;
            if(*b < 0) {
                p1 = PA + ~*b;
                x |= 1;
            } else {
                p1 = PA + *b;
            }
        } else if(r < 0) {
            if(x & 2) {
                do {
                    *a-- = *c, *c-- = *a;
                } while (*c < 0);
                x ^= 2;
            }
            *a-- = *c, *c-- = *a;
            if(c < first) {
                while (buf < b) {
                    *a-- = *b, *b-- = *a;
                }
                *a = *b, *b = t;
                break;
            }
            if(*c < 0) {
                p2 = PA + ~*c;
                x |= 2;
            } else {
                p2 = PA + *c;
            }
        } else {
            if(x & 1) {
                do {
                    *a-- = *b, *b-- = *a;
                } while (*b < 0);
                x ^= 1;
            }
            *a-- = ~*b;
            if(b <= buf) {
                *buf = t;
                break;
            }
            *b-- = *a;
            if(x & 2) {
                do {
                    *a-- = *c, *c-- = *a;
                } while (*c < 0);
                x ^= 2;
            }
            *a-- = *c, *c-- = *a;
            if(c < first) {
                while (buf < b) {
                    *a-- = *b, *b-- = *a;
                }
                *a = *b, *b = t;
                break;
            }
            if(*b < 0) {
                p1 = PA + ~*b;
                x |= 1;
            } else {
                p1 = PA + *b;
            }
            if(*c < 0) {
                p2 = PA + ~*c;
                x |= 2;
            } else {
                p2 = PA + *c;
            }
        }
    }
}

/* D&C based merge. */
static
void ss_swapmerge(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
                  saidx_t *middle, saidx_t *last, saidx_t *buf, saidx_t bufsize,
                  saidx_t depth) {
#define STACK_SIZE SS_SMERGE_STACKSIZE
#define GETIDX(a) ((0 <= (a)) ? (a) : (~(a)))
#define MERGE_CHECK(a, b, c)\
  do {\
    if(((c) & 1) ||\
       (((c) & 2) && (ss_compare(T, PA + GETIDX(*((a) - 1)), PA + *(a), depth) == 0))) {\
      *(a) = ~*(a);\
    }\
    if(((c) & 4) && ((ss_compare(T, PA + GETIDX(*((b) - 1)), PA + *(b), depth) == 0))) {\
      *(b) = ~*(b);\
    }\
  } while(0)
    struct {
        saidx_t *a, *b, *c;
        saint_t d;
    } stack[STACK_SIZE];
    saidx_t *l, *r, *lm, *rm;
    saidx_t m, len, half;
    saint_t ssize;
    saint_t check, next;

    for (check = 0, ssize = 0;;) {
        if((last - middle) <= bufsize) {
            if((first < middle) && (middle < last)) {
                ss_mergebackward(T, PA, first, middle, last, buf, depth);
            }
            MERGE_CHECK(first, last, check);
            STACK_POP(first, middle, last, check);
            continue;
        }

        if((middle - first) <= bufsize) {
            if(first < middle) {
                ss_mergeforward(T, PA, first, middle, last, buf, depth);
            }
            MERGE_CHECK(first, last, check);
            STACK_POP(first, middle, last, check);
            continue;
        }

        for (m = 0, len = MIN(middle - first, last - middle), half = len >> 1;
                0 < len; len = half, half >>= 1) {
            if(ss_compare(T, PA + GETIDX(*(middle + m + half)),
                          PA + GETIDX(*(middle - m - half - 1)), depth) < 0) {
                m += half + 1;
                half -= (len & 1) ^ 1;
            }
        }

        if(0 < m) {
            lm = middle - m, rm = middle + m;
            ss_blockswap(lm, middle, m);
            l = r = middle, next = 0;
            if(rm < last) {
                if(*rm < 0) {
                    *rm = ~*rm;
                    if(first < lm) {
                        for (; *--l < 0;) {
                        }
                        next |= 4;
                    }
                    next |= 1;
                } else if(first < lm) {
                    for (; *r < 0; ++r) {
                    }
                    next |= 2;
                }
            }

            if((l - first) <= (last - r)) {
                STACK_PUSH(r, rm, last, (next & 3) | (check & 4));
                middle = lm, last = l, check = (check & 3) | (next & 4);
            } else {
                if((next & 2) && (r == middle)) {
                    next ^= 6;
                }
                STACK_PUSH(first, lm, l, (check & 3) | (next & 4));
                first = r, middle = rm, check = (next & 3) | (check & 4);
            }
        } else {
            if(ss_compare(T, PA + GETIDX(*(middle - 1)), PA + *middle, depth)
                    == 0) {
                *middle = ~*middle;
            }
            MERGE_CHECK(first, last, check);
            STACK_POP(first, middle, last, check);
        }
    }
#undef STACK_SIZE
}

#endif /* SS_BLOCKSIZE != 0 */

/*---------------------------------------------------------------------------*/

/*- Function -*/

/* Substring sort */
void sssort(const sauchar_t *T, const saidx_t *PA, saidx_t *first,
            saidx_t *last, saidx_t *buf, saidx_t bufsize, saidx_t depth,
            saidx_t n, saint_t lastsuffix) {
    saidx_t *a;
#if SS_BLOCKSIZE != 0
    saidx_t *b, *middle, *curbuf;
    saidx_t j, k, curbufsize, limit;
#endif
    saidx_t i;

    if(lastsuffix != 0) {
        ++first;
    }

#if SS_BLOCKSIZE == 0
    ss_mintrosort(T, PA, first, last, depth);
#else
    if((bufsize < SS_BLOCKSIZE) && (bufsize < (last - first))
            && (bufsize < (limit = ss_isqrt(last - first)))) {
        if(SS_BLOCKSIZE < limit) {
            limit = SS_BLOCKSIZE;
        }
        buf = middle = last - limit, bufsize = limit;
    } else {
        middle = last, limit = 0;
    }
    for (a = first, i = 0; SS_BLOCKSIZE < (middle - a);
            a += SS_BLOCKSIZE, ++i) {
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
        ss_mintrosort(T, PA, a, a + SS_BLOCKSIZE, depth);
#elif 1 < SS_BLOCKSIZE
        ss_insertionsort(T, PA, a, a + SS_BLOCKSIZE, depth);
#endif
        curbufsize = last - (a + SS_BLOCKSIZE);
        curbuf = a + SS_BLOCKSIZE;
        if(curbufsize <= bufsize) {
            curbufsize = bufsize, curbuf = buf;
        }
        for (b = a, k = SS_BLOCKSIZE, j = i; j & 1; b -= k, k <<= 1, j >>= 1) {
            ss_swapmerge(T, PA, b - k, b, b + k, curbuf, curbufsize, depth);
        }
    }
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
    ss_mintrosort(T, PA, a, middle, depth);
#elif 1 < SS_BLOCKSIZE
    ss_insertionsort(T, PA, a, middle, depth);
#endif
    for (k = SS_BLOCKSIZE; i != 0; k <<= 1, i >>= 1) {
        if(i & 1) {
            ss_swapmerge(T, PA, a - k, a, middle, buf, bufsize, depth);
            a -= k;
        }
    }
    if(limit != 0) {
#if SS_INSERTIONSORT_THRESHOLD < SS_BLOCKSIZE
        ss_mintrosort(T, PA, middle, last, depth);
#elif 1 < SS_BLOCKSIZE
        ss_insertionsort(T, PA, middle, last, depth);
#endif
        ss_inplacemerge(T, PA, first, middle, last, depth);
    }
#endif

    if(lastsuffix != 0) {
        /* Insert last type B* suffix. */
        saidx_t PAi[2];
        PAi[0] = PA[*(first - 1)], PAi[1] = n - 2;
        for (a = first, i = *(first - 1);
                (a < last)
                        && ((*a < 0)
                                || (0 < ss_compare(T, &(PAi[0]), PA + *a, depth)));
                ++a) {
            *(a - 1) = *a;
        }
        *(a - 1) = i;
    }
}

//////////////////////////////////////////////////////////////////////
// libdivsufsort
#ifdef _OPENMP
# include <omp.h>
#endif

/*- Private Functions -*/

/* Sorts suffixes of type B*. */
static saidx_t sort_typeBstar(const sauchar_t *T, saidx_t *SA,
                              saidx_t *bucket_A, saidx_t *bucket_B, saidx_t n) {
    saidx_t *PAb, *ISAb, *buf;
#ifdef _OPENMP
    saidx_t *curbuf;
    saidx_t l;
#endif
    saidx_t i, j, k, t, m, bufsize;
    saint_t c0, c1;
#ifdef _OPENMP
    saint_t d0, d1;
    int tmp;
#endif

    /* Initialize bucket arrays. */
    for (i = 0; i < BUCKET_A_SIZE; ++i) {
        bucket_A[i] = 0;
    }
    for (i = 0; i < BUCKET_B_SIZE; ++i) {
        bucket_B[i] = 0;
    }

    /* Count the number of occurrences of the first one or two characters of each
     type A, B and B* suffix. Moreover, store the beginning position of all
     type B* suffixes into the array SA. */
    for (i = n - 1, m = n, c0 = T[n - 1]; 0 <= i;) {
        /* type A suffix. */
        do {
            ++BUCKET_A(c1 = c0);
        } while ((0 <= --i) && ((c0 = T[i]) >= c1));
        if(0 <= i) {
            /* type B* suffix. */
            ++BUCKET_BSTAR(c0, c1);
            SA[--m] = i;
            /* type B suffix. */
            for (--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1); --i, c1 = c0) {
                ++BUCKET_B(c0, c1);
            }
        }
    }
    m = n - m;
    /*
     note:
     A type B* suffix is lexicographically smaller than a type B suffix that
     begins with the same first two characters.
     */

    /* Calculate the index of start/end point of each bucket. */
    for (c0 = 0, i = 0, j = 0; c0 < ALPHABET_SIZE; ++c0) {
        t = i + BUCKET_A(c0);
        BUCKET_A(c0)= i + j; /* start point */
        i = t + BUCKET_B(c0, c0);
        for (c1 = c0 + 1; c1 < ALPHABET_SIZE; ++c1) {
            j += BUCKET_BSTAR(c0, c1);
            BUCKET_BSTAR(c0, c1) = j; /* end point */
            i += BUCKET_B(c0, c1);
        }
    }

    if(0 < m) {
        /* Sort the type B* suffixes by their first two characters. */
        PAb = SA + n - m;
        ISAb = SA + m;
        for (i = m - 2; 0 <= i; --i) {
            t = PAb[i], c0 = T[t], c1 = T[t + 1];
            SA[--BUCKET_BSTAR(c0, c1)] = i;
        }
        t = PAb[m - 1], c0 = T[t], c1 = T[t + 1];
        SA[--BUCKET_BSTAR(c0, c1)] = m - 1;

        /* Sort the type B* substrings using sssort. */
#ifdef _OPENMP
        tmp = omp_get_max_threads();
        buf = SA + m, bufsize = (n - (2 * m)) / tmp;
        c0 = ALPHABET_SIZE - 2, c1 = ALPHABET_SIZE - 1, j = m;
#pragma omp parallel default(shared) private(curbuf, k, l, d0, d1, tmp)
        {
            tmp = omp_get_thread_num();
            curbuf = buf + tmp * bufsize;
            k = 0;
            for(;;) {
#pragma omp critical(sssort_lock)
                {
                    if(0 < (l = j)) {
                        d0 = c0, d1 = c1;
                        do {
                            k = BUCKET_BSTAR(d0, d1);
                            if(--d1 <= d0) {
                                d1 = ALPHABET_SIZE - 1;
                                if(--d0 < 0) {break;}
                            }
                        }while(((l - k) <= 1) && (0 < (l = k)));
                        c0 = d0, c1 = d1, j = k;
                    }
                }
                if(l == 0) {break;}
                sssort(T, PAb, SA + k, SA + l,
                        curbuf, bufsize, 2, n, *(SA + k) == (m - 1));
            }
        }
#else
        buf = SA + m, bufsize = n - (2 * m);
        for (c0 = ALPHABET_SIZE - 2, j = m; 0 < j; --c0) {
            for (c1 = ALPHABET_SIZE - 1; c0 < c1; j = i, --c1) {
                i = BUCKET_BSTAR(c0, c1);
                if(1 < (j - i)) {
                    sssort(T, PAb, SA + i, SA + j, buf, bufsize, 2, n,
                           *(SA + i) == (m - 1));
                }
            }
        }
#endif

        /* Compute ranks of type B* substrings. */
        for (i = m - 1; 0 <= i; --i) {
            if(0 <= SA[i]) {
                j = i;
                do {
                    ISAb[SA[i]] = i;
                } while ((0 <= --i) && (0 <= SA[i]));
                SA[i + 1] = i - j;
                if(i <= 0) {
                    break;
                }
            }
            j = i;
            do {
                ISAb[SA[i] = ~SA[i]] = j;
            } while (SA[--i] < 0);
            ISAb[SA[i]] = j;
        }

        /* Construct the inverse suffix array of type B* suffixes using trsort. */
        trsort(ISAb, SA, m, 1);

        /* Set the sorted order of tyoe B* suffixes. */
        for (i = n - 1, j = m, c0 = T[n - 1]; 0 <= i;) {
            for (--i, c1 = c0; (0 <= i) && ((c0 = T[i]) >= c1); --i, c1 = c0) {
            }
            if(0 <= i) {
                t = i;
                for (--i, c1 = c0; (0 <= i) && ((c0 = T[i]) <= c1);
                        --i, c1 = c0) {
                }
                SA[ISAb[--j]] = ((t == 0) || (1 < (t - i))) ? t : ~t;
            }
        }

        /* Calculate the index of start/end point of each bucket. */
        BUCKET_B(ALPHABET_SIZE - 1, ALPHABET_SIZE - 1) = n; /* end point */
        for (c0 = ALPHABET_SIZE - 2, k = m - 1; 0 <= c0; --c0) {
            i = BUCKET_A(c0 + 1)- 1;
            for(c1 = ALPHABET_SIZE - 1; c0 < c1; --c1) {
                t = i - BUCKET_B(c0, c1);
                BUCKET_B(c0, c1) = i; /* end point */

                /* Move all type B* suffixes to the correct position. */
                for(i = t, j = BUCKET_BSTAR(c0, c1);
                        j <= k;
                        --i, --k) {SA[i] = SA[k];}
            }
            BUCKET_BSTAR(c0, c0 + 1) = i - BUCKET_B(c0, c0) + 1; /* start point */
            BUCKET_B(c0, c0) = i; /* end point */
        }
    }

    return m;
}

/* Constructs the suffix array by using the sorted order of type B* suffixes. */
static
void construct_SA(const sauchar_t *T, saidx_t *SA, saidx_t *bucket_A,
                  saidx_t *bucket_B, saidx_t n, saidx_t m) {
    saidx_t *i, *j, *k;
    saidx_t s;
    saint_t c0, c1, c2;

    if(0 < m) {
        /* Construct the sorted order of type B suffixes by using
         the sorted order of type B* suffixes. */
        for (c1 = ALPHABET_SIZE - 2; 0 <= c1; --c1) {
            /* Scan the suffix array from right to left. */
            for (i = SA + BUCKET_BSTAR(c1, c1 + 1), j = SA + BUCKET_A(c1 + 1)
                    - 1, k = NULL, c2 = -1; i <= j; --j) {
                if(0 < (s = *j)) {
                    assert(T[s] == c1);
                    assert(((s + 1) < n) && (T[s] <= T[s + 1]));
                    assert(T[s - 1] <= T[s]);
                    *j = ~s;
                    c0 = T[--s];
                    if((0 < s) && (T[s - 1] > c0)) {
                        s = ~s;
                    }
                    if(c0 != c2) {
                        if(0 <= c2) {
                            BUCKET_B(c2, c1) = k - SA;
                        }
                        k = SA + BUCKET_B(c2 = c0, c1);
                    }
                    assert(k < j);
                    *k-- = s;
                } else {
                    assert(((s == 0) && (T[s] == c1)) || (s < 0));
                    *j = ~s;
                }
            }
        }
    }

    /* Construct the suffix array by using
     the sorted order of type B suffixes. */
    k = SA + BUCKET_A(c2 = T[n - 1]);
    *k++ = (T[n - 2] < c2) ? ~(n - 1) : (n - 1);
    /* Scan the suffix array from left to right. */
    for (i = SA, j = SA + n; i < j; ++i) {
        if(0 < (s = *i)) {
            assert(T[s - 1] >= T[s]);
            c0 = T[--s];
            if((s == 0) || (T[s - 1] < c0)) {
                s = ~s;
            }
            if(c0 != c2) {
                BUCKET_A(c2)= k - SA;
                k = SA + BUCKET_A(c2 = c0);
            }
            assert(i < k);
            *k++ = s;
        } else {
            assert(s < 0);
            *i = ~s;
        }
    }
}

/* Constructs the burrows-wheeler transformed string directly
 by using the sorted order of type B* suffixes. */
static saidx_t construct_BWT(const sauchar_t *T, saidx_t *SA, saidx_t *bucket_A,
                             saidx_t *bucket_B, saidx_t n, saidx_t m) {
    saidx_t *i, *j, *k, *orig;
    saidx_t s;
    saint_t c0, c1, c2;

    if(0 < m) {
        /* Construct the sorted order of type B suffixes by using
         the sorted order of type B* suffixes. */
        for (c1 = ALPHABET_SIZE - 2; 0 <= c1; --c1) {
            /* Scan the suffix array from right to left. */
            for (i = SA + BUCKET_BSTAR(c1, c1 + 1), j = SA + BUCKET_A(c1 + 1)
                    - 1, k = NULL, c2 = -1; i <= j; --j) {
                if(0 < (s = *j)) {
                    assert(T[s] == c1);
                    assert(((s + 1) < n) && (T[s] <= T[s + 1]));
                    assert(T[s - 1] <= T[s]);
                    c0 = T[--s];
                    *j = ~((saidx_t) c0);
                    if((0 < s) && (T[s - 1] > c0)) {
                        s = ~s;
                    }
                    if(c0 != c2) {
                        if(0 <= c2) {
                            BUCKET_B(c2, c1) = k - SA;
                        }
                        k = SA + BUCKET_B(c2 = c0, c1);
                    }
                    assert(k < j);
                    *k-- = s;
                } else if(s != 0) {
                    *j = ~s;
#ifndef NDEBUG
                } else {
                    assert(T[s] == c1);
#endif
                }
            }
        }
    }

    /* Construct the BWTed string by using
     the sorted order of type B suffixes. */
    k = SA + BUCKET_A(c2 = T[n - 1]);
    *k++ = (T[n - 2] < c2) ? ~((saidx_t) T[n - 2]) : (n - 1);
    /* Scan the suffix array from left to right. */
    for (i = SA, j = SA + n, orig = SA; i < j; ++i) {
        if(0 < (s = *i)) {
            assert(T[s - 1] >= T[s]);
            c0 = T[--s];
            *i = c0;
            if((0 < s) && (T[s - 1] < c0)) {
                s = ~((saidx_t) T[s - 1]);
            }
            if(c0 != c2) {
                BUCKET_A(c2)= k - SA;
                k = SA + BUCKET_A(c2 = c0);
            }
            assert(i < k);
            *k++ = s;
        } else if(s != 0) {
            *i = ~s;
        } else {
            orig = i;
        }
    }

    return orig - SA;
}

/*---------------------------------------------------------------------------*/

/*- Function -*/

saint_t divsufsort(const sauchar_t *T, saidx_t *SA, saidx_t n) {
    saidx_t *bucket_A, *bucket_B;
    saidx_t m;
    saint_t err = 0;

    /* Check arguments. */
    if((T == NULL) || (SA == NULL) || (n < 0)) {
        return -1;
    } else if(n == 0) {
        return 0;
    } else if(n == 1) {
        SA[0] = 0;
        return 0;
    } else if(n == 2) {
        m = (T[0] < T[1]);
        SA[m ^ 1] = 0, SA[m] = 1;
        return 0;
    }

    bucket_A = (saidx_t *) malloc(BUCKET_A_SIZE * sizeof(saidx_t));
    bucket_B = (saidx_t *) malloc(BUCKET_B_SIZE * sizeof(saidx_t));

    /* Suffixsort. */
    if((bucket_A != NULL) && (bucket_B != NULL)) {
        m = sort_typeBstar(T, SA, bucket_A, bucket_B, n);
        construct_SA(T, SA, bucket_A, bucket_B, n, m);
    } else {
        err = -2;
    }

    free(bucket_B);
    free(bucket_A);

    return err;
}

saidx_t divbwt(const sauchar_t *T, sauchar_t *U, saidx_t *A, saidx_t n) {
    saidx_t *B;
    saidx_t *bucket_A, *bucket_B;
    saidx_t m, pidx, i;

    /* Check arguments. */
    if((T == NULL) || (U == NULL) || (n < 0)) {
        return -1;
    } else if(n <= 1) {
        if(n == 1) {
            U[0] = T[0];
        }
        return n;
    }

    if((B = A) == NULL) {
        B = (saidx_t *) malloc((size_t) (n + 1) * sizeof(saidx_t));
    }
    bucket_A = (saidx_t *) malloc(BUCKET_A_SIZE * sizeof(saidx_t));
    bucket_B = (saidx_t *) malloc(BUCKET_B_SIZE * sizeof(saidx_t));

    /* Burrows-Wheeler Transform. */
    if((B != NULL) && (bucket_A != NULL) && (bucket_B != NULL)) {
        m = sort_typeBstar(T, B, bucket_A, bucket_B, n);
        pidx = construct_BWT(T, B, bucket_A, bucket_B, n, m);

        /* Copy to output string. */
        U[0] = T[n - 1];
        for (i = 0; i < pidx; ++i) {
            U[i + 1] = (sauchar_t) B[i];
        }
        for (i += 1; i < n; ++i) {
            U[i] = (sauchar_t) B[i];
        }
        pidx += 1;
    } else {
        pidx = -2;
    }

    free(bucket_B);
    free(bucket_A);
    if(A == NULL) {
        free(B);
    }

    return pidx;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SuffixArray::SuffixArray(const char* text, size_t textLen)
        : text_(text),
          textLen_(textLen) {

    suffix_array_.clear();
    suffix_array_.resize(textLen);
    divsufsort((const sauchar_t *) text, &suffix_array_[0], textLen);
}
SuffixArray::SuffixArray(FILE* in, const char* text, size_t textLen) {
    uint32_t len;
    suffix_array_.clear();
    fread(&len, sizeof(len), 1, in);
    printf("%d -- %d\n", len, textLen);
    assert(len == textLen);

    text_ = text;
    textLen_ = textLen;
    suffix_array_.resize(textLen_);
    fread(&suffix_array_[0], sizeof(std::vector<saidx_t>::value_type),
          suffix_array_.size(), in);
}
SuffixArray::~SuffixArray() {
    text_ = NULL;
}

const int *SuffixArray::search(const char* pattern, int length,
                               int* numOfSolution) {
    int firstIndex;
    *numOfSolution = sa_search((const sauchar_t *) text_, textLen_,
                               (const sauchar_t *) pattern, length,
                               &suffix_array_[0], textLen_, &firstIndex);

    if(*numOfSolution == -1) {
        return numOfSolution;
    } else {
        return &suffix_array_[firstIndex];
    }
}

const int* SuffixArray::iterativeSearch(const char* pattern, int length,
                                        int startLen, int* numOfSolutions,
                                        const int solUpperLimit,
                                        const int solLowerLimit,
                                        int* finalLen) {
    int firstIndex = -1;
    int prevFirstIndex;

    int prevNumOfSolutions;
    *numOfSolutions = -1;

    int currentLen = startLen - 1;

    // TODO: +-+-  exp pa lin
    do {
        currentLen++;
        prevNumOfSolutions = *numOfSolutions;
        prevFirstIndex = firstIndex;
        *numOfSolutions = sa_search((const sauchar_t *) text_, textLen_,
                                    (const sauchar_t *) pattern, currentLen,
                                    &suffix_array_[0], textLen_, &firstIndex);

    } while (*numOfSolutions > solUpperLimit && currentLen < length);

    *finalLen = currentLen;
    if(*numOfSolutions < solLowerLimit && prevNumOfSolutions != -1) {
        *numOfSolutions = prevNumOfSolutions;
        firstIndex = prevFirstIndex;
        *finalLen = currentLen - 1;
    }

    if(*numOfSolutions == -1) {
        return numOfSolutions;
    } else {
        return &suffix_array_[firstIndex];
    }

}

const int* SuffixArray::iterativeSearchDev(const char* pattern, int length,
                                           int startLen, int* numOfSolutions,
                                           const int solUpperLimit,
                                           const int solLowerLimit,
                                           int* finalLen) {
    int firstIndex = -1;
    int prevFirstIndex;

    int prevNumOfSolutions;
    *numOfSolutions = -1;

    int currentLen = startLen - 1;

    // find upper bound
    int step;
    for (step = 1; step + startLen < length; step *= 2) {
        prevNumOfSolutions = *numOfSolutions;
        prevFirstIndex = firstIndex;
        *numOfSolutions = sa_search((const sauchar_t *) text_, textLen_,
                                    (const sauchar_t *) pattern,
                                    step + startLen, &suffix_array_[0],
                                    textLen_, &firstIndex);
        if(*numOfSolutions < solLowerLimit) break;
    }

    //
    int lo = step / 2;
    int hi = step;

    for (int i = lo + 1; i < hi; ++i) {
        prevNumOfSolutions = *numOfSolutions;
        prevFirstIndex = firstIndex;
        *numOfSolutions = sa_search((const sauchar_t *) text_, textLen_,
                                    (const sauchar_t *) pattern, i + startLen,
                                    &suffix_array_[0], textLen_, &firstIndex);
        if(*numOfSolutions < solLowerLimit) {
            currentLen = i + startLen;
            break;
        }
    }

    *finalLen = currentLen;
    if(*numOfSolutions < solLowerLimit && prevNumOfSolutions != -1) {
        *numOfSolutions = prevNumOfSolutions;
        firstIndex = prevFirstIndex;
        *finalLen = currentLen - 1;
    }

    if(*numOfSolutions == -1) {
        return numOfSolutions;
    } else {
        return &suffix_array_[firstIndex];
    }

}

void SuffixArray::saveSuffixArray(FILE* out) {
    fwrite(&textLen_, sizeof(textLen_), 1, out);
    fwrite(&suffix_array_[0], sizeof(std::vector<saidx_t>::value_type),
           suffix_array_.size(), out);
}

uint32_t SuffixArray::size() {
    return textLen_;
}
const char* SuffixArray::text() {
    return text_;
}
/*
 * lcskpp.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#include <algorithm>

using namespace std;

uint32_t LCSkpp::calcLCSkpp(uint32_t k,
                            vector<pair<uint32_t, uint32_t> > &result,
                            vector<pair<uint32_t, uint32_t> > &elements) {

    if(elements.empty()) {
        result.clear();
        return 0;
    }

    vector<event_t> events;
    events.reserve(2 * elements.size());
    uint32_t n = 0;

    for (uint32_t i = 0; i < elements.size(); ++i) {
        pair<uint32_t, uint32_t> element = elements[i];

        events.push_back(event_t(element.first, element.second, true, i));
        events.push_back(
                event_t(element.first + k, element.second + k, false, i));

        n = max(n, element.second + k);
    }
    sort(events.begin(), events.end());

    // Indexed by column, first:dp value, second:index in elements
    Fenwick<pair<uint32_t, uint32_t> > maxColDp(n);

    vector<uint32_t> dp(elements.size());
    vector<int> recon(elements.size());
    vector<int> continues(elements.size(), -1);

    // find pairs continuing each other
    if(k > 1) {
        vector<pair<uint32_t, uint32_t> >::iterator it;
        vector<pair<uint32_t, uint32_t> >::iterator prevIt;

        for (it = elements.begin(); it != elements.end(); ++it) {
            pair<uint32_t, uint32_t> prevElement = make_pair(it->first - 1,
                                                             it->second - 1);
            prevIt = lower_bound(elements.begin(), elements.end(), prevElement);
            if(*prevIt == prevElement) {
                continues[it - elements.begin()] = prevIt - elements.begin();
            }
        }
    }

    uint32_t lcskppLen = 0;
    uint32_t bestIndex = 0;

    for (auto event = events.begin(); event != events.end(); ++event) {
        int index = event->index;

        if(event->isStart) {
            pair<int, int> max = maxColDp.getMax(event->second);
            dp[index] = k;
            recon[index] = -1;

            if(max.first > 0) {
                dp[index] = max.first + k;
                recon[index] = max.second;
            }

        } else {

            if(continues[index] != -1) {
                if(dp[continues[index]] + 1 > dp[index]) {
                    dp[index] = dp[continues[index]] + 1;
                    recon[index] = continues[index];
                }
            }
            maxColDp.updateMax(event->second, make_pair(dp[index], index));

            if(dp[index] > lcskppLen) {
                lcskppLen = dp[index];
                bestIndex = index;
            }
        }
    }

    reconstructLCSkpp(elements, k, recon, bestIndex, lcskppLen, result);
    return lcskppLen;
}

uint32_t LCSkpp::calcLCSkppSlow(
        uint32_t k, vector<pair<uint32_t, uint32_t>> &reconstruction,
        vector<pair<uint32_t, uint32_t>> &matches) {

    if(matches.empty()) {
        reconstruction.clear();
        return 0;
    }

    int n = matches.size();
    vector<int> dp(n);
    vector<int> recon(n);

    int bestIndex = 0;
    int lcskLen = 0;

    for (int i = 0; i < n; ++i) {
        dp[i] = k;
        recon[i] = -1;

        int rowEnd = matches[i].first + k - 1;
        int colEnd = matches[i].second + k - 1;
        int primDiagI = rowEnd - colEnd;
        int secDiagI = (rowEnd + colEnd) / 2;

        for (int j = i - 1; j >= 0; --j) {
            if(matches[j].first + k <= matches[i].first
                    && matches[j].second + k <= matches[i].second) {
                // 1) Uzimam cijeli match interval i nastavljam neki
                // match koji je ranije vec zavrsio.
                if(dp[j] + (int) k > dp[i]) {
                    dp[i] = dp[j] + k;
                    recon[i] = j;
                }
            } else {
                // 2) Nastavak po istoj dijagonali.
                int rowEndJ = matches[j].first + k - 1;
                int colEndJ = matches[j].second + k - 1;
                int primDiagJ = rowEndJ - colEndJ;
                int secDiagJ = (rowEndJ + colEndJ) / 2;

                int extend = secDiagI - secDiagJ;
                if(primDiagI == primDiagJ && secDiagI > secDiagJ
                        && extend < (int) k) {
                    if(dp[j] + extend > dp[i]) {
                        dp[i] = dp[j] + extend;
                        recon[i] = j;
                    }
                }

            }
        }

        if(dp[i] > lcskLen) {
            bestIndex = i;
            lcskLen = dp[i];
        }
    }

    reconstructLCSkpp(matches, k, recon, bestIndex, lcskLen, reconstruction);
    return lcskLen;
}

void LCSkpp::reconstructLCSkpp(
        vector<pair<uint32_t, uint32_t>> &elements, uint32_t k,
        vector<int> &prevIndex, int lastIndex, int lcskLen,
        vector<pair<uint32_t, uint32_t>> &reconstruction) {

    reconstruction.clear();
    reconstruction.reserve(lcskLen);

    int index = lastIndex;
    while (index != -1) {
        int refEndIndex = elements[index].first + k - 1;
        int readEndIndex = elements[index].second + k - 1;

        int prev = prevIndex[index];
        uint32_t howManyElements;

        bool takeWhole = prev == -1;
        if(prev != -1 && elements[prev].first + k <= elements[index].first
                && elements[prev].second + k <= elements[index].second) {
            takeWhole = true;
        }

        if(takeWhole) {
            howManyElements = k;
        } else {
            int curr_secondary_diag = (elements[index].first
                    + elements[index].second) / 2;
            int prev_secondary_diag = (elements[prev].first
                    + elements[prev].second) / 2;
            howManyElements = curr_secondary_diag - prev_secondary_diag;
        }

        for (uint32_t j = 0; j < howManyElements; ++j) {
            reconstruction.push_back(
                    make_pair(refEndIndex - j, readEndIndex - j));
        }
        index = prevIndex[index];
    }

    reverse(reconstruction.begin(), reconstruction.end());
}

/*
 * lcskppV2.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <algorithm>

#include <cassert>

using namespace std;

bool operator<(const eventK_t &a, const eventK_t &b) {
    if(a.first != b.first) return a.first < b.first;

    if(a.second != b.second) return a.second < b.second;

    if(a.index != b.index) return a.index < b.index;

    if(a.isStart != b.isStart) return a.isStart < b.isStart;

    return a.k < b.k;

}

bool operator==(const eventK_t &a, const eventK_t &b) {

    return a.first == b.first && a.second == b.second
            && (a.isStart == b.isStart);
}

uint32_t LCSkppV2::calcLCSkpp(
        std::vector<std::pair<uint32_t, uint32_t>> &result,
        std::vector<triplet_t<uint32_t>> &elements) {

    if(elements.empty()) {
        result.clear();
        return 0;
    }

    vector<eventK_t> events;
    events.reserve(2 * elements.size());
    uint32_t n = 0;

    for (uint32_t i = 0; i < elements.size(); ++i) {
        triplet_t<uint32_t> element = elements[i];

        uint32_t k = element.third;
        events.push_back(eventK_t(element.first, element.second, k, true, i));
        events.push_back(
                eventK_t(element.first + k, element.second + k, k, false, i));

        n = max(n, element.second + k);
    }
    sort(events.begin(), events.end());

    // Indexed by column, first:dp value, second:index in elements
    Fenwick<pair<uint32_t, uint32_t>> maxColDp(n);

    vector<uint32_t> dp(elements.size());
    vector<int> recon(elements.size());
    vector<int> continues(elements.size(), -1);

    // find pairs continuing each other
    vector<triplet_t<uint32_t>>::iterator it;
    vector<triplet_t<uint32_t>>::iterator prevIt;

    for (it = elements.begin(); it != elements.end(); ++it) {
        triplet_t<uint32_t> prevElement = triplet_t<uint32_t>(it->first - 1,
                                                              it->second - 1,
                                                              0);
        prevIt = lower_bound(elements.begin(), elements.end(), prevElement);
        if(prevIt->first == prevElement.first
                && prevIt->second == prevElement.second) {
            continues[it - elements.begin()] = prevIt - elements.begin();
        }
    }

    uint32_t lcskppLen = 0;
    uint32_t bestIndex = 0;

    for (auto event = events.begin(); event != events.end(); ++event) {
        int index = event->index;

        if(event->isStart) {
            pair<int, int> max = maxColDp.getMax(event->second);
            dp[index] = event->k;
            recon[index] = -1;

            if(max.first > 0) {
                dp[index] = max.first + event->k;
                recon[index] = max.second;
            }

        } else {

            if(continues[index] != -1) {
                if(dp[continues[index]] + 1 > dp[index]) {
                    dp[index] = dp[continues[index]] + 1;
                    recon[index] = continues[index];
                }
            }
            maxColDp.updateMax(event->second, make_pair(dp[index], index));

            if(dp[index] > lcskppLen) {
                lcskppLen = dp[index];
                bestIndex = index;
            }
        }
    }
    int k = elements[bestIndex].third;
    int refEndIndex = elements[bestIndex].first + k - 1;
    int readEndIndex = elements[bestIndex].second + k - 1;

    result.push_back( { refEndIndex - lcskppLen, 0 });
    return lcskppLen;
//    reconstructLCSpp(elements, recon, bestIndex, lcskppLen, result);
//    result.erase(unique(result.begin(), result.end()), result.end());
//    if((result.size() - lcskppLen)) cout << "LCSK" << result.size() << "  " << lcskppLen  << endl;
//    return result.size();
}

void LCSkppV2::reconstructLCSpp(
        vector<triplet_t<uint32_t>> &elements, vector<int> &prevIndex,
        int lastIndex, int lcskLen,
        vector<pair<uint32_t, uint32_t>> &reconstruction) {

    reconstruction.clear();
    //reconstruction.reserve(lcskLen);

    int index = lastIndex;
    while (index != -1) {

        int k = elements[index].third;
        int refEndIndex = elements[index].first + k - 1;
        int readEndIndex = elements[index].second + k - 1;

        int prev = prevIndex[index];
        int prevK = elements[prev].third;

        uint32_t howManyElements;

        bool takeWhole = prev == -1;
        if(prev != -1 && elements[prev].first + prevK <= elements[index].first
                && elements[prev].second + prevK <= elements[index].second) {
            takeWhole = true;
        }

        if(takeWhole) {
            howManyElements = k;
        } else {
            int curr_secondary_diag = (elements[index].first
                    + elements[index].second) / 2;
            int prev_secondary_diag = (elements[prev].first
                    + elements[prev].second) / 2;
            howManyElements = curr_secondary_diag - prev_secondary_diag;
        }

        for (uint32_t j = 0; j < howManyElements; ++j) {
            reconstruction.push_back(
                    make_pair(refEndIndex - j, readEndIndex - j));
        }
        index = prevIndex[index];
    }

    reverse(reconstruction.begin(), reconstruction.end());
}

/*
 * lcsk_solver.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <algorithm>

LCSkSolver::LCSkSolver(Sequence* seq) {

    kmerK_ = KMER_K;
    windowSize_ = WINDOW_SIZE;
//  swStartOffset_ = SW_START_OFFSET;
//  swEndOffset_ = SW_END_OFFSET;
    seq_ = seq;
    sa_ = NULL;

}

LCSkSolver::~LCSkSolver() {
    delete sa_;
}

Sequence* LCSkSolver::seq() {
    return seq_;

}
SuffixArray* LCSkSolver::sa() {
    return sa_;
}

void LCSkSolver::readSuffixArrayFromFile(const char* saInPath) {
    FILE* saIn = fopen(saInPath, "rb");

    fprintf(stderr, "Reading suffix array from file\n");
    sa_ = new SuffixArray(saIn, seq_->data(), seq_->dataLen());

    fprintf(stderr, "SuffixArray read\n");
    fclose(saIn);
}

void LCSkSolver::fillPositions(Read* read) {
    std::vector<std::pair<uint32_t, uint32_t> > pos;
    assert(read->dataLen() >= kmerK_);
    for (uint32_t i = kmerK_; i < read->dataLen(); ++i) {
        getKmerPositions(read, pos, i - kmerK_);
    }
    double kmer_avg = 0;
    std::sort(pos.begin(), pos.end());
    if(pos.size() == 0) {
        return;
    }

    uint32_t startIndex = 0;
    uint32_t endIndex = 0;
    while (true) {

        // prozor dug max WINDOW SIZE * duzina reada
        uint32_t windowEnd = pos[startIndex].first
                + windowSize_ * read->dataLen();
        for (; endIndex < pos.size() && pos[endIndex].first < windowEnd;
                ++endIndex) {
            ;
        }

        runLCSkpp(startIndex, endIndex, pos, read);

        uint32_t lastPosition = pos[startIndex].first;
        // pomakni pocetak prozora za duljinu reada od proslog pocetka
        for (;
                startIndex < pos.size()
                        && pos[startIndex].first
                                < lastPosition + read->dataLen();
                ++startIndex) {
            ;
        }
        if(startIndex >= pos.size()) {
            break;
        }

    }
}

void LCSkSolver::runLCSkpp(int startIndex, int endIndex,
                           std::vector<std::pair<uint32_t, uint32_t>> &pos,
                           Read* read) {

    std::vector<std::pair<uint32_t, uint32_t> > lcsKData;

    // TODO Provjera koliko kmera se nalazi izmedu start i end
    // te ukoliko je manje od mog minimuma preskoci
    // ovisno o kvaliteti readova

    //TODO dodat reserve

    for (int i = startIndex; i < endIndex; ++i) {
        lcsKData.push_back(pos[i]);
    }

    std::vector<std::pair<uint32_t, uint32_t> > result;
    uint32_t score = LCSkpp::calcLCSkpp(kmerK_, result, lcsKData);

    //result contains pairs (refPos, readPos)
    int32_t beginPos = result[0].first - result[0].second;
    int32_t endPos = beginPos + read->dataLen() * 2l;

    int len = result.size();
    if(len > 2) {
        endPos = std::min<int32_t>(endPos, result[len - 1].first);
    }

    beginPos = std::max(beginPos, 0);
    read->addPosition(result.size(), beginPos, beginPos + read->dataLen());

}

void LCSkSolver::getKmerPositions(
        Read* read, std::vector<std::pair<uint32_t, uint32_t>> &positions,
        int kmerStart) {

    int numOfSolutions;
    const int* matches = sa_->search(read->data() + kmerStart, kmerK_,
                                     &numOfSolutions);
    if(*matches == -1) {
        return;
    }
    for (int i = 0; i < numOfSolutions; ++i) {
        positions.push_back(std::make_pair(matches[i], kmerStart));
    }

}

void LCSkSolver::findReadPosition(Read* read) {
    if(seq_->basesInt()) {
        read->allBasesToSmallInt();
    }
    Read* reverse_complement = read->getReverseComplement();
    fillPositions(read);
    pos_kmer_cnt = kmer_cnt;
    pos_tot_len = tot_len;

    // add to keep best score and skip lcsk calculations
    // that will never give best score
    Position* p = read->bestPosition(0);
    if(p) {
        reverse_complement->addPosition(p->score(), p->start(), p->end(),
                                        p->isComplement(), NULL, 0,
                                        p->secondaryScore());
    }
    fillPositions(reverse_complement);
    neg_kmer_cnt = kmer_cnt;
    neg_tot_len = tot_len;

    std::multiset<Position*, ptr_compare<Position> >::reverse_iterator it =
            reverse_complement->positions().rbegin();

    for (auto& it : reverse_complement->positions()) {
        if(p && it->score() == p->score() && it->start() == p->start()
                && it->end() == p->end()) {
            // skup position that was previously inserted to track best score
            continue;
        }

        it->setComplement(true);
        read->addPosition(it->score(), it->start(), it->end(), true, NULL, 0,
                          it->secondaryScore());
    }

    std::set<Position*, ptr_compare<Position> > tmp_set = read->positions();
    read->positions().clear();
    read->keepRatio(1.1f);

    for (it = tmp_set.rbegin(); it != tmp_set.rend(); ++it) {
        int32_t start = std::max<int32_t>(
                0, (*it)->start() - 0.545 * read->dataLen());
        uint32_t end = std::min<uint32_t>(sa()->size() - 1,
                                          (*it)->end() + 0.2 * read->dataLen());

        int score, numLocations, alignmentLength;
        int* startLocations = NULL;
        int* endLocations = NULL;
        unsigned char* alignment = NULL;

        if((*it)->isComplement()) {
            edlibCalcEditDistance(
                    (const unsigned char *) reverse_complement->data(),
                    reverse_complement->dataLen(),
                    (const unsigned char *) (sa()->text() + start),
                    end - start + 1, 5, MAX_EDIT, EDLIB_MODE_HW, FIND_STARTS,
                    ALIGN, &score, &endLocations, &startLocations,
                    &numLocations, &alignment, &alignmentLength);

        } else {
            edlibCalcEditDistance(
                    (const unsigned char *) read->data(), read->dataLen(),
                    (const unsigned char *) (sa()->text() + start),
                    end - start + 1, 5, MAX_EDIT, EDLIB_MODE_HW, FIND_STARTS,
                    ALIGN, &score, &endLocations, &startLocations,
                    &numLocations, &alignment, &alignmentLength);
        }
        //        char* cigar;
        //        edlibAlignmentToCigar(alignment, alignmentLength, EDLIB_CIGAR_EXTENDED,
        //                              &cigar);
        //if(startLocations == NULL || alignment == NULL || alignmentLength == 0) continue;
        // if(startLocations == NULL) continue;
        //        start = startLocations[0] + start;
        //        end = endLocations[0] + start;

        end = read->dataLen() + start;

        //int secondaryScore = (*it)->score();
        //int newScore = read->dataLen() - score;  // 3 * read->dataLen() - score;
        //cerr << "BLA " << newScore << endl;
        //        read->addPosition(newScore, start, end, (*it)->isComplement(), NULL, 0,
        //                          tmp_set.size());

        int myScore = score;
        //        for (int i = 0; i < alignmentLength; ++i) {
        //            if(alignment[i] == 1 || alignment[i] == 2) {
        //                myScore += 1;
        //            } else if(alignment[i] == 3) {
        //                myScore += 1;
        //            }
        //        }
        int newScore = read->dataLen() - myScore;
        read->addPosition(newScore, start, end, (*it)->isComplement(), NULL, 0,
                          (*it)->secondaryScore());

        //free(cigar);
        if(endLocations) {
            free(endLocations);
        }
        if(startLocations) {
            free(startLocations);
        }
        if(alignment) {
            free(alignment);
        }
        delete (*it);
    }

    delete reverse_complement;

}

void LCSkSolver::findReadPosition(Read* read, bool orientation) {
    if(seq_->basesInt()) {
        read->allBasesToSmallInt();
    }
    if(orientation) {
        fillPositions(read);
        pos_kmer_cnt = kmer_cnt;
        pos_tot_len = tot_len;
    } else {
        pos_kmer_cnt = 0;
        pos_tot_len = 0;

    }

    Read* reverse_complement = NULL;
    if(!orientation) reverse_complement = read->getReverseComplement();
    // add to keep best score and skip lcsk calculations
    // that will never give best score

    Position* p = read->bestPosition(0);
    if(p && !orientation) {
        reverse_complement->addPosition(p->score(), p->start(), p->end(),
                                        p->isComplement(), NULL, 0,
                                        p->secondaryScore());
    }
    if(!orientation) {
        fillPositions(reverse_complement);
        neg_kmer_cnt = kmer_cnt;
        neg_tot_len = tot_len;
        std::multiset<Position*, ptr_compare<Position> >::reverse_iterator it =
                reverse_complement->positions().rbegin();

        for (auto& it : reverse_complement->positions()) {
            if(p && it->score() == p->score() && it->start() == p->start()
                    && it->end() == p->end()) {
                // skup position that was previously inserted to track best score
                continue;
            }

            it->setComplement(true);
            read->addPosition(it->score(), it->start(), it->end(), true, NULL,
                              0, it->secondaryScore());
        }
    } else {
        neg_kmer_cnt = 0;
        neg_tot_len = 0;
    }

    std::set<Position*, ptr_compare<Position> > tmp_set = read->positions();
    read->positions().clear();
    read->keepRatio(1.1f);

    for (auto it = tmp_set.rbegin(); it != tmp_set.rend(); ++it) {
        int32_t start = std::max<int32_t>(
                0, (*it)->start() - 0.545 * read->dataLen());
        uint32_t end = std::min<uint32_t>(sa()->size() - 1,
                                          (*it)->end() + 0.2 * read->dataLen());

        int score, numLocations, alignmentLength;
        int* startLocations = NULL;
        int* endLocations = NULL;
        unsigned char* alignment = NULL;

        if((*it)->isComplement()) {
            edlibCalcEditDistance(
                    (const unsigned char *) reverse_complement->data(),
                    reverse_complement->dataLen(),
                    (const unsigned char *) (sa()->text() + start),
                    end - start + 1, 5, MAX_EDIT, EDLIB_MODE_HW, FIND_STARTS,
                    ALIGN, &score, &endLocations, &startLocations,
                    &numLocations, &alignment, &alignmentLength);

        } else {
            edlibCalcEditDistance(
                    (const unsigned char *) read->data(), read->dataLen(),
                    (const unsigned char *) (sa()->text() + start),
                    end - start + 1, 5, MAX_EDIT, EDLIB_MODE_HW, FIND_STARTS,
                    ALIGN, &score, &endLocations, &startLocations,
                    &numLocations, &alignment, &alignmentLength);
        }
        //        char* cigar;
        //        edlibAlignmentToCigar(alignment, alignmentLength, EDLIB_CIGAR_EXTENDED,
        //                              &cigar);
        //if(startLocations == NULL || alignment == NULL || alignmentLength == 0) continue;
        // if(startLocations == NULL) continue;
        //        start = startLocations[0] + start;
        //        end = endLocations[0] + start;

        end = read->dataLen() + start;

        //int secondaryScore = (*it)->score();
        //int newScore = read->dataLen() - score;  // 3 * read->dataLen() - score;
        //cerr << "BLA " << newScore << endl;
        //        read->addPosition(newScore, start, end, (*it)->isComplement(), NULL, 0,
        //                          tmp_set.size());

        int myScore = score;
        //        for (int i = 0; i < alignmentLength; ++i) {
        //            if(alignment[i] == 1 || alignment[i] == 2) {
        //                myScore += 1;
        //            } else if(alignment[i] == 3) {
        //                myScore += 1;
        //            }
        //        }
        int newScore = read->dataLen() - myScore;
        read->addPosition(newScore, start, end, (*it)->isComplement(), NULL, 0,
                          (*it)->secondaryScore());
        //free(cigar);
        if(endLocations) {
            free(endLocations);
        }
        if(startLocations) {
            free(startLocations);
        }
        if(alignment) {
            free(alignment);
        }
        delete (*it);
    }

    delete reverse_complement;

}

void LCSkSolver::printInfo() {
    fprintf(stderr, "LCSkSolver: k%d; window:%f;\n", kmerK_, windowSize_);

}

/*
 * incremental_lcsk_solver.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <stdio.h>
#include <algorithm>

IncrementalLCSkSolver::IncrementalLCSkSolver(Sequence* seq)
        : LCSkSolver(seq) {

    maxMatchNum_ = MAX_MATCH_NUM;
    minMatchNum_ = MIN_MATCH_NUM;
}
IncrementalLCSkSolver::~IncrementalLCSkSolver() {
}

void IncrementalLCSkSolver::fillPositions(Read* read) {
    std::vector<triplet_t<uint32_t>> pos;
    assert(read->dataLen() >= kmerK_);
    uint32_t prev_pos_cnt = 0;
    uint32_t len = kmerK_;

    tot_len = 0;
    kmer_cnt = 0;
    for (uint32_t i = kmerK_; i < read->dataLen(); ++i) {
        len = getKmerPositions(read, pos, i - kmerK_, len);
        len = std::max<int>(kmerK_, len);
        tot_len += len;
        kmer_cnt++;
        // nastavi s zavrsnom duzinom - 2

        if(pos.size() > prev_pos_cnt && i < read->dataLen()) {
            // naslo nesto
            uint32_t new_i = i + kmerK_ / 1.0;
            uint32_t new_len = getKmerPositions(read, pos, new_i - kmerK_,
                                                kmerK_);
            new_len = std::max<int>(kmerK_, new_len);

            if(pos.size() >= prev_pos_cnt + minMatchNum_) {
                tot_len += new_len;
                kmer_cnt++;
                // ok, skipamo
                i = new_i;
                len = new_len;  // dodano da ne nastavljamo s krivim lenom
            }
        }

        prev_pos_cnt = pos.size();
    }
    //cerr << "AVG: " << avg_kmer/kmer_cnt << "; cnt " << kmer_cnt << endl;
    std::sort(pos.begin(), pos.end());
    if(pos.size() == 0) {
        return;
    }

    uint32_t startIndex = 0;
    uint32_t endIndex = 0;
    int cnt = 0;
    while (true) {

        // prozor dug max WINDOW SIZE * duzina reada
        uint32_t windowEnd = pos[startIndex].first
                + windowSize_ * read->dataLen();
        for (; endIndex < pos.size() && pos[endIndex].first < windowEnd;
                ++endIndex) {
            ;
        }

        runLCSkpp(startIndex, endIndex, pos, read);
        cnt++;

        uint32_t lastPosition = pos[startIndex].first;
        // pomakni pocetak prozora za duljinu reada od proslog pocetka
        for (;
                startIndex < pos.size()
                        && pos[startIndex].first
                                < lastPosition + read->dataLen();  //vidit i raymislit
                ++startIndex) {
            ;
        }
        if(startIndex >= pos.size()) {
            break;
        }

    }

}

void IncrementalLCSkSolver::runLCSkpp(int startIndex, int endIndex,
                                      std::vector<triplet_t<uint32_t>> &pos,
                                      Read* read) {
    std::vector<triplet_t<uint32_t>> lcsKData;

    // TODO Provjera koliko kmera se nalazi izmedu start i end
    // te ukoliko je manje od mog minimuma preskoci
    // ovisno o kvaliteti readova

    // MOZE BEZ OVOG?????, dam iterator i velicinu
    // std::vector<int>::const_iterator a i to e to
    // jeli vridnoo

    uint32_t aproxMaxLen = 0;
    for (int i = startIndex; i < endIndex; ++i) {
        lcsKData.push_back(pos[i]);
        aproxMaxLen += pos[i].third;
    }

    if(read->bestPosition(0) != NULL) {
        uint32_t currentLen = read->bestPosition(0)->score() / KEEP_F;
        if(currentLen > aproxMaxLen) {
            // no way it could produce better position
            return;
        }
    }

    std::vector<std::pair<uint32_t, uint32_t> > result;
    uint32_t score = LCSkppV2::calcLCSkpp(result, lcsKData);

    //result contains pairs (refPos, readPos)
    int len = score;
    int32_t beginPos = result[0].first - result[0].second;
    int32_t endPos = beginPos + read->dataLen();    //result[len - 1].first;

//    if(len > 2) {
//        endPos = std::min<int32_t>(
//                endPos,
//                result[len - 1].first + read->dataLen()
//                        - result[len - 1].second);
//    }

    beginPos = std::max(beginPos, 0);
    read->addPosition(score, beginPos, endPos, false, NULL, 0, score);
}

uint32_t IncrementalLCSkSolver::getKmerPositions(
        Read* read, std::vector<triplet_t<uint32_t>> &positions, int kmerStart,
        uint32_t initialLen) {

    int numOfSolutions;
    int len;
    const int* matches = sa()->iterativeSearchDev(
            read->data() + kmerStart,
            min(MAX_KMER, read->dataLen() - kmerStart - 1), initialLen,
            &numOfSolutions, maxMatchNum_, minMatchNum_, &len);

    if(*matches == -1) {
        return kmerK_;
    }
    for (int i = 0; i < numOfSolutions; ++i) {
        positions.push_back(triplet_t<uint32_t>(matches[i], kmerStart, len));
    }

    return len;
}

void IncrementalLCSkSolver::printInfo() {
    fprintf(
    stderr,
            "IncrementalLCSkSolver: k:%d; window:%f; max_match:%d; min_match:%d\n",
            kmerK_, windowSize_, maxMatchNum_, minMatchNum_);
}
/*
 * mapper.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <algorithm>

#include <omp.h>
#include <vector>
#include <string>

using namespace std;
Mapper::Mapper(Sequence* seq, Solver* solver, uint32_t threadNum,
               float minKeepScoreRatio, uint32_t maxPositionsPerRead)
        : seq_(seq),
          solver_(solver),
          threadNum_(threadNum),
          minKeepScoreRatio_(minKeepScoreRatio),
          maxPositionsPerRead_(maxPositionsPerRead) {

}

void Mapper::copyFromTmpFileAndDelete(char* tmpFileName, FILE* src,
                                      FILE *dest) {
    char buffer[4096];

    while (fgets(buffer, sizeof buffer, src) != NULL) {
        fputs(buffer, dest);
    }

    fclose(src);
    remove(tmpFileName);
}

void Mapper::createTmpFiles(FILE* tempFiles[],
                            char tmpFileNames[][MAX_TMP_NAME_LEN],
                            uint32_t numOfFiles) {

    for (uint32_t i = 0; i < numOfFiles; ++i) {
        sprintf(tmpFileNames[i], "tmp%d.out", i);
        assert(validateOutputFile(tmpFileNames[i]));
        tempFiles[i] = fopen(tmpFileNames[i], "w");

    }
    printf("Tmp files created\n");

}

void Mapper::fillSAMHeader(FILE* out) {
//  fprintf(out, "@HD\tVN:1.4\tSQ:unsorted\n");
//
//  for (uint32_t i = 0; i < seq_->numOfSequences(); ++i) {
//    fprintf(out, "@SQ\tSN:%s\tLN:%u\n", seq_->info(i), seq_->seqLen(i));
//  }
//  fprintf(out, "@PG\tID:mapper\tPN:mapper\n");
}

void Mapper::mergeTmpFiles(char fileNames[][MAX_TMP_NAME_LEN], FILE* tmpFiles[],
                           FILE* solutionFile, int numberOfFiles) {

    for (int i = 0; i < numberOfFiles; ++i) {
        fclose(tmpFiles[i]);
        tmpFiles[i] = fopen(fileNames[i], "r");
    }
    printf("Merging tmp files\n");
    for (int i = 0; i < numberOfFiles; ++i) {
        copyFromTmpFileAndDelete(fileNames[i], tmpFiles[i], solutionFile);
    }
}
void Mapper::mapAllReads(char* readsInPath, char* solutionOutPath) {
}

}

// oni data

/**
 * Position: describe the position of a read within the genome
 */
struct Position {
    int rname;
    int from;
    int to;
    char strand;
};

/**
 * ReadResult: result of a read alignment
 */
struct ReadResult {
    double confidence;
    int r;
};

map<string, Position> truth;

/*
 * DNASequencing.h
 *
 *  Created on: Apr 17, 2016
 *      Author: marko
 */

#ifndef SRC_CORE_DNASEQUENCING_H_
#define SRC_CORE_DNASEQUENCING_H_
#include <string>
#include <vector>

using namespace std;

class DNASequencing {
 public:
    DNASequencing();
    virtual ~DNASequencing();
    int initTest(int testDifficulty);

    int preProcessing();

    int passReferenceGenome(int chromatidSequenceId,
                            const vector<string>& chromatidSequence);

    vector<string> getAlignment(int N, double normA, double normS,
                                const vector<string>& readName,
                                const vector<string>& readSequence);
    vector<uint32_t> chIds;
    vector<uint32_t> chLens;
    uint32_t totalLen = 0;
    uint32_t currentLen = 0;
    int mode;
    const int pairsDiff[4] = { 450, 450, 450, 450 };
    char* txt = NULL;

    hawker::Sequence* seq = NULL;
    hawker::SuffixArray* sa = NULL;
    hawker::Mapper* mapper = NULL;
    hawker::IncrementalLCSkSolver* solver = NULL;
    //hawker::LCSkSolver* solver;

    int threadNumber = 1;

};

#endif /* SRC_CORE_DNASEQUENCING_H_ */

#include <sstream>
#include <cstdlib>
#include <cstring>

DNASequencing::DNASequencing() {
    cerr << ">> kotor" << endl;
    txt = NULL;
}
DNASequencing::~DNASequencing() {
    cerr << ">> destr" << endl;
    delete seq;
    delete mapper;
    delete solver;
}
int DNASequencing::initTest(int testDifficulty) {
    cerr << ">> init" << endl;
    mode = testDifficulty;

    if(seq) delete seq;
    if(mapper) delete mapper;
    if(solver) delete solver;
    txt = NULL;

    // init
    seq = new hawker::Sequence;
    return 0;
}

int DNASequencing::preProcessing() {
    cerr << "preprocess" << endl;

    totalLen = currentLen;
    seq->numOfSequences_ = chIds.size();
    seq->data_ = txt;
    seq->dataLen_ = totalLen;
    seq->basesAsInt_ = false;

    uint32_t prev = 0;
    for (auto& ends : this->chLens) {
        prev += ends;
        seq->seqEndIndex_.push_back(prev);
    }
    for (auto& id : this->chIds) {
        seq->info_.push_back(id);
    }

    seq->allBasesToSmallInt();

    // TODO LOAD-SAVE SA
    string path = "/media/bio_disk/tc_dna1/demo/" + to_string(mode) + ".suffa";

    if(hawker::SAVE) {
        sa = new hawker::SuffixArray(seq->data(), seq->dataLen());
        cerr << "\tSTORE idx" << endl;
        FILE* fl = fopen(path.c_str(), "wb");
        sa->saveSuffixArray(fl);
        fclose(fl);
    } else if(hawker::LOAD) {
        cerr << "\tLOAD idx" << endl;
        FILE* fl = fopen(path.c_str(), "rb");
        sa = new hawker::SuffixArray(fl, seq->data(), seq->dataLen());
        fclose(fl);
    } else {
        sa = new hawker::SuffixArray(seq->data(), seq->dataLen());
    }

    solver = new hawker::IncrementalLCSkSolver(seq);
    //solver = new hawker::LCSkSolver(seq);
    solver->kmerK_ = hawker::KMER;
    solver->maxMatchNum_ = hawker::HI_CNT;
    solver->minMatchNum_ = hawker::LO_CNT;
    solver->sa_ = sa;

    mapper = new hawker::Mapper(seq, solver, hawker::THREADS, hawker::KEEP_F,
                                hawker::KEEP_NUM);
    return 0;
}

int DNASequencing::passReferenceGenome(
        int chromatidSequenceId, const vector<string>& chromatidSequence) {
    uint32_t len = 0;
    cerr << ">> chr" << chromatidSequenceId << endl;
    this->chIds.push_back(chromatidSequenceId);

    uint32_t vecSz = chromatidSequence.size();
    for (int i = 0; i < vecSz; ++i) {
        uint32_t ln = chromatidSequence[i].size();
        len += ln;
    }
    chLens.push_back(len);

    uint32_t newLen = currentLen + len;
    cerr << "\t total len " << newLen << endl;
    txt = (char*) realloc(txt, newLen * sizeof(char));

    cerr << "\t allocated mem" << endl;
    for (int i = 0; i < vecSz; ++i) {
        memcpy(this->txt + currentLen, chromatidSequence[i].c_str(),
               chromatidSequence[i].size());
        currentLen += chromatidSequence[i].size();
    }
    cerr << "\t copy data - done" << endl;

    return 0;
}

void fillPositions(hawker::Read* r, hawker::Sequence* s, vector<string>& ans) {
    auto mappings = r->getReadMiniSAMs(s);
    for (auto& str : mappings) {
        ans.push_back(str);
    }
}

struct ReadStats {

    uint32_t kmersLen;
    uint32_t kmersCnt;

    uint32_t revKmersLen;
    uint32_t revKmersCnt;

    ReadStats(uint32_t kLen, uint32_t kCnt, uint32_t rkLen, uint32_t rkCnt)
            : kmersLen(kLen),
              kmersCnt(kCnt),
              revKmersLen(rkLen),
              revKmersCnt(rkCnt) {
    }
};

struct BucketData {
    uint32_t id;
    double score;

    BucketData(uint32_t id_, double score_)
            : id(id_),
              score(score_) {

    }

    // desc
    bool operator <(const BucketData& b) const {
        return score > b.score;
    }
};

vector<string> DNASequencing::getAlignment(int N, double normA, double normS,
                                           const vector<string>& readName,
                                           const vector<string>& readSequence) {

    //FILE* out_file = fopen("./ans.minisam", "w");
    int notMapped = 0;
    solver->printInfo();
    vector<string> tmpOutput;

    // ideja sa bucketima
    // INIT
    const int N_BUCKETS = 10;
    const int LAST_BUCKET = 9;
    vector<BucketData> vdb;
    vector<vector<BucketData>> buckets(N_BUCKETS, vdb);

    auto addScore = [&buckets](uint32_t i, double sc, uint32_t ID) {
        BucketData bd(i, sc);
        buckets[ID].push_back(bd);
    };

    for (int i = 0; i < N; i += 2) {
        if(i > hawker::BREAK_CNT) break;
        if(i % 5000 == 0) cerr << i << endl;

        // r1
        hawker::Read* read = new hawker::Read(hawker::KEEP_F, hawker::KEEP_NUM);
        read->setData(readSequence[i], readName[i]);
        read->allBasesToSmallInt();
        solver->findReadPosition(read);
        auto p1s = read->positions();
        ReadStats readStats(hawker::pos_tot_len, hawker::pos_kmer_cnt,
                            hawker::neg_tot_len, hawker::neg_kmer_cnt);

        bool orientation = false;
        bool firstPass = true;
        bool notSame = false;
        for (auto p1 : p1s) {
            if(firstPass) {
                orientation = p1->isComplement();
                firstPass = false;
            } else {
                if(orientation != p1->isComplement()) {
                    notSame = true;
                    break;
                }

            }
        }

        // r2
        hawker::Read* read2 = new hawker::Read(hawker::KEEP_F,
                                               hawker::KEEP_NUM);
        read2->setData(readSequence[i + 1], readName[i + 1]);
        read2->allBasesToSmallInt();

        if(notSame) {
            solver->findReadPosition(read2);
        } else {
            solver->findReadPosition(read2, orientation);
            auto posss = read2->positions();
            if(posss.size() == 0) {
                solver->findReadPosition(read2, !orientation);
            }
        }
        ReadStats readStats2(hawker::pos_tot_len, hawker::pos_kmer_cnt,
                             hawker::neg_tot_len, hawker::neg_kmer_cnt);
        auto p2s = read2->positions();
        //assert(p2s.size() > 0);

        pair<hawker::Position*, hawker::Position*> ans;
        long bestDist = INT32_MAX;
        set<long> results;

        for (auto& p1 : p1s) {
            for (auto& p2 : p2s) {
                if(p1->isComplement() == p2->isComplement()) continue;
                long start1 = p1->start();
                long start2 = p2->start();
                long diff = start1 - start2;

                // NOVO: + pa - ide uvik
                diff = diff * (p1->isComplement() ? 1 : -1);
                if(diff < 0) continue;
                //diff =abs(diff);

                long offset = abs(pairsDiff[mode] - diff);
                if(diff < pairsDiff[mode] && diff < 150) {
                    offset += 450;
                }
                results.insert(offset);
                if(bestDist > offset) {
                    bestDist = offset;
                    ans = {p1, p2};
                }

            }
        }
//          REAL POSITIONS
        auto p1 = truth.find(string(read->id()));
        const Position& position1 = p1->second;
        int32_t startR1 = position1.from;

        auto p2 = truth.find(string(read2->id()));
        const Position& position2 = p2->second;
        int32_t startR2 = position2.from;

        // tocni;
        // cerr << ans.first->start() << "  "  << startR1 << endl;
//                                    long myStartR1 = seq->positionInSeq(ans.first->start());
//                                    long myStartR2 = seq->positionInSeq(ans.second->start());
//
//                                    if(abs(myStartR1 - startR1) < 300
//                                            && abs(myStartR2 - startR2) < 300) {
//                                    } else {
//                                        cerr << p1s.size() << "  " << p2s.size() << endl;
//                                        cerr << "sc " << sc << "; off " << offset << endl;
//                                        cerr << (myStartR1 - startR1) << "  "
//                                             << (myStartR2 - startR2) << endl;
//                                        cerr << ans.first->isComplement() << "    "
//                                             << ans.second->isComplement() << endl;
//
//                                        cerr << readStats.kmersCnt << ";   "
//                                             << readStats.kmersLen << "  |  ";
//                                        cerr << readStats.revKmersCnt << ";   "
//                                             << readStats.revKmersLen << endl;
//                                        cerr << readStats2.kmersCnt << ";   "
//                                             << readStats2.kmersLen << "  |  ";
//                                        cerr << readStats2.revKmersCnt << ";   "
//                                             << readStats2.revKmersLen << endl << endl;
//                                     }

        // Determine scores
        if(results.size() == 0 || *results.begin() > 2 * pairsDiff[mode]) {
            // Worst scores, ambiguous
            fillPositions(read, seq, tmpOutput);
            if(read->positionsSize() == 0) {
                addScore(i, 0, LAST_BUCKET);
                notMapped++;
            } else {
                double sc = read->bestPosition(0)->score() * 1.000;
                addScore(i, 0, LAST_BUCKET);
            }

            fillPositions(read2, seq, tmpOutput);
            if(read2->positionsSize() == 0) {
                addScore(i + 1, 0, LAST_BUCKET);
                notMapped++;
            } else {
                double sc = read2->bestPosition(0)->score() * 1.000;
                addScore(i + 1, 0, LAST_BUCKET);
            }

        } else {
            if(results.size() == 1) {
                tmpOutput.push_back(read->generateMiniSAM(ans.first, 0, seq));
                tmpOutput.push_back(read2->generateMiniSAM(ans.second, 0, seq));

                int offset = *results.begin();
                double sc = (ans.first->score() + ans.second->score()) / 2.0;

                double cntFactor;
                if(p1s.size() + p2s.size() == 2) {
                    cntFactor = 1.0;
                } else {
                    cntFactor = 1 - (p1s.size() + p2s.size()) / 200.0;  // smanji za par posto
                }

                // 100% dobri, skoro unique, rev minimalno mapiran
                // jedine pozicije
                if(ans.first->isComplement() && readStats.kmersCnt > 88
                        && readStats.revKmersCnt < 12
                        && readStats2.kmersCnt < 12
                        && !ans.second->isComplement()
                        && (p1s.size() == 1 && p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 0);
                    addScore(i + 1, sc * cntFactor, 0);

                    // bla bla TODO
                    long myStartR1 = seq->positionInSeq(ans.first->start());
                    long myStartR2 = seq->positionInSeq(ans.second->start());

                    if(abs(myStartR1 - startR1) < 300
                            && abs(myStartR2 - startR2) < 300) {
                    } else {
                        cerr << p1s.size() << "  " << p2s.size() << endl;
                        cerr << "sc " << sc << "; off " << offset << endl;
                        cerr << (myStartR1 - startR1) << "  "
                             << (myStartR2 - startR2) << endl;
                        cerr << ans.first->isComplement() << "    "
                             << ans.second->isComplement() << endl;

                        cerr << readStats.kmersCnt << ";   "
                             << readStats.kmersLen << "  |  ";
                        cerr << readStats.revKmersCnt << ";   "
                             << readStats.revKmersLen << endl;
                        cerr << readStats2.kmersCnt << ";   "
                             << readStats2.kmersLen << "  |  ";
                        cerr << readStats2.revKmersCnt << ";   "
                             << readStats2.revKmersLen << endl << endl;
                    }

                } else if(!ans.first->isComplement()
                        && readStats.revKmersCnt > 88 && readStats.kmersCnt < 12
                        && readStats2.revKmersCnt < 12
                        && ans.second->isComplement()
                        && (p1s.size() == 1 && p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 0);
                    addScore(i + 1, sc * cntFactor, 0);

                    long myStartR1 = seq->positionInSeq(ans.first->start());
                    long myStartR2 = seq->positionInSeq(ans.second->start());

                    if(abs(myStartR1 - startR1) < 300
                            && abs(myStartR2 - startR2) < 300) {
                    } else {
                        cerr << p1s.size() << "  " << p2s.size() << endl;
                        cerr << "sc " << sc << "; off " << offset << endl;
                        cerr << (myStartR1 - startR1) << "  "
                             << (myStartR2 - startR2) << endl;
                        cerr << ans.first->isComplement() << "    "
                             << ans.second->isComplement() << endl;

                        cerr << readStats.kmersCnt << ";   "
                             << readStats.kmersLen << "  |  ";
                        cerr << readStats.revKmersCnt << ";   "
                             << readStats.revKmersLen << endl;
                        cerr << readStats2.kmersCnt << ";   "
                             << readStats2.kmersLen << "  |  ";
                        cerr << readStats2.revKmersCnt << ";   "
                             << readStats2.revKmersLen << endl << endl;
                    }

                    // same but relaxed
                } else if(ans.first->isComplement() && readStats.kmersCnt > 88
                        && readStats.revKmersCnt < 47
                        && readStats2.kmersCnt < 50
                        && !ans.second->isComplement()
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 1);
                    addScore(i + 1, sc * cntFactor, 1);

                } else if(!ans.first->isComplement()
                        && readStats.revKmersCnt > 88 && readStats.kmersCnt < 47
                        && readStats2.revKmersCnt < 50
                        && ans.second->isComplement()
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 1);
                    addScore(i + 1, sc * cntFactor, 1);

                    // same but reverse
                } else if(ans.first->isComplement()
                        && readStats.revKmersCnt < 50
                        && readStats2.kmersCnt < 47
                        && readStats2.revKmersCnt > 88
                        && !ans.second->isComplement()
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 1);
                    addScore(i + 1, sc * cntFactor, 1);

                } else if(!ans.first->isComplement() && readStats.kmersCnt < 50
                        && readStats2.revKmersCnt < 47
                        && readStats2.kmersCnt > 88
                        && ans.second->isComplement()
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 1);
                    addScore(i + 1, sc * cntFactor, 1);

                    // oba sparse
                } else if(ans.first->isComplement()
                        && !ans.second->isComplement()
                        && readStats.revKmersCnt < 12
                        && readStats2.kmersCnt < 12
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 3);
                    addScore(i + 1, sc * cntFactor, 3);

                } else if(!ans.first->isComplement()
                        && ans.second->isComplement() && readStats.kmersCnt < 12
                        && readStats2.revKmersCnt < 12
                        && (p1s.size() == 1 || p2s.size() == 1)) {
                    addScore(i, sc * cntFactor, 3);
                    addScore(i + 1, sc * cntFactor, 3);

                    // negacija i drugi jako sparse
                } else if(ans.first->isComplement()
                        && !ans.second->isComplement()
                        && readStats.kmersCnt > 90 && readStats2.kmersCnt < 10
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 2);
                    addScore(i + 1, sc * cntFactor, 2);

                }

                else if(!ans.first->isComplement() && ans.second->isComplement()
                        && readStats.revKmersCnt > 90
                        && readStats2.revKmersCnt < 10
                        && (p1s.size() == 1 || p2s.size() == 1)) {

                    addScore(i, sc * cntFactor, 2);
                    addScore(i + 1, sc * cntFactor, 2);

                } else {

                    if(offset < 350 && sc > 148 && p1s.size() == 1

                    && p2s.size() == 1) {
                        addScore(i, sc * 1, 5);
                        addScore(i + 1, sc * 1, 5);

                    } else if(offset < 350) {
                        addScore(i, sc * 1.0 / p1s.size(), 6);
                        addScore(i + 1, sc * 1.0 / p1s.size(), 6);
                    } else {
                        addScore(i, sc * 1.0 / (p1s.size() + p2s.size()), 7);
                        addScore(i + 1, sc * 1.0 / (p1s.size() + p2s.size()),
                                 7);

                    }
                }
            }

            else {

                int normals = 0;
                for (auto& sc : results) {
                    if(sc > pairsDiff[mode]) break;
                    ++normals;
                }
                tmpOutput.push_back(read->generateMiniSAM(ans.first, 0, seq));
                tmpOutput.push_back(read2->generateMiniSAM(ans.second, 0, seq));

                int offset = *results.begin();
                double sc = (ans.first->score() + ans.second->score()) / 2.0;
                if(offset < 350 && sc > 148
                        && (p1s.size() == 1 || p2s.size() == 1)) {
                    addScore(i, sc * 1, 8);
                    addScore(i + 1, sc * 1, 8);

                } else {

                    addScore(i, sc * 0.0 / (p1s.size() + p2s.size()), 9);
                    addScore(i + 1, sc * 0.0 / (p1s.size() + p2s.size()), 9);
                }
            }
        }
        delete read;
    }

    for (auto&& vec : buckets) {
        sort(vec.begin(), vec.end());
    }

    vector<double> scores(N, 0);
    int cntr = 0;
    double prevScore;
    bool firstPass = true;
    int b = 0;
    for (auto&& bucket : buckets) {
        if(bucket.empty()) continue;
        int i = 0;
        if(firstPass) {
            prevScore = bucket[0].score;
            scores[bucket[0].id] = cntr;
            cntr++;
            i = 1;
            firstPass = false;
        }

        for (; i < bucket.size(); ++i) {
            //cerr << "bckt " << b <<  ";  score" << bucket[i].score << endl;
            if(abs(bucket[i].score - prevScore) < 1e-6) {
                scores[bucket[i].id] = cntr;
            } else {
                cntr++;
                scores[bucket[i].id] = cntr;
                prevScore = bucket[i].score;
            }
        }
        ++b;
    }

    int totalDiffValues = cntr;
    cerr << "Not mapped " << notMapped << endl;
    cerr << "Diff lvls " << totalDiffValues << endl;

    char score_str[30];
    vector<string> out;
    for (int i = 0; i < tmpOutput.size(); ++i) {
        double score = 1 - (1.0 * scores[i]) / totalDiffValues;
        if(score > hawker::TRESH) {
            std::snprintf(score_str, sizeof score_str, "%.6f", score);
            out.push_back(tmpOutput[i] + string(score_str));
        }

    }
    return out;

}

#ifdef LOCAL_ONLY

string PATH = string("/media/bio_disk/tc_dna1");
string TEST_NUM = to_string(5);

/**
 * Constants from the problem statement
 */
const int MAX_POSITION_DIST = 300;
const double NORM_A_SMALL = -3.392;
const double NORM_A_MEDIUM = -3.962;
const double NORM_A_LARGE = -2.710;
const double MAX_AUC = 0.999999;

/**
 * Split a comma-separated string into a vector of string
 * @param row   the string to be split
 * @return  the vector of string
 */
vector<string> tokenize(const string& row) {
    vector<string> tokens;
    for (int i = 0, pos = 0, n = row.size(); i < n; ++i) {
        if(i == n - 1 || row[i + 1] == ',') {
            string token = row.substr(pos, (i + 1) - pos);
            tokens.push_back(token);
            pos = i + 2;
        }
    }
    return tokens;
}

/**
 * Read a minisam file and build a map of ground truth
 * @param path  the path of the minisam file storing the ground truth
 * @return a map[read_name] = read_Position
 */
map<string, Position> parse_truth(const string& path) {
    map<string, Position> res;
    ifstream ifs(path);
    string s;
    while (ifs >> s) {
        vector<string> tokens = tokenize(s);
        try {
            string qname = tokens[0];
            int chromatid = stoi(tokens[1]);
            int from = stoi(tokens[2]);
            int to = stoi(tokens[3]);
            char strand = tokens[4][0];
            res[qname] = Position { chromatid, from, to, strand };
        } catch (exception& e) {
            cerr << "NOP" << endl;
        }
    }
    return res;
}

/**
 * For each string of the results vector, build a read result {confidence, r}
 * @param truth     the map of ground truth position for each read
 * @param results   the vector of results as return by getAlignment
 * @return a vector of ReadResult, that is {confidence, r}
 */
vector<ReadResult> build_read_results(const map<string, Position>& truth,
                                      const vector<string>& results) {
    vector<ReadResult> read_results;
    cerr << "Size diff " << (truth.size() - results.size()) << endl;
    ;

    int n = results.size();
    int cmpl_cnt = 0;
    int pos_cnt = 0;
    int correct = 0;
    double max_wrong = 0, min_correct = 1;
    for (int i = 0; i < n; ++i) {
        vector<string> tokens = tokenize(results[i]);
//cout << results[i] << endl;
        auto p = truth.find(tokens[0]);
        const Position& position = p->second;
        int r = 1;
        r = (stoi(tokens[1]) == position.rname) ? r : 0;
        r = (tokens[4][0] == position.strand) ? r : 0;
        int start0 = stoi(tokens[2]);
        int start1 = position.from;
//cout << start0 << " " << start1 << endl;
        r = (abs(start0 - start1) < MAX_POSITION_DIST) ? r : 0;
        double confidence = stod(tokens[5]);
        read_results.push_back(ReadResult { confidence, r });
        if(!r) {
            if(tokens[4][0] == '+') {
                pos_cnt++;
            } else {
                cmpl_cnt++;
            }
            if(abs(start0 - start1) < 1250) {
                cerr << "Read " << i << endl;

                cerr << "ja " << confidence << "; delta" << start0 - start1
                     << "; " << endl;
            }
            max_wrong = max(max_wrong, confidence);
        } else {
            min_correct = min(min_correct, confidence);
        }
        correct += r;
    }
    cerr << "Min tocni" << min_correct << endl;
    cerr << "Pos " << pos_cnt << ";  cmpl " << cmpl_cnt << endl;
    cerr << "Max krivi" << max_wrong << endl;
    cerr << "Number of correct answers: " << correct << '/' << n << " = "
         << (double) correct / (double) n << endl;
    return read_results;
}

/**
 * Compute the accuracy given the {confidence, r} pairs and the normalization facto
 * @param read_results  a vector of {confidence, r} results
 * @param norm_a        as described in the problem statement
 * @return  a double, the computed accuracy
 */
double compute_accuracy(vector<ReadResult>& read_results, double norm_a) {
    int n = read_results.size();
    sort(read_results.begin(),
         read_results.end(),
         [](const ReadResult& lhs, const ReadResult& rhs) {return (lhs.confidence>rhs.confidence);});
// merge results of equal confidence
    vector<int> cumul_si { read_results[0].r };
    vector<int> pos { 0 };
    for (int i = 1; i < n; ++i) {
        if(read_results[i].confidence == read_results[i - 1].confidence) {
            cumul_si.back() += read_results[i].r;
            pos.back() = i;
        } else {
            double cumul = cumul_si.back() + read_results[i].r;
            cumul_si.push_back(cumul);
            pos.push_back(i);
        }
    }
    cout << "Lvls: " << pos.size() << "; cmsm_el " << cumul_si.size() << endl;
// compute the AuC
    double auc = 0.0;
    double invn = 1.0 / (double) n;
    double invnp1 = 1.0 / (double) (n + 1);
    double lfmultiplier = 1.0 / log(n + 1);

    int m = cumul_si.size();
    for (int i = 0; i < m; ++i) {
        double fi = 1.0 * (2 + pos[i] - cumul_si[i]) * invnp1;
        double fi1 =
                (i == m - 1) ?
                        1.0 : 1.0 * (2 + pos[i + 1] - cumul_si[i + 1]) * invnp1;
        double lfi = lfmultiplier * log(fi);
        double lfi1 = lfmultiplier * log(fi1);
        auc += cumul_si[i] * (lfi1 - lfi) * invn;
    }
    cout << "auc = " << auc << endl;
    double tmp = log(1 - min(auc, MAX_AUC));
    cout << "log(1 - min(auc, MAX_AUC)) = " << tmp << endl;
    cout << "NormA = " << norm_a << endl;
    double accuracy = tmp / norm_a;
    cout << "accuracy = " << accuracy << endl;
    return accuracy;
}

/**
 * Perform a single test
 * @param testDifficulty    define the test type (SMALL=0, MEDIUM=1, LARGE=2)
 * @return  alignments in format specified in the problem statement
 */
vector<string> perform_test(int testDifficulty, double norm_a) {
// test data path and description
    string fa1_path, fa2_path;
    vector<int> chr_ids;
    if(testDifficulty == 0) {
        fa1_path = PATH + "/data/small" + TEST_NUM + ".fa1";
        fa2_path = PATH + "/data/small" + TEST_NUM + ".fa2";
        chr_ids = vector<int> { 20 };
    } else if(testDifficulty == 1) {
        fa1_path = PATH + "/data/medium" + TEST_NUM + ".fa1";
        fa2_path = PATH + "/data/medium" + TEST_NUM + ".fa2";
        chr_ids = vector<int> { 1, 11, 20 };
    } else if(testDifficulty == 2) {
        fa1_path = PATH + "/data/large" + TEST_NUM + ".fa1";
        fa2_path = PATH + "/data/large" + TEST_NUM + ".fa2";
        for (int i = 1; i <= 24; ++i)
            chr_ids.push_back(i);
    }
// call the MM DNASequencing methods
    DNASequencing dna_sequencing;
    dna_sequencing.initTest(testDifficulty);
// load chromatid
    for (int chromatid_seq_id : chr_ids) {
        vector<string> chromatid_seq;
        string path = PATH + "/data/chromatid" + to_string(chromatid_seq_id)
                + ".fa";
        ifstream ifs(path);
        string s;
// skip header
        getline(ifs, s);
        cerr << "Skip header: " << s << endl;
// pack all lines in chromatid_seq
        for (int i = 0; getline(ifs, s); ++i) {
            if(s.back() == '\r') s.pop_back();
            chromatid_seq.push_back(s);
        }
        dna_sequencing.passReferenceGenome(chromatid_seq_id, chromatid_seq);
    }
    dna_sequencing.preProcessing();
// load reads
    vector<string> read_id, read_seq;
    {
        ifstream ifs1(fa1_path);
        ifstream ifs2(fa2_path);
        string s1, s2;
        while (getline(ifs1, s1) && getline(ifs2, s2)) {
            if(s1.back() == '\r') s1.pop_back();
            if(s2.back() == '\r') s2.pop_back();
            read_id.push_back(s1.substr(1, s1.size() - 1));
            read_id.push_back(s2.substr(1, s2.size() - 1));
            getline(ifs1, s1);
            getline(ifs2, s2);
            if(s1.back() == '\r') s1.pop_back();
            if(s2.back() == '\r') s2.pop_back();
            read_seq.push_back(s1);
            read_seq.push_back(s2);
        }
    }
    int nreads = read_id.size();
// compute alignments
    vector<string> results = dna_sequencing.getAlignment(nreads, norm_a, 0.5,
                                                         read_id, read_seq);
    return results;
}

/**
 * Main function: read the data, perform the DNA alignments and score results
 */
int test(const int testDifficulty) {
    string minisam_path;
    double norm_a;
    if(testDifficulty == 0) {
        minisam_path = PATH + "/data/small" + TEST_NUM + ".minisam";
        norm_a = NORM_A_SMALL;
    } else if(testDifficulty == 1) {
        minisam_path = PATH + "/data/medium" + TEST_NUM + ".minisam";
        norm_a = NORM_A_MEDIUM;
    } else if(testDifficulty == 2) {
        minisam_path = PATH + "/data/large" + TEST_NUM + ".minisam";
        norm_a = NORM_A_LARGE;
    }
// load truth
    truth = parse_truth(minisam_path);

// perform test
    vector<string> results = perform_test(testDifficulty, norm_a);
    cout << "Broj rezultata:" << results.size() << endl;

    vector<ReadResult> read_results = build_read_results(truth, results);
// scoring
    double accuracy = compute_accuracy(read_results, norm_a);
    cout << "aa" << endl;
    return 0;
}

int main(int argc, char** argv) {
    int mode = 0;
    if(argc > 1) {
        mode = atoi(argv[1]);
        if(argc > 2) {
            TEST_NUM = string(argv[2]);
        }
    }
    test(mode);
}
#endif
