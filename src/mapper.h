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

#include "fenwick.h"
#include "read.h"
#include "suffix_array.h"
#include "sequence.h"

#define KMER_K 20
#define WINDOW_SIZE 2
#define MAX_NUMBER_OF_KMER_POSITIONS 1000000

#define SSW_MATCH 1
#define SSW_MISMATCH -1
#define SSW_GAP_OPEN -2
#define SSW_GAP_EXTEND -1


namespace bioutil {

class Mapper {

public:
	static void mapAllReads(char* readsInPath, char* solutionOutPath, SuffixArray* sa, Sequence* seq,
			bool generateCIGAR);

private:
	static void runLIS(int startIndex, int endIndex, std::vector<std::pair<uint32_t, uint32_t> > &pos, Read* read);
	static void runLCSk(int startIndex, int endIndex, std::vector<std::pair<uint32_t, uint32_t> > &pos, Read* read);

	static void mapReadToSuffixArray(Read* read, SuffixArray* sa, bool generateCIGAR);
	static void getKmerPositions(Read* read, SuffixArray* sa, std::vector<std::pair<uint32_t, uint32_t> > &positions,
			int kmerStart);
	static void fillMappings(Read* read, SuffixArray* sa);
	static void getPositions(Read* read, SuffixArray* sa, bool complement);

};
}
#endif
