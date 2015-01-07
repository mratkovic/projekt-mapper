/*
 * read.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SRC_CORE_READ_H_
#define SRC_CORE_READ_H_


#include <zlib.h>

#include "sequence.h"
#include "bioutil.h"
#include "external/kseq.h"
#include "suffix_array.h"
#include "mapping.h"


namespace bioutil {

KSEQ_INIT(int, read)

class Read {
private:
	char* id_;
	char* data_;
	uint32_t dataLen_;

	char* optional_identifier_;
	char* quality_;

	Mapping* mapping_;

public:
	Read();

	~Read();

	void clear();
	bool readNextFromFASTQ(kseq_t *seq);
	Read* getReverseComplement();
	void printReadSAM(FILE* outFile, Sequence* seq);

	const char* data();
	const char* id();
	uint32_t dataLen();


	void addMapping(double score, uint32_t start, uint32_t end, bool isComplement, const char* cigar, uint32_t cigarLen);
	Mapping* mapping();


};

}
/* namespace bioutil */
#endif /* SRC_CORE_READ_H_ */
