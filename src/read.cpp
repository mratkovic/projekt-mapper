/*
 * read.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */



#include "read.h"
#include "utility_functions.h"

#define LINE_SIZE 30000
#define KEEP_FACTOR 2

namespace bioutil {

Read::Read() {
	id_ = data_ = optional_identifier_ = quality_ = 0;
	dataLen_ = 0;
	mapping_ = new Mapping;
}
Read::~Read() {
	clear();
}

void Read::addMapping(double score, uint32_t start, uint32_t end, bool isComplement, const char* cigar, uint32_t cigarLen){
	if(score >= mapping_->score()) {
		delete mapping_;
		mapping_ = new Mapping(score, start, end, isComplement, cigar, cigarLen);
	}

}
Mapping* Read::mapping() {
	return mapping_;
}

const char* Read::data() {
	return data_;
}
const char* Read::id() {
	return id_;
}

uint32_t Read::dataLen() {
	return dataLen_;
}

void Read::clear() {
	if (data_) {
		delete[] data_;
		data_ = 0;
	}
	if (id_) {
		delete[] id_;
		id_ = 0;
	}
	if (optional_identifier_) {
		delete[] optional_identifier_;
		optional_identifier_ = 0;
	}
	if (quality_) {
		delete[] quality_;
		quality_ = 0;
	}
	if (mapping_) {
			delete mapping_;
			mapping_ = 0;
		}

	dataLen_ = 0;

}

bool Read::readNextFromFASTQ(kseq_t* seq) {
	if (kseq_read(seq) < 0) {
		return false;
	}

	uint32_t size = seq->name.l + seq->comment.l + 2;
	id_ = new char[size];
	memcpy(id_, seq->name.s, seq->name.l);
	id_[seq->name.l] = ' ';

	memcpy(id_ + seq->name.l + 1, seq->comment.s, seq->comment.l);
	id_[size - 1] = 0;

	dataLen_ = seq->seq.l;
	data_ = new char[dataLen_ + 1];
	memcpy(data_, seq->seq.s, dataLen_);
	data_[dataLen_] = 0;

	quality_ = new char[dataLen_ + 1];
	memcpy(quality_, seq->qual.s, dataLen_);
	quality_[dataLen_] = 0;

	return true;
}

Read* Read::getReverseComplement() {
	Read *rev = new Read();
	rev->dataLen_ = dataLen_;
	rev->data_ = new char[dataLen_ + 1];

	for (uint32_t i = 0; i < dataLen_; ++i) {
		rev->data_[i] = getACGTComplement(data_[dataLen_ - 1 - i]);
	}
	rev->data_[dataLen_] = 0;
	return rev;

}
void Read::printReadSAM(FILE* outFile, Sequence* seq) {

	if (mapping_->score() > 0) {
		fprintf(outFile, "%s\t%d\t%s\t%d\t%d\t%s\t%c\t%d\t%d\t%s\t%s\n", id_, mapping_->isComplement() ? 16 : 0,
			seq->info(), mapping_->start(), 99, mapping_->cigar(), '*', 0, 0, data_, quality_);

	//	fprintf(outFile, "%s\t%d\t%d\t%s\t%f\n", id_, mapping_->isComplement() ? 16 : 0, mapping_->start(), mapping_->cigar(),
		//		mapping_->score());

	} else {
		printf("NIJE_NASA\n");
		fprintf(outFile, "%s\t%d\n", id_, -1);

	}
}

}
