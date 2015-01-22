/*
 * sequence.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>

#include "sequence.h"
#include "external/kseq.h"
#include "utility_functions.h"
#include "bioutil.h"

namespace bioutil {

KSEQ_INIT(int, read)

Sequence::Sequence(char* data, uint32_t dataLen, char* info, char* comment) {
	data_ = data;
	dataLen_ = dataLen;
	info_ = info;
	comment_ = comment;
	basesAsInt = false;

}
Sequence::Sequence() {
	data_ = info_ = comment_ = NULL;
	dataLen_ = 0;
	basesAsInt = false;
}
Sequence::~Sequence() {
	clear();
}

void Sequence::clear() {
	if (data_) {
		delete[] data_;
		data_ = 0;
	}
	if (info_) {
		delete[] info_;
		info_ = 0;
	}
	if (comment_) {
		delete[] comment_;
		comment_ = 0;
	}
	dataLen_ = 0;
}

const char* Sequence::data() {
	return data_;
}
const char* Sequence::info() {
	return info_;
}
const char* Sequence::comment() {
	return comment_;
}
uint32_t Sequence::dataLen() {
	return dataLen_;
}

void Sequence::readSequenceFromFASTA(FILE* fastaIn) {
	clear();
	kseq_t *seq = kseq_init(fileno(fastaIn));

	assert(kseq_read(seq) >= 0);

	uint32_t size = seq->name.l + 1;
	info_ = new char[size];
	memcpy(info_, seq->name.s, seq->name.l);
	info_[seq->name.l] = 0;

	size = seq->comment.l + 1;
	comment_ = new char[size];
	memcpy(comment_, seq->comment.s, seq->comment.l);
	comment_[seq->comment.l] = 0;

	dataLen_ = seq->seq.l;
	data_ = new char[dataLen_];
	memcpy(data_, seq->seq.s, dataLen_);

	kseq_destroy(seq);
}
void Sequence::allBasesToSmallInt() {
	if (basesAsInt) {
		return;
	}

	for (uint32_t i = 0; i < dataLen_; ++i) {
		data_[i] = baseToInt(data_[i]);
	}
	basesAsInt = true;
}
void Sequence::allBasesToLetters() {
	if (!basesAsInt) {
		return;
	}

	for (uint32_t i = 0; i < dataLen_; ++i) {
		data_[i] = intToBase(data_[i]);
	}
	basesAsInt = false;
}

}

