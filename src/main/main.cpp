/*
 * main.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <utility>
#include <algorithm>

#include "../util/BioUtil.h"
#include "../util/Read.h"
#include "../util/SuffixArray.h"
#include "../util/UtilityFunctions.h"
#include "../util/Mapper.h"
#include "../util/Sequence.h"

using namespace bioutil;

#define KHMER_SIZE 18;

void solve(FILE* input, SuffixArray &sa) {
	Read tmpRead;
	while (tmpRead.readFromFASTQ(input)) {


	}
}

int main(int argc, char **argv) {
	printf("%d\n", argc);
	if (argc != 4) {
		printf("nevalja\n");
		exit(-1);
	}

	assert(isValidInputFile(argv[1]));
	assert(isValidInputFile(argv[2]));
	assert(isValidOutputFile(argv[3]));

	FILE* fastaIn = fopen(argv[1], "r");
	Sequence *seq = new Sequence;
	seq->readSequenceFromFASTA(fastaIn);

	SuffixArray *sa = new SuffixArray;
	sa->constructFromSequence(seq);
	Read r;
	r.clean();
	FILE* readsIn = fopen(argv[2], "r");
	Read read;
	int cntr = 1;
	while(read.readFromFASTQ(readsIn)) {
	printf("Ucitan readovi\n");
		int position = getPositionInSequenceFromSuffixArray(&read, sa);
		printf("%d - %d", cntr++, position);
	}
	printf("gotovooo\n");

	fclose(fastaIn);
	fclose(readsIn);
	delete seq;
	delete sa;
	return 0;
}

