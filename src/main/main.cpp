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
#include "../core/Gene.h"
#include "../core/Read.h"
#include "../core/SuffixArray.h"
#include "../util/UtilityFunctions.h"
#include "../util/Mapper.h"

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
	Gene *gen = new Gene;
	gen->readGeneFromFASTA(fastaIn);

	SuffixArray *sa = new SuffixArray;
	sa->constructFromGene(gen);

	FILE* readsIn = fopen(argv[2], "r");
	std::vector<bioutil::Read *> allReads;
	Read::getAllReadsFromFASTQ(readsIn, allReads);

	printf("Ucitani readovi\n");
	int position = getPositionInGeneFromSuffixArray(allReads[0], sa);
	printf("Pozicija\n");
	assert(allReads.size() != 0);
	printf("gotovooo %d\n", position);

	fclose(fastaIn);
	fclose(readsIn);
	delete gen;
	delete sa;
	return 0;
}

