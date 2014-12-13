/*
 * main.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "../util/Util.h"
#include "../util/BioUtil.h"
#include "../core/Gene.h"
#include "../core/SuffixArray.h"

using namespace bioutil;
int main(int argc, char **argv) {
	printf("%d\n", argc);
	if (argc != 3) {
		printf("nevalja\n");
		exit(-1);
	}

	assert(isValidInputFile(argv[1]));
	assert(isValidOutputFile(argv[2]));

	FILE* fastaIn = fopen(argv[1], "r");
	Gene *gen = new Gene;
	gen->readGeneFromFASTA(fastaIn);

	FILE* fastaOut = fopen(argv[2], "w");

	//gen->printGene(fastaOut);

	SuffixArray *sa = new SuffixArray;
	sa -> constructFromGene(gen);

	std::vector<std::pair<int, int> > positionsByKmer;
	char pattern[]= "ATTAAAACAG";

	arrayBaseToInt(pattern, strlen(pattern));
	int returnVal = sa->findStartingPositions(pattern, strlen(pattern), 1, positionsByKmer);


	printf("%lu\n", sa->getSize());


	sa->printSuffixArray(fastaOut);
	printf("gotovooo %d\n", returnVal);
	for(int i = 0; i < returnVal && i < 10; ++i) {
		printf("POS %lu\n", positionsByKmer[i].second);
	}


	fclose(fastaIn);
	fclose(fastaOut);
	delete gen;
	delete sa;
	return 0;
}

