/*
 * main.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
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

void displayInvalidCallMsg() {
	printf("Invalid call\n");
	printf("mapper construct <fastaFile> <suffixArrayOutputFile>\n");
	printf("mapper map <fastaFile> <suffixArrayOutputFile> <reads> <resultOutputFile>\n");
	exit(-1);
}

void constructSA(char* fastaInPath, char* saOutputPath) {
	assert(canReadFromFile(fastaInPath));
	assert(canWriteToFIle(saOutputPath));

	FILE* fastaIn = fopen(fastaInPath, "r");
	Sequence *seq = new Sequence;
	printf("Reading sequence file %s\n", fastaInPath);
	assert(seq->readSequenceFromFASTA(fastaIn) > 0);
	fclose(fastaIn);

	SuffixArray *sa = new SuffixArray;
	printf("Constructing suffix array started\n");
	sa->constructFromSequence(seq);
	printf("Constructing suffix array completed\n");
	printf("Size of array: %d\n", sa->getSize());

	printf("Saving suffix array to %s\n", saOutputPath);
	FILE* saOut = fopen(saOutputPath, "wb");
	assert(sa->saveSuffixArray(saOut));
	printf("Completed\n");

	delete sa;
	delete seq;

	fclose(saOut);

}

void mapReads(char* fastaInPath, char* saFile, char* readsInPath, char* outputFile) {
	assert(canReadFromFile(fastaInPath));
	assert(canReadFromFile(saFile));
	assert(canReadFromFile(readsInPath));
	assert(canWriteToFIle(outputFile));

	FILE* fastaIn = fopen(fastaInPath, "r");
	Sequence *seq = new Sequence;

	printf("Reading sequence file %s\n", fastaInPath);
	assert(seq->readSequenceFromFASTA(fastaIn) > 0);
	fclose(fastaIn);

	FILE* saIn = fopen(saFile, "r");
	printf("Reading suffix array from file\n");
	SuffixArray *sa = new SuffixArray(saIn);
	sa->setSequence(seq->data());
	printf("SuffixArray read\n");
	fclose(saIn);

	FILE* readsIn = fopen(readsInPath, "r");
	FILE* out = fopen(outputFile, "w");
	Read singleRead;

	int cntr = 0;
	while (singleRead.readFromFASTQ(readsIn)) {
		if(cntr++%1000 == 0) {
			printf("--%d\n", cntr);
		}
		mapReadToSuffixArray(&singleRead, sa);
		char mappingDetails[300];
		singleRead.getBestMapping().fillDetails(mappingDetails);
		fprintf(out, "%s - %s\n", singleRead.getId(), mappingDetails);
		singleRead.clean();
	}

	printf("Completed\n");

	fclose(readsIn);
	fclose(out);

	delete sa;
	delete seq;
}

int main(int argc, char **argv) {

	if (argc != 4 && argc != 6) {
		displayInvalidCallMsg();
	}
	if (!strcmp(argv[1], "construct") && argc == 4) {
		constructSA(argv[2], argv[3]);
	} else if (!strcmp(argv[1], "map") && argc == 6) {
		mapReads(argv[2], argv[3], argv[4], argv[5]);
	} else {
		displayInvalidCallMsg();
	}
	return 0;
}

