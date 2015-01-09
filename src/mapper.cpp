/*
 * mapper.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <algorithm>

#include "external/ssw_cpp.h"
#include "external/kseq.h"
#include "bioutil.h"
#include "mapper.h"
#include "mapping.h"
#include "lis.h"
#include "lcsk.h"
#include "validator.h"
#include "utility_functions.h"

namespace bioutil {
void Mapper::getKmerPositions(Read* read, SuffixArray* sa, std::vector<std::pair<uint32_t, uint32_t> > &positions,
		int kmerStart) {

	int numOfSolutions;
	const int* matches = sa->search(read->data() + kmerStart, KMER_K, &numOfSolutions);

	if (*matches == -1 || numOfSolutions > MAX_NUMBER_OF_KMER_POSITIONS) {
		return;
	}

	for (int i = 0; i < numOfSolutions; ++i) {
		positions.push_back(std::make_pair(matches[i], kmerStart));
	}
}
void Mapper::mapReadToSuffixArray(Read* read, SuffixArray* sa, bool generateCIGAR) {
	Read* reverse_complement = read->getReverseComplement();
	fillMappings(read, sa);
	fillMappings(reverse_complement, sa);

	std::multiset<Mapping*, ptr_compare<Mapping> >::reverse_iterator it = reverse_complement->mappings().rbegin();
	for (; it != reverse_complement->mappings().rend(); ++it) {
		(*it)->setComplement(true);
		read->addMapping((*it)->score(), (*it)->start(), (*it)->end(), true, NULL, 0);
	}

	if (!generateCIGAR) {
		delete reverse_complement;
		return;
	}

	std::multiset<Mapping*, ptr_compare<Mapping> > tmp_set = read->mappings();

	//bool multiple = false;
//	if (read->mappingsSize() > 3) {
//		multiple = true;
//		printf("%s, %d\n", read->id(), read->mappingsSize());
//		printf("%.200s\n", read->data());
//		printf("%.200s\n", reverse_complement->data());
//		for (int i = 0; i < read->mappingsSize(); ++i) {
//			printf("%d.)%f, %d %d\n", i, read->bestMapping(i)->score(), read->bestMapping(i)->start(),
//					read->bestMapping(i)->isComplement());
//
//			printf("%.200s\n", sa->text() + read->bestMapping(i)->start());
//
//		}
//	}
	read->mappings().clear();

	StripedSmithWaterman::Aligner aligner(SSW_MATCH, SSW_MISMATCH, SSW_GAP_OPEN, SSW_GAP_EXTEND);
	StripedSmithWaterman::Filter filter;

	uint32_t missmatchMaxNum = read->dataLen() * MAX_EDIT_DIST_FACTOR;

	filter.score_filter = (read->dataLen() - missmatchMaxNum) * SSW_MATCH - missmatchMaxNum * SSW_MISMATCH;

	for (it = tmp_set.rbegin(); it != tmp_set.rend(); ++it) {
		int32_t start = (*it)->start();
		uint32_t end = (*it)->end();

		StripedSmithWaterman::Alignment alignment;

		if ((*it)->isComplement()) {
			aligner.Align(reverse_complement->data(), sa->text() + start, end - start, filter, &alignment);
		} else {
			aligner.Align(read->data(), sa->text() + start, end - start, filter, &alignment);
		}

		end = start + alignment.ref_end;
		start += alignment.ref_begin;

//		if (multiple) {
//			printf("%d\n", filter.score_filter);
//			printf("SW %d\n", read->mappingsSize());
//			printf("%d.) %f, %d--miss%d %s\n", cntr++, (double) alignment.sw_score, alignment.mismatches,
//					(*it)->start(), alignment.cigar_string.c_str());
//
//		}
		if (alignment.cigar_string.size() == 0) {

			read->addMapping(alignment.sw_score, start, end, (*it)->isComplement(), "NULAAA", 6);
		} else {
			read->addMapping(alignment.sw_score, start, end, (*it)->isComplement(), alignment.cigar_string.c_str(),
					alignment.cigar_string.size());
		}
		delete (*it);
	}

	delete reverse_complement;

}

void Mapper::fillMappings(Read* read, SuffixArray* sa) {
	std::vector<std::pair<uint32_t, uint32_t> > pos;
	assert(read->dataLen() >= KMER_K);

	bool containsN = false;
	uint32_t indexOfN;

// test first K bp for N
	for (uint32_t i = 0; i < KMER_K; ++i) {
		if (read->data()[i] == 'N') {
			containsN = true;
			indexOfN = i;
		}
	}

// skip kmers with N base
	for (uint32_t i = KMER_K; i < read->dataLen(); ++i) {
		if (read->data()[i] == 'N') {
			containsN = true;
			indexOfN = i;
		}

		if (!containsN) {
			getKmerPositions(read, sa, pos, i - KMER_K);
		} else if (indexOfN < i - KMER_K) {
			containsN = false;
		}
	}

	std::sort(pos.begin(), pos.end());
	if (pos.size() == 0) {
		return;
	}

	uint32_t startIndex = 0;
	uint32_t endIndex = 0;
	while (true) {

		// prozor dug max WINDOW SIZE * duzina reada
		uint32_t windowEnd = pos[startIndex].first + WINDOW_SIZE * read->dataLen();
		for (; endIndex < pos.size() && pos[endIndex].first < windowEnd; ++endIndex) {
			;
		}

		//runLIS(startIndex, endIndex, pos, read);
		runLCSk(startIndex, endIndex, pos, read);

		uint32_t lastPosition = pos[startIndex].first;
		// pomakni pocetak prozora za duljinu reada od proslog pocetka
		for (; startIndex < pos.size() && pos[startIndex].first < lastPosition + read->dataLen(); ++startIndex) {
			;
		}
		if (startIndex >= pos.size()) {
			break;
		}

	}
}

void Mapper::runLIS(int startIndex, int endIndex, std::vector<std::pair<uint32_t, uint32_t> > &pos, Read* read) {
	std::vector<std::pair<uint32_t, uint32_t> > lisData;

	for (int i = startIndex; i < endIndex; ++i) {
		lisData.push_back(pos[i]);
	}

	std::vector<int> lisResult;
	Lis::calcLIS(&lisResult, lisData);

	int beginPos = Lis::estimateBeginingPosFromLIS(pos, lisResult);
	read->addMapping(lisResult.size(), beginPos, beginPos + read->dataLen(), false, NULL, 0);
}

void Mapper::runLCSk(int startIndex, int endIndex, std::vector<std::pair<uint32_t, uint32_t> > &pos, Read* read) {
	std::vector<std::pair<uint32_t, uint32_t> > lcsKData;

	for (int i = startIndex; i < endIndex; ++i) {
		lcsKData.push_back(pos[i]);
	}

	std::vector<std::pair<uint32_t, uint32_t> > result;
	uint32_t score = LCSk::calcLCSk(KMER_K, &result, lcsKData);

	//result contains pairs (refPos, readPos)
	int32_t beginPos = result[0].first - result[0].second;
	beginPos = std::max(beginPos, 0);

	read->addMapping(score, beginPos, beginPos + read->dataLen(), false, NULL, 0);
}

void copyFromTmpToFile(char* tmpFileName, FILE* src, FILE *dest) {
	char buffer[4096];

	while (fgets(buffer, sizeof buffer, src) != NULL) {
		fputs(buffer, dest);
	}

	fclose(src);
	remove(tmpFileName);
}

void fillSAMHeader(FILE* out, Sequence* seq) {
	fprintf(out, "@HD\tVN:1.4\tSQ:unsorted\n");
	fprintf(out, "@SQ\tSN:%s\tLN:%u\n", seq->info(), seq->dataLen());
	fprintf(out, "@PG\tID:mapper\tPN:mapper\n");
}

void Mapper::mapAllReads(char* readsInPath, char* solutionOutPath, SuffixArray* sa, Sequence* seq, bool generateCIGAR) {
	FILE* readsIn = fopen(readsInPath, "r");
	FILE* solutionOut = fopen(solutionOutPath, "w");

	fillSAMHeader(solutionOut, seq);
	kseq_t *kseq = kseq_init(fileno(readsIn));
	Read* singleRead = new Read;

	int threadNum = omp_get_num_procs();
	//threadNum = 1;

	omp_set_dynamic(0);
	omp_set_num_threads(threadNum);

	FILE* tmpOutput[threadNum];
	char tmpFilesNames[threadNum][50];

	for (int i = 0; i < threadNum; ++i) {
		sprintf(tmpFilesNames[i], "tmp%d.out", i);
		assert(validateOutputFile(tmpFilesNames[i]));
		tmpOutput[i] = fopen(tmpFilesNames[i], "w");
	}
	printf("Tmp files created\n");
	int cntr = 0;
#pragma omp parallel
	{
#pragma omp single
		{
			while (singleRead->readNextFromFASTQ(kseq)) {
				Read* read = singleRead;

				// create task solve single read
#pragma omp task firstprivate(read) shared(tmpOutput) shared(sa) shared(seq)
				{
					mapReadToSuffixArray(read, sa, generateCIGAR);
					read->printReadSAM(tmpOutput[omp_get_thread_num()], seq);
					//Validator::validateWGSIM(read);
					delete read;
				}

				singleRead = new Read();
				if (++cntr % 10000 == 0) {
					printf("--%d\n", cntr++);
				}
			}
		}
#pragma omp barrier
	}

	delete singleRead;

	for (int i = 0; i < threadNum; ++i) {
		fclose(tmpOutput[i]);
		tmpOutput[i] = fopen(tmpFilesNames[i], "r");
	}
	printf("Merging tmp files\n");
	for (int i = 0; i < threadNum; ++i) {
		copyFromTmpToFile(tmpFilesNames[i], tmpOutput[i], solutionOut);
	}

	fclose(solutionOut);
	fclose(readsIn);
	printf("Completed\n");
}
}

