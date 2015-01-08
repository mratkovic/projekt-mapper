/*
 * validator.cpp
 *
 *  Created on: Jan 4, 2015
 *      Author: marko
 */

#include "validator.h"
#include <map>
#include <vector>
#include <cstdio>
#include <stdint.h>
#include <string>

#define BUFFER_LEN 10000

namespace bioutil {

bool getReadNameAndPositionFromSAM(std::string& name, uint32_t* position, FILE* in) {
	//ime \t flag \t referenceName \t position \t cigar .....
	char* buffer = new char[BUFFER_LEN + 1];
	buffer[BUFFER_LEN] = 0;

	if (!fgets(buffer, BUFFER_LEN, in)) {
		return false;;
	}
	bool wholeLine = buffer[10000] == 0;

	char* name_c;
	name_c = strtok(buffer, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");

	char* positionStr = strtok(NULL, "\t");
	assert(sscanf(positionStr, "%u", &position) != 1);

	if (!wholeLine) {
		fscanf(in, "%*s");
	}
	name.assign(name_c);

}

void Validator::validateSAM(FILE* ref, FILE* test) {
	uint32_t totalCntr = 0;
	uint32_t wrongMapped = 0;
	uint32_t notMapped = 0;
	std::map<std::string, std::vector<uint32_t> > reads;

	// TODO ad hoc
	//skip first three lines
	fscanf(ref, "%*[^\n]\n");
	fscanf(ref, "%*[^\n]\n");
	fscanf(ref, "%*[^\n]\n");

	while (true) {
		std::string name;
		uint32_t position;
		if (!getReadNameAndPositionFromSAM(name, &position, ref)) {
			break;
		}
		std::map<std::string, std::vector<uint32_t> >::iterator it = reads.find(name);
		if (it == reads.end()) {
			std::vector<uint32_t> vector(0);
			reads.insert(std::make_pair(name, vector));
		}
		reads[name].push_back(position);
	}

	fscanf(test, "%*[^\n]\n");
	fscanf(test, "%*[^\n]\n");
	fscanf(test, "%*[^\n]\n");

	while (true) {
		std::string name;
		uint32_t position;
		getReadNameAndPositionFromSAM(name, &position, ref);

		std::map<std::string, std::vector<uint32_t> >::iterator it = reads.find(name);

		if (it == reads.end()) {
			++notMapped;
		} else if (abs(position - it->second[0]) > TOLERATED_OFFSET) {
			++wrongMapped;
		}
		scanf("%*d");
	}

}
void Validator::validateWGSIM(Read* read) {
	//wgsim header
	//readID_posStart_posEnd_.....
	// TODO
}

} /* namespace bioutil */
