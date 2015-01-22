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
		return false;
	}
	bool wholeLine = buffer[10000] == 0;

	char* name_c;
	name_c = strtok(buffer, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");

	char* positionStr = strtok(NULL, "\t");
	assert(sscanf(positionStr, "%u", position) == 1);

	if (!wholeLine) {
		fscanf(in, "%*s");
	}
	uint32_t len = strlen(name_c);
	if (name_c[len - 2] == '/') {
		// pair end in title, skip it
		name_c[len - 2] = 0;
	}
	name.assign(name_c);
	return true;

}

bool getReadNameAndPositionFromWGSIM(std::string& name, uint32_t* position, FILE* in) {
	//ime \t flag \t referenceName \t position \t cigar .....
	char* buffer = new char[BUFFER_LEN + 1];
	buffer[BUFFER_LEN] = 0;

	if (!fgets(buffer, BUFFER_LEN, in)) {
		return false;
	}
	bool wholeLine = buffer[10000] == 0;

	char* name_c;
	name_c = strtok(buffer, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");

	char* positionStr = strtok(NULL, "\t");
	assert(sscanf(positionStr, "%u", position) == 1);

	if (!wholeLine) {
		fscanf(in, "%*s");
	}
	uint32_t len = strlen(name_c);
	if (name_c[len - 2] == '/') {
		// pair end in title, skip it
		name_c[len - 2] = 0;
	}
	name.assign(name_c);
	return true;

}

void skipSAMHeader(FILE* in) {
	char* buffer = new char[BUFFER_LEN + 1];
	buffer[BUFFER_LEN] = 0;

	fpos_t prevPos;
	fgetpos(in, &prevPos);

	while (fgets(buffer, BUFFER_LEN, in)) {
		if (buffer[0] == '@') {
			fgetpos(in, &prevPos);
		} else {
			fsetpos(in, &prevPos);
			return;
		}
	}
}
void Validator::validateSAM(FILE* ref, FILE* test) {
	uint32_t totalCntr = 0;
	uint32_t wrongMapped = 0;
	uint32_t notMapped = 0;
	std::map<std::string, std::vector<uint32_t> > reads;

	skipSAMHeader(ref);

	while (true) {
		std::string name;
		uint32_t position;
		if (!getReadNameAndPositionFromSAM(name, &position, ref)) {
			break;
		}
		++totalCntr;
		std::map<std::string, std::vector<uint32_t> >::iterator it = reads.find(name);
		if (it == reads.end()) {
			std::vector<uint32_t> vector(0);
			reads.insert(std::make_pair(name, vector));
		}
		reads[name].push_back(position);
	}

	skipSAMHeader(test);

	while (true) {
		std::string name;
		uint32_t position;
		if (!getReadNameAndPositionFromSAM(name, &position, test)) {
			break;
		}

		std::map<std::string, std::vector<uint32_t> >::iterator it = reads.find(name);

		if (it == reads.end()) {
			++notMapped;
		} else if (abs(position - it->second[0]) > TOLERATED_OFFSET) {
			//printf("%s  --%u  --%u\n", name.c_str(), it->second[0], position);
			++wrongMapped;
		}
	}

	printf("Total reads in reference file: %u\n", totalCntr);
	printf("Not mapped: %u\n", notMapped);
	printf("Mapped to wrong position: %u\n", wrongMapped);
	printf("%d / %d / %d\n", notMapped, wrongMapped, totalCntr);

}
void Validator::validateWGSIM(FILE* readsIn, FILE* sam) {
	// TODO


}

} /* namespace bioutil */
