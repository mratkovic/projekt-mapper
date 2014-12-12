/*
 * Util.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <cstring>

#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <dirent.h>
#include <string>

static inline size_t trimEnd(char *data) {
	int pos = strlen(data) - 1;

	while (pos >= 0 && isspace(data[pos])) {
		data[pos] = 0;
		--pos;
	}
	return pos + 1;
}


static inline bool isValidOutputFile(const std::string& filePath) {
	FILE* f = fopen(filePath.c_str(), "w");
	if (f != NULL) {
		fclose(f);
		return true;
	}
	return false;
}

static inline bool isValidInputFile(const std::string& filePath) {
	printf("%s\n", filePath.c_str());

	FILE* f = fopen(filePath.c_str(), "r");

	if (f != NULL) {
		fclose(f);
		return true;
	}
	return false;
}

static inline bool isValidFolder(const std::string& folderPath) {
	DIR* dir = opendir(folderPath.c_str());
	bool ok = (dir != NULL);
	closedir(dir);
	return ok;
}
#endif /* UTIL_H_ */
