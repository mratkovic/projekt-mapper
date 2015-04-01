/*
 * utility_functions.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 *
 */
#ifndef SRC_UTIL_UTILITY_H_
#define SRC_UTIL_UTILITY_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <dirent.h>
#include <unistd.h>

inline bool validateOutputFile(char*filePath) {
  fprintf(stderr, "W %s\n", filePath);
  FILE* f = fopen(filePath, "w");
  if (f != NULL) {
    fclose(f);
    return true;
  }
  return false;
}

inline bool validateInputFile(char*filePath) {
  fprintf(stderr, "R %s\n", filePath);
  FILE* f = fopen(filePath, "r");

  if (f != NULL) {
    fclose(f);
    return true;
  }
  return false;
}

inline bool isFolder(char*folderPath) {
  DIR* dir = opendir(folderPath);
  if (dir != NULL) {
    closedir(dir);
    return true;
  }
  return false;
}

template<typename T, typename Pred = std::less<T> >
struct ptr_compare : Pred {
  ptr_compare(Pred const & p = Pred())
      : Pred(p) {
  }

  bool operator()(T const * p1, T const * p2) const {
    return Pred::operator()(*p1, *p2);
  }
};

#endif
