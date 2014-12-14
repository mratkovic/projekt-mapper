/*
 * Fenwick.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_FENWICK_H_
#define SRC_UTIL_FENWICK_H_

#include <vector>

template<class T>
class Fenwick {
private:
	std::vector<T> _elements;
public:
	Fenwick(int n) {
		 _elements = std::vector<T> (n+1, T());
	}
	~Fenwick() {

	}

	void updateMax(int pos, const T& val) {
		++pos;
		for (; pos < _elements.size(); pos += lobitOf(pos)) {
			_elements[pos] = std::max(_elements[pos], val);
		}
	}
	inline int lobitOf(const int& a) {
		return a & -a;
	}

	T getMax(int pos) {
		++pos;
		T ret = T();
		for (; pos > 0; pos -= lobitOf(pos)) {
			ret = std::max(ret, _elements[pos]);
		}
		return ret;
	}
};

#endif /* SRC_UTIL_FENWICK_H_ */
