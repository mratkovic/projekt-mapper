/*
 * Fenwick.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 *
 *  Template class that implements Fenwick tree data structure also known as logaritmic structure.
 */

#ifndef SRC_FENWICK_H_
#define SRC_FENWICK_H_

#include <vector>
#include <utility>

template<class T>
class Fenwick {

private:
	std::vector<T> elements_;

public:
	/**
	 * Constructor.
	 *
	 * @param n number of elements or positions
	 */
	Fenwick(int n) {
		elements_ = std::vector<T>(n + 1, T());
	}
	~Fenwick() {
	}

	/**
	 * Function that calculates value of lowest significant bit in given value
	 *
	 * @param a value
	 * @return lowest significant bit of a
	 */
	inline int lobitOf(const int& a) {
		return a & -a;
	}

	/**
	 * Method that inserts value val at position pos and updates data structure with new maximum element.
	 *
	 * @param pos position on wich the value in inserted
	 * @param val value that is inserted in data structure
	 */
	void updateMax(int pos, const T& val) {
		++pos;
		for (; pos < elements_.size(); pos += lobitOf(pos)) {
			elements_[pos] = std::max(elements_[pos], val);
		}
	}
	/**
	 * Method that returns maximum value stored in data structure in position range of [0, pos].
	 *
	 * @param pos upper bound of range
	 * @return maximum element in given range
	 */
	T getMax(int pos) {
		++pos;
		T ret = T();
		for (; pos > 0; pos -= lobitOf(pos)) {
			ret = std::max(ret, elements_[pos]);
		}
		return ret;
	}
};

#endif /* SRC_FENWICK_H_ */
