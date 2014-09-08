/*
 * utilities.h
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <vector>

using namespace std;

namespace SimpleCluster {
typedef vector<double> d_vector;
typedef vector<int> i_vector;

/**
 * Calculate the distance between two vectors
 */
double distance(d_vector, d_vector, size_t);
double distance_square(d_vector, d_vector, size_t);
d_vector mean_vector(vector<d_vector>, i_vector, size_t, d_vector);
/**
 * Get system time in milliseconds
 */
unsigned long get_millisecond_time();
void print_vector(vector<d_vector>, size_t);
/**
 * Swap two elements in an array
 * data: the array of elements
 * m,n : the indices of elements to be swapped
 * N: the size of data
 */
template<typename DataType>
void swap(DataType * data, size_t m, size_t n, size_t N);
/**
 * Separate an array into two parts: one contains only elements that < pivot,
 * another one contains only elements that >= pivot
 * data: the array of elements
 * pivot : the pivot
 * N: the size of data
 * compare: the comparator
 *
 * Return the index of the separator
 */
template<typename DataType>
int partition(DataType * data, DataType pivot, size_t N, int (*compare)(const DataType *, const DataType *));
/**
 * A bubble sorting function
 * data: the array of elements
 * N: the size of data
 * compare: the comparator
 */
template<typename DataType>
void bbsort(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *));
/**
 * Pivot selection by finding median.
 * We should be careful about the pivot,
 * choose a value that exceeds 30~70% of the members
 * data: the array of elements
 * N: the size of data
 * compare: the comparator
 *
 * This function will return a pivot with the specified type
 */
template<typename DataType>
DataType find_median(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *));
}
/**
 * Simple pivot selection
 */
template<typename DataType>
DataType * simple_find_median(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *));


#endif /* UTILITIES_H_ */
