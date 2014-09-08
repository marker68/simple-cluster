/*
 * utilities.cpp
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <sys/time.h>
#include "utilities.h"

using namespace std;

namespace SimpleCluster {
/**
 * Calculate the distance between two vectors
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
double distance(d_vector x, d_vector y, size_t d) {
	if(x.size() < d || y.size() < d) {
		cerr << "distance: Your vector have not enough dimensions!" << endl;
		cerr << "x has " << x.size() << " dimensions!" << endl;
		cerr << "y has " << y.size() << " dimensions!" << endl;
		exit(1);
	}
	size_t i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
	}

	return sqrt(dis);
}

/**
 * Calculate the square of distance between two vectors
 * @param x
 * @param y
 * @param d
 * @return the square of distance between x and y in d dimensional space
 */
double distance_square(d_vector x, d_vector y, size_t d) {
	if(x.size() < d || y.size() < d) {
		cerr << "distance_square: Your vector have not enough dimensions!" << endl;
		cerr << "x has " << x.size() << " dimensions!" << endl;
		cerr << "y has " << y.size() << " dimensions!" << endl;
		exit(1);
	}
	size_t i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
	}

	return dis;
}

/**
 * Calculate the mean of a cluster
 * @param data
 * @param index
 * @param d
 * @return the mean posize_t of a cluster
 */
d_vector mean_vector(vector<d_vector> data, i_vector index, size_t d, d_vector centroid) {
	size_t i, j, size  = index.size();
	if(size <= 0) {
		return centroid;
	}
	d_vector d_tmp;
	double tmp[d];

	for(i = 0; i < d; i++)
		tmp[i] = 0.0;

	for(i = 0; i < size; i++) {
		j = index[i];
		d_tmp = data[j];
		for(j = 0; j < d; j++) {
			tmp[j] += d_tmp[j];
		}
	}

	for(i = 0; i < d; i++)
		tmp[i] /= static_cast<double>(size);
	d_tmp.clear();
	for(i = 0; i < d; i++)
		d_tmp.push_back(tmp[i]);
	return d_tmp;
}

/**
 * Get system time in milliseconds
 */
unsigned long get_millisecond_time() {
	struct timeval tv;
	if(gettimeofday(&tv, NULL) != 0) return 0;
	return (unsigned long)((tv.tv_sec * 1000ul) + (tv.tv_usec / 1000ul));
}

/**
 * Utilities for printing vector
 */
void print_vector(vector<d_vector> data, size_t d) {
	size_t i, j, size = data.size();
	for(i = 0; i < size; i++) {
		cout << "Data point " << i << ":";
		for(j = 0; j < d; j++)
			cout << data[i][j] << " ";
		cout << endl;
	}
}

/**
 * Swap two elements in an array
 * data: the array of elements
 * m,n : the indices of elements to be swapped
 * N: the size of data
 */
template<typename DataType>
void swap(DataType * data, size_t m, size_t n, size_t N) {
	if(m < 0 || n < 0 || m >= N || n >= N) {
		perror("Out of bounds!\n");
		exit(1);
	}
	DataType c = data[m];
	data[m] = data[n];
	data[n] = c;
}

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
int partition(DataType * data, DataType pivot, size_t N, int (*compare)(const DataType *, const DataType *)) {
	size_t i=0, j=N-1;
	while(i <= j) {
		while((*compare)(&data[i],&pivot) < 0)
			i++;
		while((*compare)(&data[j],&pivot) >= 0)
			j--;
		if(i <= j) {
			swap(data,i,j,N);
			i++;
			j--;
		}
	}

	return i;
}

/**
 * A bubble sorting function
 * data: the array of elements
 * N: the size of data
 * compare: the comparator
 */
template<typename DataType>
void bbsort(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *)) {
	if(N < 2) return;
	bool swapped = true;
	size_t i, j = 0;
	while(swapped) {
		swapped = false;
		j++;
		for(i=0;i<N-j;i++) {
			if((*compare)(&data[i],&data[i+1]) > 0) {
				swap(data,i,i+1,N);
				swapped = true;
			}
		}
	}
}

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
DataType find_median(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *)) {
	if(N <= 5) {
		bbsort(data,N,*compare);
		return data[N >> 1];
	}
	size_t num_of_medians = N % 5 == 0 ? N/5 : (N/5+1);
	size_t i;
	for(i=0;i<num_of_medians;i++) {
		int size = 0;
		if(i < num_of_medians - 1) size = 5;
		else size = N - 5 * i;
		bbsort(&data[5 * i], size, *compare);
		swap(data,i,5 * i + size / 2, N);
	}
	return find_median(data,num_of_medians,*compare);
}

/**
 * Simple pivot selection
 */
template<typename DataType>
DataType * simple_find_median(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *)) {
	// First, choose an appropriate pivot
	DataType pivot = data[N >> 1];
	// Choose the pivot as the median of left, right and (left+right)/2 elements
	if(((*compare)(&pivot,&data[0]) <= 0 && (*compare)(&data[0],&data[N-1]) <= 0)
			|| ((*compare)(&pivot,&data[0]) >= 0 && (*compare)(&data[0],&data[N-1]) >= 0))
		pivot = data[0];
	if(((*compare)(&pivot,&data[N-1]) <= 0 && (*compare)(&data[N-1],&data[0]) <= 0)
			|| ((*compare)(&pivot,&data[N-1]) >= 0 && (*compare)(&data[N-1],&data[0]) >= 0))
		pivot = data[N-1];

	return &pivot;
}
}


