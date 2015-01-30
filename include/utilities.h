/*
 *  SIMPLE CLUSTERS: A simple library for clustering works.
 *  Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  utilities.h
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <vector>
#include <exception>
#include <cmath>
#include <cstring>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#define SET_THREAD_NUM omp_set_num_threads(n_thread)
#else
#define SET_THREAD_NUM 0 // disable multi-thread
#endif

using namespace std;

/**
 * The main namespace
 */
namespace SimpleCluster {

/**
 * Types of distances
 */
enum class DistanceType {
	NORM_L1,
	NORM_L2,
	HAMMING
};

void check_env();
unsigned long get_millisecond_time();
void print_vector(
		float **,
		int,
		int);
/**
 * Calculate the L1-metric distance
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType>
inline double distance_l1(
		DataType * x,
		DataType * y,
		int d) {
	int i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += fabs(tmp);
	}

	return dis;
}

/**
 * Calculate the L1-metric distance with two different data type
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType1, typename DataType2>
inline double distance_l1(
		DataType1 * x,
		DataType2 * y,
		int d) {
	int i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = static_cast<double>(x[i])
													- static_cast<double>(y[i]);
		dis += fabs(tmp);
	}

	return dis;
}

/**
 * Calculate the L1-metric distance between two vectors with multi-threading
 * @param x
 * @param y
 * @param d
 * @param n_thread number of threads
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType>
inline double distance_l1_thread(
		DataType * x,
		DataType * y,
		int d,
		int n_thread) {
#ifdef _WIN32
	int i;
#else
	int i;
#endif
	int bs = d / n_thread, st;
	double dis;
	dis = 0.0;
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for
		for(i = 0; i < n_thread; i++) {
			st = bs * i;
			dis += distance_l1<DataType>(x + st,y + st,bs);
		}
	}

	return sqrt(dis);
}

/**
 * Calculate the L2-metric distance
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType>
inline double distance_l2(
		DataType * x,
		DataType * y,
		int d) {
	int i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
	}

	return sqrt(dis);
}

/**
 * Calculate the L2-metric distance with 2 different data types
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType1, typename DataType2>
inline double distance_l2(
		DataType1 * x,
		DataType2 * y,
		int d) {
	int i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = static_cast<double>(x[i])
													- static_cast<double>(y[i]);
		dis += tmp * tmp;
	}
	return sqrt(dis);
}

/**
 * Calculate the L2-metric distance
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType>
inline double distance_l2_square(
		DataType * x,
		DataType * y,
		int d) {
	int i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
	}

	return dis;
}

/**
 * Calculate the L2-metric distance with 2 different data types
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType1, typename DataType2>
inline double distance_l2_square(
		DataType1 * x,
		DataType2 * y,
		int d) {
	int i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = static_cast<float>(x[i])
													- static_cast<float>(y[i]);
		dis += tmp * tmp;
	}

	return dis;
}

/**
 * Calculate the L2-metric distance between two vectors with multi-threading
 * @param x
 * @param y
 * @param d
 * @param n_thread number of threads
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType>
inline double distance_l2_thread(
		DataType * x,
		DataType * y,
		int d,
		int n_thread) {
#ifdef _WIN32
	int i;
#else
	int i;
#endif
	int bs = d / n_thread, st;
	double dis;
	dis = 0.0;
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for
		for(i = 0; i < n_thread; i++) {
			st = bs * i;
			dis += distance_l2_square<DataType>(x + st,y + st,bs);
		}
	}

	return sqrt(dis);
}

/**
 * Calculate the L2-metric distance between two vectors with multi-threading
 * @param x
 * @param y
 * @param d
 * @param n_thread number of threads
 * @return the distance between x and y in d dimensional space
 */
template<typename DataType>
inline double distance_l2_square_thread(
		DataType * x,
		DataType * y,
		int d,
		int n_thread) {
#ifdef _WIN32
	int i;
#else
	int i;
#endif
	int bs = d / n_thread, st;
	double dis;
	dis = 0.0;
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for
		for(i = 0; i < n_thread; i++) {
			st = bs * i;
			dis += distance_l2_square<DataType>(x + st,y + st,bs);
		}
	}

	return dis;
}

/**
 * Initialize an 1-D array.
 * @param arr the input array
 * @param N the size of the input array
 * @return true if the array was initalized successfully, otherwise return false.
 */
template<typename DataType>
inline bool init_array(
		DataType *& arr,
		size_t N) {
	if(N <= 0)
		return false;
	arr = (DataType *)::operator  new(N * sizeof(DataType));
	return true;
}

/**
 * Initialize a 2-D array.
 * @param M the size of the input array
 * @param N the size of the input array
 * @return true if the array was initialized successfully, otherwise return false.
 */
template<typename DataType>
inline bool init_array_2(
		DataType **& arr,
		int M,
		int N) {
	if(M <= 0 || N <= 0)
		return false;
	try {
		arr = (DataType **)::operator  new(M * sizeof(DataType *));
		for(int i = 0; i < M; i++) {
			::new(arr + i) DataType *;
			arr[i] = (DataType *)::operator new(N * sizeof(DataType));
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		return false;
	}

	return true;
}

/**
 * Copy an array.
 * @param from The input array
 * @param to the copy destination
 * @param N the number of the elements to be copied
 * @return true if copy was succeeded, otherwise return false.
 */
template<typename DataType>
inline bool copy_array(
		DataType * from,
		DataType *& to,
		int N) {
	if(from == nullptr || to == nullptr)
		return false;
	memcpy(to,from,N*sizeof(DataType));
	return true;
}


/**
 * Swap two elements in an array
 * @param data the array of elements
 * @param m,n the indices of elements to be swapped
 * @param N the size of data
 */
template<typename DataType>
inline void swap(
		DataType * data,
		int m,
		int n,
		int N) {
	if(m < 0 || n < 0 || m >= N || n >= N) {
		cerr << "Out of bounds!" << endl;
		exit(1);
	}
	DataType c = data[m];
	data[m] = data[n];
	data[n] = c;
}

/**
 * Separate an array into two parts: one contains only elements that < pivot,
 * another one contains only elements that >= pivot
 * @param data the array of elements
 * @param pivot the pivot
 * @param N the size of data
 * @param compare the comparator
 *
 * @return the index of the separator
 */
template<typename DataType>
inline int partition(
		DataType * data,
		DataType pivot,
		int N,
		int (*compare)(const DataType *, const DataType *)) {
	int i=0, j=N-1;
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
 * @param data the array of elements
 * @param N the size of data
 * @param compare the comparator
 */
template<typename DataType>
inline void bbsort(
		DataType * data,
		int N,
		int (*compare)(const DataType *, const DataType *)) {
	if(N < 2) return;
	bool swapped = true;
	int i, j = 0;
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
 * Select k+1-th smallest member of an array : QuickSelect by Hoare.
 * This one take O(n) time in average but O(n^2) in the worst case.
 * QuickSelect might be slower than StableSelect in worst cases,
 * but in most cases, QuickSelect outperformed StableSelect.
 * BE CAREFUL! This function will change the order of the input data.
 * Make sure you stored your input data somewhere else.
 * @param data the input data
 * @param N the size of the input
 * @param k the index
 * @param compare the comparator
 * @return the index of the k-th smallest element
 */
template<typename DataType>
inline int quick_select_k_id(
		DataType * data,
		int N,
		int k,
		int (*compare)(const DataType*, const DataType*)) {
	if(k >= N) {
		cerr << "Out of bounds!\n" << endl;
		exit(1);
	}

	if(N <= 8) {
		bbsort(data,N,*compare);
		return k;
	}

	// First, choose an appropriate pivot
	DataType pivot = data[N >> 1];
	// Choose the pivot as the median of left, right and (left+right)/2 elements
	if(((*compare)(&pivot,&data[0]) <= 0 && (*compare)(&data[0],&data[N-1]) <= 0)
			|| ((*compare)(&pivot,&data[0]) >= 0 && (*compare)(&data[0],&data[N-1]) >= 0))
		pivot = data[0];
	if(((*compare)(&pivot,&data[N-1]) <= 0 && (*compare)(&data[N-1],&data[0]) <= 0)
			|| ((*compare)(&pivot,&data[N-1]) >= 0 && (*compare)(&data[N-1],&data[0]) >= 0))
		pivot = data[N-1];

	// Then partition the array into two parts based on the pivot value.
	// The left part will contains members that are less than pivot,
	// the right part contains member that are greater than or equal pivot.
	int p = partition(data,pivot,N,*compare);
	if(p == 0) {
		bbsort(data,N,*compare);
		return k;
	}

	if(k < p) return quick_select_k_id(data,p,k,*compare);
	else if(k == p) return p;
	else return quick_select_k_id(&data[p],N-p,k-p,*compare);
}

/**
 * Select k+1-th smallest member of an array : QuickSelect by Hoare.
 * This one take O(n) time in average but O(n^2) in the worst case.
 * QuickSelect might be slower than StableSelect in worst cases,
 * but in most cases, QuickSelect outperformed StableSelect.
 * BE CAREFUL! This function will change the order of the input data.
 * Make sure you stored your input data somewhere else.
 * @param data the input data
 * @param N the size of the input
 * @param k the index
 * @param compare the comparator
 * @return the value of the k-th smallest element
 */
template<typename DataType>
inline DataType quick_select_k(
		DataType * data,
		int N,
		int k,
		int (*compare)(const DataType*, const DataType*)) {
	if(k >= N) {
		cerr << "Out of bounds!\n" << endl;
		exit(1);
	}

	if(N <= 8) {
		bbsort(data,N,*compare);
		return data[k];
	}

	// First, choose an appropriate pivot
	DataType pivot = data[N >> 1];
	// Choose the pivot as the median of left, right and (left+right)/2 elements
	if(((*compare)(&pivot,&data[0]) <= 0 && (*compare)(&data[0],&data[N-1]) <= 0)
			|| ((*compare)(&pivot,&data[0]) >= 0 && (*compare)(&data[0],&data[N-1]) >= 0))
		pivot = data[0];
	if(((*compare)(&pivot,&data[N-1]) <= 0 && (*compare)(&data[N-1],&data[0]) <= 0)
			|| ((*compare)(&pivot,&data[N-1]) >= 0 && (*compare)(&data[N-1],&data[0]) >= 0))
		pivot = data[N-1];

	// Then partition the array into two parts based on the pivot value.
	// The left part will contains members that are less than pivot,
	// the right part contains member that are greater than or equal pivot.
	int p = partition(data,pivot,N,*compare);
	if(p == 0) {
		bbsort(data,N,*compare);
		return data[k];
	}

	if(k < p) return quick_select_k(data,p,k,*compare);
	else if(k == p) return data[p];
	else return quick_select_k(&data[p],N-p,k-p,*compare);
}
}

#endif /* UTILITIES_H_ */
