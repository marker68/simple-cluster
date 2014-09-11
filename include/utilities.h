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
double distance(double *, double *, size_t);
double distance_square(double *, double *, size_t);
double * mean_vector(double **, const int *, size_t, size_t, double *);
/**
 * Get system time in milliseconds
 */
unsigned long get_millisecond_time();
void print_vector(double **, size_t, size_t);
/**
 * Swap two elements in an array
 * data: the array of elements
 * m,n : the indices of elements to be swapped
 * N: the size of data
 */
template<typename DataType>
void swap(DataType * data, int m, int n, size_t N) {
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
 * data: the array of elements
 * N: the size of data
 * compare: the comparator
 */
template<typename DataType>
void bbsort(DataType * data, size_t N, int (*compare)(const DataType *, const DataType *)) {
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
 * Select k+1-th smallest member of an array : QuickSelect by Hoare
 * This one take O(n) time in average but O(n^2) in the worst case.
 * QuickSelect might be slower than StableSelect in worst cases,
 * but in most cases, QuickSelect outperformed StableSelect.
 * BE CAREFUL! This function will change the order of the input data.
 * Make sure you stored your input data somewhere else.
 */
template<typename DataType>
size_t quick_select_k(DataType * data, size_t N, size_t k, int (*compare)(const DataType*, const DataType*)) {
  if(k >= N) {
      fprintf(stderr, "Out of bounds!\n");
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

  if(k < p) return quick_select_k(data,p,k,*compare);
  else if(k == p) return p;
  else return quick_select_k(&data[p],N-p,k-p,*compare);
}
}

#endif /* UTILITIES_H_ */
