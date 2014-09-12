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
 *  kd-tree.cpp
 *
 *  Created on: 2014/09/08
 *      Author: Nguyen Anh Tuan<t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "utilities.h"
#include "kd-tree.h"

using namespace std;

namespace SimpleCluster {

/**
 * Calculate the distances between two KDNode
 * @param _a, _b the input KDNode
 * @param N the size of the vector
 * @param verbose Just for debugging
 * @return the distance between two KDNode if no error occurs, otherwise return DBL_MAX
 */
double kd_distance(KDNode<double> * _a, KDNode<double> * _b, size_t N, bool verbose) {
	if(_a->get_size() != _b->get_size()) {
		if(verbose) {
			cerr << "Your vector has some problem" << endl;
			cerr << "a has " << _a->get_size() << " dimensions" << endl;
			cerr << "b has " << _b->get_size() << " dimensions" << endl;
		}
		return DBL_MAX;
	}
	double tmp = 0.0, res = 0.0;
	for(size_t i = 0; i < N; i++) {
		tmp = _a->get_data_at(i) - _b->get_data_at(i);
		res += tmp * tmp;
	}

	return sqrt(res);
}

/**
 * A comparator
 * @param _a, _b two float numbers
 */
int compare_double(const double * _a, const double * _b) {
	if(*_a > *_b) return 1;
	else if(*_a < *_b) return -1;
	else return 0;
}

/**
 * Find the median of the input vector data
 * @param data the input data
 * @param M,N the size of the input
 * @param id the index of the component to find median
 * @param verbose just for debugging
 * @return the index of the median
 */
size_t find_median(double ** data, size_t M, size_t N, size_t id, bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return -1;
	}
	id = id % N;
	double * arr = (double *)::operator new(M * sizeof(double));
	for(size_t i = 0; i < M; i++) {
		arr[i] = data[i][id];
	}

	size_t res = quick_select_k(arr,M,M >> 1,compare_double);
	::delete arr;
	return res;
}

/**
 * Create a balanced kd-tree
 * @param root the root node of the tree
 * @param data the data of the nodes to be inserted
 * @param M,N the size of the input
 * @param level the cut-plane level
 * @param base the base index to be added
 * @param verbose just for debugging
 */
void make_balanced_tree(KDNode<double> *& root, double ** data,
		size_t M, size_t N, size_t level, size_t base, bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return;
	}
	size_t id = find_median(data,M,N,level,verbose);

	kd_insert<double>(root,data[id],N,level,id+base,verbose);
	make_balanced_tree(root,data,id,N,(level+1)%N,base,verbose);
	if(id < M - 1) make_balanced_tree(root,&data[id+1],
			M - id - 1,N,(level+1)%N,base+id+1,verbose);
}

/**
 * Create a random kd-tree
 * @param root the root node of the tree
 * @param data the data of the nodes to be inserted
 * @param M,N the size of the input
 * @param level the cut-plane level
 * @param base the base index to be added
 * @param verbose just for debugging
 */
void make_random_tree(KDNode<double> *& root, double ** data,
		size_t M, size_t N, size_t level,size_t base, bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return;
	}
	size_t id = M >> 1;

	kd_insert<double>(root,data[id],N,level,id+base,verbose);
	make_random_tree(root,data,id,N,(level+1)%N,base,verbose);
	if(id < M - 1) make_random_tree(root,&data[id+1],
			M - id - 1,N,(level+1)%N,base+id+1,verbose);
}

/**
 * Search for the nearest neighbor in the kd-tree
 * @param root the root node of the tree
 * @param query the data of the query
 * @param result the nearest neighbor
 * @param best_dist the best distance
 * @param N the size of the input
 * @param level the cut-plane level
 * @param _id the index of the new node
 */
void nn_search(KDNode<double> * root, const double * query,
		KDNode<double> *& result,
		double& best_dist, size_t N, size_t level, bool verbose) {
	if(best_dist == 0.0) return;
	if(root == NULL || root->get_size() != N) {
		if(verbose) {
			cout << "Reached a leaf" << endl;
			if(root) cout << "ID:" << root->id << endl;
		}
		return;
	}
	if(verbose)
		cout << "Visiting node " << root->id << endl;
	KDNode<double> * tmp = ::new KDNode<double>;
	for(size_t i = 0; i < N; i++)
		tmp->add_data(query[i]);

	double d = kd_distance(root,tmp,N,verbose);
	double d1 = root->get_data_at(level) - query[level];

	if(result == NULL || d < best_dist) {
		best_dist = d;
		result = root;
		if(verbose) {
			cout << best_dist << endl;
			cout << result->id << endl;
		}
	}

	size_t l = (level + 1) % N;
	if(d1 >= 0.0) {
		nn_search(root->left,query,result,best_dist,N,l,verbose);
	} else {
		nn_search(root->right,query,result,best_dist,N,l,verbose);
	}

	if(fabs(d1) >= best_dist) return;

	if(d1 >= 0.0) {
		nn_search(root->right,query,result,best_dist,N,l,verbose);
	} else {
		nn_search(root->left,query,result,best_dist,N,l,verbose);
	}
	::delete tmp;
}
}


