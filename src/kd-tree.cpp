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
 * @param verbose Just for debugging
 * @return the distance between two KDNode if no error occurs, otherwise return DBL_MAX
 */
float kd_distance(
		KDNode<float> * _a,
		KDNode<float> * _b,
		bool verbose) {
	float * a = _a->get_data();
	float * b = _b->get_data();
	size_t N = _a->size();
	if(N != _b->size()) {
		if(verbose) {
			cerr << "Dimension are different" << endl;
			cerr << "a has " << N << " dimensions" << endl;
			cerr << "b has " << _b->size() << " dimensions" << endl;
		}
		exit(1);
	}
	return distance(a,b,N);
}

/**
 * A comparator
 * @param _a, _b two float numbers
 */
int compare_float(
		const float * _a,
		const float * _b) {
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
size_t find_median(
		float ** data,
		size_t M,
		size_t N,
		size_t id,
		bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return -1;
	}
	id = id % N;
	float * arr = (float *)::operator new(M * sizeof(float));
	for(size_t i = 0; i < M; i++) {
		arr[i] = data[i][id];
	}

	size_t res = quick_select_k_id(arr,M,M >> 1,compare_float);
	::delete arr;
	arr = nullptr;
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
void make_balanced_tree(
		KDNode<float> *& root,
		float ** data,
		size_t M,
		size_t N,
		size_t level,
		size_t base,
		bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return;
	}
	size_t id = find_median(data,M,N,level,verbose);

	kd_insert<float>(root,data[id],N,0,id+base,verbose);
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
void make_random_tree(
		KDNode<float> *& root,
		float ** data,
		size_t M,
		size_t N,
		size_t base,
		bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return;
	}
	size_t id = M >> 1;

	kd_insert<float>(root,data[id],N,0,id+base,verbose);
	make_random_tree(root,data,id,N,base,verbose);
	if(id < M - 1) make_random_tree(root,&data[id+1],
			M - id - 1,N,base+id+1,verbose);
}

/**
 * Search for the nearest neighbor in the kd-tree
 * @param root the root node of the tree
 * @param query the data of the query
 * @param result the nearest neighbor
 * @param best_dist the best distance
 * @param N the size of the input
 * @param level the cut-plane level
 * @param visited (for debugging) to detect how many nodes are visited
 * @param verbose for debugging
 */
void nn_search(
		KDNode<float> * root,
		KDNode<float> * query,
		KDNode<float> *& result,
		float& best_dist,
		size_t N,
		size_t level,
		size_t& visited,
		bool verbose) {
	if(best_dist == 0.0) return;
	if(root == nullptr || root->size() != N) {
		if(verbose) {
			cout << "Reached a leaf" << endl;
			if(root) cout << "ID:" << root->id << endl;
		}
		return;
	}
	if(verbose)
		cout << "Visting node " << root->id << " with best is " << best_dist << endl;

	float d = kd_distance(root,query,verbose);
	float d1 = root->at(level) - query->at(level);
	visited++;

	if(result == nullptr || d < best_dist) {
		best_dist = d;
		result = root;
	}

	size_t l = (level + 1) % N;
	if(d1 >= 0.0) {
		nn_search(root->left,query,result,best_dist,N,l,visited,verbose);
	} else {
		nn_search(root->right,query,result,best_dist,N,l,visited,verbose);
	}

	if(fabs(d1) >= best_dist) return;
	if(verbose)
		cout << "Right branch of node " << root->id << endl;

	if(d1 >= 0.0) {
		nn_search(root->right,query,result,best_dist,N,l,visited,verbose);
	} else {
		nn_search(root->left,query,result,best_dist,N,l,visited,verbose);
	}
}

/**
 * Search for the approximate nearest neighbor in the kd-tree
 * @param root the root node of the tree
 * @param query the data of the query
 * @param result the nearest neighbor
 * @param best_dist the best distance
 * @param alpha the parameter that set the quality of nearest neighbor
 * @param N the size of the input
 * @param level the cut-plane level
 * @param visited (for debugging) to detect how many nodes are visited
 * @param verbose for debugging
 */
void ann_search(
		KDNode<float> * root,
		KDNode<float> * query,
		KDNode<float> *& result,
		float& best_dist,
		float alpha,
		size_t N,
		size_t level,
		size_t& visited,
		bool verbose) {
	if(best_dist == 0.0) return;
	if(root == nullptr || root->size() != N) {
		if(verbose) {
			cout << "Reached a leaf" << endl;
			if(root) cout << "ID:" << root->id << endl;
		}
		return;
	}
	if(verbose)
		cout << "Visting node " << root->id << endl;

	float d = kd_distance(root,query,verbose);
	float d1 = root->at(level) - query->at(level);
	visited++;

	if(result == nullptr || d < best_dist) {
		best_dist = d;
		result = root;
	}

	level = (level + 1) % N;
	if(d1 >= 0.0) {
		nn_search(root->left,query,result,best_dist,N,level,visited,verbose);
	} else {
		nn_search(root->right,query,result,best_dist,N,level,visited,verbose);
	}

	if(fabs(d1) * alpha > best_dist) return;

	if(d1 >= 0.0) {
		nn_search(root->right,query,result,best_dist,N,level,visited,verbose);
	} else {
		nn_search(root->left,query,result,best_dist,N,level,visited,verbose);
	}
}

/**
 * A linear solution for NNS
 * @param data the database
 * @param query the input query
 * @param best the index of the NNS
 * @param best_dist the best distance
 * @param N the size of database
 * @param d the dimensions
 * @param verbose for debugging
 */
void linear_search(
		float ** data,
		float * query,
		size_t& best,
		float& best_dist,
		size_t N,
		size_t d,
		bool verbose) {
	if(N <= 0 || d <= 0) {
		if(verbose)
			cerr << "Wrong size" << endl;
		return;
	}
	best = 0;
	best_dist = SimpleCluster::distance(query,data[0],d);
	float tmp = 0.0;
	for(size_t i = 1; i < N; i++) {
		tmp = SimpleCluster::distance(query,data[i],d);
		if(tmp < best_dist) {
			best_dist = tmp;
			best = i;
		}
	}
}
}


