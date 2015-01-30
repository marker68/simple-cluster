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
 *  kd-tree.h
 *
 *  Created on: 2014/09/08
 *      Author: Nguyen Anh Tuan<t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef KD_TREE_H_
#define KD_TREE_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include "utilities.h"

using namespace std;

namespace SimpleCluster {

/**
 * KD-Tree node class
 */
template<typename DataType>
class KDNode {
private:
	DataType * data;
	int dimension;
public:
	int id;
	KDNode<DataType> * left, * right;

	/**
	 * A constructor to initialize data from another node
	 * @param other another KDNode
	 */
	KDNode(const KDNode<DataType>& other) {
		dimension = other.size();
		data = other.get_data();
		left = other.left;
		right = other.right;
		id = other.id;
	}

	/**
	 * The default constructor
	 */
	KDNode(int _d) {
		dimension = _d;
		data = (DataType *)::operator new(_d * sizeof(DataType));
		id = 0;
		left = right = nullptr;
	}

	/**
	 * Copy function
	 * @param other another KDNode
	 */
	KDNode& operator= (const KDNode<DataType>& other) {
		KDNode<DataType> tmp = ::new KDNode<DataType>(other.size());
		tmp.add_data(other.get_data());
		tmp.left = other.left;
		tmp.right = other.right;
		tmp.id = other.id;
		return tmp;
	}

	/**
	 * The destructor
	 */
	virtual ~KDNode() {
	}

	/**
	 * Insert a new data into the vector of the node
	 * @param _d The data to be inserted
	 */
	void add_data(
			DataType _d,
			int pos) {
		if(pos >= 0 && pos < dimension)
			data[pos] = _d;
	}

	/**
	 * Add an array of data
	 * @param _d data to be inserted
	 * @param N the size of the input
	 */
	void add_data(DataType * _d) {
		data = _d;
	}

	/**
	 * get a component of the vector
	 * @param _id index of the component
	 */
	DataType at(int _id) const{
		return data[_id];
	}

	/**
	 * Get all data as array
	 */
	DataType * get_data() const{
		return data;
	}

	/**
	 * Get the size or the dimensions of the vector
	 */
	int size() const{
		return dimension;
	}
};

/**
 * Insert a node into the kd-tree
 * @param root the root node of the tree
 * @param data the data of the node to be inserted
 * @param N the size of the vector
 * @param level the cut-plane level
 * @param _id the index of the new node
 * @verbose true to print verbose. Just for debugging.
 */
template<typename DataType>
inline void kd_insert(
		KDNode<DataType> *& root,
		DataType * _data,
		int N,
		int level,
		int _id,
		bool verbose) {
	// Pre-check
	if(N <= 0) return;
	level %= N;
	if(root == nullptr) {
		root = ::new KDNode<DataType>(N);
		root->add_data(_data);
		root->id = _id;
		return;
	}

	if(root->id == _id) {
		if(verbose)
			cerr << "No new node will be inserted! "
			"Constraint: Id field is an unique field! "
			"Conflicted at id:" << _id << endl;
		return;
	}

	if(root->at(level) < _data[level]) {
		kd_insert<DataType>(root->right,_data,N,(level+1)%N,_id,verbose);
	} else if(root->at(level) >= _data[level]) {
		kd_insert<DataType>(root->left,_data,N,(level+1)%N,_id,verbose);
	}
}

/**
 * Calculate the distances between two KDNode
 * @param _a, _b the input KDNode
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param verbose Just for debugging
 * @return the distance between two KDNode if no error occurs, otherwise return DBL_MAX
 */
template<typename DataType>
inline double kd_distance(
		KDNode<DataType> * _a,
		KDNode<DataType> * _b,
		DistanceType d_type,
		bool verbose) {
	DataType * a = _a->get_data();
	DataType * b = _b->get_data();
	int N = _a->size();
	if(N != _b->size()) {
		if(verbose) {
			cerr << "Dimension are different" << endl;
			cerr << "a has " << N << " dimensions" << endl;
			cerr << "b has " << _b->size() << " dimensions" << endl;
		}
		exit(1);
	}
	if(d_type == DistanceType::NORM_L1)
		return distance_l1<DataType>(a,b,N);
	return distance_l2<DataType>(a,b,N);
}

/**
 * A comparator
 * @param _a, _b two float numbers
 */
template<typename DataType>
inline int comparator(
		const DataType * _a,
		const DataType * _b) {
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
template<typename DataType>
inline int find_median(
		DataType ** data,
		int M,
		int N,
		int id,
		bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return -1;
	}
	id = id % N;
	DataType * arr = (DataType *)::operator new(M * sizeof(float));
	for(int i = 0; i < M; i++) {
		arr[i] = data[i][id];
	}

	int res = quick_select_k_id<DataType>(arr,M,M >> 1,comparator);
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
template<typename DataType>
inline void make_balanced_tree(
		KDNode<DataType> *& root,
		DataType ** data,
		int M,
		int N,
		int level,
		int base,
		bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return;
	}
	int id = find_median<DataType>(data,M,N,level,verbose);

	kd_insert<DataType>(root,data[id],N,0,id+base,verbose);
	make_balanced_tree<DataType>(root,data,id,N,(level+1)%N,base,verbose);
	if(id < M - 1) make_balanced_tree<DataType>(root,&data[id+1],
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
template<typename DataType>
inline void make_random_tree(
		KDNode<DataType> *& root,
		DataType ** data,
		int M,
		int N,
		int base,
		bool verbose) {
	if(N <= 0 || M <= 0) {
		if(verbose)
			cerr << "No data" << endl;
		return;
	}
	int id = M >> 1;

	kd_insert<DataType>(root,data[id],N,0,id+base,verbose);
	make_random_tree<DataType>(root,data,id,N,base,verbose);
	if(id < M - 1) make_random_tree<DataType>(root,&data[id+1],
			M - id - 1,N,base+id+1,verbose);
}

/**
 * Search for the nearest neighbor in the kd-tree
 * @param root the root node of the tree
 * @param query the data of the query
 * @param result the nearest neighbor
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param best_dist the best distance
 * @param N the size of the input
 * @param level the cut-plane level
 * @param visited (for debugging) to detect how many nodes are visited
 * @param verbose for debugging
 */
template<typename DataType>
inline void nn_search(
		KDNode<DataType> * root,
		KDNode<DataType> * query,
		KDNode<DataType> *& result,
		DistanceType d_type,
		double& best_dist,
		int N,
		int level,
		int& visited,
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

	double d = kd_distance<DataType>(root,query,d_type,verbose);
	double d1 = static_cast<double>(root->at(level) - query->at(level));
	visited++;

	if(result == nullptr || d < best_dist) {
		best_dist = d;
		result = root;
	}

	int l = (level + 1) % N;
	if(d1 >= 0.0) {
		nn_search<DataType>(root->left,query,result,d_type,best_dist,N,l,visited,verbose);
	} else {
		nn_search<DataType>(root->right,query,result,d_type,best_dist,N,l,visited,verbose);
	}

	if(fabs(d1) >= best_dist) return;
	if(verbose)
		cout << "Right branch of node " << root->id << endl;

	if(d1 >= 0.0) {
		nn_search<DataType>(root->right,query,result,d_type,best_dist,N,l,visited,verbose);
	} else {
		nn_search<DataType>(root->left,query,result,d_type,best_dist,N,l,visited,verbose);
	}
}

/**
 * Search for the approximate nearest neighbor in the kd-tree
 * @param root the root node of the tree
 * @param query the data of the query
 * @param result the nearest neighbor
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param best_dist the best distance
 * @param alpha the parameter that set the quality of nearest neighbor
 * @param N the size of the input
 * @param level the cut-plane level
 * @param visited (for debugging) to detect how many nodes are visited
 * @param verbose for debugging
 */
template<typename DataType>
inline void ann_search(
		KDNode<DataType> * root,
		KDNode<DataType> * query,
		KDNode<DataType> *& result,
		DistanceType d_type,
		double& best_dist,
		double alpha,
		int N,
		int level,
		int& visited,
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

	double d = kd_distance<DataType>(root,query,d_type,verbose);
	double d1 = static_cast<double>(root->at(level) - query->at(level));
	visited++;

	if(result == nullptr || d < best_dist) {
		best_dist = d;
		result = root;
	}

	level = (level + 1) % N;
	if(d1 >= 0.0) {
		ann_search<DataType>(root->left,query,result,d_type,best_dist,alpha,N,level,visited,verbose);
	} else {
		ann_search<DataType>(root->right,query,result,d_type,best_dist,alpha,N,level,visited,verbose);
	}

	if(fabs(d1) * alpha > best_dist) return;

	if(d1 >= 0.0) {
		ann_search<DataType>(root->right,query,result,d_type,best_dist,alpha,N,level,visited,verbose);
	} else {
		ann_search<DataType>(root->left,query,result,d_type,best_dist,alpha,N,level,visited,verbose);
	}
}

/**
 * A linear solution for NNS
 * @param data the database
 * @param query the input query
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param best the index of the NNS
 * @param best_dist the best distance
 * @param N the size of database
 * @param d the dimensions
 * @param verbose for debugging
 */
template<typename DataType>
inline void linear_search(
		DataType ** data,
		DataType * query,
		DistanceType d_type,
		int& best,
		double& best_dist,
		int N,
		int d,
		bool verbose) {
	if(N <= 0 || d <= 0) {
		if(verbose)
			cerr << "Wrong size" << endl;
		return;
	}
	best = 0;
	if(d_type == DistanceType::NORM_L2)
		best_dist = distance_l2<DataType>(query,data[0],d);
	else if(d_type == DistanceType::NORM_L1)
		best_dist = distance_l1<DataType>(query,data[0],d);
	float tmp = 0.0;
	for(int i = 1; i < N; i++) {
		if(d_type == DistanceType::NORM_L2)
			tmp = distance_l2<DataType>(query,data[i],d);
		else if(d_type == DistanceType::NORM_L1)
			tmp = distance_l1<DataType>(query,data[i],d);
		if(tmp < best_dist) {
			best_dist = tmp;
			best = i;
		}
	}
}

/**
 * Traveling in the kd-tree
 * @param root the root node of the tree
 * @param N the size of the vector
 * @param level the cut-plane level
 */
template<typename DataType>
inline void kd_travel(
		KDNode<DataType> * root,
		int N,
		int level) {
	if(root == nullptr) {
		cout << "Reached a leaf at level " << level << endl;
		return;
	}
	cout << "visited a node " << root->id << " at level " << level << endl;
	cout << "this node has " << root->size() << " dimensions" << endl;
	kd_travel(root->left,N,level+1);
	kd_travel(root->right,N,level+1);
	return;
}
}

#endif /* KD_TREE_H_ */
