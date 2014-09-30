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

using namespace std;

namespace SimpleCluster {

/**
 * KD-Tree node class
 */
template<typename DataType>
class KDNode {
private:
	DataType * data;
	size_t dimension;
public:
	size_t id;
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
	KDNode(size_t _d) {
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
			size_t pos) {
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
	DataType at(size_t _id) const{
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
	size_t size() const{
		return dimension;
	}
};

double kd_distance(
		KDNode<double> *,
		KDNode<double> *,
		bool);
size_t find_median(
		double **,
		size_t,
		size_t,
		size_t,
		bool);
void make_balanced_tree(
		KDNode<double> *&,
		double **,
		size_t,
		size_t,
		size_t,
		size_t,
		bool);
void make_random_tree(
		KDNode<double> *&,
		double **,
		size_t,
		size_t,
		size_t,
		bool);
void nn_search(
		KDNode<double> *,
		KDNode<double> *,
		KDNode<double> *&,
		double&,
		size_t,
		size_t,
		size_t&,
		bool);
void ann_search(
		KDNode<double> *,
		KDNode<double> *,
		KDNode<double> *&,
		double&,
		double,
		size_t,
		size_t,
		size_t&,
		bool);
void linear_search(
		double **,
		double *,
		size_t&,
		double&,
		size_t,
		size_t,
		bool);

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
void kd_insert(
		KDNode<DataType> *& root,
		DataType * _data,
		size_t N,
		size_t level,
		size_t _id,
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
		kd_insert(root->right,_data,N,(level+1)%N,_id,verbose);
	} else if(root->at(level) >= _data[level]) {
		kd_insert(root->left,_data,N,(level+1)%N,_id,verbose);
	}
}

/**
 * Traveling in the kd-tree
 * @param root the root node of the tree
 * @param N the size of the vector
 * @param level the cut-plane level
 */
template<typename DataType>
void kd_travel(
		KDNode<DataType> * root,
		size_t N,
		size_t level) {
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
