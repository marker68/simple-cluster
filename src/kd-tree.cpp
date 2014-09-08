/*
 * kd-tree.cpp
 *
 *  Created on: 2014/09/08
 *      Author: Nguyen Anh Tuan<t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <stdlib.h>
#include "utilities.h"
#include "kd-tree.h"

using namespace std;

namespace SimpleCluster {
d_vector KDNode::get_data() const {
	return data;
}

void KDNode::set_data(d_vector _data) {
	data = _data;
}

size_t KDNode::get_compare_dim() const {
	return compare_dim;
}

void KDNode::set_compare_dim(size_t _cd) {
	compare_dim = _cd;
}

double distance(KDNode * _a, KDNode * _b) {
	d_vector _av = _a->get_data();
	d_vector _bv = _b->get_data();
	if(_av.size() != _bv.size()) {
		cerr << "Your vectors have some problems." << endl;
		exit(1);
	}
	size_t d = _av.size();
	return distance(_av,_bv,d);
}

int node_compare(const KDNode * _a, const KDNode * _b) {
	d_vector _av = _a->get_data();
	d_vector _bv = _b->get_data();
	size_t id = _a->get_compare_dim();
	if(_av.size() != _bv.size() || _av.size() < id || _bv.size() < id) {
		cerr << "Your vectors have some problems." << endl;
		exit(1);
	}
	if(_av[id] < _bv[id]) return -1;
	else if(_av[id] > _bv[id]) return 1;
	else return 0;
}

KDNode * make_tree(KDNode * data, size_t N, size_t id, size_t d) {
	KDNode * root;
	if(N <= 0)
		return(0);

	size_t i;
	for(i = 0; i < N; i++) {
		data[i].set_compare_dim(id);
	}

	if((root = simple_find_median<KDNode>(data,N,node_compare))) {
		id = (id + 1) % d;
		root->left = make_tree(data, root - data, id, d);
		root->right = make_tree(root + 1, data + N - (root + 1), id, d);
	}

	return root;
}

void find_nearest(KDNode * root, KDNode * node, KDNode ** best,
		double * best_dist, size_t id, size_t d) {
	double tmp, tmp2, tmp3;
	if(!root) return;
	tmp = distance(root,node);
	d_vector _rv = root->get_data();
	d_vector _nv = node->get_data();
	tmp2 = _rv[id] - _nv[id];
	tmp3 = tmp2 * tmp2;

	if (!*best || tmp < *best_dist) {
		*best_dist = tmp;
		*best = root;
	}

	if(!*best_dist) return;
	if(++id >= d) id = 0;
	find_nearest(tmp2 > 0?root->left:root->right,
			node,best,best_dist,id,d);
	if(tmp3 >= *best_dist) return;
	find_nearest(tmp2 > 0?root->right:root->left,
			node,best,best_dist,id,d);
}
}


