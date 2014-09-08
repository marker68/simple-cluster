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
//d_vector KDNode::get_data() const {
//	return data;
//}
//
//void KDNode::set_data(d_vector _data) {
//	data = _data;
//}
//
//size_t KDNode::get_id() const {
//	return id;
//}
//
//void KDNode::set_id(size_t _id) {
//	id = _id;
//}

double distance(KDNode * _a, KDNode * _b) {
	d_vector _av = _a->data;
	d_vector _bv = _b->data;
	if(_av.size() != _bv.size()) {
		cerr << "Your vectors have some problems." << endl;
		exit(1);
	}
	size_t d = _av.size();
	return distance(_av,_bv,d);
}

KDNode * make_tree(KDNode * data, size_t N, size_t id, size_t d) {
	KDNode * root;
	if(N <= 1)
		return data;

	if((root = &data[N >> 1])) {
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
	d_vector _rv = root->data;
	d_vector _nv = node->data;
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


