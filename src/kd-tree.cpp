/*
 * kd-tree.cpp
 *
 *  Created on: 2014/09/08
 *      Author: Nguyen Anh Tuan<t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "utilities.h"
#include "kd-tree.h"

using namespace std;

namespace SimpleCluster {

double distance(KDNode * _a, KDNode * _b) {
	double * _av = _a->data;
	double * _bv = _b->data;
	if(_a->dim != _b->dim) {
		cerr << "Your vectors have some problems." << endl;
		cerr << "a has " << _a->dim << " dimensions" << endl;
		cerr << "b has " << _b->dim << " dimensions" << endl;
		exit(1);
	}
	size_t d = _a->dim;
	double tmp = 0.0, res = 0.0;
	for(size_t i = 0; i < d; i++) {
		tmp = _av[i] - _bv[i];
		res += tmp * tmp;
	}

	return sqrt(res);
}

KDNode * make_tree(KDNode * data, size_t N, size_t id, size_t d) {
	KDNode * root;
	if(N == 1)
		return data;
	if(N < 1)
		return NULL;

	if((root = &data[N >> 1])) {
		id = (id + 1) % d;
		root->left = make_tree(data, root - data, id, d);
		root->right = make_tree(root + 1, data + N - (root + 1), id, d);
	}

	return root;
}

void find_nearest(KDNode * root, KDNode * node, KDNode ** best,
		double * best_dist, size_t id) {
	double tmp, tmp2, tmp3;
	if(root == NULL) return;
	tmp = distance(root,node);
	double * _rv = root->data;
	double * _nv = node->data;
	size_t d = node->dim;

	tmp2 = _rv[id] - _nv[id];
	tmp3 = tmp2 * tmp2;

	if (!*best || tmp < *best_dist) {
		*best_dist = tmp;
		*best = root;
	}

	if(!*best_dist) return;
	if(++id >= d) id = 0;
	find_nearest(tmp2 > 0?root->left:root->right,
			node,best,best_dist,id);
	if(tmp3 >= *best_dist) return;
	find_nearest(tmp2 > 0?root->right:root->left,
			node,best,best_dist,id);
}
}


