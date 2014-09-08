/*
 * kd-tree.h
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
typedef vector<double> d_vector;

typedef struct kd_node {
	d_vector data;
	struct kd_node * left, * right;
} KDNode;

double distance(KDNode *, KDNode *);
KDNode * make_tree(KDNode *, size_t, size_t, size_t);
void find_nearest(KDNode *, KDNode *, KDNode **, double *, size_t, size_t);
}

#endif /* KD_TREE_H_ */
