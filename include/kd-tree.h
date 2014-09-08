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
//class KDNode {
//private:
//	d_vector data;
//	size_t id;
//public:
//	KDNode *left, *right;
//	KDNode() {
//		left = right = NULL;
//		id = 0;
//	}
//	virtual ~KDNode() {}
//	d_vector get_data() const;
//	void set_data(d_vector);
//	size_t get_id() const;
//	void set_id(size_t);
//};

double distance(KDNode *, KDNode *);
KDNode * make_tree(KDNode *, size_t, size_t, size_t);
void find_nearest(KDNode *, KDNode *, KDNode **, double *, size_t, size_t);
}

#endif /* KD_TREE_H_ */
