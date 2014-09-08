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

class KDNode {
private:
	d_vector data;
	size_t compare_dim;
public:
	KDNode *left, *right;
	KDNode() {
		left = right = NULL;
		compare_dim = 0;
	}
	virtual ~KDNode() {}
	d_vector get_data() const;
	void set_data(d_vector);
	size_t get_compare_dim() const;
	void set_compare_dim(size_t);
};

double distance(KDNode *, KDNode *);
KDNode * make_tree(KDNode *, size_t, size_t, size_t);
void find_nearest(KDNode *, KDNode *, KDNode **, double *, size_t, size_t);
}

#endif /* KD_TREE_H_ */
