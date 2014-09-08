/*
 * test_kmeans.cpp
 *
 *  Created on: 2014/09/06
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <random>
#include <gtest/gtest.h>
#include <math.h>
#include <stdlib.h>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "k-means.h"
#include "utilities.h"

using namespace std;
using namespace cv;
using namespace SimpleCluster;

/**
 * A converter that converts vector to Mat
 */
void convert_dvector_mat(d_vector _f, Mat& _t) {
	_t.release();
	d_vector::iterator it = _f.begin();
	d_vector::iterator ie = _f.end();
	while(it != ie) {
		_t.push_back(*it);
		++it;
	}
	_t.reshape(0,1);
}

/**
 * Another converter
 */
void convert_vector_mat(vector<d_vector> _f, Mat& _t) {
	_t.release();
	vector<d_vector>::iterator it = _f.begin();
	vector<d_vector>::iterator ie = _f.end();
	while(it != ie) {
		d_vector::iterator idt = (*it).begin();
		d_vector::iterator ide = (*it).end();
		while(idt != ide) {
			_t.push_back(*idt);
			++idt;
		}
		++it;
	}
	_t = _t.reshape(0,_f.size());
}

/**
 * Customized test case for testing
 */
class KmeansTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 1000;
		d = 128;
		k = 256;
		int i, j, scale;
		d_vector tmp;
		double x;

		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> real_dis(0.0, static_cast<double>(N));

		for(i = 0; i < N; i++) {
			for(j = 0; j < d; j++) {
				scale = real_dis(gen);
				x = sqrt(scale);
				tmp.push_back(x);
			}
			data.push_back(tmp);
			tmp.clear();
		}
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		data.clear();
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static vector<d_vector> data;
	static int N;
	static int d;
	static int k;
};

vector<d_vector> KmeansTest::data;
int KmeansTest::N;
int KmeansTest::d;
int KmeansTest::k;

TEST_F(KmeansTest, test1) {
	vector<d_vector> centroids, seeds;
	vector<i_vector> clusters;
	KmeansCriteria criteria = {1.0,10};
	cout << "Our method:" << endl;
	unsigned long int t1, t2;
	t1 = get_millisecond_time();
	simple_k_means(KmeansType::KMEANS_PLUS_SEEDS,N,k,criteria,d,data,centroids,clusters,seeds);
	t2 = get_millisecond_time();
	cout << "Finished in " << t2-t1 << "[ms]" << endl;
	cout << "The distortion is " << distortion(d,N,k,data,centroids,clusters) << endl;
	cout << "=======================" << endl;
	cout << "OpenCV method:" << endl;
	Mat _data;
	convert_vector_mat(data,_data);
	_data.convertTo(_data,CV_32F);
	Mat _labels, _centers(k, 1, _data.type());
	TermCriteria opencv_criteria {CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1.0};
	t1 = get_millisecond_time();
	kmeans(_data, k, _labels, opencv_criteria, 3, KMEANS_PP_CENTERS, _centers);
	t2 = get_millisecond_time();
	cout << "Finished in " << t2-t1 << "[ms]" << endl;
}

/*int main(int argc, char * argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}*/
