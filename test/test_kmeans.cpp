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
 * Another converter
 */
void convert_vector_mat(vector<double *> _f, Mat& _t, size_t d) {
	size_t i = 0;
	_t.release();
	size_t rows = 0;
	vector<double *>::iterator it = _f.begin();
	vector<double *>::iterator ie = _f.end();
	while(it != ie) {
		i = 0;
		while(i < d) {
			_t.push_back((*it)[i]);
			++i;
		}
		++it;
		++rows;
	}
	cout << rows << endl;
	cout << "(" << _t.rows << "," << _t.cols << ")" << endl;
	_t = _t.reshape(0,rows);
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
		double x;

		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> real_dis(0.0, static_cast<double>(N));

		for(i = 0; i < N; i++) {
			double * tmp = (double *)malloc(d * sizeof(double));
			if(tmp == NULL) {
				cerr << "Cannot allocate memory at loop " << i << endl;
				exit(1);
			}
			for(j = 0; j < d; j++) {
				scale = real_dis(gen);
				x = sqrt(scale);
				tmp[j] = x;
			}
			data.push_back(tmp);
		}

//		print_vector(data,d);
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
	static vector<double *> data;
	static int N;
	static int d;
	static int k;
};

vector<double *> KmeansTest::data;
int KmeansTest::N;
int KmeansTest::d;
int KmeansTest::k;

TEST_F(KmeansTest, test1) {
	vector<double *> centroids, seeds;
	vector<int *> clusters;
	KmeansCriteria criteria = {1.0,10};
	cout << "Our method:" << endl;
	unsigned long int t1, t2;
	t1 = get_millisecond_time();
	simple_k_means(KmeansType::RANDOM_SEEDS,N,k,criteria,d,data,centroids,clusters,seeds);
	t2 = get_millisecond_time();
	cout << "Finished in " << t2-t1 << "[ms]" << endl;
	cout << "The distortion is " << distortion(d,N,k,data,centroids,clusters) << endl;
	cout << "=======================" << endl;
	cout << "OpenCV method:" << endl;
	Mat _data;
	convert_vector_mat(data,_data,d);
	_data.convertTo(_data,CV_32F);
	Mat _labels, _centers(k, 1, _data.type());
	TermCriteria opencv_criteria {CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1.0};
	t1 = get_millisecond_time();
	kmeans(_data, k, _labels, opencv_criteria, 3, KMEANS_PP_CENTERS, _centers);
	t2 = get_millisecond_time();
	cout << "Finished in " << t2-t1 << "[ms]" << endl;
}
