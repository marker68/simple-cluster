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
 *  test_kmeans.cpp
 *
 *  Created on: 2014/09/12
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */


#include <iostream>
#include <vector>
#include <random>
#include <float.h>
#include <gtest/gtest.h>
#include <math.h>
#include <stdlib.h>
/*#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"*/
#include "k-means.h"
#include "utilities.h"

#ifdef _OPENMP
#include <omp.h>
#define SET_THREAD_NUM omp_set_num_threads(n_thread)
#else
#define SET_THREAD_NUM 0 // disable multi-thread
#endif

using namespace std;
using namespace SimpleCluster;

/**
 * Customized test case for testing
 */
class KmeansTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 10000;
		d = 128;
		k = 256;
		int i, j;

		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<float> real_dis(0.0, 255.0);

		if(!init_array<float>(data,N*d)) {
			cerr << "Cannot allocate memory for test data!" << endl;
			exit(1);
		}
		int base = 0;
		for(i = 0; i < N; i++) {
			for(j = 0; j < d; j++) {
				data[base++] = real_dis(gen);
			}
		}

		if(!init_array<float>(seeds,k*d)) {
			cerr << "Cannot allocate memory for seeds data!" << endl;
			exit(1);
		}
		if(!init_array<float>(centers,k*d)) {
			cerr << "Cannot allocate memory for centers data!" << endl;
			exit(1);
		}
		label = (int *)calloc(N,sizeof(int));
		if(label == NULL) {
			cerr << "Cannot allocate memory for label data!" << endl;
			exit(1);
		}
		//		print_vector(data,d,N);
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		::delete data;
		data = nullptr;
		::delete seeds;
		seeds = nullptr;
		::delete centers;
		centers = nullptr;
		::delete label;
		label = nullptr;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static float * data;
	static float * seeds;
	static float * centers;
	static int * label;
	static int N, d, k;
};

float * KmeansTest::data;
float * KmeansTest::seeds;
float * KmeansTest::centers;
int * KmeansTest::label;
int KmeansTest::N;
int KmeansTest::d;
int KmeansTest::k;

TEST_F(KmeansTest, test1) {
	random_seeds<float>(data,seeds,d,N,k,8,false);
}

TEST_F(KmeansTest, test2) {
	kmeans_pp_seeds<float>(data,seeds,DistanceType::NORM_L2,d,N,k,8,false);
}

TEST_F(KmeansTest, test3) {
	KmeansCriteria criteria = {2.0,1.0,100};
	simple_kmeans<float>(
			data,centers,label,seeds,
			KmeansType::KMEANS_PLUS_SEEDS,
			KmeansAssignType::LINEAR,
			criteria,
			DistanceType::NORM_L2,
			EmptyActs::SINGLETON,
			N,k,d,8,
			false);
	cout << "ANN: Distortion is " << distortion<float>(data,centers,label,DistanceType::NORM_L2,d,N,k,false) << endl;
}

TEST_F(KmeansTest, test4) {
	KmeansCriteria criteria = {2.0,1.0,100};
	simple_kmeans<float>(
			data,centers,label,seeds,
			KmeansType::RANDOM_SEEDS,
			KmeansAssignType::LINEAR,
			criteria,
			DistanceType::NORM_L2,
			EmptyActs::SINGLETON,
			N,k,d,8,
			false);
	cout << "ANN: Distortion is " << distortion<float>(data,centers,label,DistanceType::NORM_L2,d,N,k,false) << endl;
}

TEST_F(KmeansTest, test5) {
	KmeansCriteria criteria = {2.0,1.0,100};
	greg_kmeans<float>(
			data,centers,label,seeds,
			KmeansType::RANDOM_SEEDS,
			criteria,
			DistanceType::NORM_L2,
			EmptyActs::SINGLETON,
			N,k,d,4,
			false);
	cout << "LINEAR: Distortion is " << distortion<float>(data,centers,label,DistanceType::NORM_L2,d,N,k,false) << endl;
}

TEST_F(KmeansTest, test6) {
	KmeansCriteria criteria = {2.0,1.0,100};
	greg_kmeans<float>(
			data,centers,label,seeds,
			KmeansType::KMEANS_PLUS_SEEDS,
			criteria,
			DistanceType::NORM_L2,
			EmptyActs::SINGLETON,
			N,k,d,4,
			false);
	cout << "LINEAR: Distortion is " << distortion<float>(data,centers,label,DistanceType::NORM_L2,d,N,k,false) << endl;
}

TEST_F(KmeansTest, test7) {
	KmeansCriteria criteria = {2.0,1.0,100};
	float * _data, * _centers, * _seeds;
	int * _labels;
	init_array(_data,6);
	_data[0] = _data[5] = 0.2f;
	_data[1] = _data[4] = 0.1f;
	_data[2] = _data[3] = 0.0f;
	init_array(_centers,12);
	init_array(_seeds,12);
	init_array(_labels,6);

	greg_kmeans<float>(
			_data,_centers,_labels,_seeds,
			KmeansType::KMEANS_PLUS_SEEDS,
			criteria,
			DistanceType::NORM_L2,
			EmptyActs::SINGLETON,
			3,6,2,4,
			false);
	cout << "LINEAR: Distortion is " << distortion<float>(_data,_centers,_labels,DistanceType::NORM_L2,2,3,6,false) << endl;
	cout << "The content of _centers:" << endl;
	for(int i = 0; i < 12; i++) {
		cout << _centers[i] << " ";
	}
	cout << endl;
}

/*TEST_F(KmeansTest, test6) {
	Mat _data;
	convert_array_to_mat(data,_data,N,d);
	_data.convertTo(_data,CV_32F);
	Mat _labels, _centers(k, 1, _data.type());
	TermCriteria opencv_criteria {CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 1.0};
	kmeans(_data, k, _labels, opencv_criteria, 3, KMEANS_PP_CENTERS, _centers);
	// Find the distortion
	float distortion = 0.0;
	int i;
	int n_thread = 8;
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < _labels.rows; i++) {
			int id = _labels.at<int>(i,0);
			Mat tmp = _data.row(i);
			Mat c = _centers.row(id);
			float d_t = norm(tmp,c,NORM_L2);
			distortion += d_t * d_t;
		}
	}
	cout << "Distortion is " << sqrt(distortion) << endl;
}*/

int main(int argc, char * argv[])
{
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}

