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
 *  test_utilities.cpp
 *
 *  Created on: 2014/09/09
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <random>
#include <gtest/gtest.h>
#include <math.h>
#include <stdlib.h>
#include "utilities.h"

using namespace std;
using namespace SimpleCluster;

/**
 * A comparator
 */
int compare_double(const double * _a, const double * _b) {
	if(*_a > *_b) return 1;
	else if(*_a < *_b) return -1;
	else return 0;
}

/**
 * Customized test case for testing
 */
class UtilTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 10000;
		d = 128;
		int i, j;

		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> real_dis(0.0, static_cast<double>(N));

		data = (double **)::operator new(N * sizeof(double *));
		for(i = 0; i < N; i++) {
			::new(data +i) double *;
			data[i] = (double *)::operator new(d * sizeof(double));
			for(j = 0; j < d; j++) {
				::new(data[i] + j) double;
				data[i][j] = real_dis(gen);
			}
		}
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		size_t i;
		for(i = 0; i < N; i++)
			::delete data[i];
		::delete data;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static double ** data;
	static int N;
	static int d;
};

double ** UtilTest::data;
int UtilTest::N;
int UtilTest::d;

// Let's start with some testcases
TEST_F(UtilTest, test1) {
	EXPECT_EQ(0.0, distance(data[0],data[0],d));
}

TEST_F(UtilTest, test2) {
	EXPECT_EQ(0.0, distance_square(data[0],data[0],d));
}

TEST_F(UtilTest, test3) {
	EXPECT_LT(0.0, distance(data[0],data[1],d));
}

TEST_F(UtilTest, test4) {
	double d1 = distance(data[0],data[4],d);
	double d2 = distance_square(data[0],data[4],d);
	EXPECT_LT(abs(d2 - d1 * d1), 1e-04);
}

TEST_F(UtilTest, test5) {
	int index[] = {0,1,2,3,4};
	int * index2 = (int *)malloc(5 * sizeof(int));
	for(int i = 0; i < 5; i++)
		index2[i] = i;

	double * tmp = (double *)calloc(d,sizeof(double));
	double * v1 = mean_vector(data,index,d,5,tmp);
	double * v2 = mean_vector(data,index2,d,5,tmp);
	EXPECT_EQ(0.0, distance(v1,v2,d));

	free(index2);
	free(tmp);
	free(v1);
	free(v2);
}

TEST_F(UtilTest, test6) {
	double arr[] = {1.0, 3.0, 5.0, 7.0, 9.0};
	int m = quick_select_k(arr,5,3,compare_double);
	EXPECT_EQ(3,m);
}

TEST_F(UtilTest, test7) {
	double arr[10000];
	for(int i = 0; i < 10000; i++)
		arr[i] = i;
	int m = quick_select_k(arr,10000,5000,compare_double);
	EXPECT_EQ(5000,m);
}

TEST_F(UtilTest, test8) {
	double * t;
	EXPECT_TRUE(init_array<double>(t,1000) && t!=nullptr);
}

TEST_F(UtilTest, test9) {
	double * t;
	EXPECT_TRUE(init_array<double>(t,d) && t != nullptr);
	EXPECT_TRUE(copy_array<double>(data[0],t,d));
	for(int i = 0; i < d; i++) {
		EXPECT_TRUE(data[0][i] == t[i]);
	}
}

TEST_F(UtilTest, test10) {
	double ** t;
	EXPECT_TRUE(init_array_2<double>(t,1000,200) && t!=nullptr);
}

TEST_F(UtilTest, test11) {
	double ** t;
	EXPECT_TRUE(init_array_2<double>(t,N,d) && t != nullptr);
	EXPECT_TRUE(copy_array_2<double>(data,t,N,d));
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < d; j++) {
			EXPECT_TRUE(data[i][j] == t[i][j]);
		}
	}
}

TEST_F(UtilTest, test12) {
	double ** t;
	EXPECT_TRUE(init_array_2<double>(t,N,d) && t != nullptr);
	EXPECT_TRUE(dealloc_array_2<double>(t,N));
}
