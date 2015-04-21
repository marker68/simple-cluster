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
 *  test_okmeans.cpp
 *
 *  Created on: 2015/03/30
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */


#include <iostream>
#include <vector>
#include <random>
#include <cfloat>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>

#include "okmeans.h"

using namespace std;
using namespace SimpleCluster;

class OkmeansTest : public ::testing::Test {
protected:
	static void SetUpTestCase() {
		X = (float *)::operator new (9 * sizeof(float));
		X_mu = (float *)::operator new (9 * sizeof(float));
		X[0] = 0.9f;
		X[1] = 1.2f;
		X[2] = 2.3f;
		X[3] = 0.1f;
		X[4] = 1.3f;
		X[5] = 4.3f;
		X[6] = 0.4f;
		X[7] = 0.3f;
		X[8] = 0.5f;

		mu = (float *)::operator new(3 * sizeof(float));
		R = (float *)::operator new(4 * sizeof(float));
		R_pc = (float *)::operator new(6 * sizeof(float));
	}

	static void TearDownTestCase() {

	}

	virtual void SetUp() {}
	virtual void TearDown() {}
public:
	static float * X, * X_mu, * mu;
	static float * R, * R_pc;
	static float * C;
};

float * OkmeansTest::X;
float * OkmeansTest::X_mu;
float * OkmeansTest::mu;
float * OkmeansTest::R;
float * OkmeansTest::R_pc;
float * OkmeansTest::C;

TEST_F(OkmeansTest, test1) {
	ok_init(
			X,
			2,
			3,
			3,
			1,
			mu,
			X_mu,
			C,
			R,
			R_pc,
			true);

	cout << "mu:" << endl << "[";
	for(int i = 0; i < 3; i++) {
		cout << mu[i] << " ";
	}
	cout << "]" << endl;

	cout << "X_mu:" << endl << "[";
	for(int i = 0; i < 9; i++) {
		cout << X_mu[i] << " ";
	}
	cout << "]" << endl;

	cout << "C:" << endl << "[";
	for(int i = 0; i < 9; i++) {
		cout << C[i] << " ";
	}
	cout << "]" << endl;

	cout << "R_pc:" << endl << "[";
	for(int i = 0; i < 6; i++) {
		cout << R_pc[i] << " ";
	}
	cout << "]" << endl;
}

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
