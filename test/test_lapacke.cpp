/*
 *  Copyright (C) 2015 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
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
 *  test_lapacke.cpp
 *
 *  Created on: 2015/03/29
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

#include "rand.h"

using namespace std;
using namespace SimpleCluster;

/**
 * Test class
 */
class LapackeTest : public ::testing::Test {
protected:
	static void SetUpTestCase() {
		C = (float *)::operator new(9 * sizeof(float));
		C[0] = C[4] = C[5] = C[7] = C[8] = 0.5f;
		C[1] = C[2] = C[3] = C[6] = -0.5f;

		R = (float *)::operator new(4 * sizeof(float));
		R_pc = (float *)::operator new(6 * sizeof(float));

		w = (float *)::operator new(3 * sizeof(float)); // eigenvalues
		z = (float *)::operator new(9 * sizeof(float)); // eigenvectors ASC
		isuppz = (int *)::operator new(6 * sizeof(float));
	}

	static void TearDownTestCase() {
//		::delete C;
//		::delete w;
//		::delete z;
//		::delete z_pc;
//		::delete isuppz;
	}

	virtual void SetUp() {}
	virtual void TearDown() {}

public:
	static float * C, * R, * R_pc;
	static float * w, * z, * z_pc;
	static int * isuppz;
};

float * LapackeTest::C;
float * LapackeTest::R;
float * LapackeTest::R_pc;
float * LapackeTest::w;
float * LapackeTest::z;
float * LapackeTest::z_pc;
int * LapackeTest::isuppz;

TEST_F(LapackeTest, test1) {
	int found;
	// Calculate m largest magnitude eigenvalues and corresponding eigenvectors
	LAPACKE_ssyevr(
			LAPACK_COL_MAJOR,
			'V',
			'A',
			'U',
			3,
			C,
			3,
			0.0f,
			0.0f,
			0,
			0,
			LAPACKE_slamch('S'),
			&found,
			w, // l
			z, // pc
			3,
			isuppz);

	cout << "w:" << endl << "[";
	for(int i = 0; i < 3; i++) {
		cout << w[i] << " ";
	}
	cout << "]" << endl;

	cout << "z:" << endl;
	z_pc = z;
	for(int i = 0; i < 3; i++) {
		cout << "[";
		for(int j = 0; j < 3; j++) {
			cout << *(z_pc++) << " ";
		}
		cout << "]" << endl;
	}
	EXPECT_LT(fabs(w[0]),1e-06);
	EXPECT_LT(fabs(w[1]),1e-06);
	EXPECT_LT(fabs(w[2] - 1.5f),1e-06);
	EXPECT_EQ(3, found);
}

TEST_F(LapackeTest, test2) {
	z_pc = (float *)::operator new(6 * sizeof(float));
	cblas_scopy(3,z+6,1,z_pc,1);
	cblas_scopy(3,z+3,1,z_pc+3,1);
	EXPECT_LT(fabs(0.408248f - z_pc[3]), 1e-06);
}

TEST_F(LapackeTest, test3) {
	// Generate a random m * m matrix and decompose it by SVD
	float * rm = (float *)::operator new(4 * sizeof(float));
	rm[0] = 9.0f;
	rm[1] = 3.5f;
	rm[2] = 7.6f;
	rm[3] = 5.6f;
	float * S = (float *)::operator new(4 * sizeof(float)); // S
	float * V = (float *)::operator new(4 * sizeof(float)); // V
	float * superb = (float *)::operator new(2 * sizeof(float));
	LAPACKE_sgesvd(
			LAPACK_ROW_MAJOR,
			'A',
			'A',
			2,2,
			rm,2,
			S,
			R,2,
			V,2,
			superb);
	cout << "R:" << endl;
	float * tmp = R;
	for(int i = 0; i < 2; i++) {
		cout << "[";
		for(int j = 0; j < 2; j++) {
			cout << *(tmp++) << " ";
		}
		cout << "]" << endl;
	}

//	EXPECT_EQ(fabs(-0.715353f - R[0]), 1e-06);

	cout << "S:" << endl;
	tmp = S;
	for(int i = 0; i < 2; i++) {
		cout << "[";
		for(int j = 0; j < 2; j++) {
			cout << *(tmp++) << " ";
		}
		cout << "]" << endl;
	}

	cout << "V:" << endl;
	tmp = V;
	for(int i = 0; i < 2; i++) {
		cout << "[";
		for(int j = 0; j < 2; j++) {
			cout << *(tmp++) << " ";
		}
		cout << "]" << endl;
	}
}

TEST_F(LapackeTest, test4) {
	// R_pc = z_pc * R
	cblas_sgemm(
			CblasColMajor,
			CblasNoTrans,
			CblasTrans,
			3, 2, 2,
			-1.0f,
			z_pc, 3,
			R, 2,
			0.0f,
			R_pc, 3);
	cout << "R_pc:" << endl;
	float * tmp = R_pc;
	for(int i = 0; i < 2; i++) {
		cout << "[";
		for(int j = 0; j < 3; j++) {
			cout << *(tmp++) << " ";
		}
		cout << "]" << endl;
	}

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
