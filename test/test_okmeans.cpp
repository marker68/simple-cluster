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

#include "rand.h"
#include "utilities.h"
#include "mat_utilities.h"
#include "okmeans.h"

using namespace std;
using namespace SimpleCluster;

class OkmeansTest : public ::testing::Test {
protected:
	static void SetUpTestCase() {

	}

	static void TearDownTestCase() {

	}

	virtual void SetUp() {}
	virtual void TearDown() {}
public:
	static float * X;
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
