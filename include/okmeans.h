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
 *  okmeans.h
 *
 *  Created on: 2015/03/12
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef INCLUDE_OKMEANS_H_
#define INCLUDE_OKMEANS_H_

#include <iostream>
#include <exception>
#include <algorithm>
#include <vector>
#include <random>
#include <cstring>
#include <climits>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include "utilities.h"
#include "rand_utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace SimpleCluster {

typedef struct {
	int m;
	int p;
	int type;
} CKModel;

void ok_init(
		float * X,
		int m,
		int n,
		int p,
		int nthread,
		CKModel& model,
		float *& mu,
		float *& X_mu,
		float *& R,
		float *& R_pc,
		bool verbose);

void ok_loop(
		int p,
		int m,
		int n,
		float * X,
		float *& R, // = R in ok_init
		float *& R_pc,
		float *& X_mu, // = X_mu in ok_init
		float *& B,
		float *& D,
		float *& DB,
		float *& RX,
		float *& RDB,
		float *& XDB,
		float *& S,
		float *& V,
		float *& U,
		float *& mu,
		int nthread,
		bool verbose
		);
}

#endif /* INCLUDE_OKMEANS_H_ */
