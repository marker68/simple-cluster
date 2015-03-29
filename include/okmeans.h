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
#include "rand.h"

// Use BLAS/LAPACK packages
#include <lapacke.h>
#include <cblas.h>

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
		bool verbose) {
	if(nthread <= 0) nthread = 1;
	if(m <= 0 || n <= 0 || p <= 0)
		return;

	// Init mu
	mean(X,n,p,1,mu,verbose);

	// Init X_mu
	cblas_scopy(n * p,X,1,X_mu,1);
	int i0, i, j;
	int blk = n / nthread;
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#pragma omp parallel
	{
#pragma omp for private(i, i0)
#endif
		for(i0 = 0; i0 < nthread; i0++) {
			int start = i0 * blk;
			int end = start + blk;
			if(end > n) end = n;
			float * tmp1 = X_mu + start;
			for(i = start; i < end; i++) {
				cblas_saxpy(p,-1.0f,mu,1,tmp1,1);
			}
		}
#ifdef _OPENMP
	}
#endif

	// Compute the co-variance matrix of X_mu
	// C = 1/p * X_mu' * X_mu (note: E[X_mu] = 0)
	float * C = (float *)::operator new (p * p * sizeof(float));
	cblas_sgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			p, p, n,
			1.0f / p,
			X_mu, p,
			X_mu, n,
			0.0f,
			C, p);

	// Generate a random m * m matrix and decompose it by SVD
	float * rm;
	gen_rand_vector_float(rm,m * m,0.0f,1.0f * m,true);
	R = (float *)::operator new(m * m * sizeof(float)); // R
	float * S = (float *)::operator new(m * m * sizeof(float)); // S
	float * V = (float *)::operator new(m * m * sizeof(float)); // V
	float * superb = (float *)::operator new(m * sizeof(float));
	LAPACKE_sgesvd(
			LAPACK_ROW_MAJOR,
			'A',
			'A',
			m,m,
			rm,m,
			S,
			R,m,
			V,m,
			superb);

	// Calculate eigenvalues of C and init R.
	// Note: C is a symmetric matrix.
	int found, * isuppz;
	float * w, * z, * z_pc;
	if(p > m) {
		R_pc = (float *)::operator new(p * m * sizeof(float)); // R
		w = (float *)::operator new(p * sizeof(float)); // eigenvalues
		z = (float *)::operator new(p * p * sizeof(float)); // eigenvectors ASC
		z_pc = (float *)::operator new(p * m * sizeof(float)); // eigenvectors DESC
		isuppz = (int *)::operator new((m << 1) * sizeof(float));
		// Calculate m largest magnitude eigenvalues and corresponding eigenvectors
		LAPACKE_ssyevr(
				LAPACK_COL_MAJOR,
				'V',
				'A',
				'U',
				p,
				C,
				p,
				0.0f,
				0.0f,
				1,
				p,
				LAPACKE_slamch('S'),
				&found,
				w, // l
				z, // pc
				p,
				isuppz);

		// convert into DESC order
		for(i = 0; i < m; i++) {
			cblas_scopy(p,z+(p-i-1) * p,1,z_pc+i * p,1);
		}

		// R_pc = z_pc * R
		cblas_sgemm(
				CblasRowMajor,
				CblasTrans,
				CblasNoTrans,
				p, m, m,
				1.0f,
				z_pc, p,
				R, m,
				0.0f,
				R_pc, p);
	}
}

void ok_loop(
		) {

}
}

#endif /* INCLUDE_OKMEANS_H_ */
