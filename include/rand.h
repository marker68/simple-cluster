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
 *  rand.h
 *
 *  Created on: 2015/03/16
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef INCLUDE_RAND_H_
#define INCLUDE_RAND_H_

#include <iostream>
#include <random>


using namespace std;

namespace SimpleCluster {

void gen_rand_vector_float(float *& seed, int size, float low, float high, bool verbose) {
	if(size <= 0) return;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> real_dis(low, high);

	seed = (float *)::operator new(size * sizeof(float));
	for(int i = 0; i  < size; i++) {
		seed[i] = real_dis(gen);
	}
}
}



#endif /* INCLUDE_RAND_H_ */
