/*
 * high_kmeans.h
 *
 *  Created on: 2014/12/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */
#ifndef EXT_K_MEANS_H_
#define EXT_K_MEANS_H_
#include <iostream>
#include <k-means.h>

using namespace std;

namespace SimpleCluster {
/**
 * Update centers in 2*d dimensional space
 */
void update_center_2d(
		float *& centers,
		float *& means,
		int k,
		int d,
		bool verbose) {
	int i, j, base = 0;
	for(i = 0; i < k * d; i++) {
		if(i < d) means[i] = 0.0;
		means[i % d] += centers[i];
	}
	for(i = 0; i < d; i++)
		means[i] /= (1.0f * k);
	for(i = 0; i < k; i++) {
		for(j = 0; j < d; j++) {
			centers[base] += means[j];
			centers[base] /= 2.0f;
			base++;
		}
	}
}

/**
 * A kmeans method to extend kmeans to higher dimensional space
 */
template<typename DataType>
void extended_kmeans(
		DataType * data,
		float *& centers,
		int *& labels,
		float *& seeds,
		KmeansType type,
		KmeansAssignType assign,
		KmeansCriteria criteria,
		DistanceType d_type,
		EmptyActs ea,
		int N,
		int k,
		int d,
		int n_thread,
		bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	if(seeds == nullptr) {
		init_array<float>(seeds,k*d);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds<DataType>(data,seeds,d,N,k,n_thread,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds<DataType>(data,seeds,d_type,d,N,k,n_thread,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	int iters = criteria.iterations, it = 0, count = 0;
	float error = criteria.accuracy, e = error, e_prev,
			d_tmp, dfst;
	int i, j, f_tmp, fst,s_max, l_tmp;

	// Initialize the centers
	copy_array<float>(seeds,centers,k*d);
	init_array<int>(labels,N);
	for(i = 0; i < N; i++) labels[i] = -1;
	int * size;
	init_array<int>(size,k);
	float * sum;
	float * c_tmp, * means;
	init_array<float>(c_tmp,d);
	init_array<float>(means,d);
	init_array<float>(sum,k*d);
	int base = 0, base1, base2;
	for(i = 0; i < k; i++) {
		for(j = 0; j < d; j++)
			sum[base++] = 0.0;
		size[i] = 0;
	}
	if(verbose)
		cout << "Finished initialization" << endl;

	while (1) {
		// Update centers in 2*d dim space
		update_center_2d(centers,means,k,d,verbose);
		// Assigning
		linear_assign<DataType>(data,centers,labels,size,sum,
				d_type,d,N,k,n_thread,verbose);
		// Check for empty clusters
		if(ea != EmptyActs::NONE) {
			for(i = 0; i < k; i++) {
				if(size[i] <= 0) {
					if(ea == EmptyActs::SINGLETON_2) {
						l_tmp = 0;
						for(j = 0; j < k; j++) {
							if(l_tmp < size[j]) {
								l_tmp = size[j];
								s_max = j;
							}
						}
					}
					// Move the centers
					base = i * d;
					if(ea == EmptyActs::SINGLETON)
						find_lonely<DataType>(data,centers,labels,d_type,
								dfst,fst,N,k,d,verbose);
					else if(ea == EmptyActs::SINGLETON_2)
						find_farthest<DataType>(data,centers + base,labels,d_type,
								s_max,dfst,fst,N,k,d,verbose);
					base1 = fst * d;
					base2 = labels[fst] * d;
					for(j = 0; j < d; j++) {
						centers[base] = static_cast<float>(data[base1++]);
						sum[base] = centers[base];
						sum[base2++] -= centers[base++];
					}
					size[i] = 1;
					size[labels[fst]]--;
					labels[fst] = i;
				}
			}
		}

		if(verbose) {
			cout << "Spec " << it << ":";
			for(i = 0; i < k; i++) {
				cout << size[i] << " ";
			}
			cout << endl;
		}

		// Update centers in 2d dimensional space
		e_prev = e;
		e = 0.0;
		base = 0;
		for(i = 0; i < k; i++) {
			for(j = 0; j < d; j++) {
				c_tmp[j] = centers[base];
				centers[base] = sum[base++] / size[i];
			}
			e += distance_l2_square<float>(c_tmp,centers + (base - d),d);
		}
		e = sqrt(e);
		count += (fabs(e-e_prev) < error? 1 : 0);

		if(verbose)
			cout << "Iterator " << it
			<< "-th with error = " << e
			<< " and distortion = " <<
			distortion(data,centers,labels,d_type,
					d,N,k,false)
					<< endl;
		it++;

		if(
				it >= iters ||
				fabs(e-e_prev) < error ||
				count >= 10
		) break;
	}

	if(verbose)
		cout << "Finished clustering with error is " <<
		e << " after " << it << " iterations." << endl;
}

/**
 * Extended version of greg_kmeans
 */
template<typename DataType>
void extended_greg_kmeans(
		DataType * data,
		float *& centers,
		int *& label,
		float *& seeds,
		KmeansType type,
		KmeansCriteria criteria,
		DistanceType d_type,
		EmptyActs ea,
		int N,
		int k,
		int d,
		int n_thread,
		bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	if(seeds == nullptr) {
		init_array<float>(seeds,k * d);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds<DataType>(data,seeds,d,N,k,n_thread,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds<DataType>(data,seeds,d_type,d,N,k,n_thread,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	int iters = criteria.iterations, it = 0, count = 0;
	float error = criteria.accuracy, e = error, e_prev;

	// Variables for Greg's method
	float * c_sum;
	float * means;
	float * moved;
	float * closest;
	float * upper;
	float * lower;
	int * size;

	init_array<float>(c_sum,k * d);
	init_array<float>(means,d);
	init_array<float>(moved,k);
	init_array<float>(closest,k);
	init_array<float>(upper,N);
	init_array<float>(lower,N);
	init_array<int>(size,k);

	int i, j, s_max, l_tmp, fst, base, base0, base1, base2;
	float * fpt1, * fpt2;
	DataType * dpt = data;
	int tmp = 0;
	float min, min2, min_tmp = 0.0,
			d_tmp = 0.0, m, dfst;

	// Initialize the centers
	copy_array<float>(seeds,centers,k * d);
	greg_initialize<DataType>(data,centers,c_sum,upper,lower,
			label,size,d_type,ea,N,k,d,verbose);
	if(verbose)
		cout << "Finished initialization" << endl;

	while (1) {
		// Update centers in 2*d dim space
		update_center_2d(centers,means,k,d,verbose);
		// Update the closest distances
		fpt1 = centers;
		for(i = 0; i < k; i++) {
			min2 = min = FLT_MAX;
			fpt2 = centers;
			for(j = 0; j < k; j++) {
				if(j != i) {
					if(d_type == DistanceType::NORM_L2)
						min_tmp = distance_l2<float>(fpt1,fpt2,d);
					else if(d_type == DistanceType::NORM_L1)
						min_tmp = distance_l1<float>(fpt1,fpt2,d);
					if(min > min_tmp) min = min_tmp;
				}
				fpt2 += d;
			}
			closest[i] = min;
			fpt1 += d;
		}

#ifdef _OPENMP
		omp_set_num_threads(n_thread);
#pragma omp parallel
		{
#pragma omp for private(i,j,d_tmp,m,min,min2,tmp)
#endif
			for(i = 0; i < N; i++) {
				// Update m for bound test
				d_tmp = closest[label[i]]/2.0;
				m = std::max(d_tmp,lower[i]);
				// First bound test
				if(upper[i] > m) {
					// We need to tighten the upper bound
					if(d_type == DistanceType::NORM_L2)
						upper[i] = distance_l2<DataType,float>(data + i * d,centers + label[i] * d,d);
					else if(d_type == DistanceType::NORM_L1)
						upper[i] = distance_l1<DataType,float>(data + i * d,centers + label[i] * d,d);
					// Second bound test
					if(upper[i] > m) {
						int l = label[i];
						min2 = min = FLT_MAX;
						tmp = -1;
						// Assign the data to clusters
						fpt1 = centers;
						for(j = 0; j < k; j++) {
							if(d_type == DistanceType::NORM_L2)
								d_tmp = distance_l2<float,DataType>(fpt1,data + i * d,d);
							else if(d_type == DistanceType::NORM_L1)
								d_tmp = distance_l1<float,DataType>(fpt1,data + i * d,d);
							if(min >= d_tmp) {
								min2 = min;
								min = d_tmp;
								tmp = j;
							} else {
								if(min2 > d_tmp) min2 = d_tmp;
							}
							fpt1 += d;
						}

						// Assign the data[i] into cluster tmp
						label[i] = tmp; // Update the label
						upper[i] = min; // Update the upper bound on this distance
						lower[i] = min2; // Update the lower bound on this distance

						if(l != tmp) {
							size[tmp]++;
							size[l]--;
							if(size[l] == 0) {
								if(verbose)
									cout << "An empty cluster was found!"
									" label = " << l << endl;
							}
							base = i * d;
							base0 = tmp * d;
							base1 = l * d;
							for(j = 0; j < d; j++) {
								c_sum[base0++] += static_cast<float>(data[base]);
								c_sum[base1++] -= static_cast<float>(data[base++]);
							}
						}
					}
				}
			}
#ifdef _OPENMP
		}
#endif
		// Check for empty clusters
		if(ea != EmptyActs::NONE) {
			for(i = 0; i < k; i++) {
				if(size[i] <= 0) {
					if(ea == EmptyActs::SINGLETON_2) {
						l_tmp = 0;
						for(j = 0; j < k; j++) {
							if(l_tmp < size[j]) {
								l_tmp = size[j];
								s_max = j;
							}
						}
					}
					// Move the centers
					base = i * d;
					if(ea == EmptyActs::SINGLETON)
						find_lonely<DataType>(data,centers,label,d_type,
								dfst,fst,N,k,d,verbose);
					else if(ea == EmptyActs::SINGLETON_2)
						find_farthest<DataType>(data,centers + base,label,d_type,
								s_max,dfst,fst,N,k,d,verbose);
					base1 = fst * d;
					base2 = label[fst] * d;
					for(j = 0; j < d; j++) {
						centers[base] = static_cast<float>(data[base1++]);
						c_sum[base] = centers[base];
						c_sum[base2++] -= centers[base++];
					}
					size[i] = 1;
					size[label[fst]]--;
					label[fst] = i;
				}
			}
		}
		// Move the centers
		update_center(c_sum,size,centers,moved,d_type,k,d,n_thread);
		// Update the bounds
		update_bounds(moved,label,upper,lower,N,k,n_thread);

		if(verbose) {
			cout << "Spec " << it << ":";
			for(i = 0; i < k; i++) {
				cout << size[i] << " ";
			}
			cout << endl;
		}

		// Calculate the distortion
		e_prev = e;
		e = 0.0;
		for(i = 0; i < k; i++) {
			e += moved[i];
		}
		e = sqrt(e);
		count += (fabs(e-e_prev) < error? 1 : 0);
		if(verbose)
			cout << "Iterator " << it
			<< "-th with error = " << e
			<< " and distortion = "
			<< distortion(data,centers,label,d_type,d,N,k,false)
			<< endl;
		it++;
		if(it >= iters || e < error || count >= 10) break;
	}

	if(verbose)
		cout << "Finished clustering with error is " <<
		e << " after " << i << " iterations." << endl;
}
}

#endif
