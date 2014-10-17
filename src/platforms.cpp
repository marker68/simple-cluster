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
 * platforms.cpp
 *
 *  Created on: 2014/10/16
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <platforms.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>

#ifdef _OPENMP
#include <omp.h>
#define THREAD_NUM omp_get_thread_num()
#else
#define THREAD_NUM 0 // disable multi-threading
#endif

using namespace std;

/*! exectutes a set of tasks in parallel using a thread pool
 *
 * @param n            number of tasks to execute
 * @param nthread      number of threads that will run the tasks
 * @param task_fun     this callback will be called with
 *              -  arg = task_arg
 *              -  tid = identifier of the thread in 0..nthread-1
 *              -  i = call number in 0..n-1
 */
void compute_tasks (int n, int nthread,
		void (*task_fun) (void *arg, int tid, int i),
		void *task_arg) {
	int i;
#pragma omp parallel for schedule(dynamic) num_threads(nthread)
	for(i = 0; i < n; i++)
		(*task_fun)(task_arg, THREAD_NUM, i);
}

/**
 * Count the number of CPUs in the machine
 */
void cpu_count(int& n_cpu, int& n_cpu_max) {
	n_cpu = -1;
	n_cpu_max = -1;

#ifdef _WIN32
#ifndef _SC_NPROCESSORS_ONLN
	SYSTEM_INFO info;
	GetSystemInfo(&info);
#define sysconf(a) info.dwNumberOfProcessors
#define _SC_NPROCESSORS_ONLN
#endif
#endif

#ifdef _SC_NPROCESSORS_ONLN
	n_cpu = sysconf(_SC_NPROCESSORS_ONLN);
	if (n_cpu < 1)
	{
		cerr << "Could not determine number of CPUs online:" << strerror (errno) << endl;
		exit (1);
	}
	n_cpu_max = sysconf(_SC_NPROCESSORS_CONF);
	if (n_cpu_max < 1)
	{
		cerr << "Could not determine number of CPUs configured:" << strerror (errno) << endl;
		exit (1);
	}
	cout << n_cpu << " of " << n_cpu_max << " processors online" << endl;
	exit (0);
#else
	cerr << "Could not determine number of CPUs" << endl;
	exit (1);
#endif
}
