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
 * platforms.h
 *
 *  Created on: 2014/10/16
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef PLATFORMS_H_
#define PLATFORMS_H_

void compute_tasks (int, int,
                    void (*) (void *, int, int),
                    void *);

void cpu_count(int&, int&);

#endif /* PLATFORMS_H_ */
