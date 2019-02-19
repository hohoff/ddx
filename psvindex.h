/*
 * Copyright 2019 Holger Hoffmann, Ted Moldenhawer, Thomas Münch, Jürgen Schmidt, Michael Seiler, Frank Spahn
 * 
 * This file is part of DDX.
 * 
 * DDX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * DDX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DDX.  If not, see <https://www.gnu.org/licenses/>.
 * 
 */


#ifndef dyne_psv_index
#define dyne_psv_index

#include <stdbool.h>

typedef struct {
    int nvars;          // number of phase space variables
    int neq;            // number of equations to be integrated
    int num_moons;
    int potential_base_idx;
    int moons_base_idx;
    int sun_base_idx;
    int time_base_idx;
}
psvidx_t;               // metadata for phase space vector indices

void idx_initialise(psvidx_t *idx, int num_moons, int integrate_grain_potential);
int grain_pos_idx(psvidx_t *idx, int component);
int grain_vel_idx(psvidx_t *idx, int component);
int grain_pot_idx(psvidx_t *idx);
int time_idx(psvidx_t *idx);
int sun_pos_idx(psvidx_t *idx, int component);
int sun_vel_idx(psvidx_t *idx, int component);
int moon_pos_idx(psvidx_t *idx, int moon_num, int component);
int moon_vel_idx(psvidx_t *idx, int moon_num, int component);
int moon_start_idx(void);
int moon_end_idx(psvidx_t *idx);
bool exist_moon(psvidx_t *idx, int j);
int moon_array_idx(psvidx_t *idx, int moon_num);
#endif
