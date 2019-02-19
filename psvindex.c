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


#include "psvindex.h"


void idx_initialise(psvidx_t *idx, int num_moons, int integrate_grain_potential)
{
    idx->num_moons = num_moons;
    idx->moons_base_idx = 6;
    idx->sun_base_idx = 6 + num_moons*6;
    // idx->nvars = 6 * (num_moons + 1) + 2;        // without the sun
    idx->nvars = 6 * (num_moons + 2) + 2;           // with the sun

    idx->neq = 6;                   // dust grain
    idx->neq += 6 * num_moons;      // moons
    idx->neq += 6;                  // sun
    if (integrate_grain_potential)
        idx->neq += 1;

    
    if (integrate_grain_potential)
        idx->potential_base_idx = idx->neq - 1;
    else
        idx->potential_base_idx = idx->nvars - 2;

    idx->time_base_idx = idx->nvars - 1;
}

int grain_pos_idx(psvidx_t *idx, int component)
{
    return component;
}

int grain_vel_idx(psvidx_t *idx, int component)
{
    return component + 3;
}

int moon_pos_idx(psvidx_t *idx, int moon_num, int component)
{
    return idx->moons_base_idx + 6*(moon_num - 1) + component;
}

int moon_vel_idx(psvidx_t *idx, int moon_num, int component)
{
    return idx->moons_base_idx + 6*(moon_num - 1) + 3 + component;
}

int sun_pos_idx(psvidx_t *idx, int component)
{
    return idx->sun_base_idx + component;
}

int sun_vel_idx(psvidx_t *idx, int component)
{
    return idx->sun_base_idx + 3 + component;
}

int grain_pot_idx(psvidx_t *idx)
{
    return idx->potential_base_idx;
}

int time_idx(psvidx_t *idx)
{
    return idx->time_base_idx;
}

int moon_start_idx(void)
{
    return 1;
}

int moon_end_idx(psvidx_t *idx)
{
    return idx->num_moons;
}

bool exist_moon(psvidx_t *idx, int j)
{
    return (j >= 1 && j <= idx->num_moons);
}

int moon_array_idx(psvidx_t *idx, int moon_num)
{
    return moon_num - 1;
}

