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


#ifndef MOON_NAME_LENGTH
#define MOON_NAME_LENGTH 50
#endif

#ifndef dyne_impact_h
#define dyne_impact_h

#include "psvindex.h"
#include "moons.h"

typedef enum {
    NO_IMPACT,
    IMPACT
} impact_status_t;

typedef struct {
    impact_status_t status;
    char sink_name[MOON_NAME_LENGTH];
    double time;
} impact_data_t;

void check4impacts(double GMplanet,
                   type_moons *moons,
                   psvidx_t *idx,
                   double oldpsv[],
                   double newpsv[],
                   impact_data_t *impact);
#endif

