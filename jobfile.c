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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include "jobfile.h"


static bool is_comment(char *line)
{
    if (isalpha(line[0]))
    {
        return false;
    }
    return true;
}


static void skip_comment_lines(FILE *stream, char *line)
{
    do
    {   
        fgets(line, LINE_LENGTH, stream);
    }   
    while (is_comment(line));
}


static void parse_option_line(char *line, option_t type, char *name, void *storage)
{   
    int n = 0;
    char varname[LINE_LENGTH + 1];

    switch (type)
    {   
        case INTEGER:
            n = sscanf(line, "%s " INPUTSEP " %d", varname, (int *) storage);
            break;
        case LONG_INTEGER:
            n = sscanf(line, "%s " INPUTSEP " %ld", varname, (long int *) storage);
            break;
        case REAL:
            n = sscanf(line, "%s " INPUTSEP " %lf", varname, (double *) storage);
            break;
        case STRING:
            n = sscanf(line, "%s " INPUTSEP " %s", varname, (char *) storage);
            break;
    }

    // Did we succesfully parsed the input line?
    if (n != 2)
    {
        printf("%s:%d: Read too few tokens.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    // Are we storing the right line into our storage?
    if (strncmp(varname, name, sizeof(varname)) != 0)
    {
        printf("%s:%d: Read variable '%s', but expected variable '%s'.\n", __FILE__, __LINE__, varname, name);
        exit(EXIT_FAILURE);
    }
}


void read_option(FILE *stream, option_t type, char *name, void *storage)
{
    char line[LINE_LENGTH + 1];
    skip_comment_lines(stream, line);
    parse_option_line(line, type, name, storage);
}

