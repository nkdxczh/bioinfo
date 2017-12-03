#ifndef __CONF__H__
#define __CONF__H__

#include <stdio.h>

#include "defs.h"

BEGIN_C_DECL

/** 
 * Reads in a configuration file.  Reads are additive and copy over.
 *
 * @param fp a file pointer to the file.
 *
 * @return the number of variables read or -1 on an error.
 */

int read_conf(FILE *fp);

/**
 * Writes a configuration file from the current state of the machine.
 *
 * @param fp a file pointer to the file.
 *
 * @return the number of variables written or -1 on an error.
 */

int write_conf(FILE *fp);

/**
 * Checks if a variable is set.
 *
 * @param name the name of the variable.
 *
 * @return TRUE if the variable is set and false otherwise.
 */

boolean_t_new check_variable_conf(string_t name);

/**
 * Read double variable.
 *
 * @param name the name of the variable.
 * @param def a default value for the variable.
 *
 * @return the value if set or the default value if not set.
 */

double read_double_conf(string_t name, double def);

/**
 * Read string variable.
 *
 * @param name the name of the variable.
 * @param def a default value for the variable.
 *
 * @return the value if set or the default value if not set.
 */

string_t read_string_conf(string_t name, string_t def);


/**
 * Read int variable.
 *
 * @param name the name of the variable.
 * @param def a default value for the variable.
 *
 * @return the value if set or the default value if not set.
 */

int read_int_conf(string_t name, int def);

/**
 * Read data variable.
 *
 * @param name the name of the variable.
 * @param length set to the length of the data.
 * @param def a filename to read the data from.
 *
 * @return a pointer to the data or NULL on a failure.
 */

ptr_t read_data_conf(string_t name, int *length, string_t def);

/**
 * Read int range variable.
 * a range string can be of the following forms:
 * a:b:c ==> rmin=a, rstep=b, rmax=c
 * a:b   ==> rmin=a, rstep=1, rmax=b
 * a     ==> rmin=a, rstep=1, rmax=a
 *
 * @param name the name of the variable.
 * @param def a default value for the variable.
 *
 * @return the int array [rmin,rstep,rmax] if the value is set or the default
 * value if not set. 
 */

int* read_int_range_conf(string_t name, int* def);

END_C_DECL

#endif



