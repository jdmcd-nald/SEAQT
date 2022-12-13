#ifndef GENERAL_12_9_22_H
#define GENERAL_12_9_22_H
#include "Constants.h"

void init_poly_cores(int); // init the polymer coordinates

void init_en_and_dist_array(int); // initialize the energy and the distance array

void poly_coord_decompose(int, int); // decomompes a given polymer index into cartesian coordinates
double point_distance(int, int); //gives the distance between two points
int int_point_distance(int, int);
double point_distance(int, int,int(&arr)[3]); //gives the distance between two points
double poly_point_distance(int*,int, int); //gives the distance between two points

double poly_point_distance(int*,int, int,int(&arr)[3]); //gives the distance between two points
int general_coord(int(&arr)[3], int); // output the x y and z compenent for any index
int general_index(int(&arr)[3]); // output index from any x , y , z components
int general_index_op(int(&arr)[3]); // output index from any x , y , z components
int poly_latticepoint_clear(int*, int, int); // clears a set group of lattice points (to prevent duplication of the lattice to perform a move
int poly_latticepoint_clear(int);
int latticepoint_assign(int,int); // assigns lattice point based on index values taken from argument array (also checks for occupation will output 0 if assignment overlap will occur
int latticepoint_assign2(int,int);

int reduced_lattice_optimization(int,int);
int anint(int, int);
int anint(float);

void lattice_polymer_reset(int* p_chain, int* c_chain, int start, int end);
int lattice_poly_index_reset(int* p_chain, int* c_chain, int start, int end);
int poly_coord_reset(int* p_chain, int* c_chain, int start, int end);
void solvent_reset(int*,int*);

int init_solvent(int number);

void optomize();

void eye();



#endif
