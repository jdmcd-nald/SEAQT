#ifndef MOV_12_9_22_H
#define MOV_12_9_22_H

#include "General.h"
#include "Energy.h"

int propose_update();    // returns energy change if a spin would be updated
int poly_mov();
int mov_crankshaft();
int mov_pivot();
int mov_inversion();
int mov_planar_reflection();
int mov_translation();

int solvent_translation();


int reflection_wrt_plane(double(&s_pos)[3], double(&e_pos)[3], int(&init_p)[3], double(&refl_p)[3], int axis);
int rotation_wrt_vector(double(&s_pos)[3], double(&e_pos)[3], int(&init_p)[3], double(&rot_p)[3]);
#endif
