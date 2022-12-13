#ifndef ENERGY_12_9_22_H
#define ENERGY_12_9_22_H

int total_energy() ;
int total_energy(int* p_l_c,int* p_l_i, int* s_l) ;
int local_energy(int* p_l_c,int* p_l_i, int* s_l) ;
int local_energy_op(int* p_l_i,int* s_l_i);
void exclusion(int (&arr)[120],int index);
int local_energy_op_2(int* p_l_i,int* s_l_i,int reset_ind, int s, int e);
#endif
