#ifndef MOV_PIVOT_12_9_22
#define MOV_PIVOT_12_9_22


int init_pivot() ;// instances the vector that will hold the pivot indexes
int pivot_rotations(int place);
int pivot_assignment(int place, int chain_y, int chain_pivot);
int pivot_extract(int reference_index, int chain_reference, int type) ;// will extract the relative indexes for the pivot operation reference_index is the relative pivot index and chain reference is the poly chain referecnce

int pivot_movement() ;// will perform the main pivot operation


#endif
