#ifndef MOV_REBRIDGE_12_9_22
#define MOV_REBRIDGE_12_9_22

int circle(int index, int next_index);
int h_t_wiggle() ;// this funtion is for rebridging of bonds around the head or tail (Flipping of the bonds)
int rebridging_h_t() ;// this funtion is for rebridging of bonds around the head or tail (Flipping of the bonds)
int rebridging_2(int index, int tagged) ;// Reconnects the segment which has formed a loop
int rebridging();

#endif
