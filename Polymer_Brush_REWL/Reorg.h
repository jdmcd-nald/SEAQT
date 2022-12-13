#ifndef REORG_12_9_22
#define REORG_12_9_22

int lp_reorg_non_rebrid(int* p_chain, int* c_chain, int start, int end); // reorganizes the lattice_polymer for non rebridging moves takes the input of the former chain and the updated chain
int lattice_reorg_non_rebrid(int* p_chain, int* c_chain, int fill, int start, int end); // reorganizes the latticepoints more efficient than previous scrubbing over entire lattice
int reorg_wl(int start, int end); // reorganizes the lattice_polymer for non rebridging moves takes the input of the former chain and the updated chain
int lp_reorg_rebrid(int index, int h_t); //The reorganization program for the rebridgining operations would like for the tqo peograms to be joined but the rebridging is too different to be generalized

#endif
