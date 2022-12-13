#ifndef MOV_PULL_12_10_22
#define MOV_PULL_12_10_22

int pull_2(int p_index, int c_index, int q, int c_chainlocation, int ht);
int pull_1(int q, int ht); //the chain forms a kink along an empty plane along random points crucial to avoid bottlenecking
int pull_relaxation(int k); //this will accomate both relaxation and repulling of the polymer chain

#endif
