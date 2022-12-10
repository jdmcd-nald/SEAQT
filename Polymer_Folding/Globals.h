#ifndef GLOBALS_12_9
#define GLOBALS_12_9


#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <float.h>


extern int switched;

extern int iterr; // Apart of an added debug output. Effictively just a counter for the number of repeated histogram values when flatness is checked
extern int repv; // Holds the value of the (potentially) repeated minvalue of the histogram
extern int phrepv; // Holds the value of the previous minimum histogram value to compare with the prior variable

extern int localgroup;
extern int local_dup_group;
extern int num_dup_proc;
extern int dup_headproc;

const int L1dim = 60;                            // linear dimension of the system
const int Dimension = 3;                        // number of dimensions in the system
const int numberspins = pow(L1dim, Dimension);  // total number of spins
const int numberneighbors = 6;                  // for 3D square lattice (others not implemented right now)
const int Eglobalmin = -48;          // minimum energy for 3D square lattice Potts model
const int Eglobalmax = 0;// -188/2;                       // maximum energy of Potts model
const int Eglobalwidth = abs(Eglobalmin - Eglobalmax);
const int bctype = 0;                           // type of boundary condition: 0 - periodic; 1 - Braskamp Kunz
const int q = 2;                               // number of different possible spin states (q-state Potts model)
const int chain_length = 58;                  // Length of the chain

extern int* latticepoint;                      // list containing values of all spins
extern int* neighbor;                          // list containing indices of neighbors for all spins

//Value buffers for recombining function
extern double* HE;                             // energy histogram
extern double* lngE;                           // ln g(E) estimator
extern double* lngE_buf;                       // ln g(E) estimator buffer for exchange
extern double* pseudolngE;
extern double* real_lngE;
extern double* real_lngE_buf;
extern double* microT;
extern double* microT_buf;
const int hist_size = (-Eglobalmin) + 1;        // histogram size

//Polymer chain constants
extern int* poly_chain;                        // Polymer chain indexholder
extern int* pseudo_poly_chain;                 // Pseudo polymer to hold previous indexes if change is rejected

extern int* doublepseudo_poly_chain;                 // Pseudo polymer to hold previous indexes if change is rejected

extern int* pull_poly_chain;                 // Pseudo polymer to hold previous indexes if change is rejected
extern int* wl_pseudo_chain;                   // records an image of the polymer chain in case the movement is rejected
extern int* chain_sequence;                    // stores the polymer sequence so it can be reset for all operations specifically rebridging
extern int* pivot_indexes;



extern int* digital_eye;  // holds the pivot indexes
extern int rebridge_index;

extern double* centeromass;
extern int* s_vector;
extern int* s_vector1;
extern int* s_vector2;
extern double* tortuosity;
extern double* rog;
extern long long* visits;
extern double* tortuosity_buf;
extern double* rog_buf;
extern long long* visits_buf;
extern double* s_avg;
extern double*  endtoend;
extern double* endtoend_buf;



//polymer rebridging accomadations
extern int* latticepolymer;                          // list containing indices of neighbors for all spins

extern int rseed;                              // seed for random number generator
extern int energy, localenergy;

extern double Emin, Emax;                       // doubles as calculation of boundaries in Energy is not in 'int'
extern int Eminindex, Emaxindex, Estartindex;  // local boundaries as index 

// MPI; set up with local communicators for replica exchange (RE)
// Needed so that two processes can talk w/o depending on / bothering the others
extern int numprocs, myid, multiple, comm_id;
// each process belongs to two local groups,
// one to communicate to left neighbor and one to communicate to right neighbor
// each process has different loca IDs in different communicatore, in general 
extern int mylocalid[2];                       // id in local communicators 
extern MPI_Comm* mpi_local_comm;
extern MPI_Group* mpi_local_group;
extern MPI_Group world;
extern MPI_Status status;
extern int merge_hists;                    // flag whether or not to merge histograms between iterations for multiple walkers on the same energy range

// to keep track of exchange statistics
extern int tryleft, tryright, exchangeleft, exchangeright;

// File handlers for I/O
extern FILE* file;
extern FILE* stdoutlog;
extern FILE* wanderlog;
extern char filename[50];
extern char resetbuffer[50];
extern char stdoutlogname[128];
extern char wanderlogname[128];
extern int ret_status;

extern double flatratio;
extern double flatmin;

extern int happened;
extern int happened2;
extern int happened3;

extern int cf_value;

extern int MoveProposal;

extern int elocvor;
extern int elocnach;
extern int r_picker;
extern int r_picke;
extern int compl_check;
//int buffer=0;
extern int hundred;

extern int axis;
extern int cos_value;
extern int sin_value;

extern int* pivot_indexes_r;
extern int xr;
extern int yr;
extern int zr;
extern int seq_plhdr;
extern int check_plhdr;
extern int pivotflag;

extern int xx;
extern int yy;
extern int zz;

extern int search_limit;
extern int search_checks;
extern int* configuration_search_r;
extern int* configuration_search_t;
// to track execution times
extern time_t timenow, timestart, timeend;

extern int pull_loc;
extern int pull_hort;

extern int pull_relax_loc;
extern int pull_relax_hort;

extern int wiggle_ht;

extern int pivot_loc;

extern int re_ht_si;
extern int re_ht_ht;

extern double* enindex;// hold values for the radius of gyration at specfic energy levels
extern double* enindex2;// hold values for the tortuosity at specfic energy levels
extern double* enindex3;// hold values for the tortuosity at specfic energy levels
extern double* lowchecks;
extern double* highchecks;
extern int* checks;
extern double* checks_countdown;
extern long double unused1;
extern long double unused2;

extern int working;

#endif
