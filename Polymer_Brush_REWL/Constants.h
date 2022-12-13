#ifndef CONSTANTS_12_9_22_H
#define CONSTANTS_12_9_22_H

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <random>

extern int switched;

extern int cccounter;

extern int iterr; // Apart of an added debug output. Effictively just a counter for the number of repeated histogram values when flatness is checked
extern int repv; // Holds the value of the (potentially) repeated minvalue of the histogram
extern int phrepv; // Holds the value of the previous minimum histogram value to compare with the prior variable

// The lattice is effectively going to mimic a BCC Lattice
// I have code to track neighbors and second nearest neighbors for a BCC Lattice but the code models the unit cells exactly
// meaning that the end of the unit cell is cut off to make the modeled space periodic
// So given I would have to redo the math for neighboring sites regardless for this system
// I am implementing a different structure that utilizes a outer lattice (outer atoms of a BCC lattice, denoted as L because the size it's the largest lattice, and an inner lattice (interior atom of BCC (lattice), denoted S because it is one length in the three cardinal directions smaller than the larger lattice


extern int minmaxid;

const int L1dim_L = 34;                            // linear dimension of the
const int L1dim_S = L1dim_L;                            // linear dimension

const int L1dim_L_xy = 34;                            // linear dimension of the
const int L1dim_S_xy = 34;                            // linear dimension

const int L1dim_L_z = 80;                            // linear dimension of the
const int L1dim_S_z = 80;                            // linear dimension


const int Dimension = 3;                        // not sure at this moment if other than 2D is actually implemented
const int numberspins_L = L1dim_L_xy * L1dim_L_xy * L1dim_L_z;  // total number of spins
const int numberspins_S = L1dim_S_xy * L1dim_S_xy * L1dim_S_z;  // total number of spins
const int numberneighbors_L = 8;                  // for 2D square lattice (others not implemented right now)
const int numberneighbors_S = 6;                  // for 2D square lattice (others not implemented right now)
const int numbercores = 8;
const int lengthpercore = 40;
const int pertrusionpercore = 1; // number of arms per core


extern int* poly_lattice_indexes;
extern int* lattice_polymer;
extern int* poly_lattice_coordinates;
extern int* poly_lattice_connection;

const int Eglobalmin = -200000;          // minimum energy for 2D square lattice Potts model
const int Eglobalmax = 500000;// -188/2;                       // maximum energy of Potts model
const int Eglobalwidth = -(Eglobalmin - Eglobalmax);
const int bctype = 0;                           // type of boundary condition: 0 - periodic; 1 - Braskamp Kunz
extern int q;                               // number of different possible spin states (q-state Potts model)

extern int* latticepoint_L;                      // list containing values of all spins in the lattice (Checks occupation)
extern int* latticepoint_S;                      // list containing values of all spins (checks occupation)
extern int* latticepoint_S;                      // list containing values of all spins (checks occupation)

extern int* neighbor_L;                          // list containing indices of neighbors for all spins (Contains Large lattice indexes for the smaller lattice
extern int* neighbor_S;                          // list containing indices of neighbors for all spins
                            // energy histogram
extern double* HE;
                        // ln g(E) estimator
extern double* lngE;

extern double* lngE_buf;                       // ln g(E) estimator buffer for exchange
extern double* pseudolngE;
extern double* real_lngE;
extern double* real_lngE_buf;
extern double* microT;
extern double* microT_buf;
const int hist_size = (-Eglobalmin)+Eglobalmax;        // histogram size


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

extern FILE* DOSEstimate;//holds values for the dos estimate all values for the rog and tort are placed here to later be deemed unique or not
extern char dosestimate[50];

extern char filename[50];
extern char resetbuffer[50];
extern char stdoutlogname[128];
extern char wanderlogname[128];
extern int ret_status;

extern double flatratio;
extern double flatmin;

extern int en_choice;
extern int en_array[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int en_array_pp[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int en_array_ps[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int en_array_ss[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int en_array_bond[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

extern int int_en_array[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int int_en_array_pp[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int int_en_array_ps[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int int_en_array_ss[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern int int_en_array_bond[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

extern int en_array_wall[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
extern double dist_array[4][4][4]; // plan to hold the available distances in this array so the square root does not need to be called when values are needed
extern int int_dist_array[4][4][4]; // plan to hold the available distances in this array so the square root does not need to be called when values are neededxxxxxx
// to track execution times
extern time_t timenow, timestart, timeend;

extern int* solvent_loc;
extern int* ph_solvent_loc;
extern int* wl_solvent_loc;

extern int reset_solvent_index;
const int no_solvent_sites = 5800;

extern int* ph_poly_lattice_indexes;
extern int* ph_poly_lattice_coordinates;
extern int* wl_pseudo_chain;
extern int* wl_pseudo_chain_coordinates;

const double max_dist = sqrt(11);//maximum allowable distance any nieghboring monomers can be from each other
const int int_max_dist = 11;//maximum allowable distance any nieghboring monomers can be from each other
const int plane_z_max = L1dim_S_z-2;// axis upon which brush plane exists. Is a set value to check for any chains that try to go beneath it
const int plane_z_min =3;

extern int MoveProposal; // tracks if a move has been performed

//checks to calculate local energy
extern int poly_solvent_mov; // says wether a chosen mov moved solvent or a monomer
extern int chosen_mov ; // tellw which move was last chosen
extern int poly_mov_start; // tells which monomwrs were involved with a given move the starting position (on thr chain)
extern int poly_mov_end; // tells which monomer is the end considered loop

extern int reset_indexes[2];

const double max_en_dist = sqrt(10.5);//maximum checked distance for any monomer interaction
const int int_max_en_dist = 10;//maximum checked distance for any monomer interaction

//std::mt19937 gen(74); // Standard mersenne_twister_engine seeded with rd()
//std::uniform_real_distribution<> polymermovment(0, 99);

extern int* lattice_dist_optomize;
extern int* lattice_dist_optomize_2;
extern int* lattice_dist_optomize_2_lattice;

extern int* index_optomize_x;
extern int* index_optomize_y;
extern int* index_optomize_z;

extern int* calc_index_optomize_x;
extern int* calc_index_optomize_y;
extern int* calc_index_optomize_z;

extern int* polymer_optomize;

extern int what_to_exclude[120*numberneighbors_L];
extern int* local_en_indexes_op;
extern int local_en_dist[120];
extern int local_en_indexes_op_x[120];
extern int local_en_indexes_op_y[120];
extern int local_en_indexes_op_z[120];

extern int* xyz_index;
extern int* index_x;
extern int* index_y;
extern int* index_z;

const int neg_xy = L1dim_S_xy*6;
const int neg_z = L1dim_S_z*3;

const int lower_bounds = L1dim_S_xy*L1dim_S_xy*(plane_z_max-1);
const int upper_bounds = L1dim_S_xy*L1dim_S_xy*(plane_z_min+1);

// ------------System Parameter vectors---------------- Initialized in General.cpp
extern double* density_profile_buf;
extern double* density_profile;
extern double* density_profile_solvent_buf;
extern double* density_profile_solvent;

extern double* indiv_rog; // radius of gyration per chain
extern double* indiv_rog_buf;

extern double* tortuosity_buf;
extern double* tortuosity;

extern double* rog;
extern double* rog_buf;

extern double* rog_xy;
extern double* rog_xy_buf;

extern double* rog_z;
extern double* rog_z_buf;

extern long long* visits;
extern long long* visitdp;// unused
extern long long* visits_buf;

extern double*  endtoend;
extern  double* endtoend_buf;

// ------------values used for output maximum and minimum energy lattices found---------------- Utilized in Main.cpp
extern int minimum_en_config;
extern int maximum_en_config;
extern int printed;
extern int max_printed;


extern int loc_en_before;

//const int solvent_max = numberspins_S-( L1dim_S_xy*L1dim_S_xy*(plane_z_min+1)) -(L1dim_S_xy*L1dim_S_xy*(L1dim_S_z-plane_z_max));
//const int solvent_offset = L1dim_S_xy*L1dim_S_xy*(plane_z_min+1);

// ------------vectors for outputing representative microstructures----------------
extern double* output_table;
extern int table_length;
#endif
