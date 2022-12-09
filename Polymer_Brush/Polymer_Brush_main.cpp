// mpirun --oversubscribe -np 9 ./MPITest 0.8 1 8000 0789152 2 5 0
// Fully adapted to periodic boundary condtions utilized Dr. Bai's lecture notes in order to get everthing worked out correctly
// The only changes that should be left are verifying current equations
// Adapting the remaining function (should require no change just adaptations)
// Also running a small test case
// work on a local energy function this time
/*
  Replica Exchange Wang Landau demo code for simulating the 2D Potts model
  (c) Thomas Vogel and Ying Wai Li (2013-2019)

  License: Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
           https://creativecommons.org/licenses/by-sa/4.0/lecgalcode

       You are free to:

       - Share, copy and redistribute the material in any medium or format
       - Adapt, transform, and build upon the material for any purpose, even commercially

       (The licensor cannot revoke these freedoms as long as you
       follow the following license terms.)

       Under the following terms:

       1) Attribution: You must give appropriate credit and must
       not remove this preface from any file containing parts of
       this code. You may give credit in any reasonable manner,
       but not in any way that suggests the licensor endorses you
       or your use.

       If you publish data that was created using this code, or
       parts of it, we ask to cite the original papers:
         - T. Vogel et al., Phys. Rev. Lett. 110 (2013) 210603
             - T. Vogel et al., Phys. Rev. E 90 (2014) 023302
       Accompanying publications include:
         - T. Vogel et al., J. Phys.: Conf. Ser. 487 (2014) 014001
             - Y.W. Li et al., J. Phys.: Conf. Ser. 510 (2014) 014014

       2) ShareAlike: If you modify, transform, or build upon the
       material, you must distribute your contributions under the
       same license as the original.

       3) No additional restrictions: You may not apply legal terms
       or technological measures that legally restrict others from
       doing anything the license permits.

       Notices:

       You do not have to comply with the license for elements of
       the material in the public domain or where your use is
       permitted by an applicable exception or limitation.  No
       warranties are given. The license may not give you all of
       the permissions necessary for your intended use. For
       example, other rights such as publicity, privacy, or moral
       rights may limit how you use the material.

 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS# endif
#endif

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <iostream>
#include <float.h>
#include <random>

int confused = 0;
int confused2 = 0;

int switched = INT_MAX;

int uni_counter = 0;

int iterr;    // Apart of an added debug output. Effictively just a counter for the number of repeated histogram values when flatness is checked
int repv = 0;    // Holds the value of the (potentially) repeated minvalue of the histogram
int phrepv;    // Holds the value of the previous minimum histogram value to compare with the prior variable

// The lattice is effectively going to mimic a BCC Lattice
// I have code to track neighbors and second nearest neighbors for a BCC Lattice but the code models the unit cells exactly
// meaning that the end of the unit cell is cut off to make the modeled space periodic
// So given I would have to redo the math for neighboring sites regardless for this system
// I am implementing a different structure that utilizes a outer lattice (outer atoms of a BCC lattice, denoted as L because the size it's the largest lattice, and an inner lattice (interior atom of BCC (lattice), denoted S because it is one length in the three cardinal directions smaller than the larger lattice

int localgroup;
int local_dup_group;
int num_dup_proc;
int dup_headproc;



const int L1dim_L_xy = 30;    // linear dimension of the
const int L1dim_S_xy = 30;    // linear dimension5

const int L1dim_L_z = 30;    // linear dimension of the
const int L1dim_S_z = 30;    // linear dimension

const int L1dim_L = L1dim_L_xy;    // linear dimension of the
const int L1dim_S = L1dim_S_xy;    // linear dimension

const int Dimension = 3;    // not sure at this moment if other than 2D is actually implemented
const int numberspins_L = L1dim_L_xy *L1dim_L_xy * L1dim_L_z;    // total number of spins
const int numberspins_S = L1dim_S_xy *L1dim_S_xy * L1dim_S_z;    // total number of spins
const int numberneighbors_L = 8;    // for 2D square lattice (others not implemented right now)
const int numberneighbors_S = 6;    // for 2D square lattice (others not implemented right now)
const int numbercores = 39;
const int lengthpercore = 2;
const int pertrusionpercore = 1;    // number of arms per core

const int backbone_length = 45;
const int numberofchains=3;
const int lengthperchain=15;
const int startofchan[numberofchains] ={0,15,30};

const int min_max_id=3;

int reset_arm = -1;
int reset_ab = 0;

int *poly_lattice_indexes;

int *poly_spin_reference;

int *lattice_polymer;
int *poly_lattice_coordinates;
int *poly_lattice_connection;

const int Eglobalmin = 164701604;    // minimum energy for 2D square lattice Potts model
const int Eglobalmax = 173701604;                        // maximum energy of Potts model
const int Eglobalwidth = -(Eglobalmin - Eglobalmax);
const int bctype = 0;    // type of boundary condition: 0 - periodic; 1 - Braskamp Kunz
int q = 2;    // number of different possible spin states (q-state Potts model)


int *latticepoint_L;    // list containing values of all spins in the lattice (Checks occupation)
int *latticepoint_S;    // list containing values of all spins (checks occupation)
int *neighbor_L;    // list containing indices of neighbors for all spins (Contains Large lattice indexes for the smaller lattice
int *neighbor_S;    // list containing indices of neighbors for all spins

double *HE;    // energy histogram
double *lngE;    // ln g(E) estimator
double *lngE_buf;    // ln g(E) estimator buffer for exchange
double *pseudolngE;
double *real_lngE;
double *real_lngE_buf;
double *microT;
double *microT_buf;
int hist_size = (-Eglobalmin) + Eglobalmax;    // histogram size

int rseed;    // seed for random number generator
int energy, localenergy;

double Emin, Emax;    // doubles as calculation of boundaries in Energy is not in 'int'
int Eminindex, Emaxindex, Estartindex;    // local boundaries as index

// MPI; set up with local communicators for replica exchange (RE)
// Needed so that two processes can talk w/o depending on / bothering the others
int numprocs, myid, multiple, comm_id;
// each process belongs to two local groups,
// one to communicate to left neighbor and one to communicate to right neighbor
// each process has different loca IDs in different communicatore, in general
int mylocalid[2];    // id in local communicators
MPI_Comm * mpi_local_comm;
MPI_Group * mpi_local_group;
MPI_Group world;
MPI_Status status;
int merge_hists = 1;    // flag whether or not to merge histograms between iterations for multiple walkers on the same energy range

// to keep track of exchange statistics
int tryleft, tryright, exchangeleft, exchangeright;

// File handlers for I/O
FILE * file;
FILE * stdoutlog;
FILE * wanderlog;
char filename[50];
char filename2[50];
char resetbuffer[50];
char stdoutlogname[128];
char wanderlogname[128];
int ret_status;

double flatratio;
double flatmin;

int en_array[4][4][4];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

int en_array_coloumb_attr[4][4][4];
int en_array_coloumb_repul[4][4][4];
int long_range_ion_en_attr;
int long_range_ion_en_repul;

int en_array_bond[4][4][4];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

int int_en_array[14];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

int int_en_array_coloumb_attr[14];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 fo
int int_en_array_coloumb_repul[14];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 fo

int int_en_array_bond[14];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

int en_array_wall[4][4][4];    // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1


double dist_array[4][4][4];    // plan to hold the available distances in this array so the square root does not need to be called when values are needed
int int_dist_array[4][4][4];    // plan to hold the available distances in this array so the square root does not need to be called when values are needed


int no_pos_electrostatic_inter[2];
int no_neg_electrostatic_inter[2];

int en_array_PMPM[4][4][4];
int en_array_PMPO[4][4][4];
int en_array_PMO[4][4][4];
int en_array_PMH[4][4][4];
int en_array_PMEu[4][4][4];
int int_en_array_PMPM[14];
int int_en_array_PMPO[14];
int int_en_array_PMO[14];
int int_en_array_PMH[14];
int int_en_array_PMEu[14];

int en_array_POPO[4][4][4];
int en_array_POO[4][4][4];
int en_array_POH[4][4][4];
int en_array_POEu[4][4][4];
int int_en_array_POPO[14];
int int_en_array_POO[14];
int int_en_array_POH[14];
int int_en_array_POEu[14];


int en_array_OO[4][4][4];
int en_array_OH[4][4][4];
int en_array_OEu[4][4][4];
int int_en_array_OO[14];
int int_en_array_OH[14];
int int_en_array_OEu[14];


int en_array_HH[4][4][4];
int en_array_HEu[4][4][4];
int int_en_array_HH[14];
int int_en_array_HEu[14];


int en_array_EuEu[4][4][4];
int int_en_array_EuEu[14];

int en_array_OH_Att[4][4][4];
int en_array_OEu_Att[4][4][4];
int int_en_array_OH_Att[14];
int int_en_array_OEu_Att[14];

int en_array_OO_Rep[4][4][4];
int en_array_HH_Rep[4][4][4];
int int_en_array_OO_Rep[14];
int int_en_array_HH_Rep[14];

int en_array_EuEu_Rep[4][4][4];
int en_array_HEu_Rep[4][4][4];
int int_en_array_EuEu_Rep[14];
int int_en_array_HEu_Rep[14];

int long_range_OH_en_attr;
int long_range_OEu_en_attr;
int long_range_OO_en_repul;
int long_range_HH_en_repul;
int long_range_EuEu_en_repul;
int long_range_HEu_en_repul;
int no_OH_electrostatic_inter[2];
int no_OEu_electrostatic_inter[2];
int no_OO_electrostatic_inter[2];
int no_HH_electrostatic_inter[2];
int no_EuEu_electrostatic_inter[2];
int no_HEu_electrostatic_inter[2];
long long iqq = 0;
// to track execution times
time_t timenow, timestart, timeend;

void keypressed()    // just for developing / manual debugging
{
    for (;;)
        if (getchar() == 27) break;
}

int anint(int, int);
int anint(float);

//general set up
void init_poly_cores();    // init the polymer coordinates
void init_en_and_dist_array();    // initialize the energy and the distance array
void poly_coord_decompose(int, int);    // decomompes a given polymer index into cartesian coordinates
double point_distance(int, int);    //gives the distance between two points
int int_point_distance(int, int);
double point_distance(int, int, int(&arr)[3]);    //gives the distance between two points
double poly_point_distance(int *, int, int);    //gives the distance between two points

double poly_point_distance(int *, int, int, int(&arr)[3]);    //gives the distance between two points
int general_coord(int(&arr)[3], int);    // output the x y and z compenent for any index
int general_index(int(&arr)[3]);    // output index from any x, y, z components
int general_index_op(int(&arr)[3]);    // output index from any x, y, z components
int poly_latticepoint_clear(int *, int, int);    // clears a set group of lattice points (to prevent duplication of the lattice to perform a move
int poly_latticepoint_clear(int);
int latticepoint_assign(int, int);    // assigns lattice point based on index values taken from argument array (also checks for occupation will output 0 if assignment overlap will occur

int latticepoint_assign2(int);    // assigns lattice point based on index values taken from argument array (also checks for occupation will output 0 if assignment overlap will occur

int mov_crankshaft();
int mov_crankshaft_backbone();
int mov_pivot();
int mov_inversion();
int mov_planar_reflection();

int mov_translation();
int mov_translation_backbone();

int poly_mov();    // Selects movement options for polymer

int *solvent_loc;
int *ph_solvent_loc;
int *wl_solvent_loc;

int *ion_loc;
int *ph_ion_loc;
int *wl_ion_loc;

int init_solvent(int *, int *, int *, int, int);
int solvent_translation(int *, int *, int, int);
int solvent_translation2(int *, int *, int, int);
int no_solvent_sites =900;

int no_ion_sites = 15;
int attachment(int, int);
int reattachment(int, int);
int detachment(int);
int *attachment_ref;
int *attachment_ion_ref;

int *ph_poly_lattice_indexes;
int *ph_poly_lattice_coordinates;
int *wl_pseudo_chain;
int *wl_pseudo_chain_coordinates;

void lattice_polymer_reset(int *p_chain, int *c_chain, int start, int end);
int lattice_poly_index_reset(int *p_chain, int *c_chain, int start, int end);
int poly_coord_reset(int *p_chain, int *c_chain, int start, int end);
void solvent_reset(int *, int *, int);
void solvent_reset2(int *, int *, int);

void accessiblelevels();
void sysparam(int);

const double max_dist = sqrt(11);    //maximum allowable distance any nieghboring monomers can be from each other
const int int_max_dist = 11;    //maximum allowable distance any nieghboring monomers can be from each other
const int plane_z_max = L1dim_S_z - 1;    // axis upon which brush plane exists. Is a set value to check for any chains that try to go beneath it
const int plane_z_min = 3;
int MoveProposal = 0;    // tracks if a move has been performed

//checks to calculate local energy
int poly_solvent_mov = 0;    // says wether a chosen mov moved solvent or a monomer
int poly_ion_mov = 0;    // says wether a chosen mov moved solvent or a monomer
int chosen_mov = 0;    // tellw which move was last chosen
int poly_mov_start = 0;    // tells which monomwrs were involved with a given move the starting position (on thr chain)
int poly_mov_end = 0;    // tells which monomer is the end considered loop

void optomize();

int reset_indexes[7] = { 0, 0, 0, 0, 0, 0, 0 };

const double max_en_dist = sqrt(13.5);    //maximum checked distance for any monomer interaction
const int int_max_en_dist = 13;    //maximum checked distance for any monomer interaction

const int int_max_attach_dist = 8;

std::mt19937 gen(74);    // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution <> polymermovment(0, 99);

int *lattice_dist_optomize;
int* lattice_dist_optomize_2;
int* lattice_dist_optomize_2_lattice;

int *index_optomize_x;
int *index_optomize_y;
int *index_optomize_z;

int *calc_index_optomize_x;
int *calc_index_optomize_y;
int *calc_index_optomize_z;

int *polymer_optomize;

int what_to_exclude[176*numberneighbors_L]={176};
int* local_en_indexes_op;
int local_en_dist[176]={0};
int local_en_indexes_op_x[176]={0};
int local_en_indexes_op_y[176]={0};
int local_en_indexes_op_z[176]={0};

int *xyz_index;
int *index_x;
int *index_y;
int *index_z;

const int neg_xy = L1dim_S_xy * 6;
const int neg_z = L1dim_S_z * 3;

const int lower_bounds = L1dim_S_xy *L1dim_S_xy *(plane_z_max - 1);
const int upper_bounds = L1dim_S_xy *L1dim_S_xy *(plane_z_min + 1);


int current_attach=0;
int current_attachment[numbercores*pertrusionpercore*lengthpercore] = {0};
int vicinity_check=0;

double *tortuosity_buf;
double *tortuosity;

double *rog;
double *rog_buf;

double *tortuositybb_buf;
double *tortuositybb;

double *rogbb;
double *rogbb_buf;


float* radial_dist_EuPO;
float* radial_dist_EuPO_buf;
float* radial_dist_EuEu;
float* radial_dist_EuEu_buf;
float* radial_dist_EuO;
float* radial_dist_EuO_buf;
float* radial_dist_PMPO;
float* radial_dist_PMPO_buf;
float* radial_dist_POPO;
float* radial_dist_POPO_buf;

float* radial_visits_ion;
float* radial_visits_poly;
float* radial_visits_ion_buf;
float* radial_visits_poly_buf;

int* part_seq;
        float* func_per_seq;
            int* perf_seq;
int* part_seq_buf;
float* func_per_seq_buf;
    int* perf_seq_buf;

int index_track=-1;

long long * visits;
long long * visits_buf;

double *endtoend;
double *endtoend_buf;

int minimum_en_config =  2147483000;
int maximum_en_config = 0;
int printed = 0;
int max_printed = 0;

int loc_en_before = 0;
int loc_en_before_2 = 0;

int core_movement = 0;
int random_index;    // chooses a random insex
//int random_index = (rand()%solvent_max) + solvent_offset;    // chooses a random insex
int random_solvent_site;
//const int solvent_max = numberspins_S-(L1dim_S_xy*L1dim_S_xy*(plane_z_min+1)) -(L1dim_S_xy*L1dim_S_xy*(L1dim_S_z-plane_z_max));
//const int solvent_offset = L1dim_S_xy*L1dim_S_xy*(plane_z_min+1);
int local_energy(int *p_l_c, int *p_l_i, int *s_l, int *i_l,int inc) ;

void init_neighbors()    // create neighbor list first for smaller lattice
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to init neighbors()\n");
    fclose(stdoutlog);
    // neighbor contains the index of the neighboring spin
    // for each spin there are four neighbors in this order: above, right, below, left
    neighbor_S = (int*) malloc(numberspins_S *numberneighbors_S* sizeof(int));

    for (int i = 0; i < numberspins_S; i++)    // in general
    {
        neighbor_S[numberneighbors_S *i] = i - L1dim_S_xy * L1dim_S_xy;    // out of plane
        neighbor_S[numberneighbors_S *i + 1] = i - L1dim_S_xy;    // above
        neighbor_S[numberneighbors_S *i + 2] = i + 1;    // right
        neighbor_S[numberneighbors_S *i + 3] = i + L1dim_S_xy * L1dim_S_xy;    // into plane
        neighbor_S[numberneighbors_S *i + 4] = i + L1dim_S_xy;    // below
        neighbor_S[numberneighbors_S *i + 5] = i - 1;    // left
    }

    if (bctype == 0)    // periodic BC
        for (int i = 0; i < numberspins_S; i++)    // now treat boundaries separately
    {
        if (i < L1dim_S_xy *L1dim_S_xy)    // highest plane
            neighbor_S[numberneighbors_S *i] = i + (L1dim_S_xy *L1dim_S_xy *(L1dim_S_z - 1));
        //neighbor_S[numberneighbors_S *i] = numberspins_S;

        if (i % (L1dim_S_xy *L1dim_S_xy) < L1dim_S_xy)    // top row
            neighbor_S[numberneighbors_S *i + 1] = i + (L1dim_S_xy *(L1dim_S_xy - 1));
        //neighbor_S[numberneighbors_S *i] = numberspins_S;

        if ((i + 1) % L1dim_S_xy == 0)    // rightmost column
            neighbor_S[numberneighbors_S *i + 2] = i - (L1dim_S_xy - 1);
        //neighbor_S[numberneighbors_S *i] = numberspins_S;

        if (i >= L1dim_S_xy *L1dim_S_xy *(L1dim_S_z - 1))    // deepest row
            neighbor_S[numberneighbors_S *i + 3] = i - (L1dim_S_xy *L1dim_S_xy *(L1dim_S_z - 1));
        //neighbor_S[numberneighbors_S *i] = numberspins_S;

        if ((i % (L1dim_S_xy *L1dim_S_xy)) > (L1dim_S_xy *(L1dim_S_xy - 1) - 1))    // bottom row
            neighbor_S[numberneighbors_S *i + 4] = i - (L1dim_S_xy *(L1dim_S_xy - 1));
        //neighbor_S[numberneighbors_S *i] = numberspins_S;

        if (i % L1dim_S_xy == 0)    // leftmost column
            neighbor_S[numberneighbors_S *i + 5] = i + (L1dim_S_xy - 1);
        //neighbor_S[numberneighbors_S *i] = numberspins_S;
    }

    //print neighbor list (for manual debugging)
    //  for (int i=0; i<numberspins*numberneighbors; i++)
    //    {
    //      printf("%4d", neighbor[i]);
    //      if ((i+1)%numberneighbors == 0) printf("\n");
    //    };

    // contains indexes of larger lattice in indexes for the smaller lattice
    neighbor_L = (int*) malloc(numberspins_S *numberneighbors_L* sizeof(int));

    int reference = 0;

    for (int i = 0; i < numberspins_L; i++)    // in general
    {
        reference = i;

        int next_plane_ref = ((reference + L1dim_L_xy *L1dim_L_xy) % (L1dim_L_xy *L1dim_L_xy *L1dim_L_z));

        neighbor_L[numberneighbors_L *i] = reference;    // upper topleft
        neighbor_L[numberneighbors_L *i + 1] = reference + ((reference + 1) % L1dim_L_xy - ((reference) % L1dim_L_xy));;    // upper topright
        neighbor_L[numberneighbors_L *i + 2] = (reference + L1dim_L_xy) % (L1dim_L_xy *L1dim_L_xy) + (reference / (L1dim_L_xy *L1dim_L_xy)) *L1dim_L_xy * L1dim_L_xy;    // upper bottom left
        neighbor_L[numberneighbors_L *i + 3] = ((reference + L1dim_L_xy) % (L1dim_L_xy *L1dim_L_xy)) + ((reference + 1) % L1dim_L_xy - ((reference) % L1dim_L_xy)) + (reference / (L1dim_L_xy *L1dim_L_xy)) *L1dim_L_xy * L1dim_L_xy;    // upper bottom right
        neighbor_L[numberneighbors_L *i + 4] = next_plane_ref;    //lower top right
        neighbor_L[numberneighbors_L *i + 5] = next_plane_ref + ((next_plane_ref + 1) % L1dim_L_xy - ((next_plane_ref) % L1dim_L_xy));    // lower topright
        neighbor_L[numberneighbors_L *i + 6] = (next_plane_ref + L1dim_L_xy) % (L1dim_L_xy *L1dim_L_xy) + (next_plane_ref / (L1dim_L_xy *L1dim_L_xy)) *L1dim_L_xy * L1dim_L_xy;    // lower bottom left
        neighbor_L[numberneighbors_L *i + 7] = ((next_plane_ref + L1dim_L_xy) % (L1dim_L_xy *L1dim_L_xy)) + ((next_plane_ref + 1) % L1dim_L_xy - ((next_plane_ref) % L1dim_L_xy)) + (next_plane_ref / (L1dim_L_xy *L1dim_L_xy)) *L1dim_L_xy * L1dim_L_xy;    // lower bottom right
    }

    //print neighbor list (for manual debugging)
    stdoutlog = fopen(stdoutlogname, "a");
    int q = 0;
    fprintf(stdoutlog, "\ns%i: ", q);
    for (int i = 0; i < numberspins_S * numberneighbors_S; i++)
    {
        fprintf(stdoutlog, " %4d ", neighbor_S[i]);
        if ((i + 1) % numberneighbors_S == 0)
        {
            q++;
            fprintf(stdoutlog, "\n%i: ", q);
        };

    };
    
    stdoutlog = fopen(stdoutlogname, "a");
    q = 0;
    fprintf(stdoutlog, "\nl%i: ", q);
    for (int i = 0; i < numberspins_L * numberneighbors_L; i++)
    {
        fprintf(stdoutlog, " %4d ", neighbor_L[i]);
        if ((i + 1) % numberneighbors_L == 0)
        {
            q++;
            fprintf(stdoutlog, "\n%i: ", q);
        };

    };

    fclose(stdoutlog);
}

int init_solvent_sub2(int index, int spin,int &dipole_index,int &dipole_index2)
{
    int check=0;
    
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        if (latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] != 0)    // checks to ensure lattice is assignable
        {
            return 0;
        }
    }
    
    dipole_index = neighbor_S[numberneighbors_S*index+2];
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        if (latticepoint_L[neighbor_L[numberneighbors_L *dipole_index + iii]] != 0)    // checks to ensure lattice is assignable
        {
            return 0;
        }
    }
    
    dipole_index2 = neighbor_S[numberneighbors_S*index+1];
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        if (latticepoint_L[neighbor_L[numberneighbors_L *dipole_index2 + iii]] != 0)    // checks to ensure lattice is assignable
        {
            return 0;
        }
    }

    return 1;
}

int init_solvent2(int *s_l, int *ph_s_l, int *wl_s_l, int number, int spin)    // randomly instance solvent atoms
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "start solvent 2\n");
    fclose(stdoutlog);
        int counter = 0;    // number of solvent placed
        for (int index = 0; index < numberspins_S; index++)    // sequences through lattice
        {
            if (counter < number)    //&& index>=((upper_bounds)) && index<(L1dim_S_xy*L1dim_S_xy*L1dim_S_z-L1dim_S_xy*L1dim_S_xy))
            {
                // while all solvent has not been placed and while the solvent is within a given range

                if (latticepoint_S[index] == 0)    // checks for if assignment is possible
                {
                    int neighbor_check = 0;    // checks if neighbor assignment is possible
                    int dipole_index = 0;
                    int dipole_index2 = 0;
                    neighbor_check = init_solvent_sub2(index, spin,dipole_index,dipole_index2);    // calls sub function to reduce size of function

                    if (neighbor_check == 1)    // if assignment was possible
                    {
                        latticepoint_S[index] = spin;    // assigns solvent spin
                        latticepoint_S[dipole_index] = -spin;
                        latticepoint_S[dipole_index2] = -spin;
                        
                        stdoutlog = fopen(stdoutlogname, "a");
                        fprintf(stdoutlog, "\t%i\t%i\t%i\n",counter,index,dipole_index);
                        fclose(stdoutlog);
                        
                        for (int iii = 0; iii < numberneighbors_L; iii++)
                        {
                            latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] = spin;
                            latticepoint_L[neighbor_L[numberneighbors_L *dipole_index + iii]] = -spin;
                            latticepoint_L[neighbor_L[numberneighbors_L *dipole_index2 + iii]] = -spin;
                        }
                        
                        s_l[3*counter] = index;    // assigns the spin to solvent_loc and increments
                        ph_s_l[3*counter] = index;
                        wl_s_l[3*counter] = index;
                        
                        s_l[3*counter+1] = dipole_index;    // assigns the spin to solvent_loc and increments
                        ph_s_l[3*counter+1] = dipole_index;
                        wl_s_l[3*counter+1] = dipole_index;
                        
                        s_l[3*counter+2] = dipole_index2;    // assigns the spin to solvent_loc and increments
                        ph_s_l[3*counter+2] = dipole_index2;
                        wl_s_l[3*counter+2] = dipole_index2;
                        counter++;
                    }
                }
            }
        }
    

    
    reset_indexes[0]=-1;
    reset_indexes[1]=3*no_solvent_sites;
    solvent_reset2(s_l, ph_s_l, spin);
    solvent_reset2(wl_s_l, s_l, spin);

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "exit solvent 2\n");
    fclose(stdoutlog);
    return 0;
}

int init_solvent_sub(int index, int spin)
{
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        if (latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] == 0)    // checks to ensure lattice is assignable
        {
            latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] = spin;
        }
        else    // resets the points switched to occupaton in the case where occupation is not allowed. Just performs cleanup itself so no other function has to
        {
            latticepoint_S[index] = 0;
            for (int j = 0; j < iii; j++)
            {
                latticepoint_L[neighbor_L[numberneighbors_L *index + j]] = 0;
            }

            return 0;
        }
    }

    return 1;
}

int init_solvent(int *s_l, int *ph_s_l, int *wl_s_l, int number, int spin)    // randomly instance solvent atoms
{
    char solventbuffer[50];
    if (myid < min_max_id)
    {

        int counter = 0;    // number of solvent placed
        for (int index = 0; index < numberspins_S; index++)    // sequences through lattice
        {
            if (counter < number)    //&& index>=((upper_bounds)) && index<(L1dim_S_xy*L1dim_S_xy*L1dim_S_z-L1dim_S_xy*L1dim_S_xy))
            {
                // while all solvent has not been placed and while the solvent is within a given range

                if (latticepoint_S[index] == 0)    // checks for if assignment is possible
                {
                    latticepoint_S[index] = spin;    // assigns solvent spin

                    int neighbor_check = 0;    // checks if neighbor assignment is possible
                    neighbor_check = init_solvent_sub(index, spin);    // calls sub function to reduce size of function

                    if (neighbor_check == 1)    // if assignment was possible
                    {
                        stdoutlog = fopen(stdoutlogname, "a");
                        fprintf(stdoutlog, "\t%i\t%i\n",counter,index);
                        fclose(stdoutlog);
                        s_l[counter] = index;    // assigns the spin to solvent_loc and increments
                        ph_s_l[counter] = index;
                        wl_s_l[counter] = index;
                        counter++;
                    }
                }
            }
        }

    }

    if (myid >= min_max_id)
    {
        int counter = 0;    // number of solvent placed
        for (int index = 0; index < numberspins_S; index++)    // sequences through lattice
        {
            if (counter < number)    //&& index>=((upper_bounds)) && index<(L1dim_S_xy*L1dim_S_xy*L1dim_S_z-L1dim_S_xy*L1dim_S_xy))
            {
                // while all solvent has not been placed and while the solvent is within a given range

                if (latticepoint_S[index] == 0)    // checks for if assignment is possible
                {
                    latticepoint_S[index] = spin;    // assigns solvent spin

                    int neighbor_check = 0;    // checks if neighbor assignment is possible
                    neighbor_check = init_solvent_sub(index, spin);    // calls sub function to reduce size of function

                    if (neighbor_check == 1)    // if assignment was possible
                    {
                        stdoutlog = fopen(stdoutlogname, "a");
                        fprintf(stdoutlog, "\t%i\t%i\n",counter,index);
                        fclose(stdoutlog);
                        s_l[counter] = index;    // assigns the spin to solvent_loc and increments
                        ph_s_l[counter] = index;
                        wl_s_l[counter] = index;
                        counter++;
                    }
                }
            }
        }
    }
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "reset ion 1\n");
    fclose(stdoutlog);
    reset_indexes[0]=-1;
    reset_indexes[1]=no_ion_sites;
    solvent_reset(s_l, ph_s_l, spin);
    solvent_reset(wl_s_l, s_l, spin);
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "freset ion 1\n");
    fclose(stdoutlog);
    return 0;
}


int attachment(int ion_position, int ion_index)
{
    int random = rand() % numbercores;
    for (int i = 0; i < numbercores; i++)
    {
        //random=random%numbercores;

        //            stdoutlog = fopen(stdoutlogname, "a");
        //            fprintf(stdoutlog, "\n\t\tdistance_calc: %i \t poly_pos: %i ",lattice_dist_optomize[ion_index*numberspins_S+poly_lattice_indexes[((i*pertrusionpercore*lengthpercore)+(lengthpercore-1))]],((i*pertrusionpercore*lengthpercore)+(lengthpercore-1)));
        //            fclose(stdoutlog);

        if (int_point_distance(ion_index, poly_lattice_indexes[((i *pertrusionpercore *lengthpercore) + (lengthpercore - 1))]) <= int_max_attach_dist && attachment_ion_ref[2 *ion_position] == 0 && attachment_ref[2 *((i *pertrusionpercore *lengthpercore) + (lengthpercore - 1))] == 0)
        {
            attachment_ion_ref[2 *ion_position] = 1;
            attachment_ion_ref[2 *ion_position + 1] = i *lengthpercore *pertrusionpercore + (lengthpercore - 1);

            attachment_ref[2 *((i *pertrusionpercore *lengthpercore) + (lengthpercore - 1))] = 1;
            attachment_ref[2 *((i *pertrusionpercore *lengthpercore) + (lengthpercore - 1)) + 1] = ion_position;

            return 1;
        }

        //random++;
    }

    return 0;
}

int attachment_poly(int poly_position, int poly_index)
{
    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "aapttach\n");
    //    fclose(stdoutlog);
    int random = rand() % numbercores;
    for (int i = 0; i < no_ion_sites; i++)
    {
        //random=random%numbercores;
        //attach  1  2  575  575     8  10 f:2
        //            stdoutlog = fopen(stdoutlogname, "a");
        //            fprintf(stdoutlog, "\n\tpolytdistance_calc: %i \t ion_pos: %i \t%i %i\t %i %i\n attachment_ion_ref[2*i]:%i", lattice_dist_optomize[ion_loc[i]*numberspins_S+ph_poly_lattice_indexes[poly_position]],i,poly_position,poly_index,ion_loc[i],ph_ion_loc[i], attachment_ion_ref[2*i]);
        //            fclose(stdoutlog);

        if (int_point_distance(poly_index, ion_loc[i]) <= int_max_attach_dist && attachment_ref[2 *poly_position] == 0 && attachment_ion_ref[2 *i] == 0)
        {
            random_solvent_site = i;
            random_index = solvent_loc[i];

            return 1;
        }

        //random++;
    }

    return 0;
}

int detachment(int ion_position)
{
    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "ddrettach\n");
    //    fclose(stdoutlog);
    int i = attachment_ion_ref[2 *ion_position + 1];

    attachment_ref[2 *i] = 0;
    attachment_ref[2 *i + 1] = 0;

    attachment_ion_ref[2 *ion_position] = 0;
    attachment_ion_ref[2 *ion_position + 1] = 0;

    return 0;
}

int reattachment(int ion_position, int poly_position)
{
    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "rettach\n");
    //    fclose(stdoutlog);

    attachment_ion_ref[2 *ion_position] = 1;
    attachment_ion_ref[2 *ion_position + 1] = poly_position;

    //int i = attachment_ion_ref[2*ion_position+1];

    attachment_ref[2 *poly_position] = 1;
    attachment_ref[2 *poly_position + 1] = ion_position;

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "\nreattach  ion_position:%i attachment_ref[2*poly_position] %i\n",ion_position,attachment_ref[2*poly_position+1]);
    //    fclose(stdoutlog);

    return 0;
}

int solvent_translation2(int *s_l, int *ph_s_l, int no_sites, int spin)    // attempts to move a solvent molecule to a different location
{
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Made it to solvent translation  %i   %i\n",spin,latticepoint_L[0]);
//        fclose(stdoutlog);

    random_index = rand() % numberspins_S;    // chooses a random insex
    //int random_index = (rand()%solvent_max) + solvent_offset;    // chooses a random insex
    random_solvent_site = rand() % no_sites;    // chooses random_solvent molecule

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "random_index:%i   random_solvent_site:%i   solvent_loc[random_solvent_site]: %i \n\n",random_index,random_solvent_site,solvent_loc[random_solvent_site]);
    //    fclose(stdoutlog);

    reset_indexes[0] = 3*random_solvent_site - 1;
    reset_indexes[1] = 3*random_solvent_site + 3;
    
    poly_latticepoint_clear(s_l[3*random_solvent_site]);    // clears the solvent in the previous localtion
    poly_latticepoint_clear(s_l[3*random_solvent_site+1]);
    
    poly_latticepoint_clear(s_l[3*random_solvent_site+2]);
    
    int rand_neigh = rand()%numberneighbors_S;
    int rand_neigh2 = (rand_neigh+(rand()%(numberneighbors_S-1)))%numberneighbors_S;
    while(rand_neigh2==(rand_neigh+3)%numberneighbors_S || rand_neigh2==rand_neigh)
    {
        rand_neigh2 = rand()%numberneighbors_S;
    }
    int new_dipole_index = neighbor_S[numberneighbors_S*random_index+rand_neigh];
    int new_dipole_index2 = neighbor_S[numberneighbors_S*random_index+rand_neigh2];
    
    
    int checker = 1;
    for(int i=0;i<numberneighbors_L;i++)
    {
        if(latticepoint_L[neighbor_L[numberneighbors_L*new_dipole_index+i]]!=0)
        {
            checker = 0;
        }
        if(latticepoint_L[neighbor_L[numberneighbors_L*new_dipole_index2+i]]!=0)
        {
            checker = 0;
        }
    }
    
    int assignment_possible;
    if(checker ==1)
    {
         assignment_possible = latticepoint_assign2(random_index);    // checks to see if assignment was possible
    }
    else{
        assignment_possible = 0;
    }

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "assignment possible %i\n",assignment_possible);
    //    fclose(stdoutlog);
    reset_indexes[4] = 0;
    if (assignment_possible)    //if assignment is possible
    {
        loc_en_before = local_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc,0);
        
        ph_s_l[3*random_solvent_site] = random_index;    // new index is assigned to solvent_loc
        ph_s_l[3*random_solvent_site+1] = new_dipole_index;
        ph_s_l[3*random_solvent_site+2] = new_dipole_index2;
        
        latticepoint_S[new_dipole_index]=-2;
        for(int i=0;i<numberneighbors_L;i++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L*new_dipole_index+i]]=-2;
        }
        latticepoint_S[new_dipole_index2]=-2;
        for(int i=0;i<numberneighbors_L;i++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L*new_dipole_index2+i]]=-2;
        }

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "random_index:%i   random_solvent_site:%i   solvent_loc[random_solvent_site]: %i   new: %i  %i  %i \n\n",random_index,random_solvent_site,ph_solvent_loc[random_solvent_site],new_dipole_index,ph_s_l[2*random_solvent_site],ph_s_l[2*random_solvent_site+1]);
//                fclose(stdoutlog);

        return 1;
    }
    else
    {
        reset_indexes[0] = 3*random_solvent_site - 1;
        reset_indexes[1] = 3*random_solvent_site + 3;
        solvent_reset2(ph_s_l, s_l, spin);
    }

    return 0;
}

int solvent_translation(int *s_l, int *ph_s_l, int no_sites, int spin)    // attempts to move a solvent molecule to a different location
{
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Made it to solvent translation  %i   %i\n",spin,latticepoint_L[0]);
//        fclose(stdoutlog);

    random_index = rand() % numberspins_S;    // chooses a random insex
    //int random_index = (rand()%solvent_max) + solvent_offset;    // chooses a random insex
    random_solvent_site = rand() % no_sites;    // chooses random_solvent molecule

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "random_index:%i   random_solvent_site:%i   solvent_loc[random_solvent_site]: %i \n\n",random_index,random_solvent_site,solvent_loc[random_solvent_site]);
    //    fclose(stdoutlog);

    reset_indexes[0] = random_solvent_site - 1;
    reset_indexes[1] = random_solvent_site + 1;

    poly_latticepoint_clear(s_l[random_solvent_site]);    // clears the solvent in the previous localtion
    int assignment_possible = latticepoint_assign2(random_index);    // checks to see if assignment was possible

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "assignment possible %i\n",assignment_possible);
    //    fclose(stdoutlog);
    reset_indexes[4] = 0;
    if (assignment_possible)    //if assignment is possible
    {
        loc_en_before = local_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc,0);
        
        ph_s_l[random_solvent_site] = random_index;    // new index is assigned to solvent_loc
        return 1;
    }
    else
    {
//        reset_indexes[0] = random_solvent_site - 1;
//        reset_indexes[1] = random_solvent_site + 1;
        solvent_reset(ph_s_l, s_l, spin);
    }

    return 0;
}

void solvent_reset2(int *p_index, int *c_index, int spin)
{
    for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++)    // clears the lattice polymer in the range set
    {
        latticepoint_S[p_index[i]] = 0;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *p_index[i] + iii]] = 0;
        }
    }
    
    for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++)    // clears the lattice polymer in the range set
    {
        if(i%3==0) spin =2;
        if(i%3!=0) spin = -2;
        
        latticepoint_S[c_index[i]] = spin;
        
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *c_index[i] + iii]] = spin;
        }

        p_index[i] = c_index[i];
    }
}

void solvent_reset(int *p_index, int *c_index, int spin)
{
    for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++)    // clears the lattice polymer in the range set
    {
        latticepoint_S[p_index[i]] = 0;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *p_index[i] + iii]] = 0;
        }
    }

    for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++)    // clears the lattice polymer in the range set
    {
        latticepoint_S[c_index[i]] = spin;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *c_index[i] + iii]] = spin;
        }

        p_index[i] = c_index[i];
    }
}

int anint(int a, int b)
{
    float c = (float)(a) / ((float) b);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to anint()  c:%f\n",c);
    //fclose(stdoutlog);

    if (c >= 0)
    {
        return ((int)(c + .5));
    }

    if (c < 0)
    {
        return ((int)(c - .5));
    }

    if (c < -0.5 && c > -1.5)
    {
        return -1;
    }

    if (c <= 0.5 && c >= -0.5)
    {
        return 0;
    }

    if (c < 1.5 && c > 0.5)
    {
        return 1;
    }

    return 0;
}

int anint(float c)
{
    if (c >= 0)
    {
        return ((int)(c + .5));
    }

    if (c < 0)
    {
        return ((int)(c - .5));
    }

    if (c < -0.5 && c > -1.5)
    {
        return -1;
    }

    if (c <= 0.5 && c >= -0.5)
    {
        return 0;
    }

    if (c < 1.5 && c > 0.5)
    {
        return 0;
    }

    return 0;
}

void init_optomize_energy_index()
{
    int index = 3*3+4+4*3+4;
        int c_x;
        int c_y;
        int c_z;
        int counter_2=0;
        c_x = index_x[index];
        c_y = index_y[index];
        c_z = index_z[index];
        
        int hold_indexes[176]={0};
        
    
        for(int a =-3;a<=3;a++)
        {
            for(int b = -3;b<=3;b++)
            {
                for(int c = -3 ; c<=3;c++)
                {
                    int checked_index = xyz_index[(a+c_z+5)*(L1dim_S_xy+10)*(L1dim_S_xy+10)+(b+c_y+5)*(L1dim_S_xy+10)+(c+c_x+5)];
                    
                    if(point_distance(checked_index,index)<max_en_dist && point_distance(checked_index,index)>=2)
                    {
                        local_en_indexes_op_x[counter_2]=c;
                        local_en_indexes_op_y[counter_2]=b;
                        local_en_indexes_op_z[counter_2]=a;
                        hold_indexes[counter_2]=checked_index;
                        counter_2++;
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "ioei %i: {%i,%i,%i} {%i,%i,%i}\n\t",checked_index,c,b,a,local_en_indexes_op_x[counter_2-1],local_en_indexes_op_y[counter_2-1],local_en_indexes_op_z[counter_2-1]);
//                        fclose(stdoutlog);
                        
                        for(int i=0;i<numberneighbors_L;i++)
                        {
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "%i  ",neighbor_L[checked_index*numberneighbors_L+i]);
//                            fclose(stdoutlog);
                        }
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "\n");
//                        fclose(stdoutlog);
                    }
                }
            }
        }
        
    for(int q =0;q<176;q++)
    {
        for(int r=1;r<numberneighbors_L;r++)
        {
            what_to_exclude[q*numberneighbors_L+r]=176;
            for(int g = 0;g<176;g++)
            {
                if(hold_indexes[g] == neighbor_L[hold_indexes[q]*numberneighbors_L+r])
                {
                    what_to_exclude[q*numberneighbors_L+r]=g;
                }
            }
        }
    }
    
    for(int q =0;q<176;q++)
    {
        stdoutlog = fopen(stdoutlogname, "a");
                               fprintf(stdoutlog, "%i:  %i\n\t",q,hold_indexes[q]);
                               fclose(stdoutlog);
        for(int r=1;r<numberneighbors_L;r++)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "%i  ",what_to_exclude[q*numberneighbors_L+r]);
            fclose(stdoutlog);
        }
        
        stdoutlog = fopen(stdoutlogname, "a");
                               fprintf(stdoutlog, "\n");
                               fclose(stdoutlog);
    }
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\n\n\n%i\n",counter_2);
        fclose(stdoutlog);
    
   //MPI_Abort(MPI_COMM_WORLD,1);
}


void optomize_energy_index()
{
    local_en_indexes_op = (int*)malloc(numberspins_S*176 * sizeof(int));
    
    for(int y = 0;y<numberspins_S;y++)
    {
        int index = y;
        int c_x;
        int c_y;
        int c_z;
        int counter_2=0;
        c_x = index_x[index];
        c_y = index_y[index];
        c_z = index_z[index];
    
        for(int a =-3;a<=3;a++)
        {
            for(int b = -3;b<=3;b++)
            {
                for(int c = -3 ; c<=3;c++)
                {
                    int checked_index = xyz_index[(a+c_z+5)*(L1dim_S_xy+10)*(L1dim_S_xy+10)+(b+c_y+5)*(L1dim_S_xy+10)+(c+c_x+5)];
                    
                    if(point_distance(checked_index,index)<max_en_dist && point_distance(checked_index,index)>=2)
                    {
                        local_en_indexes_op[index*176+counter_2]=checked_index;
                        local_en_dist[counter_2]=int_point_distance(checked_index,index);
                        counter_2++;
                    }
                }
            }
        }
    }
}


void optomize()
{
    //lattice_dist_optomize = (int*)malloc(L1dim_S_xy *L1dim_S_xy *L1dim_S_z *L1dim_S_xy *L1dim_S_xy *L1dim_S_z* sizeof(int));

    index_optomize_x = (int*) malloc((L1dim_S_xy *12) *sizeof(int));
    index_optomize_y = (int*) malloc((L1dim_S_xy *12) *sizeof(int));
    index_optomize_z = (int*) malloc((L1dim_S_z *6) *sizeof(int));

    calc_index_optomize_x = (int*) malloc((L1dim_S_xy *12) *sizeof(int));
    calc_index_optomize_y = (int*) malloc((L1dim_S_xy *12) *sizeof(int));
    calc_index_optomize_z = (int*) malloc((L1dim_S_z *6) *sizeof(int));

    //polymer_optomize = (int*)malloc(L1dim_S_xy *L1dim_S_xy *L1dim_S_z* sizeof(int));

    //    for(int i=0;i < numberspins_S;i++)
    //    {
    //        polymer_optomize[i] = i*numberspins_S;
    //        for(int j=0;j < numberspins_S;j++)
    //        {
    //            lattice_dist_optomize[i*numberspins_S +j]=int_point_distance(i,j);
    //        }

    //    }

    lattice_dist_optomize_2 = (int*)malloc(L1dim_S_xy * L1dim_S_xy * L1dim_S_z * sizeof(int));
    for(int j=0;j<numberspins_S;j++)
    {
        float x = j % L1dim_S_xy; // column
        float y = (j % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
        
        lattice_dist_optomize_2[j]=(x+y*L1dim_S_xy)+L1dim_S_xy*L1dim_S_xy*3;
    }
    
    lattice_dist_optomize_2_lattice = (int*)malloc((L1dim_S_xy*L1dim_S_xy)*(L1dim_S_xy * L1dim_S_xy*7) * sizeof(int));
    for(int i=0+(L1dim_S_xy*L1dim_S_xy*3);i<(L1dim_S_xy*L1dim_S_xy*4);i++)
    {
        for(int j=0;j<(L1dim_S_xy * L1dim_S_xy*7);j++)
        {
            lattice_dist_optomize_2_lattice[(i-(L1dim_S_xy*L1dim_S_xy*3))*(L1dim_S_xy * L1dim_S_xy*7) + j ]=int_point_distance(i,j);
        }
    }
    
    for (int i = 0; i < L1dim_S_xy * 12; i++)
    {
        index_optomize_x[i] = i % L1dim_S_xy;
        index_optomize_y[i] = i % L1dim_S_xy;

        calc_index_optomize_x[i] = (i % L1dim_S_xy);
        calc_index_optomize_y[i] = (i % L1dim_S_xy) *L1dim_S_xy;
    }

    for (int i = 0; i < L1dim_S_z * 6; i++)
    {
        index_optomize_z[i] = i % L1dim_S_z;
        calc_index_optomize_z[i] = (i % L1dim_S_z) *L1dim_S_xy * L1dim_S_xy;
    }

    xyz_index = (int*) malloc((L1dim_S_z + 10) *(L1dim_S_xy + 10) *(L1dim_S_xy + 10) *sizeof(int));
    int xyz_array[3] = { 0 };
    int x = -5;
    int y = -5;
    int z = -5;
    for (int i = 0; i < (L1dim_S_xy + 10); i++)
    {
        y = -5;
        for (int j = 0; j < L1dim_S_xy + 10; j++)
        {
            z = -5;
            for (int k = 0; k < (L1dim_S_z + 10); k++)
            {
                xyz_array[0] = i - 5;
                xyz_array[1] = j - 5;
                xyz_array[2] = k - 5;

                xyz_index[k *(L1dim_S_xy + 10) *(L1dim_S_xy + 10) + j *(L1dim_S_xy + 10) + i] = general_index(xyz_array);
                z++;
            }

            y++;
        }

        x++;
    }

    index_x = (int*) malloc(numberspins_S* sizeof(int));
    index_y = (int*) malloc(numberspins_S* sizeof(int));
    index_z = (int*) malloc(numberspins_S* sizeof(int));
    int ph_array[3] = { 0 };

    for (int i = 0; i < numberspins_S; i++)
    {
        general_coord(ph_array, i);
        index_x[i] = ph_array[0];
        index_y[i] = ph_array[1];
        index_z[i] = ph_array[2];
    }
    
    init_optomize_energy_index();
    optomize_energy_index();
}

int reduced_lattice_optimization(int index_1,int index_2)
{
    const int reduced_optomize_value =  (L1dim_S_xy*L1dim_S_xy*3);
    int shifted_index_1 = lattice_dist_optomize_2[index_1];
    int index_1_difference = index_1 - shifted_index_1;
    int shifted_index_2 = index_2 - index_1_difference;
    
    if(shifted_index_2>=0 && shifted_index_2<7*L1dim_S_xy*L1dim_S_xy)
    {
        return lattice_dist_optomize_2_lattice[(shifted_index_1-(L1dim_S_xy*L1dim_S_xy*3))*(L1dim_S_xy * L1dim_S_xy*7)+shifted_index_2];
    }
        
    return INT_MAX;
}

void init_poly_cores()
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to init poly_core()\n");
    fclose(stdoutlog);
    //holds the indexes of the chains
    poly_lattice_indexes = (int*) malloc((numbercores *lengthpercore *pertrusionpercore + backbone_length) *sizeof(int));
    wl_pseudo_chain = (int*) malloc((numbercores *lengthpercore *pertrusionpercore + backbone_length) *sizeof(int));

    poly_spin_reference = (int*) malloc((numbercores *lengthpercore *pertrusionpercore + backbone_length) *sizeof(int));

    //may not need forgot why this is here in the first place aa it will be tracked by the lattice polymer
    poly_lattice_connection = (int*) malloc((numbercores *lengthpercore *pertrusionpercore + backbone_length) *sizeof(int));

    //hold the physical coordinates of a given monomer
    poly_lattice_coordinates = (int*) malloc((numbercores *lengthpercore *pertrusionpercore *3 + 3 *backbone_length) *sizeof(int));
    wl_pseudo_chain_coordinates = (int*) malloc((numbercores *lengthpercore *pertrusionpercore *3 + 3 *backbone_length) *sizeof(int));

    ph_poly_lattice_indexes = (int*) malloc((numbercores *lengthpercore *pertrusionpercore + backbone_length) *sizeof(int));
    ph_poly_lattice_coordinates = (int*) malloc((numbercores *lengthpercore *pertrusionpercore *3 + 3 *backbone_length) *sizeof(int));

    int starter = L1dim_S_xy + 1;
    int array_inc = 0;
    int sign = -1;
    int column_inc = -1;
    int core_backbone_array[backbone_length]  ={0,1,1,1,1,1,1,1,1,1,1,1,0,1,1, 1,1,1,1,1,1,1,1,1,1,0,1,1,1,0, 1,0,1,1,1,1,1,1,1,1,1,1,1,1,0}; // needa to be ordeered least to greatest

//    for(int i=0;i<numbercores;i++)
//    {
//        core_backbone_array[i] = i;
//    }
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "init_poly_cores %i\n",poly_lattice_indexes[36]);
    fclose(stdoutlog);

    int displacement_index = (((L1dim_S_xy / 2) - 1) % (L1dim_S_xy / 2)) *2 *sign * L1dim_S_xy;
    for (int i = 0; i < backbone_length; i++)
    {
        if (i % (L1dim_S_xy / 2) == 0)
        {
            column_inc++;
            sign = sign *-1;

            displacement_index += (((L1dim_S_xy / 2) - 1) % (L1dim_S_xy / 2)) *2 *sign * L1dim_S_xy;
        }

        poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i] = starter + (i % (L1dim_S_xy / 2)) *2 *sign *L1dim_S_xy + 2 *column_inc - displacement_index;

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "%i %i\n", numbercores *lengthpercore *pertrusionpercore + i,poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i]);
        fclose(stdoutlog);

        if (core_backbone_array[i] == 1 && array_inc<numbercores)
        {
            for (int h = 0; h < lengthpercore; h++)
            {
                poly_lattice_indexes[array_inc *lengthpercore+h] = starter + (i % (L1dim_S_xy / 2)) *2 *sign *L1dim_S_xy + 2 *column_inc - displacement_index + L1dim_S_xy *L1dim_S_xy *h * 2;

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "\t%i %i\n",array_inc *lengthpercore+h, poly_lattice_indexes[array_inc *lengthpercore+h]);
                fclose(stdoutlog);
                

            }

            array_inc++;

        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "\n%i\n", poly_lattice_indexes[36]);
    fclose(stdoutlog);
    
    if (myid < min_max_id)
    {
        if (fopen("Minimum_Config_Poly.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;

            fptr = fopen("Minimum_Config_Poly.txt", "r");

            while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
            {
                poly_lattice_indexes[value] = valuee;
            }

            fclose(fptr);
        }
    }

    if (myid >= min_max_id)
    {
        if (fopen("Maximum_Config_Poly.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;

            fptr = fopen("Maximum_Config_Poly.txt", "r");

            while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
            {
                poly_lattice_indexes[value] = valuee;
            }

            fclose(fptr);
        }
    }

    //    if (fopen("Minimum_Config_Poly.txt","r") == NULL)
    //    {
    //        stdoutlog = fopen(stdoutlogname, "a");
    //        fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
    //        fclose(stdoutlog);
    //    }

    //    else    // rewrite of prior file reading system to be more portable
    //    {
    //        FILE* fptr;
    //        int value;
    //        int valuee;
    //
    //        fptr = fopen("Minimum_Config_Poly.txt", "r");
    //
    //        while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
    //        {
    //            poly_lattice_indexes[value] = valuee;
    //        }

    //
    //        fclose(fptr);
    //    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Assigned index init poly_core()%i\n",poly_lattice_indexes[36]);
    fclose(stdoutlog);

    // writes the chain occupation to lattice point
    int poly_index = -1;
    int spin_inc=0;
    for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
    {
        int value =1;
        
        if((i+1)%lengthpercore==0) value = 5+numbercores*pertrusionpercore;
        if((i)%lengthpercore==0){ value = 5+spin_inc;spin_inc++;}
        
        latticepoint_S[poly_lattice_indexes[i]] = value;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = value;
        }

        poly_spin_reference[i] = value;
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "psr:%i  %i\n",i,value);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    int spin_increment = 0;
    for (int i = numbercores *lengthpercore * pertrusionpercore; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int spin = 1;

        for (int j = 0; j < numbercores *lengthpercore * pertrusionpercore; j++)
        {
            if (poly_lattice_indexes[i] == poly_lattice_indexes[j])
            {
                spin = 5 + spin_increment;
                spin_increment++;
            }
        }

        latticepoint_S[poly_lattice_indexes[i]] = spin;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = spin;
        }

        poly_spin_reference[i] = spin;

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "psr2:%i  %i %i\n",i,spin,spin_increment);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "\n");
    fclose(stdoutlog);

    
    array_inc = 0;
    for (int i = 0; i < backbone_length; i++)
    {
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "%i  %i\n", poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i],latticepoint_S[poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i]]);
        fclose(stdoutlog);

        if (core_backbone_array[array_inc] == i)
        {
            for (int h = 0; h < lengthpercore; h++)
            {
                

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "\t%i   %i\n", poly_lattice_indexes[array_inc *lengthpercore+h],latticepoint_S[poly_lattice_indexes[array_inc *lengthpercore+h]]);
                fclose(stdoutlog);

            }

            array_inc++;

        }
    }
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice point init poly_core()%i\n",poly_lattice_indexes[36]);
    fclose(stdoutlog);

    // writes the lattice polymer data
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int chain_position = i % lengthpercore;

        lattice_polymer[3 *poly_lattice_indexes[i]] = poly_lattice_indexes[i];

        if (chain_position == 0)    // core definition needs to be removed or ignored if there are multiple petrusions from a single core
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = -100;
        }
        else if (chain_position == (lengthpercore - 1))    // head of the poly chain
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = -200;
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
        else
        {
            // aqny intermediate monomer
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice polymer init poly_core()\n");
    fclose(stdoutlog);

    // fills poly_lattice_coordinates
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        poly_coord_decompose(poly_lattice_indexes[i], i);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned coord init poly_core()\n");
    fclose(stdoutlog);

    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
        ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "created place holders init poly_core()\n");
    fclose(stdoutlog);

    init_en_and_dist_array();

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Called distance array init poly_core()%i\n",poly_lattice_indexes[0]);
    fclose(stdoutlog);

    solvent_loc = (int*) malloc(3*no_solvent_sites* sizeof(int));    // initializes main function to track solvent
    ph_solvent_loc = (int*) malloc(3*no_solvent_sites* sizeof(int));    // initializes main function to track solvent
    wl_solvent_loc = (int*) malloc(3*no_solvent_sites* sizeof(int));    // initializes wl function to track solvent vector for replacemt if movement fails
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Init solvent vectos\n");
    fclose(stdoutlog);
    init_solvent2(solvent_loc, ph_solvent_loc, wl_solvent_loc, no_solvent_sites, 2);

    attachment_ref = (int*) malloc(2 *((pertrusionpercore *lengthpercore *numbercores) + backbone_length) *sizeof(int));
    attachment_ion_ref = (int*) malloc(2 *no_ion_sites* sizeof(int));

    for (int i = 0; i < (pertrusionpercore *lengthpercore *numbercores) + backbone_length; i++)
    {
        attachment_ref[2 *i] = 0;
        attachment_ref[2 *i + 1] = 0;
    }

    for (int i = 0; i < no_ion_sites; i++)
    {
        attachment_ion_ref[2 *i] = 0;
        attachment_ion_ref[2 *i + 1] = 0;
    }

    ion_loc = (int*) malloc(no_ion_sites* sizeof(int));    // initializes main function to track solvent
    ph_ion_loc = (int*) malloc(no_ion_sites* sizeof(int));    // initializes main function to track solvent
    wl_ion_loc = (int*) malloc(no_ion_sites* sizeof(int));    // initializes wl function to track solvent vector for replacemt if movement fails
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "start ion\n");
    fclose(stdoutlog);
    init_solvent(ion_loc, ph_ion_loc, wl_ion_loc, no_ion_sites, 3);
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "exit ion\n");
    fclose(stdoutlog);
    
    if (myid < min_max_id)
    {
        if (fopen("Minimum_Config_Poly_Attach.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            int valueee;
            fptr = fopen("Minimum_Config_Poly_Attach.txt", "r");

            while (fscanf(fptr, "%i\t%i\t%i\n", &value, &valuee, &valueee) > 0)
            {
                attachment_ref[2 *value] = valuee;
                attachment_ref[2 *value + 1] = valueee;
            }

            fclose(fptr);
        }

        if (fopen("Minimum_Config_Ion_Attach.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            int valueee;
            fptr = fopen("Minimum_Config_Ion_Attach.txt", "r");

            while (fscanf(fptr, "%i\t%i\t%i\n", &value, &valuee, &valueee) > 0)
            {
                attachment_ion_ref[2 *value] = valuee;
                attachment_ion_ref[2 *value + 1] = valueee;
            }

            fclose(fptr);
        }

        if (fopen("Minimum_Config_Solvent.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            int valueee;
            int valueeee;
            fptr = fopen("Minimum_Config_Solvent.txt", "r");

            while (fscanf(fptr, "%i\t%i\t%i\t%i\n", &value, &valuee,&valueee,&valueeee) > 0)
            {
                solvent_loc[3*value] = valuee;
                solvent_loc[3*value+1] = valueee;
                solvent_loc[3*value+2] = valueeee;
            }

            fclose(fptr);
        }
        reset_indexes[0]=-1;
        reset_indexes[1]=3*no_solvent_sites;
        solvent_reset2(ph_solvent_loc, solvent_loc, 2);
        solvent_reset2(wl_solvent_loc, solvent_loc, 2);

        if (fopen("Minimum_Config_Ion.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            fptr = fopen("Minimum_Config_Ion.txt", "r");

            while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
            {
                ion_loc[value] = valuee;
            }

            fclose(fptr);
        }
        reset_indexes[0]=-1;
        reset_indexes[1]=no_ion_sites;
        solvent_reset(ph_ion_loc, ion_loc, 3);
        solvent_reset(wl_ion_loc, ion_loc, 3);
    }

    if (myid >= min_max_id)
    {
        if (fopen("Maximum_Config_Poly_Attach.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            int valueee;
            fptr = fopen("Maximum_Config_Poly_Attach.txt", "r");

            while (fscanf(fptr, "%i\t%i\t%i\n", &value, &valuee, &valueee) > 0)
            {
                attachment_ref[2 *value] = valuee;
                attachment_ref[2 *value + 1] = valueee;
            }

            fclose(fptr);
        }

        if (fopen("Maximum_Config_Ion_Attach.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            int valueee;
            fptr = fopen("Maximum_Config_Ion_Attach.txt", "r");

            while (fscanf(fptr, "%i\t%i\t%i\n", &value, &valuee, &valueee) > 0)
            {
                attachment_ion_ref[2 *value] = valuee;
                attachment_ion_ref[2 *value + 1] = valueee;
            }

            fclose(fptr);
        }

        if (fopen("Maximum_Config_Solvent.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
                        FILE * fptr;
                int value;
                int valuee;
                int valueee;
                int valueeee;
                fptr = fopen("Maximum_Config_Solvent.txt", "r");

                while (fscanf(fptr, "%i\t%i\t%i\t%i\n", &value, &valuee,&valueee,&valueeee) > 0)
                {
                    solvent_loc[3*value] = valuee;
                    solvent_loc[3*value+1] = valueee;
                    solvent_loc[3*value+2] = valueeee;
                }

                fclose(fptr);
            }

        reset_indexes[0]=-1;
        reset_indexes[1]=3*no_solvent_sites;
        
            solvent_reset2(ph_solvent_loc, solvent_loc, 2);
            solvent_reset2(wl_solvent_loc, solvent_loc, 2);
        
        if (fopen("Maximum_Config_Ion.txt", "r") == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
            fclose(stdoutlog);
        }
        else    // rewrite of prior file reading system to be more portable
        {
            FILE * fptr;
            int value;
            int valuee;
            fptr = fopen("Maximum_Config_Ion.txt", "r");

            while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
            {
                ion_loc[value] = valuee;
            }

            fclose(fptr);
        }

        reset_indexes[0]=-1;
        reset_indexes[1]=no_ion_sites;
        
        solvent_reset(ph_ion_loc, ion_loc, 3);
        solvent_reset(wl_ion_loc, ion_loc, 3);
    }

    optomize();

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Called init solvent()%i\n",poly_lattice_indexes[36]);
    fclose(stdoutlog);

    if (myid == 0)
    {
        printed = 1;

        sprintf(filename, "Minimum_Config_Poly2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
            }
        }

        fclose(file);
        sprintf(filename, "Minimum_Config_Poly_Attach2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ref[2 *i], attachment_ref[2 *i + 1]);
            }
        }

        fclose(file);
        sprintf(filename, "Minimum_Config_Solvent2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_solvent_sites; i++)
            {
                fprintf(file, "%i\t%i\t%i\t%i\n", i, solvent_loc[3*i],solvent_loc[3*i+1],solvent_loc[3*i+2]);
            }
        }

        fclose(file);

        sprintf(filename, "Minimum_Config_Ion2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\n", i, ion_loc[i]);
            }
        }

        fclose(file);
        
        sprintf(filename, "Minimum_Config_Ion_Attach2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
            }
        }

        fclose(file);
    }

    if (myid == numprocs - 1)
    {
        max_printed = 1;

        sprintf(filename, "Maximum_Config_Poly2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
            }
        }

        fclose(file);

        sprintf(filename, "Maximum_Config_Poly_Attach2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ref[2 *i], attachment_ref[2 *i + 1]);
            }
        }

        fclose(file);

        sprintf(filename, "Maximum_Config_Solvent2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_solvent_sites; i++)
            {
               fprintf(file, "%i\t%i\t%i\t%i\n", i, solvent_loc[3*i],solvent_loc[3*i+1],solvent_loc[3*i+2]);
            }
        }

        fclose(file);
        sprintf(filename, "Maximum_Config_Ion2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\n", i, ion_loc[i]);
            }
        }

        fclose(file);
        sprintf(filename, "Maximum_Config_Ion_Attach2.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
            }
        }

        fclose(file);

    }

    for (int i = 0; i < numberspins_S; i++)
    {
        latticepoint_S[i] = 0;
        lattice_polymer[i *3] = i;
        lattice_polymer[i *3 + 1] = -1;
        lattice_polymer[i *3 + 2] = -1;
    }

    for (int i = 0; i < numberspins_L; i++)
    {
        latticepoint_L[i] = 0;
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Assigned index init poly_core()\n");
    fclose(stdoutlog);

    // writes the chain occupation to lattice point
    poly_index = -1;
    spin_inc=0;
    for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
    {
        int value =1;
        
        if((i+1)%lengthpercore==0) value = 5+numbercores*pertrusionpercore;
        if((i)%lengthpercore==0){ value = 5+spin_inc;spin_inc++;}
        
        latticepoint_S[poly_lattice_indexes[i]] = value;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = value;
        }

        poly_spin_reference[i] = value;
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "2psr:%i  %i\n",i,value);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    spin_increment = 0;
    for (int i = numbercores *lengthpercore * pertrusionpercore; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int spin = 1;

        for (int j = 0; j < numbercores *lengthpercore * pertrusionpercore; j++)
        {
            if (poly_lattice_indexes[i] == poly_lattice_indexes[j])
            {
                spin = 5 + spin_increment;
                spin_increment++;
            }
        }

        latticepoint_S[poly_lattice_indexes[i]] = spin;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = spin;
        }

        poly_spin_reference[i] = spin;

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "2psr2:%i  %i %i\n",i,spin,spin_increment);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice point init poly_core()\n");
    fclose(stdoutlog);

    // writes the lattice polymer data
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int chain_position = i % lengthpercore;

        lattice_polymer[3 *poly_lattice_indexes[i]] = poly_lattice_indexes[i];

        if (chain_position == 0)    // core definition needs to be removed or ignored if there are multiple petrusions from a single core
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = -100;
        }
        else if (chain_position == (lengthpercore - 1))    // head of the poly chain
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = -200;
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
        else
        {
            // aqny intermediate monomer
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice polymer init poly_core()\n");
    fclose(stdoutlog);

    // fills poly_lattice_coordinates
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        poly_coord_decompose(poly_lattice_indexes[i], i);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned coord init poly_core()2\n");
    fclose(stdoutlog);

    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
        ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
        wl_pseudo_chain[i] = poly_lattice_indexes[i];
        wl_pseudo_chain_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        wl_pseudo_chain_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        wl_pseudo_chain_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];

    }

    for (int i = 0; i < no_solvent_sites; i++)
    {
        wl_solvent_loc[3*i] = solvent_loc[3*i];
        wl_solvent_loc[3*i+1] = solvent_loc[3*i+1];
        wl_solvent_loc[3*i+2] = solvent_loc[3*i+2];
        ph_solvent_loc[3*i] = solvent_loc[3*i];
        ph_solvent_loc[3*i+1] = solvent_loc[3*i+1];
        ph_solvent_loc[3*i+2] = solvent_loc[3*i+2];

        latticepoint_S[solvent_loc[3*i]] = 2;
        latticepoint_S[solvent_loc[3*i+1]] = -2;
        latticepoint_S[solvent_loc[3*i+2]] = -2;
        
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *solvent_loc[3*i] + iii]] = 2;
            latticepoint_L[neighbor_L[numberneighbors_L *solvent_loc[3*i+1] + iii]] = -2;
            latticepoint_L[neighbor_L[numberneighbors_L *solvent_loc[3*i+2] + iii]] = -2;
        }
    }

    for (int i = 0; i < no_ion_sites; i++)
    {
        ph_ion_loc[i] = ion_loc[i];
        wl_ion_loc[i] = ion_loc[i];

        latticepoint_S[ion_loc[i]] = 3;
        init_solvent_sub(ion_loc[i], 3);
    }

    //MPI_Abort(MPI_COMM_WORLD,1);
    tortuosity = (double*) malloc(hist_size* sizeof(double));
    tortuosity_buf = (double*) malloc(hist_size* sizeof(double));
    rog = (double*) malloc(hist_size* sizeof(double));
    rog_buf = (double*) malloc(hist_size* sizeof(double));

    tortuositybb = (double*) malloc(hist_size* sizeof(double));
    tortuositybb_buf = (double*) malloc(hist_size* sizeof(double));
    rogbb = (double*) malloc(hist_size* sizeof(double));
    rogbb_buf = (double*) malloc(hist_size* sizeof(double));

    visits = (long long *) malloc(hist_size* sizeof(long long));
    visits_buf = (long long *) malloc(hist_size* sizeof(long long));
    endtoend = (double*) malloc(hist_size* sizeof(double));
    endtoend_buf = (double*) malloc(hist_size* sizeof(double));

    radial_dist_EuPO = (float*) malloc(hist_size*14* sizeof(float));
    radial_dist_EuPO_buf =(float*) malloc(hist_size*14* sizeof(float));
	radial_dist_EuEu = (float*) malloc(hist_size*14* sizeof(float));
    radial_dist_EuEu_buf =(float*) malloc(hist_size*14* sizeof(float));
    radial_dist_EuO = (float*) malloc(hist_size*14*sizeof(float));
    radial_dist_EuO_buf = (float*) malloc(hist_size*14* sizeof(float));
    radial_dist_PMPO = (float*) malloc(hist_size*14* sizeof(float));
    radial_dist_PMPO_buf = (float*) malloc(hist_size*14* sizeof(float));
    radial_dist_POPO = (float*) malloc(hist_size*14* sizeof(float));
    radial_dist_POPO_buf = (float*) malloc(hist_size*14* sizeof(float));
    
    radial_visits_ion = (float*) malloc(hist_size* sizeof(float));
    radial_visits_poly = (float*) malloc(hist_size* sizeof(float));
    radial_visits_ion_buf = (float*) malloc(hist_size* sizeof(float));
    radial_visits_poly_buf = (float*) malloc(hist_size* sizeof(float));
    
    part_seq = (int*) malloc(hist_size* sizeof(int));
    func_per_seq = (float*) malloc(hist_size* sizeof(float));
    perf_seq = (int*) malloc(hist_size* sizeof(int));
    part_seq_buf = (int*) malloc(hist_size* sizeof(int));
       func_per_seq_buf = (float*) malloc(hist_size* sizeof(float));
       perf_seq_buf = (int*) malloc(hist_size* sizeof(int));
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned coord init poly_core()2\n");
    fclose(stdoutlog);
    
    for (int i = 0; i < hist_size; i++)
    {
        tortuosity[i] = 0.0;
        tortuosity_buf[i] = 0.0;
        rog[i] = 0.0;
        rog_buf[i] = 0.0;
        tortuositybb[i] = 0.0;
        tortuositybb_buf[i] = 0.0;
        rogbb[i] = 0.0;
        rogbb_buf[i] = 0.0;
        visits[i] = 0.0;
        visits_buf[i] = 0.0;
        endtoend[i] = 0.0;
        endtoend_buf[i] = 0.0;
        
            for(int k =0;k<14;k++)
            {
                radial_dist_EuEu[i*14+k] = 0;
                radial_dist_EuEu_buf[i*14+k] =0;
				radial_dist_EuPO[i*14+k] = 0;
                radial_dist_EuPO_buf[i*14+k] =0;
                radial_dist_EuO[i*14+k] = 0;
                radial_dist_EuO_buf[i*14+k] = 0;
            }
            radial_visits_ion[i]=0;
            radial_visits_ion_buf[i]=0;
        
            for(int k =0;k<14;k++)
            {
                radial_dist_PMPO[i*14+k] = 0;
                radial_dist_PMPO_buf[i*14+k] = 0;
                radial_dist_POPO[i*14+k] = 0;
                radial_dist_POPO_buf[i*14+k] = 0;
            }
            radial_visits_poly[i]=0;
            radial_visits_poly_buf[i]=0;
        
        part_seq[i]=0;
        func_per_seq[i]=0;
        perf_seq[i]=0;
        
        part_seq_buf[i]=0;
        func_per_seq_buf[i]=0;
        perf_seq_buf[i]=0;
    }
    
    array_inc = 0;
    for (int i = 0; i < backbone_length; i++)
    {
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "%i  %i\n", poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i],latticepoint_S[poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i]]);
        fclose(stdoutlog);

        if (core_backbone_array[array_inc] == i)
        {
            for (int h = 0; h < lengthpercore; h++)
            {
                

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "\t%i   %i\n", poly_lattice_indexes[array_inc *lengthpercore+h],latticepoint_S[poly_lattice_indexes[array_inc *lengthpercore+h]]);
                fclose(stdoutlog);

            }

            array_inc++;

        }
    }
    

}

void poly_coord_decompose(int reference_index, int polymer_position)
{
    int x = reference_index % L1dim_S_xy;    // column
    int y = (reference_index % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    int z = (reference_index / (L1dim_S_xy *L1dim_S_xy));    // depth

    poly_lattice_coordinates[3 *polymer_position] = x;
    poly_lattice_coordinates[3 *polymer_position + 1] = y;
    poly_lattice_coordinates[3 *polymer_position + 2] = z;
}

double point_distance(int reference_index_1, int reference_index_2)
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to point_distance() reference_index_1:%i reference_index_2:%i\n",reference_index_1,reference_index_2);
    //fclose(stdoutlog);

    float x_1 = reference_index_1 % L1dim_S_xy;    // column
    float y_1 = (reference_index_1 % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    float z_1 = (reference_index_1 / (L1dim_S_xy *L1dim_S_xy));    // depth

    float x_2 = reference_index_2 % L1dim_S_xy;    // column
    float y_2 = (reference_index_2 % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    float z_2 = (reference_index_2 / (L1dim_S_xy *L1dim_S_xy));    // depth

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i \n",(x_1 - x_2),(y_1 - y_2),(z_1 - z_2));

    int dx = abs((x_1 - x_2) - anint((x_1 - x_2), L1dim_S_xy) *L1dim_S_xy);
    int dy = abs((y_1 - y_2) - anint((y_1 - y_2), L1dim_S_xy) *L1dim_S_xy);
    int dz = abs((z_1 - z_2) - anint((z_1 - z_2), L1dim_S_z) *L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,sqrt(dx *dx + dy *dy + dz *dz));
    //fclose(stdoutlog);

    return (abs(sqrt(dx *dx + dy *dy + dz *dz)));
}

int int_point_distance(int reference_index_1, int reference_index_2)
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to point_distance() reference_index_1:%i reference_index_2:%i\n",reference_index_1,reference_index_2);
    //fclose(stdoutlog);

    int x_1 = reference_index_1 % L1dim_S_xy;    // column
    int y_1 = (reference_index_1 % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    int z_1 = (reference_index_1 / (L1dim_S_xy *L1dim_S_xy));    // depth

    int x_2 = reference_index_2 % L1dim_S_xy;    // column
    int y_2 = (reference_index_2 % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    int z_2 = (reference_index_2 / (L1dim_S_xy *L1dim_S_xy));    // depth

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i \n",(x_1 - x_2),(y_1 - y_2),(z_1 - z_2));

    int dx = abs((x_1 - x_2) - anint((x_1 - x_2), L1dim_S_xy) *L1dim_S_xy);
    int dy = abs((y_1 - y_2) - anint((y_1 - y_2), L1dim_S_xy) *L1dim_S_xy);
    int dz = abs((z_1 - z_2) - anint((z_1 - z_2), L1dim_S_z) *L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,sqrt(dx *dx + dy *dy + dz *dz));
    //fclose(stdoutlog);

    return (dx *dx + dy *dy + dz *dz);
}

double point_distance(int reference_index_1, int reference_index_2, int(&arr)[3])
{
    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "\t\t\tMade it to point_distance() reference_index_1:%i reference_index_2:%i\n",reference_index_1,reference_index_2);
    //    fclose(stdoutlog);

    int x_1 = reference_index_1 % L1dim_S_xy;    // column
    int y_1 = (reference_index_1 % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    int z_1 = (reference_index_1 / (L1dim_S_xy *L1dim_S_xy));    // depth

    int x_2 = reference_index_2 % L1dim_S_xy;    // column
    int y_2 = (reference_index_2 % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    int z_2 = (reference_index_2 / (L1dim_S_xy *L1dim_S_xy));    // depth

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "\t\t\tafter dz:%i \t dy: %i \t dz: %i\n",(x_1 - x_2),(y_1 - y_2),(z_1 - z_2));
    //    fclose(stdoutlog);

    int dx = abs((x_1 - x_2) - anint((x_1 - x_2), L1dim_S_xy) *L1dim_S_xy);
    int dy = abs((y_1 - y_2) - anint((y_1 - y_2), L1dim_S_xy) *L1dim_S_xy);
    int dz = abs((z_1 - z_2) - anint((z_1 - z_2), L1dim_S_z) *L1dim_S_z);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "\t\t\tafter dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,sqrt(dx *dx + dy *dy + dz *dz));
    //    fclose(stdoutlog);

    arr[0] = 0;
    arr[1] = 0;
    arr[2] = 0;

    if (dx < 4 && dy < 4 && dz < 4)
    {
        arr[0] = abs(dx);
        arr[1] = abs(dy);
        arr[2] = abs(dz);
    }

    return 0;
}

double poly_point_distance(int *poly_lattice_coordinates, int polymer_position_1, int polymer_position_2)
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to poly_point_distance()\n");
    //fclose(stdoutlog);

    int dx = poly_lattice_coordinates[3 *polymer_position_1] - poly_lattice_coordinates[3 *polymer_position_2];

    int dy = poly_lattice_coordinates[3 *polymer_position_1 + 1] - poly_lattice_coordinates[3 *polymer_position_2 + 1];

    int dz = poly_lattice_coordinates[3 *polymer_position_1 + 2] - poly_lattice_coordinates[3 *polymer_position_2 + 2];

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "befor dz:%i \t dy: %i \t dz: %i\n",dx,dy,dz);
    //fclose(stdoutlog);

    dx = abs(dx - anint(dx, L1dim_S_xy) *L1dim_S_xy);
    dy = abs(dy - anint(dy, L1dim_S_xy) *L1dim_S_xy);
    dz = abs(dz - anint(dz, L1dim_S_z) *L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,dist_array[dx][dy][dz]);
    //fclose(stdoutlog);

    if (dx < 4 && dy < 4 && dz < 4)
    {
        return dist_array[dx][dy][dz];
    }

    return dist_array[0][0][0];
}

double poly_point_distance(int *poly_lattice_coordinates, int polymer_position_1, int polymer_position_2, int(&arr)[3])
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to poly_point_distance()\n");
    //fclose(stdoutlog);

    int dx = poly_lattice_coordinates[3 *polymer_position_1] - poly_lattice_coordinates[3 *polymer_position_2];

    int dy = poly_lattice_coordinates[3 *polymer_position_1 + 1] - poly_lattice_coordinates[3 *polymer_position_2 + 1];

    int dz = poly_lattice_coordinates[3 *polymer_position_1 + 2] - poly_lattice_coordinates[3 *polymer_position_2 + 2];

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "before dz:%i \t dy: %i \t dz: %i \n",dx,dy,dz);
    //fclose(stdoutlog);

    dx = abs(dx - anint(dx, L1dim_S_xy) *L1dim_S_xy);
    dy = abs(dy - anint(dy, L1dim_S_xy) *L1dim_S_xy);
    dz = abs(dz - anint(dz, L1dim_S_z) *L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,dist_array[dx][dy][dz]);
    //fclose(stdoutlog);

    if (dx < 4 && dy < 4 && dz < 4)
    {
        arr[0] = dx;
        arr[1] = dy;
        arr[2] = dz;
        return dist_array[dx][dy][dz];
    }

    return dist_array[0][0][0];
}

void init_en_and_dist_array()
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to en array and distance()\n");
    fclose(stdoutlog);

    long_range_ion_en_attr = -20;
    long_range_ion_en_repul = 10;

    long_range_OH_en_attr = (-170.-45.)/2.;
    long_range_OEu_en_attr=(-1173.-313.)/2.;
    long_range_OO_en_repul = (603.+161.)/2.;
    long_range_HH_en_repul =(164.+44.)/2.;
    long_range_EuEu_en_repul = (7805.+2081.)/2.;
    long_range_HEu_en_repul = (1132.+302.)/2.;
    
    
    int PMPM[9]={-8, -6, -4, -2, -1, -1, -1, -1, 0};
    int PMPO[9]={-7, -5, -3, -2, -1, -1, -1, 0};
    int PMO[9]={-8, -6, -4, -2, -1, -1, -1, -1, 0};
    int PMH[9]={-6, -4, -3, -1, -1, -1, -1, 0, 0};
    int PMEu[9]={-3, -2, -2, -1, -1, 0, 0, 0, 0};
    
    int POPO[9]={-6, -4, -3, -1, -1, -1, -1, 0, 0};
    int POO[9]={-7, -5, -3, -2, -1, -1, -1, 0, 0};
    int POH[9]={-5, -4, -2, -1, -1, -1, 0, 0, 0};
    int POEu[9]={-3, -2, -1, -1, 0, 0, 0, 0, 0};
    
    int OO[9]={-8, -6, -4, -2, -1, -1, -1, -1, 0};
    int OH[9]={-6, -5, -3, -1, -1, -1, -1, 0, 0};
    int OEu[9]={-3, -2, -2, -1, -1, 0, 0, 0, 0};

    int HH[9]={-4, -3, -2, -1, -1, -1, 0, 0, 0};
    int HEu[9]={-2, -2, -1, -1, 0, 0, 0, 0, 0};
    
    int EuEu[9] = {-1, -1, -1, 0, 0, 0, 0, 0, 0};
    
    int AttOH[9]={-340, -304, -278, -240, -227, -215, -205, -196, -189};
    int AttOEu[9]={-2346, -2098, -1915, -1659, -1564, -1483, -1414, -1354, -1301};
    
    int RepOO[9]={1206, 1079, 985, 853, 804, 763, 727, 696, 669};
    int RepHH[9]={328, 294, 268, 232, 219, 208, 198, 189, 182};
    int RepEuEu[9]={15609, 13961, 12745, 11037, 10406, 9872, 9413, 9012, 8658};
    int RepHEu[9]={2263, 2024, 1848, 1600, 1509, 1431, 1365, 1307, 1255};
    
    int Bond[9] = {613, 1615, 2932, 6231, 8128, 10155,0,0,0};
    
    for (int i = 0; i < 14; i++)
    {
        int_en_array[i] = 0;
        int_en_array_bond[i] = 0;
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                en_array[i][j][k] = 0;
                en_array_bond[i][j][k] = 0;
                en_array_wall[i][j][k] = 0;

                if (i *i + j *j + k *k == 4)
                {
                    en_array[i][j][k] = (int) 303.974;
                    en_array_bond[i][j][k] = (int) 2028.01;
                    en_array_wall[i][j][k] = (int) - 95.7;

                    en_array[i][j][k] = (int) - 90;
                    en_array_coloumb_attr[i][j][k] = (int) - 140;
                    en_array_coloumb_repul[i][j][k] = (int) - 40;
                    en_array_bond[i][j][k] = (int) Bond[0];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[0];
                    en_array_PMPO[i][j][k] = (int) PMPO[0];
                    en_array_PMO[i][j][k] = (int) PMO[0];
                    en_array_PMH[i][j][k] = (int) PMH[0];
                    en_array_PMEu[i][j][k] = (int) PMEu[0];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[0];
                    en_array_POO[i][j][k] = (int) POO[0];
                    en_array_POH[i][j][k] = (int) POH[0];
                    en_array_POEu[i][j][k] = (int) POEu[0];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[0];
                    en_array_OH[i][j][k] = (int) OH[0];
                    en_array_OEu[i][j][k] = (int) OEu[0];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[0];
                    en_array_HEu[i][j][k] = (int) HEu[0];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[0];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[0];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[0];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[0];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[0];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[0];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[0];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }

                if (i *i + j *j + k *k == 5)
                {
                    en_array[i][j][k] = (int) 16.3753;
                    en_array_bond[i][j][k] = (int) 2355.75;
                    en_array_wall[i][j][k] = (int) - 75.9574;

                    en_array[i][j][k] = (int) - 60;
                    en_array_coloumb_attr[i][j][k] = (int) - 110;
                    en_array_coloumb_repul[i][j][k] = (int) - 40;
                    en_array_bond[i][j][k] = (int) Bond[1];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[1];
                    en_array_PMPO[i][j][k] = (int) PMPO[1];
                    en_array_PMO[i][j][k] = (int) PMO[1];
                    en_array_PMH[i][j][k] = (int) PMH[1];
                    en_array_PMEu[i][j][k] = (int) PMEu[1];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[1];
                    en_array_POO[i][j][k] = (int) POO[1];
                    en_array_POH[i][j][k] = (int) POH[1];
                    en_array_POEu[i][j][k] = (int) POEu[1];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[1];
                    en_array_OH[i][j][k] = (int) OH[1];
                    en_array_OEu[i][j][k] = (int) OEu[1];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[1];
                    en_array_HEu[i][j][k] = (int) HEu[1];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[1];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[1];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[1];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[1];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[1];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[1];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[1];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }

                if (i *i + j *j + k *k == 6)
                {
                    en_array[i][j][k] = (int) 3.50123;
                    en_array_bond[i][j][k] = (int) 3095.98;
                    en_array_wall[i][j][k] = (int) - 60.2979;

                    en_array[i][j][k] = (int) - 40;
                    en_array_coloumb_attr[i][j][k] = (int) - 80;
                    en_array_coloumb_repul[i][j][k] = (int) - 20;
                    en_array_bond[i][j][k] = (int) Bond[1];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[2];
                    en_array_PMPO[i][j][k] = (int) PMPO[2];
                    en_array_PMO[i][j][k] = (int) PMO[2];
                    en_array_PMH[i][j][k] = (int) PMH[2];
                    en_array_PMEu[i][j][k] = (int) PMEu[2];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[2];
                    en_array_POO[i][j][k] = (int) POO[2];
                    en_array_POH[i][j][k] = (int) POH[2];
                    en_array_POEu[i][j][k] = (int) POEu[2];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[2];
                    en_array_OH[i][j][k] = (int) OH[2];
                    en_array_OEu[i][j][k] = (int) OEu[2];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[2];
                    en_array_HEu[i][j][k] = (int) HEu[2];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[2];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[2];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[2];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[2];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[2];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[2];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[2];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }

                if (i *i + j *j + k *k == 8)
                {
                    en_array[i][j][k] = (int) 43.1734;
                    en_array_bond[i][j][k] = (int) 1175;
                    en_array_wall[i][j][k] = (int) - 40.462;

                    en_array[i][j][k] = (int) - 10;
                    en_array_coloumb_attr[i][j][k] = (int) - 50;
                    en_array_coloumb_repul[i][j][k] = (int) 10;
                    en_array_bond[i][j][k] = (int) Bond[3];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[3];
                    en_array_PMPO[i][j][k] = (int) PMPO[3];
                    en_array_PMO[i][j][k] = (int) PMO[3];
                    en_array_PMH[i][j][k] = (int) PMH[3];
                    en_array_PMEu[i][j][k] = (int) PMEu[3];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[3];
                    en_array_POO[i][j][k] = (int) POO[3];
                    en_array_POH[i][j][k] = (int) POH[3];
                    en_array_POEu[i][j][k] = (int) POEu[3];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[3];
                    en_array_OH[i][j][k] = (int) OH[3];
                    en_array_OEu[i][j][k] = (int) OEu[3];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[3];
                    en_array_HEu[i][j][k] = (int) HEu[3];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[3];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[3];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[3];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[3];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[3];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[3];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[3];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }

                if (i *i + j *j + k *k == 9)
                {
                    en_array[i][j][k] = (int) 57.6302;
                    en_array_bond[i][j][k] = (int) 7828.85;
                    en_array_wall[i][j][k] = (int) - 34.1454;

                    en_array[i][j][k] = (int) - 10;
                    en_array_coloumb_attr[i][j][k] = (int) - 40;
                    en_array_coloumb_repul[i][j][k] = (int) 30;
                    en_array_bond[i][j][k] = (int) Bond[4];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[4];
                    en_array_PMPO[i][j][k] = (int) PMPO[4];
                    en_array_PMO[i][j][k] = (int) PMO[4];
                    en_array_PMH[i][j][k] = (int) PMH[4];
                    en_array_PMEu[i][j][k] = (int) PMEu[4];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[4];
                    en_array_POO[i][j][k] = (int) POO[4];
                    en_array_POH[i][j][k] = (int) POH[4];
                    en_array_POEu[i][j][k] = (int) POEu[4];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[4];
                    en_array_OH[i][j][k] = (int) OH[4];
                    en_array_OEu[i][j][k] = (int) OEu[4];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[4];
                    en_array_HEu[i][j][k] = (int) HEu[4];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[4];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[4];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[4];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[4];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[4];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[4];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[4];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }
                if (i *i + j *j + k *k == 10)
                {
                    en_array[i][j][k] = (int) 57.6302;
                    en_array_bond[i][j][k] = (int) 7828.85;
                    en_array_wall[i][j][k] = (int) - 34.1454;

                    en_array[i][j][k] = (int) - 10;
                    en_array_coloumb_attr[i][j][k] = (int) - 40;
                    en_array_coloumb_repul[i][j][k] = (int) 30;
                    en_array_bond[i][j][k] = (int) Bond[5];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[5];
                    en_array_PMPO[i][j][k] = (int) PMPO[5];
                    en_array_PMO[i][j][k] = (int) PMO[5];
                    en_array_PMH[i][j][k] = (int) PMH[5];
                    en_array_PMEu[i][j][k] = (int) PMEu[5];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[5];
                    en_array_POO[i][j][k] = (int) POO[5];
                    en_array_POH[i][j][k] = (int) POH[5];
                    en_array_POEu[i][j][k] = (int) POEu[5];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[5];
                    en_array_OH[i][j][k] = (int) OH[5];
                    en_array_OEu[i][j][k] = (int) OEu[5];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[5];
                    en_array_HEu[i][j][k] = (int) HEu[5];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[5];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[5];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[5];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[5];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[5];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[5];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[5];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }
                if (i *i + j *j + k *k == 11)
                               {
                                   en_array[i][j][k] = (int) 57.6302;
                                   en_array_bond[i][j][k] = (int) 7828.85;
                                   en_array_wall[i][j][k] = (int) - 34.1454;

                                   en_array[i][j][k] = (int) - 10;
                                   en_array_coloumb_attr[i][j][k] = (int) - 40;
                                   en_array_coloumb_repul[i][j][k] = (int) 30;
                                   en_array_bond[i][j][k] = (int) Bond[6];
                                   en_array_wall[i][j][k] = (int) - 100;

                                   int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                                   int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                                   int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                                   int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                                   // en_array[i][j][k] =(int) 300;
                                   //en_array_bond[i][j][k] =(int) 2000;
                                   //en_array_wall[i][j][k] =(int) -100;
                                   en_array_PMPM[i][j][k] = (int) PMPM[6];
                                   en_array_PMPO[i][j][k] = (int) PMPO[6];
                                   en_array_PMO[i][j][k] = (int) PMO[6];
                                   en_array_PMH[i][j][k] = (int) PMH[6];
                                   en_array_PMEu[i][j][k] = (int) PMEu[6];
                                   int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                                   int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                                   int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                                   int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                                   int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                                   
                                   
                                   en_array_POPO[i][j][k] = (int) POPO[6];
                                   en_array_POO[i][j][k] = (int) POO[6];
                                   en_array_POH[i][j][k] = (int) POH[6];
                                   en_array_POEu[i][j][k] = (int) POEu[6];
                                   int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                                   int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                                   int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                                   int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                                   
                                   
                                   en_array_OO[i][j][k] = (int) OO[6];
                                   en_array_OH[i][j][k] = (int) OH[6];
                                   en_array_OEu[i][j][k] = (int) OEu[6];
                                   int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                                   int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                                   int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                                   
                                   en_array_HH[i][j][k] = (int) HH[6];
                                   en_array_HEu[i][j][k] = (int) HEu[6];
                                   int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                                   int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                                   
                                   en_array_EuEu[i][j][k] = (int) EuEu[6];
                                   int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                                   
                                   en_array_OH_Att[i][j][k]  = (int) AttOH[6];
                                   en_array_OEu_Att[i][j][k] = (int) AttOEu[6];
                                   int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                                   int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                                   
                                   en_array_OO_Rep[i][j][k]  = (int) RepOO[6];
                                   en_array_HH_Rep[i][j][k] = (int) RepHH[6];
                                   int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                                   int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                                   
                                   en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[6];
                                   en_array_HEu_Rep[i][j][k] = (int) RepHEu[6];
                                   int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                                   int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                               }
                if (i *i + j *j + k *k == 12)
                {
                    en_array[i][j][k] = (int) 57.6302;
                    en_array_bond[i][j][k] = (int) 7828.85;
                    en_array_wall[i][j][k] = (int) - 34.1454;

                    en_array[i][j][k] = (int) - 10;
                    en_array_coloumb_attr[i][j][k] = (int) - 40;
                    en_array_coloumb_repul[i][j][k] = (int) 30;
                    en_array_bond[i][j][k] = (int) Bond[7];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[7];
                    en_array_PMPO[i][j][k] = (int) PMPO[7];
                    en_array_PMO[i][j][k] = (int) PMO[7];
                    en_array_PMH[i][j][k] = (int) PMH[7];
                    en_array_PMEu[i][j][k] = (int) PMEu[7];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[7];
                    en_array_POO[i][j][k] = (int) POO[7];
                    en_array_POH[i][j][k] = (int) POH[7];
                    en_array_POEu[i][j][k] = (int) POEu[7];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[7];
                    en_array_OH[i][j][k] = (int) OH[7];
                    en_array_OEu[i][j][k] = (int) OEu[7];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[7];
                    en_array_HEu[i][j][k] = (int) HEu[7];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[7];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[7];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[7];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[7];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[7];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[7];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[7];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }
                if (i *i + j *j + k *k == 13)
                {
                    en_array[i][j][k] = (int) 57.6302;
                    en_array_bond[i][j][k] = (int) 7828.85;
                    en_array_wall[i][j][k] = (int) - 34.1454;

                    en_array[i][j][k] = (int) - 10;
                    en_array_coloumb_attr[i][j][k] = (int) - 40;
                    en_array_coloumb_repul[i][j][k] = (int) 30;
                    en_array_bond[i][j][k] = (int) Bond[8];
                    en_array_wall[i][j][k] = (int) - 100;

                    int_en_array[i *i + j *j + k *k] = en_array[i][j][k];
                    int_en_array_coloumb_attr[i *i + j *j + k *k] = en_array_coloumb_attr[i][j][k];
                    int_en_array_coloumb_repul[i *i + j *j + k *k] = en_array_coloumb_repul[i][j][k];
                    int_en_array_bond[i *i + j *j + k *k] = en_array_bond[i][j][k];

                    // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                    en_array_PMPM[i][j][k] = (int) PMPM[8];
                    en_array_PMPO[i][j][k] = (int) PMPO[8];
                    en_array_PMO[i][j][k] = (int) PMO[8];
                    en_array_PMH[i][j][k] = (int) PMH[8];
                    en_array_PMEu[i][j][k] = (int) PMEu[8];
                    int_en_array_PMPM[i *i + j *j + k *k]=en_array_PMPM[i][j][k];
                    int_en_array_PMPO[i *i + j *j + k *k]=en_array_PMPO[i][j][k];
                    int_en_array_PMO[i *i + j *j + k *k]=en_array_PMO[i][j][k];
                    int_en_array_PMH[i *i + j *j + k *k]=en_array_PMH[i][j][k];
                    int_en_array_PMEu[i *i + j *j + k *k]=en_array_PMEu[i][j][k];
                    
                    
                    en_array_POPO[i][j][k] = (int) POPO[8];
                    en_array_POO[i][j][k] = (int) POO[8];
                    en_array_POH[i][j][k] = (int) POH[8];
                    en_array_POEu[i][j][k] = (int) POEu[8];
                    int_en_array_POPO[i *i + j *j + k *k]=en_array_POPO[i][j][k];
                    int_en_array_POO[i *i + j *j + k *k]=en_array_POO[i][j][k];
                    int_en_array_POH[i *i + j *j + k *k]=en_array_POH[i][j][k];
                    int_en_array_POEu[i *i + j *j + k *k]=en_array_POEu[i][j][k];
                    
                    
                    en_array_OO[i][j][k] = (int) OO[8];
                    en_array_OH[i][j][k] = (int) OH[8];
                    en_array_OEu[i][j][k] = (int) OEu[8];
                    int_en_array_OO[i *i + j *j + k *k]=en_array_OO[i][j][k];
                    int_en_array_OH[i *i + j *j + k *k]=en_array_OH[i][j][k];
                    int_en_array_OEu[i *i + j *j + k *k]=en_array_OEu[i][j][k];
                    
                    en_array_HH[i][j][k] = (int) HH[8];
                    en_array_HEu[i][j][k] = (int) HEu[8];
                    int_en_array_HH[i *i + j *j + k *k]=en_array_HH[i][j][k];
                    int_en_array_HEu[i *i + j *j + k *k]=en_array_HEu[i][j][k];
                    
                    en_array_EuEu[i][j][k] = (int) EuEu[8];
                    int_en_array_EuEu[i *i + j *j + k *k]=en_array_EuEu[i][j][k];
                    
                    en_array_OH_Att[i][j][k]  = (int) AttOH[8];
                    en_array_OEu_Att[i][j][k] = (int) AttOEu[8];
                    int_en_array_OH_Att[i *i + j *j + k *k]=en_array_OH_Att[i][j][k];
                    int_en_array_OEu_Att[i *i + j *j + k *k]=en_array_OEu_Att[i][j][k];
                    
                    en_array_OO_Rep[i][j][k]  = (int) RepOO[8];
                    en_array_HH_Rep[i][j][k] = (int) RepHH[8];
                    int_en_array_OO_Rep[i *i + j *j + k *k]=en_array_OO_Rep[i][j][k];
                    int_en_array_HH_Rep[i *i + j *j + k *k]=en_array_HH_Rep[i][j][k];
                    
                    en_array_EuEu_Rep[i][j][k]  = (int) RepEuEu[8];
                    en_array_HEu_Rep[i][j][k] = (int) RepHEu[8];
                    int_en_array_EuEu_Rep[i *i + j *j + k *k]=en_array_EuEu_Rep[i][j][k];
                    int_en_array_HEu_Rep[i *i + j *j + k *k]=en_array_HEu_Rep[i][j][k];
                }
            }
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "en array created()\n");
    fclose(stdoutlog);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                dist_array[i][j][k] = INT_MAX;
                dist_array[0][0][0] = INT_MAX;
                if (i *i + j *j + k * k < 14 && i *i + j *j + k * k >= 4)
                {
                    dist_array[i][j][k] = sqrt(i *i + j *j + k *k);
                }
            }
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "dist array created()\n");
    fclose(stdoutlog);
}

inline int general_coord(int(&arr)[3], int reference_index)    // output the x y and z compenent for any index
{
    float x = reference_index % L1dim_S_xy;    // column
    float y = (reference_index % (L1dim_S_xy *L1dim_S_xy)) / L1dim_S_xy;    // row
    float z = (reference_index / (L1dim_S_xy *L1dim_S_xy));    // depth

    arr[0] = x;
    arr[1] = y;
    arr[2] = z;

    return 0;
}

int general_index(int(&arr)[3])    // output index from any x, y, z components
{
    while (arr[0] < 0)
    {
        arr[0] = arr[0] + L1dim_S_xy;
    }

    while (arr[0] >= L1dim_S_xy)
    {
        arr[0] = arr[0] - L1dim_S_xy;
    }

    while (arr[1] < 0)
    {
        arr[1] = arr[1] + L1dim_S_xy;
    }

    while (arr[1] >= L1dim_S_xy)
    {
        arr[1] = arr[1] - L1dim_S_xy;
    }

    while (arr[2] < 0)
    {
        arr[2] = arr[2] + L1dim_S_z;
    }

    while (arr[2] >= L1dim_S_z)
    {
        arr[2] = arr[2] - L1dim_S_z;
    }

    return (arr[0] + arr[1] *L1dim_S_xy + arr[2] *L1dim_S_xy *L1dim_S_xy);
}

inline int general_index_op(int(&arr)[3])    // output index from any x, y, z components
{
    arr[0] = index_optomize_x[arr[0] + neg_xy];
    arr[1] = index_optomize_y[arr[1] + neg_xy];
    arr[2] = index_optomize_z[arr[2] + neg_z];
    return (calc_index_optomize_x[arr[0]] + calc_index_optomize_y[arr[1]] + calc_index_optomize_z[arr[2]]);
}

inline int poly_latticepoint_clear(int *c_chain, int start, int end)    // set a range of monomers between two to zero for unimpeded movement
{
    int spin_holder = 0;

    for (int i = start + 1; i < end; i++)
    {
//        if(iqq>=29800 && myid == 2)
//        {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "Made it to poly_clear()  start: %i  end: %i\n",start,end);
//                fprintf(stdoutlog, "latticepoint_S[c_chain[i]]: %i\t",latticepoint_S[c_chain[i]]);
//                fprintf(stdoutlog, "neighbor_L[numberneighbors_L *c_chain[i] + iii]:");
//                fclose(stdoutlog);
//        }
        spin_holder = latticepoint_S[c_chain[i]];
        latticepoint_S[c_chain[i]] = 0;
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
//            if(iqq>=29800 && myid == 2)
//            {
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "%i\t",latticepoint_L[neighbor_L[numberneighbors_L *c_chain[i] + iii]]);
//                        fclose(stdoutlog);
//            }

            latticepoint_L[neighbor_L[numberneighbors_L *c_chain[i] + iii]] = 0;
        }
//        if(iqq>=29800 && myid == 2)
//        {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "\n");
//                fclose(stdoutlog);
//        }
    }

    return spin_holder;
}

int poly_latticepoint_clear(int index)    // set a range of monomers between two to zero for unimpeded movement
{
    latticepoint_S[index] = 0;

    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] = 0;
    }

    return 0;
}

inline int latticepoint_assign(int index, int spin)
{
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "checking assign index: %i    %i\n",index,latticepoint_S[index]);
//        fclose(stdoutlog);

    if (latticepoint_S[index] == 0)
    {
        latticepoint_S[index] = spin;
    }
    else
    {
        return 0;
    }

    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        //int x= latticepoint_L[neighbor_L[numberneighbors_L *index + iii]];
        
        if (latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] == 0)
        {
            

            latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] = spin;

//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "ni:%i\t%i\t%i\n",neighbor_L[numberneighbors_L *index + iii],x,latticepoint_L[neighbor_L[numberneighbors_L *index + iii]]);
//                        fclose(stdoutlog);
        }
        else    // resets the points switched to occupaton in the case where occupation is not allowed. Just performs cleanup itself so no other function has to
        {
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "checking assign index: %i  iii: %i\n",index,iii);
//                        fprintf(stdoutlog, "b latticepoint_S[index]: %i\n",latticepoint_S[index]);
//                        fclose(stdoutlog);

            latticepoint_S[index] = 0;

//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "a latticepoint_S[index]: %i\n\n",latticepoint_S[index]);
//                        fclose(stdoutlog);

            for (int j = 0; j < iii; j++)
            {
                //int x= latticepoint_L[neighbor_L[numberneighbors_L *index + j]];

                latticepoint_L[neighbor_L[numberneighbors_L *index + j]] = 0;

//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "%i\t%i\n",x,latticepoint_L[neighbor_L[numberneighbors_L *index + j]]);
//                                fclose(stdoutlog);
            }

//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "\n");
//                        fclose(stdoutlog);
            return 0;
        }
    }

    return 1;
}

inline int latticepoint_assign2(int index)
{
    if (latticepoint_S[index] != 0)
    {
        return 0;
    }

    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        if (latticepoint_L[neighbor_L[numberneighbors_L *index + iii]] != 0)
        {
            return 0;
        }
    }

    return 1;
}

int reflection_wrt_plane(double(&s_pos)[3], double(&e_pos)[3], int(&init_p)[3], double(&refl_p)[3], int axis)    // output index from any x, y, z components
{
    if (axis == 0)    // axis_end_pos[0]+1
    {
        if ((e_pos[1] - s_pos[1]) == 0 && (e_pos[2] - s_pos[2]) == 0)
        {
            refl_p[0] = init_p[0];
            refl_p[1] = init_p[1];
            refl_p[2] = init_p[2];
            return 0;
        }

        refl_p[0] = init_p[0];    //Good

        refl_p[1] = ((e_pos[1] *e_pos[1]) *init_p[1] - (e_pos[2] *e_pos[2]) *(init_p[1] - 2 *s_pos[1]) + init_p[1] *(s_pos[1] *s_pos[1]) + 2 *e_pos[2] *init_p[1] *s_pos[2] + 2 *init_p[2] *s_pos[1] *s_pos[2] - init_p[1] *(s_pos[2] *s_pos[2]) - 2 *e_pos[2] *s_pos[1] *(init_p[2] + s_pos[2]) + 2 *e_pos[1] *(e_pos[2] *init_p[2] - init_p[1] *s_pos[1] - (e_pos[2] + init_p[2]) *s_pos[2] + (s_pos[2] *s_pos[2]))) / ((e_pos[1] - s_pos[1]) *(e_pos[1] - s_pos[1]) + (e_pos[2] - s_pos[2]) *(e_pos[2] - s_pos[2]));

        refl_p[2] = init_p[2] + (2 *(e_pos[1] - s_pos[1]) *(e_pos[2] *(init_p[1] - s_pos[1]) + init_p[2] *s_pos[1] - init_p[1] *s_pos[2] + e_pos[1] *(-init_p[2] + s_pos[2]))) / (((e_pos[1] - s_pos[1]) *(e_pos[1] - s_pos[1])) + ((e_pos[2] - s_pos[2]) *(e_pos[2] - s_pos[2])));
    }

    if (axis == 1)
    {
        if ((e_pos[0] - s_pos[0]) == 0 && (e_pos[2] - s_pos[2]) == 0)
        {
            refl_p[0] = init_p[0];
            refl_p[1] = init_p[1];
            refl_p[2] = init_p[2];
            return 0;
        }

        refl_p[0] = ((e_pos[0] *e_pos[0]) *
            init_p[0] - (e_pos[2] *e_pos[2]) *(init_p[0] - 2 *s_pos[0]) +
            init_p[0] *(s_pos[0] *s_pos[0]) + 2 *e_pos[2] *init_p[0] *s_pos[2] +
            2 *init_p[2] *s_pos[0] *s_pos[2] - init_p[0] *(s_pos[2] *s_pos[2]) -
            2 *e_pos[2] *s_pos[0] *(init_p[2] + s_pos[2]) +
            2 *e_pos[0] *(e_pos[2] *init_p[2] -
                init_p[0] *s_pos[0] - (e_pos[2] + init_p[2]) *
                s_pos[2] + (s_pos[2] *s_pos[2])));

        refl_p[0] /= (((e_pos[0] - s_pos[0]) *(e_pos[0] - s_pos[0])) + ((e_pos[2] -
            s_pos[2]) *(e_pos[2] - s_pos[2])));

        refl_p[1] = init_p[1];

        refl_p[2] = init_p[2] + (2 *(e_pos[0] - s_pos[0]) *(e_pos[2] *(init_p[0] - s_pos[0]) + init_p[2] *s_pos[0] - init_p[0] *s_pos[2] + e_pos[0] *(-init_p[2] + s_pos[2]))) / ((e_pos[0] - s_pos[0]) *(e_pos[0] - s_pos[0]) + (e_pos[2] - s_pos[2]) *(e_pos[2] - s_pos[2]));

    }

    if (axis == 2)
    {
        if ((e_pos[0] - s_pos[0]) == 0 && (e_pos[1] - s_pos[1]) == 0)
        {
            refl_p[0] = init_p[0];
            refl_p[1] = init_p[1];
            refl_p[2] = init_p[2];
            return 0;
        }

        refl_p[0] = ((e_pos[0] *e_pos[0]) *
            init_p[0] - (e_pos[1] *e_pos[1]) *(init_p[0] - 2 *s_pos[0]) +
            init_p[0] *(s_pos[0] *s_pos[0]) + 2 *e_pos[1] *init_p[0] *s_pos[1] +
            2 *init_p[1] *s_pos[0] *s_pos[1] - init_p[0] *(s_pos[1] *s_pos[1]) -
            2 *e_pos[1] *s_pos[0] *(init_p[1] + s_pos[1]) +
            2 *e_pos[0] *(e_pos[1] *init_p[1] -
                init_p[0] *s_pos[0] - (e_pos[1] + init_p[1]) *
                s_pos[1] + (s_pos[1] *s_pos[1])));

        refl_p[0] /= (((e_pos[0] - s_pos[0]) *(e_pos[0] - s_pos[0])) + ((e_pos[1] -
            s_pos[1]) *(e_pos[1] - s_pos[1])));

        refl_p[1] = init_p[1] + (2 *(e_pos[0] - s_pos[0]) *(e_pos[1] *(init_p[0] - s_pos[0]) + init_p[1] *s_pos[0] - init_p[0] *s_pos[1] + e_pos[0] *(-init_p[1] + s_pos[1]))) / ((e_pos[0] - s_pos[0]) *(e_pos[0] - s_pos[0]) + (e_pos[1] - s_pos[1]) *(e_pos[1] - s_pos[1]));

        refl_p[2] = init_p[2];
    }

    //return (arr[0]+arr[1]*L1dim_S+arr[2]*L1dim_S*L1dim_S)
    return 0;
}

int rotation_wrt_vector(double(&s_pos)[3], double(&e_pos)[3], int(&init_p)[3], double(&rot_p)[3])    // output index from any x, y, z components
{
    double angle = (((rand() % 23 + 1) *15) *3.14159265) / 180.0;

    double ph =
        sqrt((s_pos[0] - e_pos[0]) *(s_pos[0] - e_pos[0]) + (s_pos[1] -
            e_pos[1]) *(s_pos[1] - e_pos[1]) + (s_pos[2] -
            e_pos[2]) *(s_pos[2] - e_pos[2]));

    double phh = (s_pos[0] - e_pos[0]) *(s_pos[0] -
        e_pos[0]) + (s_pos[1] - e_pos[1]) *(s_pos[1] -
        e_pos[1]) + (s_pos[2] - e_pos[2]) *(s_pos[2] - e_pos[2]);
    double sin_sub = sin(angle);
    double cos_sub = cos(angle);

    rot_p[0] = ((s_pos[0] *s_pos[0]) *init_p[0] + (e_pos[0] *e_pos[0]) *init_p[0] +
        cos_sub *
        init_p[0] *((s_pos[1] - e_pos[1]) *(s_pos[1] -
            e_pos[1]) + (s_pos[2] - e_pos[2]) *(s_pos[2] - e_pos[2])) -
        sin_sub *e_pos[1] *s_pos[2] *ph +
        sin_sub *init_p[1] *s_pos[2] *ph +
        sin_sub *s_pos[1] *ph *e_pos[2] -
        sin_sub *init_p[1] *ph *e_pos[2] +
        s_pos[0] *(-2 *e_pos[0] *
            init_p[0] + (-1 +
                cos_sub) *((s_pos[1] - e_pos[1]) *(e_pos[1] -
                init_p[1]) + (s_pos[2] - e_pos[2]) *(e_pos[2] -
                init_p[2]))) - (-1 + cos_sub) *
        e_pos[0] *((s_pos[1] - e_pos[1]) *(s_pos[1] -
            init_p[1]) + (s_pos[2] - e_pos[2]) *(s_pos[2] - init_p[2])) -
        sin_sub *s_pos[1] *ph *init_p[2] +
        sin_sub *e_pos[1] *ph *init_p[2]) / (phh);

    rot_p[1] = init_p[1] + (sin_sub *e_pos[0] *s_pos[2]) /
        ph - (sin_sub *init_p[0] *s_pos[2]) /
        ph - (sin_sub *s_pos[0] *e_pos[2]) /
        ph + (sin_sub *init_p[0] *e_pos[2]) / ph + (sin_sub *s_pos[0] *init_p[2]) / ph - ((sin_sub *e_pos[0] *init_p[2]) / ph) + ((-1 + cos_sub) *(e_pos[0] *init_p[0] *(s_pos[1] - e_pos[1]) +
            s_pos[0] *(init_p[0] *(-s_pos[1] + e_pos[1]) +
                e_pos[0] *(s_pos[1] + e_pos[1] - 2 *init_p[1])) +
            e_pos[0] *e_pos[0] *(-s_pos[1] + init_p[1]) +
            s_pos[0] *
            s_pos[0] *(-e_pos[1] + init_p[1]) - (s_pos[2] -
                e_pos[2]) *(init_p[1] *(-s_pos[2] + e_pos[2]) +
                e_pos[1] *(s_pos[2] - init_p[2]) +
                s_pos[1] *(-e_pos[2] + init_p[2])))) / phh;

    rot_p[2] = -(sin_sub *e_pos[0] *s_pos[1]) / ph + (sin_sub *init_p[0] *s_pos[1]) /
        ph + (sin_sub *s_pos[0] *e_pos[1]) /
        ph - (sin_sub *init_p[0] *e_pos[1]) /
        ph - (sin_sub *s_pos[0] *init_p[1]) /
        ph + (sin_sub *e_pos[0] *init_p[1]) +
        e_pos[2] + ((-1 + cos_sub) *(s_pos[2] -
                e_pos[2]) *((s_pos[0] - e_pos[0]) *(e_pos[0] -
                init_p[0]) + (s_pos[1] - e_pos[1]) *(e_pos[1] -
                init_p[1]) + (s_pos[2] - e_pos[2]) *(e_pos[2] - init_p[2])) /
            phh) + cos_sub *(-e_pos[2] + init_p[2]);

    return 0;
}

int mov_crankshaft()
{
    int arm = rand() % (numbercores * pertrusionpercore);// choose the arm the crankshaft will be performed on
    int poly_index_p1 = arm * lengthpercore + (rand() % (lengthpercore - 2));
    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
    general_coord(point1, poly_lattice_indexes[poly_index_p1]);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov crankshaft()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);

    uni_counter++;
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "k%i %i\n",uni_counter,latticepoint_S[80]);
//    fclose(stdoutlog);
    
    int poly_index_p2 = arm * lengthpercore + (rand() % lengthpercore);
    while (poly_index_p2 <= poly_index_p1 + 1)
    {
        poly_index_p2 = arm * lengthpercore + (rand() % lengthpercore);
    }
    int point2[3]; // choose the end point on the arm to perform the movement (must be farther along the arm and 2 indexes from point 1
    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "2nd point defined\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
//    fclose(stdoutlog);
    
    int gap = poly_index_p2 - poly_index_p1 - 1; // defines the length of the gap between the monomers
    
    reset_indexes[0] = poly_index_p1;
    reset_indexes[1] = poly_index_p2;
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc,0);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1, poly_index_p2); // clears all of the lattice points between the two axis points

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i  latticepoint_S[poly_lattice_indexes[poly_index_p1+1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p2]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1]],latticepoint_S[poly_lattice_indexes[poly_index_p1+1]],latticepoint_S[poly_lattice_indexes[poly_index_p2]]);
//    fclose(stdoutlog);
    
    double rotated_vector[3]; // holds the rotated point values
    int int_rotated_vector[3]; // holds the rounnded interger point values

    double starting_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    double end_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p2],(double)poly_lattice_coordinates[3 * poly_index_p2 + 1],(double)poly_lattice_coordinates[3 * poly_index_p2 + 2] };
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
//    fclose(stdoutlog);
    
    int assignment_possible = 0; // checks to ensure assignment is possible

    uni_counter++;
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "o%i %i\n",uni_counter,latticepoint_S[80]);
//    fclose(stdoutlog);
    
    for (int i = 0; i < gap; i++)
    {
        int initial_point[3]={ poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i)],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 1],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 2] };


//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
//        fclose(stdoutlog);
        
        rotation_wrt_vector(starting_point, end_point, initial_point, rotated_vector);

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "rotated={%f,%f,%f}\n",rotated_vector[0],rotated_vector[1],rotated_vector[2]);
//        fclose(stdoutlog);
//
        for (int iii = 0; iii < 3; iii++) // rounds the outputted values from the rotated vector
        {
            if (rotated_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(rotated_vector[iii] + .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] + .5);
            }
            if (rotated_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(rotated_vector[iii] - .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] - .5);
            }
        }
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "int_rotated_vector={%i,%i,%i}\n",int_rotated_vector[0],int_rotated_vector[1],int_rotated_vector[2]);
//        fclose(stdoutlog);
        
        int new_index = general_index(int_rotated_vector);
                
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//        fclose(stdoutlog);
        
        assignment_possible = latticepoint_assign(new_index,poly_spin_reference[poly_index_p1 + 1 + i]); // checks too see if assignment of the new lattice point wwas possible
        
        uni_counter++;
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "p%i %i\n",uni_counter,latticepoint_S[80]);
//        fclose(stdoutlog);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//        fclose(stdoutlog);
        
        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p1 + 1 + i] = new_index; // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+2] = index_z[new_index];
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1 + 1 + i], ph_poly_lattice_indexes[poly_index_p1 + i]) >= int_max_dist && assignment_possible == 1) //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit2\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2], ph_poly_lattice_indexes[poly_index_p2 - 1]) >= int_max_dist)  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "early exit3\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
//        uni_counter++;
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "r%i %i\n",uni_counter,latticepoint_S[80]);
//        fclose(stdoutlog);
//
        return 1;
    }

    return 1;
}

int mov_crankshaft_backbone()
{
    int chain_choice  =rand()%3;
    int arm = rand() % (backbone_length-2);// choose the arm the crankshaft will be performed on
    int place = (rand()%(lengthperchain-4));
    int poly_index_p1 = numbercores *pertrusionpercore *lengthpercore + chain_choice*lengthperchain + place;
    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
    general_coord(point1, poly_lattice_indexes[poly_index_p1]);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov crankshaft_backbone()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);

    int poly_index_p2 = poly_index_p1+ (rand()%(lengthperchain-1 - place));

    //while (poly_index_p2 <= poly_index_p1 + 1)
    //{
    //    poly_index_p2 = numbercores *pertrusionpercore *lengthpercore + rand() % (backbone_length);
    //}
    int point2[3]; // choose the end point on the arm to perform the movement (must be farther along the arm and 2 indexes from point 1
    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "2nd point defined\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
//    fclose(stdoutlog);

    int gap = poly_index_p2 - poly_index_p1 - 1; // defines the length of the gap between the monomers
    
    reset_indexes[0] = poly_index_p1;
    reset_indexes[1] = poly_index_p2;
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc,0);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1, poly_index_p2); // clears all of the lattice points between the two axis points

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i  latticepoint_S[poly_lattice_indexes[poly_index_p1+1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p2]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1]],latticepoint_S[poly_lattice_indexes[poly_index_p1+1]],latticepoint_S[poly_lattice_indexes[poly_index_p2]]);
//    fclose(stdoutlog);
    
    double rotated_vector[3]; // holds the rotated point values
    int int_rotated_vector[3]; // holds the rounnded interger point values

    double starting_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    double end_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p2],(double)poly_lattice_coordinates[3 * poly_index_p2 + 1],(double)poly_lattice_coordinates[3 * poly_index_p2 + 2] };
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
//    fclose(stdoutlog);
    
    int assignment_possible = 0; // checks to ensure assignment is possible

    
    for (int i = 0; i < gap; i++)
    {
        int spin_check = poly_spin_reference[poly_index_p1 + 1 + i];
        
        int initial_point[3]={ poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i)],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 1],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 2] };


//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
//        fclose(stdoutlog);
        
        rotation_wrt_vector(starting_point, end_point, initial_point, rotated_vector);

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "rotated={%f,%f,%f}\n",rotated_vector[0],rotated_vector[1],rotated_vector[2]);
//        fclose(stdoutlog);
//
        for (int iii = 0; iii < 3; iii++) // rounds the outputted values from the rotated vector
        {
            if (rotated_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(rotated_vector[iii] + .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] + .5);
            }
            if (rotated_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(rotated_vector[iii] - .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] - .5);
            }
        }
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "int_rotated_vector={%i,%i,%i}\n",int_rotated_vector[0],int_rotated_vector[1],int_rotated_vector[2]);
//        fclose(stdoutlog);
        
        int new_index = general_index(int_rotated_vector);
                
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//        fclose(stdoutlog);
        
        assignment_possible = latticepoint_assign(new_index, poly_spin_reference[poly_index_p1 + 1 + i]); // checks too see if assignment of the new lattice point wwas possible
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i  \n\n",assignment_possible);
//        fclose(stdoutlog);
        
        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p1 + 1 + i] = new_index; // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+2] = index_z[new_index];
            
            if (spin_check != 1)
            {
                ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore] = new_index;    // generates an index from the rotated vector
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3] = index_x[new_index];
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 1] = index_y[new_index];
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 2] = index_z[new_index];

                core_movement = 1;
            }
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1 + 1 + i], ph_poly_lattice_indexes[poly_index_p1 + i]) >= int_max_dist) //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit2\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

        if (spin_check != 1 && int_point_distance(ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore], ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore + 1]) >= int_max_dist)    //checks new distance between created monomer and previous monomer for every change
        {
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore+backbone_length; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit3\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore+backbone_length; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }
    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2], ph_poly_lattice_indexes[poly_index_p2 - 1]) >= int_max_dist)  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "early exit4\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 1;
    }

    return 1;
}

int mov_pivot()
{
    int chain_choice = rand()%3;
    int arm = rand() % (backbone_length-2);// choose the arm the crankshaft will be performed on
    int place = (rand()%(lengthperchain-4));
    int poly_index_p1 = numbercores *pertrusionpercore *lengthpercore + chain_choice*lengthperchain + place;
    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
    general_coord(point1, poly_lattice_indexes[poly_index_p1]);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov crankshaft_backbone()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);

    //int poly_index_p2 = poly_index_p1+ (rand()%(lengthperchain-1 - place));
    int poly_index_p2 = poly_index_p1 + 1;
    //    int point2[3];    // initializes end point that is one away from previous to create the pivot axis
    //    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "2nd point defined\n");
    //    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
    //    fclose(stdoutlog);

    int poly_index_arm_end = numbercores*pertrusionpercore*lengthpercore + chain_choice*lengthperchain + lengthperchain;

    int gap = poly_index_arm_end - poly_index_p2;

    //int gap = poly_index_p2 - poly_index_p1 - 1;    // defines the length of the gap between the monomers
    reset_indexes[0] = poly_index_p2;
    reset_indexes[1] = numbercores*pertrusionpercore*lengthpercore + chain_choice*lengthperchain + lengthperchain;
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc,0);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p2, arm *lengthpercore + lengthpercore);    // clears all of the lattice points between the two axis points

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
    //    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i\n\n",latticepoint_S[poly_lattice_indexes[arm *lengthpercore + lengthpercore-1]],latticepoint_S[poly_lattice_indexes[arm *lengthpercore + lengthpercore-1]]);
    //    fclose(stdoutlog);

    double rotated_vector[3];    // holds the rotated point values
    int int_rotated_vector[3];    // holds the rounnded interger point values

    double starting_point[3] = { (double) poly_lattice_coordinates[3 *poly_index_p1], (double) poly_lattice_coordinates[3 *poly_index_p1 + 1], (double) poly_lattice_coordinates[3 *poly_index_p1 + 2]
    };
    double end_point[3] = { (double) poly_lattice_coordinates[3 *poly_index_p2], (double) poly_lattice_coordinates[3 *poly_index_p2 + 1], (double) poly_lattice_coordinates[3 *poly_index_p2 + 2]
    };

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
    //    fclose(stdoutlog);

    int assignment_possible = 0;    // checks to ensure assignment is possible


    for (int i = 0; i < gap; i++)
    {
        int initial_point[3] = { poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i)], poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + 1], poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + 2]
        };

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
        //        fclose(stdoutlog);

        rotation_wrt_vector(starting_point, end_point, initial_point, rotated_vector);

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "rotated={%f,%f,%f}\n",rotated_vector[0],rotated_vector[1],rotated_vector[2]);
        //        fclose(stdoutlog);

        for (int iii = 0; iii < 3; iii++)    // rounds the outputted values from the rotated vector
        {
            if (rotated_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + iii] = (int)(rotated_vector[iii] + .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] + .5);
            }

            if (rotated_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + iii] = (int)(rotated_vector[iii] - .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] - .5);
            }
        }

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "int_rotated_vector={%i,%i,%i}\n",int_rotated_vector[0],int_rotated_vector[1],int_rotated_vector[2]);
        //        fclose(stdoutlog);

        int new_index = general_index(int_rotated_vector);

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
        //        fclose(stdoutlog);

        assignment_possible = latticepoint_assign(new_index, poly_spin_reference[poly_index_p2 + 1 + i]);    // checks too see if assignment of the new lattice point wwas possible

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
        //        fclose(stdoutlog);

        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p2 + 1 + i] = new_index;    // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i) *3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i) *3 + 1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i) *3 + 2] = index_z[new_index];
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

            return 0;
        }

        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2 + 1 + i], ph_poly_lattice_indexes[poly_index_p2 + i]) >= int_max_dist && assignment_possible == 1)    //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

            //            stdoutlog = fopen(stdoutlogname, "a");
            //            fprintf(stdoutlog, "early exit\n\n");
            //            fclose(stdoutlog);

            //            for (int j = 0; j < gap; j++)
            //            {
            //                stdoutlog = fopen(stdoutlogname, "a");
            //                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
            //                fclose(stdoutlog);
            //            }

            //            stdoutlog = fopen(stdoutlogname, "a");
            //            for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
            //            {
            //                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
            //            }

            //            fclose(stdoutlog);

            return 0;
        }
    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_arm_end - 2], ph_poly_lattice_indexes[poly_index_arm_end - 1]) >= int_max_dist)    //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "early exit 2\n\n");
        //        fclose(stdoutlog);

        //        for (int j = 0; j < gap; j++)
        //        {
        //            stdoutlog = fopen(stdoutlogname, "a");
        //            fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
        //            fclose(stdoutlog);
        //        }

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
        //        {
        //            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
        //        }

        //        fclose(stdoutlog);
        return 0;
    }
    if ( (poly_index_p1 + 1) % lengthpercore == 0 && attachment_ref[2 *poly_index_p1] == 0)
    {
        int success = attachment_poly(poly_index_p1, ph_poly_lattice_indexes[poly_index_p1]);
        if (success) reset_indexes[4] = 2;
    }

    return 1;

}

int mov_pivot_backbone()
{
    int chain_choice = rand()%numberofchains;
    int arm = 0;    // choose the arm the crankshaft will be performed on
    int poly_index_p1 = numbercores*pertrusionpercore*lengthpercore + arm *lengthpercore + chain_choice*lengthperchain +(rand() % (lengthperchain - 4));
    //    int point1[3];    // choose the initial point to perform the movement (must be 2 less than length of arm)
    //    general_coord(point1, poly_lattice_indexes[poly_index_p1]);

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Made it to mov pivot()\n");
//        fprintf(stdoutlog, "arm = %i poly_index_p1=%i  poly_lattice_indexes[poly_index_p1] = %i   %i     %i       %i\n\n",arm,poly_index_p1,poly_lattice_indexes[poly_index_p1],chain_choice,numbercores*pertrusionpercore*lengthpercore + arm *lengthpercore + chain_choice*lengthperchain,numbercores*pertrusionpercore*lengthpercore);
//        fclose(stdoutlog);

    int poly_index_p2 = poly_index_p1 + 1;
    //    int point2[3];    // initializes end point that is one away from previous to create the pivot axis
    //    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "2nd point defined\n");
    //    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
    //    fclose(stdoutlog);

    int poly_index_arm_end = numbercores*pertrusionpercore*lengthpercore + arm *lengthpercore + (chain_choice+1)*lengthperchain;

    int gap = (poly_index_arm_end - poly_index_p2)-1;

    //int gap = poly_index_p2 - poly_index_p1 - 1;    // defines the length of the gap between the monomers
    reset_indexes[0] = poly_index_p2;
    reset_indexes[1] = poly_index_arm_end;

//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc,0);

    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p2, arm *lengthpercore + lengthpercore);    // clears all of the lattice points between the two axis points

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//        fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i\n\n",latticepoint_S[poly_lattice_indexes[arm *lengthpercore + lengthpercore-1]],latticepoint_S[poly_lattice_indexes[arm *lengthpercore + lengthpercore-1]]);
//        fclose(stdoutlog);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "77 check  %i    %i\n\n",ph_poly_lattice_indexes[77],ph_poly_lattice_indexes[77]);
//    fclose(stdoutlog);
    double rotated_vector[3];    // holds the rotated point values
    int int_rotated_vector[3];    // holds the rounnded interger point values

    double starting_point[3] = { (double) poly_lattice_coordinates[3 *poly_index_p1], (double) poly_lattice_coordinates[3 *poly_index_p1 + 1], (double) poly_lattice_coordinates[3 *poly_index_p1 + 2]
    };
    double end_point[3] = { (double) poly_lattice_coordinates[3 *poly_index_p2], (double) poly_lattice_coordinates[3 *poly_index_p2 + 1], (double) poly_lattice_coordinates[3 *poly_index_p2 + 2]
    };

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
    //    fclose(stdoutlog);

    int assignment_possible = 0;    // checks to ensure assignment is possible

    for (int i = 0; i < gap; i++)
    {
        int spin_check = poly_spin_reference[poly_index_p2 + 1 + i];
        
        int initial_point[3] = { poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i)], poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + 1], poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + 2]
        };

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
//                fclose(stdoutlog);

        rotation_wrt_vector(starting_point, end_point, initial_point, rotated_vector);

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "rotated={%f,%f,%f}\n",rotated_vector[0],rotated_vector[1],rotated_vector[2]);
//                fclose(stdoutlog);

        for (int iii = 0; iii < 3; iii++)    // rounds the outputted values from the rotated vector
        {
            if (rotated_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + iii] = (int)(rotated_vector[iii] + .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] + .5);
            }

            if (rotated_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p2 + 1 + i) + iii] = (int)(rotated_vector[iii] - .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] - .5);
            }
        }

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "int_rotated_vector={%i,%i,%i}\n",int_rotated_vector[0],int_rotated_vector[1],int_rotated_vector[2]);
//                fclose(stdoutlog);

        int new_index = general_index(int_rotated_vector);

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//                fclose(stdoutlog);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "5Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
//        fclose(stdoutlog);
        assignment_possible = latticepoint_assign(new_index, poly_spin_reference[poly_index_p2 + 1 + i]);    // checks too see if assignment of the new lattice point wwas possible
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "6Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
//        fclose(stdoutlog);
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//                fclose(stdoutlog);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//        fclose(stdoutlog);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//        fclose(stdoutlog);

        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p2 + 1 + i] = new_index;    // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i) *3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i) *3 + 1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i) *3 + 2] = index_z[new_index];
            
            if (spin_check != 1)
            {
                ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore] = new_index;    // generates an index from the rotated vector
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3] = index_x[new_index];
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 1] = index_y[new_index];
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 2] = index_z[new_index];

                core_movement = 1;
            }
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

            return 0;
        }

        if ((((poly_index_p2 + 1 + i)-numbercores*pertrusionpercore*lengthpercore)%lengthperchain)!=0 &&  reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2 + 1 + i], ph_poly_lattice_indexes[poly_index_p2 + i]) >= int_max_dist && assignment_possible == 1)    //checks new distance between created monomer and previous monomer for every change
        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit1  %i    %i\n\n",ph_poly_lattice_indexes[77],ph_poly_lattice_indexes[77]);
//            fclose(stdoutlog);
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "early exit14  %i    %i\n\n",ph_poly_lattice_indexes[77],ph_poly_lattice_indexes[77]);
//                        fclose(stdoutlog);
//
//                        for (int j = 0; j < gap; j++)
//                        {
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "%i ph_poly:%i  poly:%i\n",j,ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                            fclose(stdoutlog);
//                        }
//
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
//                        {
//                            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
//                        }
//
//                        fclose(stdoutlog);

            return 0;
        }
        if (spin_check != 1 && int_point_distance(ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore], ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore + 1]) >= int_max_dist)    //checks new distance between created monomer and previous monomer for every change
        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit3\n\n");
//            fclose(stdoutlog);

            
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "j: %i  %i  ph_poly:%i  poly:%i\n",j,(poly_index_p2 + 1 + j),ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore+backbone_length; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit31\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore+backbone_length; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }
    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_arm_end - 2], ph_poly_lattice_indexes[poly_index_arm_end - 1]) >= int_max_dist)    //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "early exit 2\n\n");
//                fclose(stdoutlog);
//
//                for (int j = 0; j < gap; j++)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                    fclose(stdoutlog);
//                }
//
//                stdoutlog = fopen(stdoutlogname, "a");
//                for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
//                {
//                    fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
//                }
//
//                fclose(stdoutlog);
        return 0;
    }
    else    //checks new distance between last created monomer and considered end monomer for last change
    {
        // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "pivot success\n\n");
//                fclose(stdoutlog);
//
//                for (int j = 0; j < gap; j++)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                    fclose(stdoutlog);
//                }
//
//                stdoutlog = fopen(stdoutlogname, "a");
//                for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
//                {
//                    fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
//                }
//
//                fclose(stdoutlog);

        return 1;
    }

    return 1;

}

int mov_inversion()
{
    return 0;
}

int mov_planar_reflection()
{
    int arm = rand() % (numbercores *pertrusionpercore);    // choose the arm the crankshaft will be performed on
    int poly_index_p1 = arm *lengthpercore + (rand() % (lengthpercore - 2));
    //    int point1[3];    // choose the initial point to perform the movement (must be 2 less than length of arm)
    //    general_coord(point1, poly_lattice_indexes[poly_index_p1]);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Made it to mov planar reflection()\n");
    //    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
    //    fclose(stdoutlog);

    int poly_index_p2 = arm *lengthpercore + (rand() % lengthpercore);
    while (poly_index_p2 <= poly_index_p1 + 1)
    {
        poly_index_p2 = arm *lengthpercore + (rand() % lengthpercore);
    }

    //    int point2[3];    // choose the end point on the arm to perform the movement (must be farther along the arm and 2 indexes from point 1
    //    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "2nd point defined\n");
    //    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
    //    fclose(stdoutlog);

    int gap = poly_index_p2 - poly_index_p1 - 1;    // defines the length of the gap between the monomers
    reset_indexes[0] = poly_index_p1;
    reset_indexes[1] = poly_index_p2;
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc,0);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1, poly_index_p2);    // clears all of the lattice points between the two axis points

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
    //    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i poly_lattice_indexes[poly_index_p1+1] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1]],latticepoint_S[poly_lattice_indexes[poly_index_p1+1]] ,latticepoint_S[poly_lattice_indexes[poly_index_p2]]);
    //    fclose(stdoutlog);

    int axis = rand() % 3;

    double reflected_vector[3];    // holds the rotated point values
    int int_reflected_vector[3];    // holds the rounnded interger point values

    double starting_point[3] = { (double) ph_poly_lattice_coordinates[3 *poly_index_p1], (double) ph_poly_lattice_coordinates[3 *poly_index_p1 + 1], (double) ph_poly_lattice_coordinates[3 *poly_index_p1 + 2]
    };
    double end_point[3] = { (double) ph_poly_lattice_coordinates[3 *poly_index_p2], (double) ph_poly_lattice_coordinates[3 *poly_index_p2 + 1], (double) ph_poly_lattice_coordinates[3 *poly_index_p2 + 2]
    };

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "sp={%f,%f,%f} axis = %i ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2], axis,end_point[0],end_point[1],end_point[2]);
    //    fclose(stdoutlog);

    int assignment_possible = 0;    // checks to ensure assignment is possible

    for (int i = 0; i < gap; i++)
    {
        int initial_point[3] = { poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i)], poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + 1], poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + 2]
        };

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
        //        fclose(stdoutlog);

        reflection_wrt_plane(starting_point, end_point, initial_point, reflected_vector, axis);

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "reflected={%f,%f,%f}\n",reflected_vector[0],reflected_vector[1],reflected_vector[2]);
        //        fclose(stdoutlog);

        for (int iii = 0; iii < 3; iii++)    // rounds the outputted values from the rotated vector
        {
            if (reflected_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + iii] = (int)(reflected_vector[iii] + .5);
                int_reflected_vector[iii] = (int)(reflected_vector[iii] + .5);
            }

            if (reflected_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + iii] = (int)(reflected_vector[iii] - .5);
                int_reflected_vector[iii] = (int)(reflected_vector[iii] - .5);
            }
        }

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "int_reflected_vector={%i,%i,%i}\n",int_reflected_vector[0],int_reflected_vector[1],int_reflected_vector[2]);
        //
        //        fclose(stdoutlog);

        int new_index = general_index(int_reflected_vector);

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
        //        fclose(stdoutlog);

        assignment_possible = latticepoint_assign(new_index, poly_spin_reference[poly_index_p1 + 1 + i]);    // checks too see if assignment of the new lattice point wwas possible

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
        //        fclose(stdoutlog);

        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p1 + 1 + i] = new_index;    // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i) *3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i) *3 + 1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i) *3 + 2] = index_z[new_index];
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            return 0;
        }

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog,"ph_poly_lattice_indexes[poly_index_p1 + 1 + i]:%i   ph_poly_lattice_indexes[poly_index_p1 + i]:%i   max_dist: %f\n\n",poly_index_p1 + 1 + i,poly_index_p1 + i,max_dist);
        //        fclose(stdoutlog);

        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1 + 1 + i], ph_poly_lattice_indexes[poly_index_p1 + i]) >= int_max_dist && assignment_possible == 1)    //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            return 0;
        }
    }

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "max_dist: %f\n\n",max_dist);
    //    fclose(stdoutlog);

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2], ph_poly_lattice_indexes[poly_index_p2 - 1]) >= int_max_dist)    //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }
    else    //checks new distance between last created monomer and considered end monomer for last change
    {
        // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

        return 1;
    }

    return 1;
}

int mov_planar_reflection_backbone()
{
    int chain_choice  =rand()%3;
    int arm = rand() % (backbone_length-2);// choose the arm the crankshaft will be performed on
    int place = (rand()%(lengthperchain-4));
    int poly_index_p1 = numbercores *pertrusionpercore *lengthpercore + chain_choice*lengthperchain + place;
    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
    general_coord(point1, poly_lattice_indexes[poly_index_p1]);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov crankshaft_backbone()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);

    int poly_index_p2 = poly_index_p1+ (rand()%(lengthperchain-1 - place));

    //while (poly_index_p2 <= poly_index_p1 + 1)
    //{
    //    poly_index_p2 = numbercores *pertrusionpercore *lengthpercore + rand() % (backbone_length);
    //}
    int point2[3]; // choose the end point on the arm to perform the movement (must be farther along the arm and 2 indexes from point 1
    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "2nd point defined\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
//    fclose(stdoutlog);

    int gap = poly_index_p2 - poly_index_p1 - 1; // defines the length of the gap between the monomers
    
    reset_indexes[0] = poly_index_p1;
    reset_indexes[1] = poly_index_p2;
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc,0);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1, poly_index_p2); // clears all of the lattice points between the two axis points

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i  latticepoint_S[poly_lattice_indexes[poly_index_p1+1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p2]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1]],latticepoint_S[poly_lattice_indexes[poly_index_p1+1]],latticepoint_S[poly_lattice_indexes[poly_index_p2]]);
//    fclose(stdoutlog);
    
    double rotated_vector[3]; // holds the rotated point values
    int int_rotated_vector[3]; // holds the rounnded interger point values

    double starting_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    double end_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p2],(double)poly_lattice_coordinates[3 * poly_index_p2 + 1],(double)poly_lattice_coordinates[3 * poly_index_p2 + 2] };
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
//    fclose(stdoutlog);
    
    int assignment_possible = 0; // checks to ensure assignment is possible

    uni_counter++;
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "o%i %i\n",uni_counter,latticepoint_S[80]);
//    fclose(stdoutlog);
    
    int axis = rand() % 3;

    double reflected_vector[3];    // holds the rotated point values
    int int_reflected_vector[3];    // holds the rounnded interger point values
    
    for (int i = 0; i < gap; i++)
    {
        int spin_check = poly_spin_reference[poly_index_p1 + 1 + i];
        
        int initial_point[3] = { poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i)], poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + 1], poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + 2]
        };

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
        //        fclose(stdoutlog);

        reflection_wrt_plane(starting_point, end_point, initial_point, reflected_vector, axis);

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "reflected={%f,%f,%f}\n",reflected_vector[0],reflected_vector[1],reflected_vector[2]);
        //        fclose(stdoutlog);

        for (int iii = 0; iii < 3; iii++)    // rounds the outputted values from the rotated vector
        {
            if (reflected_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + iii] = (int)(reflected_vector[iii] + .5);
                int_reflected_vector[iii] = (int)(reflected_vector[iii] + .5);
            }

            if (reflected_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 *(poly_index_p1 + 1 + i) + iii] = (int)(reflected_vector[iii] - .5);
                int_reflected_vector[iii] = (int)(reflected_vector[iii] - .5);
            }
        }

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "int_reflected_vector={%i,%i,%i}\n",int_reflected_vector[0],int_reflected_vector[1],int_reflected_vector[2]);
        //
        //        fclose(stdoutlog);

        int new_index = general_index(int_rotated_vector);
                
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//        fclose(stdoutlog);
        
        assignment_possible = latticepoint_assign(new_index, poly_spin_reference[poly_index_p1 + 1 + i]); // checks too see if assignment of the new lattice point wwas possible
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i  \n\n",assignment_possible);
//        fclose(stdoutlog);
        
        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p1 + 1 + i] = new_index; // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+2] = index_z[new_index];
            
            if (spin_check != 1)
            {
                ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore] = new_index;    // generates an index from the rotated vector
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3] = index_x[new_index];
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 1] = index_y[new_index];
                ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 2] = index_z[new_index];

                core_movement = 1;
            }
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1 + 1 + i], ph_poly_lattice_indexes[poly_index_p1 + i]) >= int_max_dist) //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit2\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

        if (spin_check != 1 && int_point_distance(ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore], ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore + 1]) >= int_max_dist)    //checks new distance between created monomer and previous monomer for every change
        {
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore+backbone_length; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit3\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore+backbone_length; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }
    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2], ph_poly_lattice_indexes[poly_index_p2 - 1]) >= int_max_dist)  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "early exit4\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 1;
    }

    return 1;
}

inline int mov_translation_backbone()
{
    int chain_choice  =rand()%3;
    int arm = rand() % (lengthperchain);    // choose the arm the translation mov will be performed on
    int poly_index_p1 = numbercores *pertrusionpercore *lengthpercore + chain_choice*lengthperchain + arm;
    //int point1[3];    // choose the initial point to perform the movement

    //general_coord(point1, poly_lattice_indexes[poly_index_p1]);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Made it to mov translation()\n");
    //    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
    //    fclose(stdoutlog);

    int spin_check = poly_spin_reference[poly_index_p1];
    
    reset_indexes[0] = poly_index_p1 - 1;
    reset_indexes[1] = poly_index_p1 + 1;

    //reset_arm = -1;

    if(spin_check != 1)
    {
        core_movement = 1;
    }
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
    
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1 - 1, poly_index_p1 + 1);    // clears all of the lattice points between the two axis points

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called poly_clear()\n");
    //    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1-1]],latticepoint_S[poly_lattice_indexes[poly_index_p1]]);
    //    fclose(stdoutlog);

    int int_translated_vector[3] = { poly_lattice_coordinates[3 *poly_index_p1], poly_lattice_coordinates[3 *poly_index_p1 + 1], poly_lattice_coordinates[3 *poly_index_p1 + 2]
    };
    //double end_point[3] = { poly_lattice_coordinates[3 *poly_index_p2],poly_lattice_coordinates[3 *poly_index_p2 + 1],poly_lattice_coordinates[3 *poly_index_p2 + 2] };

    //int int_translated_vector[3] = { starting_point[0],starting_point[1],starting_point[2] };    // holds the translated interger point values

    int assignment_possible = 0;    // checks to ensure assignment is possible

    int random_comp = (rand() % 3);
    int value;
    if ((rand() % 2) == 0)
    {
        value = -1;
    }
    else
    {
        value = 1;
    }

    int multi = (rand()%2)+1;
    int_translated_vector[random_comp] += value*multi;

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "sp={%i,%i,%i} comp = %i tv={%i,%i,%i}\n",starting_point[0],starting_point[1],starting_point[2], random_comp,int_translated_vector[0],int_translated_vector[1],int_translated_vector[2]);
    //    fclose(stdoutlog);

    
    int new_index = general_index_op(int_translated_vector);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
    //    fclose(stdoutlog);
    
    assignment_possible = latticepoint_assign2(new_index);    // checks too see if assignment of the new lattice point wwas possible

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
    //    fclose(stdoutlog);
    
    if (assignment_possible == 1)
    {
        
        
        ph_poly_lattice_indexes[poly_index_p1] = new_index;    // generates an index from the rotated vector
        ph_poly_lattice_coordinates[poly_index_p1 * 3] = int_translated_vector[0];
        ph_poly_lattice_coordinates[poly_index_p1 * 3 + 1] = int_translated_vector[1];
        ph_poly_lattice_coordinates[poly_index_p1 * 3 + 2] = int_translated_vector[2];

        if (spin_check != 1)
        {
            ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore] = new_index;    // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3] = int_translated_vector[0];
            ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 1] = int_translated_vector[1];
            ph_poly_lattice_coordinates[(spin_check - 5) *lengthpercore *3 + 2] = int_translated_vector[2];

            core_movement = 1;
        }
    }
    else
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

        return 0;
    }
    
    //assignment_possible = latticepoint_assign(new_index,1);    // checks too see if assignment of the new lattice point wwas possible

    if (poly_index_p1 != chain_choice*lengthperchain+lengthperchain - 1 && int_point_distance(ph_poly_lattice_indexes[poly_index_p1], ph_poly_lattice_indexes[poly_index_p1 + 1]) >= int_max_dist)    //checks new distance between created monomer and previous monomer for every change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

        return 0;
    }

    if (spin_check != 1 && int_point_distance(ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore], ph_poly_lattice_indexes[(spin_check - 5) *lengthpercore *pertrusionpercore + 1]) >= int_max_dist)    //checks new distance between created monomer and previous monomer for every change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }

    if (poly_index_p1 != chain_choice*lengthperchain && int_point_distance(ph_poly_lattice_indexes[poly_index_p1], ph_poly_lattice_indexes[poly_index_p1 - 1]) >= int_max_dist)    //checks new distance between last created monomer and considered end monomer for last change
    {

        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

        return 0;
    }
    else    //checks new distance between last created monomer and considered end monomer for last change
    {
        // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        loc_en_before = local_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc,0);
        return 1;
    }

    return 1;
}

inline int mov_translation()
{
    int arm = rand() % (numbercores *pertrusionpercore);    // choose the arm the translation mov will be performed on
    int poly_index_p1 = arm *lengthpercore + (rand() % (lengthpercore - 1)) + 1;    // cant eqaul 0 otherwise the core would move
    //int point1[3];    // choose the initial point to perform the movement

    //general_coord(point1, poly_lattice_indexes[pventoly_index_p1]);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Made it to mov translation()\n");
    //    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
    //    fclose(stdoutlog);

    reset_indexes[0] = poly_index_p1 - 1;
    reset_indexes[1] = poly_index_p1 + 1;
    
//    no_neg_electrostatic_inter[0]=0;
//    no_pos_electrostatic_inter[0]=0;
   
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1 - 1, poly_index_p1 + 1);    // clears all of the lattice points between the two axis points

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called poly_clear()\n");
    //    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1-1]],latticepoint_S[poly_lattice_indexes[poly_index_p1]]);
    //    fclose(stdoutlog);

    int int_translated_vector[3] = { poly_lattice_coordinates[3 *poly_index_p1], poly_lattice_coordinates[3 *poly_index_p1 + 1], poly_lattice_coordinates[3 *poly_index_p1 + 2]
    };
    //double end_point[3] = { poly_lattice_coordinates[3 *poly_index_p2],poly_lattice_coordinates[3 *poly_index_p2 + 1],poly_lattice_coordinates[3 *poly_index_p2 + 2] };

    //int int_translated_vector[3] = { starting_point[0],starting_point[1],starting_point[2] };    // holds the translated interger point values

    int assignment_possible = 0;    // checks to ensure assignment is possible

    int random_comp = (rand() % 3);
    int value;
    if ((rand() % 2) == 0)
    {
        value = -1;
    }
    else
    {
        value = 1;
    }

    int multi = (rand()%2)+1;
    int_translated_vector[random_comp] += value*multi;

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "sp={%i,%i,%i} comp = %i tv={%i,%i,%i}\n",starting_point[0],starting_point[1],starting_point[2], random_comp,int_translated_vector[0],int_translated_vector[1],int_translated_vector[2]);
    //    fclose(stdoutlog);

//    reset_indexes[0] = poly_index_p1 - 1;
//    reset_indexes[1] = poly_index_p1 + 1;


    int new_index = general_index_op(int_translated_vector);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
    //    fclose(stdoutlog);

    assignment_possible = latticepoint_assign2(new_index);    // checks too see if assignment of the new lattice point wwas possible

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
    //    fclose(stdoutlog);

    if (assignment_possible == 1)
    {
        
        
        ph_poly_lattice_indexes[poly_index_p1] = new_index;    // generates an index from the rotated vector
        ph_poly_lattice_coordinates[poly_index_p1 * 3] = int_translated_vector[0];
        ph_poly_lattice_coordinates[poly_index_p1 * 3 + 1] = int_translated_vector[1];
        ph_poly_lattice_coordinates[poly_index_p1 * 3 + 2] = int_translated_vector[2];
    }
    else
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }

    //assignment_possible = latticepoint_assign(new_index,1);    // checks too see if assignment of the new lattice point wwas possible
//
//    if (attachment_ref[2 *poly_index_p1] == 1 && int_point_distance(ph_poly_lattice_indexes[poly_index_p1], attachment_ref[2 *poly_index_p1]) > int_max_attach_dist)    //checks new distance between created monomer and previous monomer for every change
//    {
//        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//        return 0;
//    }

    if ((poly_index_p1 + 1) % lengthpercore != 0 && int_point_distance(ph_poly_lattice_indexes[poly_index_p1], ph_poly_lattice_indexes[poly_index_p1 + 1]) >= int_max_dist)    //checks new distance between created monomer and previous monomer for every change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }

    if (int_point_distance(ph_poly_lattice_indexes[poly_index_p1], ph_poly_lattice_indexes[poly_index_p1 - 1]) >= int_max_dist)    //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }
    else
    {
         loc_en_before = local_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc,0);
        return 1;
    }
//    if (attachment_ref[2 *poly_index_p1] == 0 && (poly_index_p1 + 1) % lengthpercore == 0)
//    {
//        int success = attachment_poly(poly_index_p1, ph_poly_lattice_indexes[poly_index_p1]);
//        if (success) reset_indexes[4] = 2;
//    }

    return 1;
}

inline int poly_mov()
{
    //confused2 =0;
    
    reset_arm = -1;
    int mov_choice;
    //std::random_device rd;     // Will be used to obtain a seed for the random number engine

    core_movement = 0;

    mov_choice = (rand() % (9999)) + 1;
    //mov_choice = 1;
    //int mov_completion = 0;

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Made it to poly move() Move choice:%i  %i\n",mov_choice,latticepoint_L[0]);
//        fclose(stdoutlog);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
//    fclose(stdoutlog);

    if (mov_choice < 1500)
    {
        poly_solvent_mov = 0;    //monomer mov
        poly_ion_mov = 0;
        chosen_mov = 0;    //translation
        MoveProposal = mov_translation();
        return MoveProposal;
    }

    if (mov_choice >= 1500 && mov_choice < 3000)
    {
        chosen_mov = 6;    //solvent
        poly_ion_mov = 1;    //solvent mov
        poly_solvent_mov = 0;    //monomer mov
        MoveProposal = solvent_translation(ion_loc, ph_ion_loc, no_ion_sites, 3);
        return MoveProposal;
    }
    if (mov_choice >= 3000 && mov_choice < 8500)
    {
        chosen_mov = 4;    //solvent
        poly_solvent_mov = 1;    //solvent mov

        poly_ion_mov = 0;
        MoveProposal = solvent_translation2(solvent_loc, ph_solvent_loc, no_solvent_sites, 2);
        
        return MoveProposal;
    }
    if (mov_choice >= 8500 && mov_choice < 9999)
    {
        chosen_mov = 5;    //solvent
        poly_solvent_mov = 0;    //monomer mov
        poly_ion_mov = 0;
        chosen_mov = 5;    //translation
        loc_en_before_2=0;
        
        MoveProposal = mov_translation_backbone();
        
        return MoveProposal;
    }
    else
    {
        int sec_choice = (rand() % (5)) + 5;
        //int sec_choice = 4;

        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "Had to make a second choice poly move() sec choice:%i\n",sec_choice);
        //        fclose(stdoutlog);

        if (sec_choice == 1)
        {
            chosen_mov = 1;    //crankshaft
            poly_solvent_mov = 0;    //monomer mov
            poly_ion_mov = 0;
            //MoveProposal = mov_pivot();
            return MoveProposal;
        }

        if (sec_choice == 2)
        {
            chosen_mov = 2;    //pivot
            poly_solvent_mov = 0;    //monomer mov
            poly_ion_mov = 0;
            //MoveProposal = mov_planar_reflection();
            return MoveProposal;
        }

        if (sec_choice == 3)
        {
            chosen_mov = 3;    //reflection
            poly_solvent_mov = 0;    //monomer mov
            poly_ion_mov = 0;
            //MoveProposal = mov_crankshaft();
            return MoveProposal;
        }

        if (sec_choice == 4)
        {
            chosen_mov = 4;    //solvent
            poly_solvent_mov = 1;    //solvent mov

            poly_ion_mov = 0;
            MoveProposal = solvent_translation2(solvent_loc, ph_solvent_loc, no_solvent_sites, 2);
            
            return MoveProposal;
        }

        if (sec_choice == 5)
        {
            chosen_mov = 5;    //solvent
            poly_solvent_mov = 0;    //monomer mov
            poly_ion_mov = 0;
            chosen_mov = 5;    //translation
            loc_en_before_2=0;
            
            MoveProposal = mov_translation_backbone();
            
            return MoveProposal;
        }

        if (sec_choice == 6)
        {
            chosen_mov = 6;    //solvent
            poly_ion_mov = 1;    //solvent mov
            poly_solvent_mov = 0;    //monomer mov
            MoveProposal = solvent_translation(ion_loc, ph_ion_loc, no_ion_sites, 3);
            return MoveProposal;
        }
        if (sec_choice == 7)
        {
            chosen_mov = 7;    //solvent
            poly_ion_mov = 0;    //solvent mov
            poly_solvent_mov = 0;    //monomer mov
            MoveProposal = mov_crankshaft_backbone();
            return MoveProposal;
        }
        if (sec_choice == 8)
        {
            chosen_mov = 8;    //pivot
            poly_solvent_mov = 0;    //monomer mov
            poly_ion_mov = 0;
            MoveProposal = mov_planar_reflection_backbone();
            return MoveProposal;
        }
        if (sec_choice == 9)
        {
            chosen_mov = 9;    //crankshaft
            poly_solvent_mov = 0;    //monomer mov
            poly_ion_mov = 0;
            MoveProposal = mov_pivot_backbone();
            return MoveProposal;
        }
    }

    return 0;
}

void lattice_polymer_reset(int *p_chain, int *c_chain, int start, int end)    //retrofit of reorganize function from previous code saved some time rewriting parts of the function
{
    for (int i = start + 1; i < end; i++)    // clears the lattice polymer in the range set
    {
        lattice_polymer[3 *p_chain[i]] = p_chain[i];
        lattice_polymer[3 *p_chain[i] + 1] = -1;    // next lattice point on chain (if not applicable -1 assigned)
        lattice_polymer[3 *p_chain[i] + 2] = -1;
    }

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "qq\n");
    //fclose(stdoutlog);
    for (int i = start + 1; i < end; i++)
    {
        lattice_polymer[3 *p_chain[i]] = p_chain[i];

        if (i % lengthpercore == 0)    //core definition
        {
            lattice_polymer[3 *c_chain[i] + 1] = c_chain[i + 1];    // next lattice point on chain (if not applicable -1 assigned)
            lattice_polymer[3 *c_chain[i] + 2] = -100;    //previous lattice point on chain (if applicable -1 assigned)
        }

        if (i % lengthpercore != 0 && i % (lengthpercore) != lengthpercore - 1)    // intermediate monomer
        {
            lattice_polymer[3 *c_chain[i] + 1] = c_chain[i + 1];    // next lattice point on chain (if not applicable -1 assigned)
            lattice_polymer[3 *c_chain[i] + 2] = c_chain[i - 1];    //previous lattice point on chain (if applicable -1 assigned)
        }

        if (i % (lengthpercore) == lengthpercore - 1)    // length of an arm
        {
            lattice_polymer[3 *c_chain[i] + 1] = -200;    // next lattice point on chain (if not applicable -1 assigned)
            lattice_polymer[3 *c_chain[i] + 2] = c_chain[i - 1];    //previous lattice point on chain (if applicable -1 assigned)
        }
    }
}

int lattice_poly_index_reset(int *p_chain, int *c_chain, int start, int end)    // reorganizes the latticepoints more efficient than previous scrubbing over entire lattice
{
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//    fclose(stdoutlog);
    for (int i = start + 1; i < end; i++)    // clears the lattice polymer in the range set
    {
        latticepoint_S[p_chain[i]] = 0;
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        if(poly_spin_reference[i]!=1 && poly_spin_reference[i]!=5+numbercores*pertrusionpercore)
        {
            latticepoint_S[p_chain[(poly_spin_reference[i]-5)*lengthpercore]] =0;
        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *p_chain[i] + iii]] = 0;
            if(poly_spin_reference[i]!=1 && poly_spin_reference[i]!=5+numbercores*pertrusionpercore)
            {
                
                latticepoint_L[neighbor_L[numberneighbors_L *p_chain[(poly_spin_reference[i]-5)*lengthpercore] + iii]] =0;
            }
        }
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "qjj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
    }

    for (int i = start + 1; i < end; i++)    // clears the lattice polymer in the range set
    {
//         stdoutlog = fopen(stdoutlogname, "a");
//         fprintf(stdoutlog, "353f%i %i\n",uni_counter,latticepoint_L[0]);
//         fclose(stdoutlog);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        
        latticepoint_S[c_chain[i]] = poly_spin_reference[i];
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "353f%i %i\n",uni_counter,latticepoint_L[0]);
//        fclose(stdoutlog);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        if(poly_spin_reference[i]!=1 && poly_spin_reference[i]!=5+numbercores*pertrusionpercore)
        {
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "jhj%i %i  %i\n",uni_counter,c_chain[(poly_spin_reference[i]-5)*lengthpercore],c_chain[i]);
//            fclose(stdoutlog);
            
            latticepoint_S[c_chain[(poly_spin_reference[i]-5)*lengthpercore]] = poly_spin_reference[i];
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "3534f%i %i\n",uni_counter,latticepoint_L[0]);
//            fclose(stdoutlog);
        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *c_chain[i] + iii]] = poly_spin_reference[i];
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "35351f%i %i\n",uni_counter,latticepoint_L[0]);
//            fclose(stdoutlog);
            if(poly_spin_reference[i]!=1 && poly_spin_reference[i]!=5+numbercores*pertrusionpercore)
            {
                latticepoint_L[neighbor_L[numberneighbors_L *c_chain[(poly_spin_reference[i]-5)*lengthpercore] + iii]] =poly_spin_reference[i];
                
            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "3535f%i %i\n",uni_counter,latticepoint_L[0]);
//            fclose(stdoutlog);
        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        p_chain[i] = c_chain[i];
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
        if(poly_spin_reference[i]!=1 && poly_spin_reference[i]!=5+numbercores*pertrusionpercore)
        {
            p_chain[(poly_spin_reference[i]-5)*lengthpercore] = c_chain[(poly_spin_reference[i]-5)*lengthpercore];
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "3536f%i %i\n",uni_counter,latticepoint_L[0]);
//            fclose(stdoutlog);
        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
//        fclose(stdoutlog);
    }

    return 0;
}

int poly_coord_reset(int *p_chain, int *c_chain, int start, int end)
{
    for (int i = start + 1; i < end; i++)    // clears the lattice polymer in the range set
    {
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "253f%i %i\n",uni_counter,latticepoint_L[0]);
//        fclose(stdoutlog);
        p_chain[i *3 + 0] = c_chain[i *3 + 0];
        p_chain[i *3 + 1] = c_chain[i *3 + 1];
        p_chain[i *3 + 2] = c_chain[i *3 + 2];
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "253f%i %i\n",uni_counter,latticepoint_L[0]);
//        fclose(stdoutlog);
        if(poly_spin_reference[i]!=1 && poly_spin_reference[i]!=5+numbercores*pertrusionpercore)
        {
            p_chain[((poly_spin_reference[i]-5)*lengthpercore) *3 + 0] = c_chain[((poly_spin_reference[i]-5)*lengthpercore) *3 + 0];
            p_chain[((poly_spin_reference[i]-5)*lengthpercore) *3 + 1] = c_chain[((poly_spin_reference[i]-5)*lengthpercore)*3 + 1];
            p_chain[((poly_spin_reference[i]-5)*lengthpercore) *3 + 2] = c_chain[((poly_spin_reference[i]-5)*lengthpercore)*3 + 2];
        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "253f%i %i\n",uni_counter,latticepoint_L[0]);
//        fclose(stdoutlog);
    }

    return 0;
}



int total_energy(int *p_l_c, int *p_l_i, int *s_l, int *i_l)    // returns total energy of system
{
    
    int e = 0;
    int array[3];
    //
    for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
    {
        for (int j = 0; j < no_ion_sites; j++)    // calculatw=es energy for solvent interactions
        {
            point_distance(p_l_i[i], i_l[j], array);
            if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
            {
                if ((i) % lengthpercore == 0)
                {
                   // e -= 2 *en_array_PMEu[array[0]][array[1]][array[2]];
                    e += 2 *en_array_HEu_Rep[array[0]][array[1]][array[2]];
                    
                    
                }
            }

            if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
            {
                if (i % lengthpercore == 0)
                {
                    e += 2 * long_range_HEu_en_repul;
                }
            }
        }

        if (i % lengthpercore != 0)    // excludes the monomer in the backbone
        {
            if (poly_point_distance(p_l_c, i - 1, i, array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
            {
                e += 2 *en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
            }
            for (int j = 0; j < 3*no_solvent_sites; j++)    // calculatw=es energy for solvent interactions
            {
                
                point_distance(p_l_i[i], s_l[j], array);
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    if(j%3==0)
                    {
                    e += 2 *en_array_POO[array[0]][array[1]][array[2]];
                    }
                    if(j%3!=0)
                    {
                    e += 2 *en_array_POH[array[0]][array[1]][array[2]];
                    }
                }
                if(poly_spin_reference[i]==5+numbercores*pertrusionpercore)
                {
                    if(j%3==0)
                    {
                        if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                        {
                            //e -= 2 *en_array_POO[array[0]][array[1]][array[2]];
                            e += 2 *en_array_OO_Rep[array[0]][array[1]][array[2]];
                            
//                            if(j==14 || j==15)
//                            {
//                                confused2+=2 *en_array_coloumb_repul[array[0]][array[1]][array[2]];
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "te PO 2*rep    %i\n",confused2);
//                                fclose(stdoutlog);
//                            }
                        }
                    
                        if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                        {
                            e += 2 *long_range_OO_en_repul;
                            
                            
                        }
                    }
                    if(j%3!=0)
                    {
                        if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                        {
                            //e -= 2 *en_array_POH[array[0]][array[1]][array[2]];
                            e += 2 *en_array_OH_Att[array[0]][array[1]][array[2]];
                            
//                            if(j==14 || j==15)
//                            {
//                                confused2+=2 *en_array_coloumb_attr[array[0]][array[1]][array[2]];
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "te PO 2*attr  %i\n",confused2);
//                            fclose(stdoutlog);
//                            }
                        }
                        if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                        {
                            e += 2 *long_range_OH_en_attr;
                            
                            
                        }
                    }
                }
            }

            for (int j = 0; j < no_ion_sites; j++)    // calculatw=es energy for solvent interactions
            {
                point_distance(p_l_i[i], i_l[j], array);
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    e += 2 *en_array_POEu[array[0]][array[1]][array[2]];

                    if ((i + 1) % lengthpercore == 0)
                    {
                        //e -= 2 *en_array_POEu[array[0]][array[1]][array[2]];
                        e += 2 *en_array_OEu_Att[array[0]][array[1]][array[2]];
                        
//                        if(j==14)
//                        {
//                            confused2 += 2 *en_array_OEu_Att[array[0]][array[1]][array[2]];
//                            confused2 += 2 *en_array_POEu[array[0]][array[1]][array[2]];
//                            if(iqq == 22981)
//                            {
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "\tPO and PO OEu  %i\n",confused2);
//                                fclose(stdoutlog);
//                            }
//                        }
                    }
                }

                if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                {
                    if ((i + 1) % lengthpercore == 0)
                    {
                        e += 2 * long_range_OEu_en_attr;
                    }
                }
            }
            for (int j = 0; j < numbercores *lengthpercore *pertrusionpercore + backbone_length; j++)
            {
                if (j % lengthpercore != 0 || j >= numbercores *lengthpercore *pertrusionpercore)    // excludes the monomer in the backbone
                {
                    point_distance(p_l_i[i], p_l_i[j], array);

                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        if(i< numbercores *lengthpercore *pertrusionpercore && j >= numbercores *lengthpercore *pertrusionpercore && (i+1)%lengthpercore==0)
                        {
                            e += en_array_PMPO[array[0]][array[1]][array[2]];
                            
                            if(j==48||j==2)
                            {
                                confused2 += en_array_PMPO[array[0]][array[1]][array[2]];
                                if(iqq == 22981)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "\t PMPO %i\n",confused2);
                                    fclose(stdoutlog);
                                }
                            }
                        }
                        else
                        {
                            e += en_array_POPO[array[0]][array[1]][array[2]];
                            
                            if(j==48||j==2)
                            {
                                confused2 += en_array_POPO[array[0]][array[1]][array[2]];
                                if(iqq == 22981)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "\t PMPO %i\n",confused2);
                                    fclose(stdoutlog);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = numbercores *lengthpercore *pertrusionpercore + 1; i < numbercores *lengthpercore *pertrusionpercore + backbone_length-1; i++)
    {
            if(((i-(numbercores *lengthpercore *pertrusionpercore)))%lengthperchain==0)
            {
                poly_point_distance(p_l_c, i-1, i, array);
                
                e -= 2 *en_array_bond[array[0]][array[1]][array[2]];
                
            }
    }
    
    for (int i = numbercores *lengthpercore *pertrusionpercore; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        if(i!=numbercores *lengthpercore *pertrusionpercore)
        {
            if (poly_point_distance(p_l_c, i - 1, i, array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
            {
                e += 2 *en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
            }
        }
        for (int j = 0; j < 3*no_solvent_sites; j++)    // calculatw=es energy for solvent interactions
        {
            point_distance(p_l_i[i], s_l[j], array);
            if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
            {
                if(j%3==0)
                {
                    e += 2 *en_array_PMO[array[0]][array[1]][array[2]];
                    
                   
                }
                else{
                    if(j%3!=0)
                    {
                        e += 2 *en_array_PMH[array[0]][array[1]][array[2]];
                        
                        
                    }
                    
                }
            }
            if(poly_spin_reference[i]>=5)
            {
                if(j%3==0)
                {
                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        //e -= 2 *en_array_PMO[array[0]][array[1]][array[2]];
                        e += 2 *en_array_OH_Att[array[0]][array[1]][array[2]];
                        
                        if(i==48||i==2)
                        {
                            confused2 += 2 *en_array_OH_Att[array[0]][array[1]][array[2]];
                            if(iqq == 22981)
                            {
                                stdoutlog = fopen(stdoutlogname, "a");
                                fprintf(stdoutlog, "\t PMO Att  %i\n",confused2);
                                fclose(stdoutlog);
                            }
                        }
                        
//                        if(j==14 || j==15)
//                        {
//                            confused2 += 2 *en_array_coloumb_attr[array[0]][array[1]][array[2]];
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "te PM 2*attr  %i\n",confused2);
//                        fclose(stdoutlog);
//                        }
                    }
                    if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                    {
                        e += 2 *long_range_OH_en_attr;
                    }
                }
                if(j%3!=0)
                {
                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        //e -= 2 *en_array_PMH[array[0]][array[1]][array[2]];
                        e += 2 *en_array_HH_Rep[array[0]][array[1]][array[2]];
                        
                        if(i==48||i==2)
                        {
                            confused2 += 2 *en_array_HH_Rep[array[0]][array[1]][array[2]];
                            if(iqq == 22981)
                            {
                                stdoutlog = fopen(stdoutlogname, "a");
                                fprintf(stdoutlog, "\t PMH Rep  %i\n",confused2);
                                fclose(stdoutlog);
                            }
                        }
//                        if(j==14 || j==15)
//                        {
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "te PM 2*rep   %i\n",confused2);
//                        fclose(stdoutlog);
//                        }
                    }
                    if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                    {
                        e += 2 *long_range_HH_en_repul;
                        
                        
                    }
                }
            }
            
        }
        for (int j = 0; j < no_ion_sites; j++)    // calculatw=es energy for solvent interactions
        {
            point_distance(p_l_i[i], i_l[j], array);
            if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
            {
                e += 2 *en_array_PMEu[array[0]][array[1]][array[2]];
                
            }
        }
        for (int j = 0; j < numbercores *lengthpercore *pertrusionpercore + backbone_length; j++)
        {
            if (j % lengthpercore != 0 || j >= numbercores *lengthpercore *pertrusionpercore)    // excludes the monomer in the backbone
            {
                if(i!=j)
                {
                    point_distance(p_l_i[i], p_l_i[j], array);

                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        if(j<numbercores *lengthpercore *pertrusionpercore && i>=numbercores *lengthpercore *pertrusionpercore && (j+1)%lengthpercore==0)
                        {
                            e += en_array_PMPO[array[0]][array[1]][array[2]];
                            
                            if((i==48||i==2) || (j==48||j==2))
                            {
                                confused2 += en_array_PMPO[array[0]][array[1]][array[2]];
                                if(iqq == 22981)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "\t PMPO  %i\n",confused2);
                                    fclose(stdoutlog);
                                }
                            }
                        }
                        else{
                        
                            e += en_array_PMPM[array[0]][array[1]][array[2]];
                            
                            if((i==48||i==2) || (j==48||j==2))
                            {
                                confused2 += en_array_PMPM[array[0]][array[1]][array[2]];
                                if(iqq == 22981)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "\t PMPM  %i\n",confused2);
                                    fclose(stdoutlog);
                                }
                            }
                            
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < no_ion_sites; i++)    // calculatw=es energy for solvent interactions
    {
        for (int j = 0; j < no_ion_sites; j++)    // calculatw=es energy for solvent interactions
        {
            if (i != j)
            {
                point_distance(i_l[i], i_l[j], array);
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    e += en_array_EuEu[array[0]][array[1]][array[2]];
                    //e -= en_array_EuEu[array[0]][array[1]][array[2]];
                    e += en_array_EuEu_Rep[array[0]][array[1]][array[2]];
                    
//                    if(j==14 || i==14)
//                    {
//                        confused2 += en_array_EuEu_Rep[array[0]][array[1]][array[2]];
//                        confused2 += en_array_EuEu[array[0]][array[1]][array[2]];
//                        if(iqq == 22981)
//                        {
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "\tEuEu x1 %i %i  %i\n",i,j,confused2);
//                            fclose(stdoutlog);
//                        }
//                    }
                }

                if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                {
                    e += long_range_EuEu_en_repul;
                }
            }
        }
    }

    
    for (int j = 0; j < 3*no_solvent_sites; j++)    // calculatw=es energy for solvent interactions
    {
        for (int k = 0; k < 3*no_solvent_sites; k++)    // calculatw=es energy for solvent interactions
        {
            if(j!=k)
            {
                point_distance(s_l[j], s_l[k], array);
                
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    //e += en_array[array[0]][array[1]][array[2]];
                }
                if(j%3==0 && k%3!=0)
                {
                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        e += en_array_OH[array[0]][array[1]][array[2]];
                        //e -= 1 *en_array_OH[array[0]][array[1]][array[2]];
                        e += 1 *en_array_OH_Att[array[0]][array[1]][array[2]];
                        
//                        if(j==14 || j==15 || k==14 || k==15)
//                        {
//                            confused2+= 1 *en_array_coloumb_attr[array[0]][array[1]][array[2]];
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "te OH attr %i\n",confused2);
//                        fclose(stdoutlog);
//                        }
                    }
                    if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                    {
                        e += long_range_OH_en_attr;
                       
                    }
                }
                if(j%3!=0 && k%3==0)
                {
                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        e += en_array_OH[array[0]][array[1]][array[2]];
                        //e -= 1 *en_array_OH[array[0]][array[1]][array[2]];
                        e += 1 *en_array_OH_Att[array[0]][array[1]][array[2]];
                        
//                        if(j==14 || j==15 || k==14 || k==15)
//                        {
//                            confused2 += 1 *en_array_coloumb_attr[array[0]][array[1]][array[2]];
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "te HO attr  %i\n",confused2);
//                        fclose(stdoutlog);
//                        }
                    }
                    if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                    {
                        e += long_range_OH_en_attr;
                        
                        
                        
                    }
                }
                if(j%3==0 && k%3==0)
                {
                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        e += en_array_OO[array[0]][array[1]][array[2]];
                        //e -= en_array_OO[array[0]][array[1]][array[2]];
                        e += en_array_OO_Rep[array[0]][array[1]][array[2]];
                        
//                        if(j==14 || j==15 || k==14 || k==15)
//                        {
//                            confused2 += en_array_coloumb_repul[array[0]][array[1]][array[2]];
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "te OO rep  %i\n",confused2);
//                        fclose(stdoutlog);
//                        }
                    }
                    if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                    {
                        e += long_range_OO_en_repul;
                        
                        
                    }
                }
                if(j%3!=0 && k%3!=0)
                {
                    if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                    {
                        e += en_array_HH[array[0]][array[1]][array[2]];
                        //e -= 1 *en_array_HH[array[0]][array[1]][array[2]];
                        e += 1 *en_array_HH_Rep[array[0]][array[1]][array[2]];
                        
//                        if(j==14 || j==15 || k==14 || k==15)
//                        {
//                            confused2 += 1 *en_array_coloumb_repul[array[0]][array[1]][array[2]];
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "te HH rep   %i\n",confused2);
//                            fclose(stdoutlog);
//                        }
                    }
                    if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                    {
                        e += long_range_HH_en_repul;
                    }
                }
            }
        }
    }

    for (int j = 0; j < 3*no_solvent_sites; j++)    // calculatw=es energy for solvent interactions
    {
        for (int k = 0; k < no_ion_sites; k++)    // calculatw=es energy for solvent interactions
        {
            point_distance(s_l[j], i_l[k], array);

            if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
            {
                
            }
            if(j%3==0) //oxygen
            {
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    e += 2 *en_array_OEu[array[0]][array[1]][array[2]];
                    //e -= 2 *en_array_OEu[array[0]][array[1]][array[2]];
                    e += 2 *en_array_OEu_Att[array[0]][array[1]][array[2]];
                    
                    
                    
//                    if(j==14 || j==15)
//                    {
//                        confused2 += 2 *en_array_coloumb_attr[array[0]][array[1]][array[2]];
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "te OEu attr   %i\n",confused2);
//                    fclose(stdoutlog);
//                    }
                }
                if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                {
                    e += 2 *long_range_OEu_en_attr;
                    
                    
                }
            }
            if(j%3!=0) // hydrogen
            {
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    e += 2 *en_array_HEu[array[0]][array[1]][array[2]];
                    //e -= 2 *en_array_HEu[array[0]][array[1]][array[2]];
                    e += 2 *en_array_HEu_Rep[array[0]][array[1]][array[2]];
                    
//                    if(j==14 || j==15)
//                    {
//                        confused2 += 2 *en_array_coloumb_repul[array[0]][array[1]][array[2]];
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "te HEu rep   %i\n",confused2);
//                    fclose(stdoutlog);
//                    }
                    
                   
                }
                if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                {
                    e += 2 *long_range_HEu_en_repul;
                }
            }
        }
    }

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "status q1\n");
//        fclose(stdoutlog);
        
    for (int i = 0; i < numbercores * pertrusionpercore; i++)
    {
        for (int j = 0; j < numbercores * pertrusionpercore; j++)
        {
            if (i != j)
            {
                point_distance(p_l_i[i *lengthpercore], p_l_i[j *lengthpercore], array);

                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
                    //e -= en_array_PMPM[array[0]][array[1]][array[2]];
                    e += en_array_HH_Rep[array[0]][array[1]][array[2]];
                    
                    if(i==0||j==0)
                    {
                        confused2 += en_array_HH_Rep[array[0]][array[1]][array[2]];
                        if(iqq == 22981)
                        {
                            stdoutlog = fopen(stdoutlogname, "a");
                            fprintf(stdoutlog, "\t PMPM Rep 1x %i\n",confused2);
                            fclose(stdoutlog);
                        }
                    }
                }

                if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                {
                    e += long_range_HH_en_repul;
                    
                    if(i==0||j==0)
                    {
                        confused2 += long_range_HH_en_repul;
                        if(iqq == 22981)
                        {
                            stdoutlog = fopen(stdoutlogname, "a");
                            fprintf(stdoutlog, "\t PMPM HH Lonog Rep 1x %i\n",confused2);
                            fclose(stdoutlog);
                        }
                    }
                }

                point_distance(p_l_i[i *lengthpercore + lengthpercore - 1], p_l_i[j *lengthpercore + lengthpercore - 1], array);

                
                if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "status q3 %i  %i  %i\n",i,j,en_array_OO_Rep[array[0]][array[1]][array[2]]);
//                    fclose(stdoutlog);
                    //e -= en_array_POPO[array[0]][array[1]][array[2]];
                    e += en_array_OO_Rep[array[0]][array[1]][array[2]];
                }

                if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
                {
                    e += long_range_OO_en_repul;
                }
            }

            point_distance(p_l_i[i *lengthpercore], p_l_i[j *lengthpercore + lengthpercore - 1], array);

            if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist)
            {
                //e -= 2 *en_array_PMPO[array[0]][array[1]][array[2]];
                e += 2 *en_array_OH_Att[array[0]][array[1]][array[2]];
                
                if(i *lengthpercore==0)
                {
                    confused2+= 2 *en_array_OH_Att[array[0]][array[1]][array[2]];
                    if(iqq == 22981)
                    {
                        stdoutlog = fopen(stdoutlogname, "a");
                        fprintf(stdoutlog, "\t PMPO Att x2  %i\n",confused2);
                        fclose(stdoutlog);
                    }
                }
            }

            if (dist_array[array[0]][array[1]][array[2]] > max_en_dist)
            {
                e += 2 * long_range_OH_en_attr;
                
                if(i *lengthpercore==0)
                {
                    confused2+= 2 * long_range_OH_en_attr;
                    if(iqq == 22981)
                    {
                        stdoutlog = fopen(stdoutlogname, "a");
                        fprintf(stdoutlog, "\t PMPO OH Long Att x2  %i\n",confused2);
                        fclose(stdoutlog);
                    }
                }
            }
        }
    }

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "status q2  %i\n",e/2);
//        fclose(stdoutlog);
//
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Leave it to totalenergy\n");
//        fclose(stdoutlog);
    return (e / 2);
}

int local_energy_op(int *p_l_i, int *s_l_i, int *i_l_i);
int local_energy_op_2(int *p_l_i, int *s_l_i, int *i_l_i,int reset_ind,int inc);
int local_energy(int *p_l_c, int *p_l_i, int *s_l, int *i_l,int inc)    // returns energy of a single spin
{
    int eloc = 0;
    int array[3];

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "\nMade it to local_energy\n");
//        fprintf(stdoutlog, "chosen_mov:%i \t poly_solv:%i \t poly_ion %i \n",chosen_mov,poly_solvent_mov,poly_ion_mov);
//        fprintf(stdoutlog, "reset_indexes[0]+1:%i \t reset_indexes[1]:%i\n",reset_indexes[0]+1,reset_indexes[1]);
//        fprintf(stdoutlog, "reset_indexes[2]:%i \t reset_indexes[3]:%i  \t reset_indexes[4]:%i \t reset_indexes[5]:%i  \t reset_indexes[6]:%i\n",reset_indexes[2],reset_indexes[3],reset_indexes[4],reset_indexes[5],reset_indexes[6]);
//        fclose(stdoutlog);
    

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
//    fclose(stdoutlog);
    
    if (poly_solvent_mov == 0 && poly_ion_mov == 0)    // for polymwe movements only increments between moved monomer
    {
        
        for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++)
        {
            //if(i%lengthpercore !=0)
            //{
            
            if(i>=0 && i < numbercores *lengthpercore * pertrusionpercore)
            {
                if (i % lengthpercore != 0)    // excludes the monomer in the backbone
                {
                    if (poly_point_distance(p_l_c, i - 1, i, array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
                    {
                        if(i==reset_indexes[0]+1)
                        {
                            eloc += 2 *en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                        }
                        else
                        {
                            eloc += en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                        }
                    }
                }
            
                if ((i+1) % lengthpercore != 0)    // excludes the monomer in the backbone
                {
                    if (poly_point_distance(p_l_c, i + 1, i, array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
                    {
                        if(i==reset_indexes[1]-1)
                        {
                            eloc += 2 *en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                        }
                        else
                        {
                            eloc += en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                        }
                    }
                }
            }
            
            if(i>=numbercores *lengthpercore *pertrusionpercore+1 && i < numbercores *lengthpercore *pertrusionpercore + backbone_length)
            {
                if (poly_point_distance(p_l_c, i - 1, i, array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
                {
                    if(i==reset_indexes[0]+1)
                    {
                        eloc += 2 *en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                    }
                    else
                    {
                        eloc += en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                    }
                    
                    if(((i)-(numbercores *lengthpercore *pertrusionpercore))%lengthperchain==0)
                    {
                        eloc -= 2 *en_array_bond[array[0]][array[1]][array[2]];
                    }
                }
                
                if(poly_spin_reference[i]>1 && poly_spin_reference[i] !=5+numbercores*pertrusionpercore)
                {
                    if (poly_point_distance(p_l_c, lengthpercore*(poly_spin_reference[i]-5) +  1, lengthpercore*(poly_spin_reference[i]-5), array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
                    {
                        eloc += 2 *en_array_bond[array[0]][array[1]][array[2]];
                    }
                    
                }
                
                    
                
            }
            if(i>=numbercores *lengthpercore *pertrusionpercore && i < numbercores *lengthpercore *pertrusionpercore + backbone_length-1)
            {
                if (poly_point_distance(p_l_c, i + 1, i, array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
                {
                    if(i==reset_indexes[1]-1)
                    {
                        eloc += 2 *en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                    }
                    else
                    {
                        eloc += en_array_bond[array[0]][array[1]][array[2]];    // abs are taken because array[x] may not be positive
                    }
                    
                    if(((i+1)-(numbercores *lengthpercore *pertrusionpercore))%lengthperchain==0)
                    {
                        eloc -= 2 *en_array_bond[array[0]][array[1]][array[2]];
                    }
                }
                
                
                if(i==numbercores*pertrusionpercore*lengthpercore)
                {
                    if(poly_spin_reference[i]>1 && poly_spin_reference[i] !=5+numbercores*pertrusionpercore)
                    {
                        if (poly_point_distance(p_l_c, lengthpercore*(poly_spin_reference[i]-5) +  1, lengthpercore*(poly_spin_reference[i]-5), array) < max_dist)    // checks the distance between the current monomer and the one below for bond energy alculation
                        {
                            eloc += 2 *en_array_bond[array[0]][array[1]][array[2]];
                        }
                        
                    }
                }
            }
            
            eloc += local_energy_op_2(p_l_i, s_l, i_l,i,inc);
        }
    }

    if (poly_solvent_mov == 1)
    {
        for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++)
        {
            eloc += local_energy_op_2(p_l_i, s_l, i_l,i,inc);
        }
    }

    if (poly_ion_mov == 1)
    {
        eloc += local_energy_op_2(p_l_i, s_l, i_l,reset_indexes[0]+1,inc);
    }

    if(inc==1)
    {
        int rep_diff=no_neg_electrostatic_inter[0]-no_neg_electrostatic_inter[1];
        
        if(rep_diff!=0)
        {
            if(rep_diff>0)
            {
                while(rep_diff>0)
                {
                    eloc +=long_range_ion_en_repul;
                    rep_diff--;
                }
            }
            if(rep_diff<0)
            {
                while(rep_diff<0)
                {
                    eloc -= long_range_ion_en_repul;
                    rep_diff++;
                }
            }
        }
        
        int att_diff=no_pos_electrostatic_inter[0]-no_pos_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_ion_en_attr;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_ion_en_attr;
                att_diff++;
            }
        }
        
        att_diff=no_OH_electrostatic_inter[0]-no_OH_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_OH_en_attr;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_OH_en_attr;
                att_diff++;
            }
        }
        
        att_diff=no_OEu_electrostatic_inter[0]-no_OEu_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_OEu_en_attr;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_OEu_en_attr;
                att_diff++;
            }
        }
        
        att_diff=no_OO_electrostatic_inter[0]-no_OO_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_OO_en_repul;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_OO_en_repul;
                att_diff++;
            }
        }
        
        att_diff=no_HH_electrostatic_inter[0]-no_HH_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_HH_en_repul;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_HH_en_repul;
                att_diff++;
            }
        }
        
        att_diff=no_EuEu_electrostatic_inter[0]-no_EuEu_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_EuEu_en_repul;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_EuEu_en_repul;
                att_diff++;
            }
        }
        
        att_diff=no_HEu_electrostatic_inter[0]-no_HEu_electrostatic_inter[1];
        
        if(att_diff>0)
        {
            while(att_diff>0)
            {
                eloc += long_range_HEu_en_repul;
                att_diff--;
            }
        }
        if(att_diff<0)
        {
            while(att_diff<0)
            {
                eloc -= long_range_HEu_en_repul;
                att_diff++;
            }
        }
    }

    return (eloc / 2);
}

void exclusion(int(&arr)[121], int index)
{
    int ph = index * numberneighbors_L;
    for (int i = 1; i < numberneighbors_L; i++)
    {
        arr[what_to_exclude[ph + i]] = 1;
    }
}

int local_energy_op_2(int *p_l_i, int *s_l_i, int* i_l_i,int reset_ind,int inc)
{
    int eloc = 0;
    //int exclude_indexes[121] = { 0 };
    //int considered_spin = poly_solvent_mov;
//    int i = 2;
//    int j = 0;
//    int k = 0;
//    int x;
//    int y;
//    int z;

    //int phh;
    //int phhh;
    //int phhhh;
    //int phhhhh;

    int index_of_concern;
    int spin_of_concerned_index;
    int o_or_h;
    int spin_ref;
    
    int balance_check =0;
    int balance_check2 =0;
    
    if (poly_solvent_mov == 0 && poly_ion_mov == 0)
    {
        //phh = p_l_i[reset_ind];
        //phhh = p_l_i[(reset_ind) - 1];
        //phhhh = p_l_i[(reset_ind) + 1];
        //phhhhh = (reset_ind + 1) % lengthpercore;
        index_of_concern = p_l_i[reset_ind];
        spin_of_concerned_index =latticepoint_S[p_l_i[reset_ind]];
        
        spin_ref=poly_spin_reference[reset_ind];
    }

    if (poly_solvent_mov == 1)
    {
        index_of_concern = s_l_i[reset_ind];
        spin_of_concerned_index =latticepoint_S[s_l_i[reset_ind]];
        o_or_h = reset_ind%3;
        if(o_or_h==0)
        {
            spin_ref=2;
        }
        else{
            spin_ref=-2;
        }
    }
    
    if (poly_ion_mov == 1)
    {
        index_of_concern = i_l_i[reset_indexes[0] + 1];
        spin_of_concerned_index =latticepoint_S[i_l_i[reset_ind]];
    }

    for (int i = 0; i < 176; i++)
    {

        //if (exclude_indexes[i] == 0)
        //{
            int checked_index = local_en_indexes_op[index_of_concern *176 + i];

//            x = local_en_indexes_op_x[i];
//            y = local_en_indexes_op_y[i];
//            z = local_en_indexes_op_z[i];
            
            if (poly_solvent_mov == 0 && poly_ion_mov == 0)
            {
                int sum = local_en_dist[i];
                int holder = latticepoint_S[checked_index];
                if (holder == 1 || holder >= 5)
                {
                    for(int i = reset_indexes[0]+1;i<reset_indexes[1];i++)
                    {
                        if(checked_index == p_l_i[i])
                        {
                            if((reset_indexes[0]+1)>=numbercores*pertrusionpercore*lengthpercore)
                            {
                                eloc-=int_en_array_PMPM[sum];
                            }
                            else
                            {
                                eloc-=int_en_array_PMPO[sum];
                            }
                            if(spin_ref>=5 && spin_ref < 5+numbercores*pertrusionpercore)
                            {
                                if(holder>=5 && holder<5+numbercores*pertrusionpercore)
                                {
                                    eloc -= int_en_array_HH_Rep[sum];
                                    no_HH_electrostatic_inter[inc]--;
                                }
                                if( holder==5+numbercores*pertrusionpercore)
                                {
                                    eloc -= int_en_array_OH_Att[sum];
                                    no_OH_electrostatic_inter[inc]--;
                                }
                            }
                            if(spin_ref == 5+numbercores*pertrusionpercore)
                            {
                                if(holder>=5 && holder<5+numbercores*pertrusionpercore)
                                {
                                    eloc -= int_en_array_OH_Att[sum];
                                    no_OH_electrostatic_inter[inc]--;
                                }
                                if(holder==5+numbercores*pertrusionpercore)
                                {
                                    eloc -= int_en_array_OO_Rep[sum];
                                    no_OO_electrostatic_inter[inc]--;
                                }
                            }
                        }
                    }
                    if(spin_ref>=5 && spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        if(holder>=5 && holder<5+numbercores*pertrusionpercore)
                        {
                            eloc += 2 *int_en_array_HH_Rep[sum];
                            no_HH_electrostatic_inter[inc]+=2;
                        }
                        if(holder==5+numbercores*pertrusionpercore)
                        {
                            eloc += 2 *int_en_array_OH_Att[sum];
                            no_OH_electrostatic_inter[inc]+=2;
                        }
                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        if(holder>=5 && holder<5+numbercores*pertrusionpercore)
                        {
                            eloc += 2 *int_en_array_OH_Att[sum];
                            no_OH_electrostatic_inter[inc]+=2;
                            
                            if(index_track!=-1 && inc == 0)
                            {
                                radial_dist_PMPO[index_track*14+sum]++;
                                
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "\n rdPMPO  %f: %i  %i  %i\n", radial_dist_PMPO[index_track*numbercores+reset_ind*14+sum],index_track,reset_ind,sum);
//                                fclose(stdoutlog);
                                
                                 balance_check = 3;
                                
                            }
                        }
                        if(holder==5+numbercores*pertrusionpercore)
                        {
                            eloc += 2 *int_en_array_OO_Rep[sum];
                            no_OO_electrostatic_inter[inc]+=2;
                        }
                    }
                    
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        if(holder ==1 || (holder>=5 && holder<5+numbercores*pertrusionpercore))
                        {
                            eloc += 2 *int_en_array_PMPO[sum];
                        }
                        if(holder==5+numbercores*pertrusionpercore)
                        {
                            eloc += 2 *int_en_array_POPO[sum];
                        }
                        
                    }
                    else{
                        if(spin_ref ==1 || (spin_ref>=5 && spin_ref < 5+numbercores*pertrusionpercore))
                        {
                            if( holder==5+numbercores*pertrusionpercore)
                            {
                                eloc += 2 *int_en_array_PMPO[sum];
                            }
                            if(holder ==1 || (holder>=5 && holder<5+numbercores*pertrusionpercore))
                            {
                                eloc += 2 *int_en_array_PMPM[sum];
                            }
                            
                        }
                    }
                }
                if (holder == 2)
                {
                    if(spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_PMO[sum];
                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_POO[sum];
                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_OO_Rep[sum];
                        no_OO_electrostatic_inter[inc]+=2;
                    }
                    if(spin_ref >= 5 && spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_OH_Att[sum];
                        no_OH_electrostatic_inter[inc]+=2;
                    }
                }

                if (holder == -2)
                {
                    if(spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_PMH[sum];
                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_POH[sum];
                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_OH_Att[sum];
                        no_OH_electrostatic_inter[inc]+=2;
                    }
                    if(spin_ref >= 5 && spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_HH_Rep[sum];
                        no_HH_electrostatic_inter[inc]+=2;
                    }
                }
                
                if (holder == 3)
                {
                    if(spin_ref>=5 && spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_HEu_Rep[sum];
                        no_HEu_electrostatic_inter[inc]+=2;

                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_OEu_Att[sum];
                        no_OEu_electrostatic_inter[inc]+=2;
                    }
                    
                    if(spin_ref < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_PMEu[sum];
                    }
                    if(spin_ref == 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_POEu[sum];
                    }
                }
            }

            if (poly_solvent_mov == 1)
            {
                int sum = local_en_dist[i];
                int holder = latticepoint_S[checked_index];
                if (holder == 1 || holder >= 5)
                {
                    if(spin_ref==2 && holder < 5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_PMO[sum];
                    }
                    else{
                        if(spin_ref==-2 && holder < 5+numbercores*pertrusionpercore)
                        {
                            eloc += 2 *int_en_array_PMH[sum];
                        }
                        else
                        {
                            if(spin_ref==2 && holder == 5+numbercores*pertrusionpercore)
                            {
                                eloc += 2 *int_en_array_POO[sum];
                            }
                            if(spin_ref==-2 && holder == 5+numbercores*pertrusionpercore)
                            {
                                eloc += 2 *int_en_array_POH[sum];
                            }
                        }
                    }
                    if(holder == 5+numbercores*pertrusionpercore && spin_ref==2)
                    {
                        eloc += 2 *int_en_array_OO_Rep[sum];
                        no_OO_electrostatic_inter[inc]+=2;
                    }
                    if(holder == 5+numbercores*pertrusionpercore && spin_ref==-2)
                    {
                        eloc += 2 *int_en_array_OH_Att[sum];
                        no_OH_electrostatic_inter[inc]+=2;
                    }
                    if(holder >= 5 && holder < 5+numbercores*pertrusionpercore && spin_ref==2)
                    {
                        eloc += 2 *int_en_array_OH_Att[sum];
                        no_OH_electrostatic_inter[inc]+=2;
                    }
                    if(holder >= 5 && holder < 5+numbercores*pertrusionpercore && spin_ref==-2)
                    {
                        eloc += 2 *int_en_array_HH_Rep[sum];
                        no_HH_electrostatic_inter[inc]+=2;
                    }
                }

                if (holder == 2 && checked_index != ph_solvent_loc[reset_ind])
                {
                    if(spin_ref==-2)
                    {
                        eloc += 2 *int_en_array_OH[sum];
                        eloc += 2 *int_en_array_OH_Att[sum];
                        no_OH_electrostatic_inter[inc]+=2;
                    }
                    if(spin_ref==2)
                    {
                            
                        eloc += 2 *int_en_array_OO[sum];
                        eloc += 2 *int_en_array_OO_Rep[sum];
                        no_OO_electrostatic_inter[inc]+=2;
                    }

                }
                if (holder == -2 && checked_index != ph_solvent_loc[reset_ind])
                {
                    if(spin_ref==2)
                    {
                        eloc += 2 *int_en_array_OH[sum];
                        eloc += 2 *int_en_array_OH_Att[sum];
                        no_OH_electrostatic_inter[inc]+=2;
                    }
                    if(spin_ref==-2)
                    {
                        eloc += 2 *int_en_array_HH[sum];
                        eloc += 2 *int_en_array_HH_Rep[sum];
                        no_HH_electrostatic_inter[inc]+=2;
                    }
                }
                if (holder == 3)
                {
                    if(spin_ref==2)
                    {
                        eloc += 2 *int_en_array_OEu[sum];
                        eloc += 2 *int_en_array_OEu_Att[sum];
                        no_OEu_electrostatic_inter[inc]+=2;
                    }
                    if(spin_ref==-2)
                    {
                        eloc += 2 *int_en_array_HEu[sum];
                        eloc += 2 *int_en_array_HEu_Rep[sum];
                        no_HEu_electrostatic_inter[inc]+=2;
                    }
                    
                }
            }

            if (poly_ion_mov == 1)
            {
                int sum = local_en_dist[i];
                int holder = latticepoint_S[checked_index];
                if (holder == 1 || holder >= 5)
                {
                    if(holder>=5 && holder<5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_HEu_Rep[sum];
                        no_HEu_electrostatic_inter[inc]+=2;
                    }
                    if( holder==5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_OEu_Att[sum];
                        no_OEu_electrostatic_inter[inc]+=2;
                    }
                    if(holder<5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_PMEu[sum];
                    }
                    if(holder==5+numbercores*pertrusionpercore)
                    {
                        eloc += 2 *int_en_array_POEu[sum];
                    }
                }

                if (holder == 2)
                {
                    eloc += 2 *int_en_array_OEu[sum];
                    eloc += 2 *int_en_array_OEu_Att[sum];
                    no_OEu_electrostatic_inter[inc]+=2;
                }
                if (holder == -2)
                {
                    eloc += 2 *int_en_array_HEu[sum];
                    eloc += 2 *int_en_array_HEu_Rep[sum];
                    no_HEu_electrostatic_inter[inc]+=2;
                }
                if (holder == 3)
                {
                    eloc += 2 *int_en_array_EuEu[sum];
                    eloc += 2 *int_en_array_EuEu_Rep[sum];
                    no_EuEu_electrostatic_inter[inc]+=2;
                }
            }
        //}
        
    }

 
    return eloc;
}

void radial_counter(int *p_l_i, int *s_l_i, int* i_l_i,int index_tr)
{
    int seq_counter=0;
    int counter_cores = 0;
    int counter_ions = 0;
    int reset_ind;
    for(int i =0;i<numbercores+no_ion_sites;i++)
    {
        if(i<numbercores)
        {
            reset_ind=i*2+1;
            poly_solvent_mov = 0;
            poly_ion_mov = 0;
        }
        if(i>=numbercores)
        {
            reset_ind=(i-numbercores);
            poly_ion_mov = 1;
        }
        
        int eloc = 0;
        int index_of_concern;
        int spin_of_concerned_index;
        int o_or_h;
        int spin_ref;
        
        int balance_check =0;
        int balance_check2 =0;
        
        if (poly_solvent_mov == 0 && poly_ion_mov == 0)
        {
            index_of_concern = p_l_i[reset_ind];
            spin_of_concerned_index =latticepoint_S[p_l_i[reset_ind]];
            
            spin_ref=poly_spin_reference[reset_ind];
            radial_visits_poly[index_tr]++;
            balance_check2 = -3;
        }
        if (poly_ion_mov == 1)
        {
            index_of_concern = i_l_i[reset_ind];
            spin_of_concerned_index =latticepoint_S[i_l_i[reset_ind]];
            
            radial_visits_ion[index_tr]++;
            balance_check2 = -2;
        }

        int no_func_int=0;
        int no_solvent_int=0;
        
        for (int i = 0; i < 176; i++)
        {

            //if (exclude_indexes[i] == 0)
            //{
                int checked_index = local_en_indexes_op[index_of_concern *176 + i];

    //            x = local_en_indexes_op_x[i];
    //            y = local_en_indexes_op_y[i];
    //            z = local_en_indexes_op_z[i];
                
                if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                {
                    int sum = local_en_dist[i];
                    int holder = latticepoint_S[checked_index];
                    if (holder == 1 || holder >= 5)
                    {
                        if(spin_ref == 5+numbercores*pertrusionpercore)
                        {
                            if(holder>=5 && holder<5+numbercores*pertrusionpercore)
                            {
                                radial_dist_PMPO[index_tr*14+sum]++;
                                balance_check = 3;
                            }
                            if(holder== 5+numbercores*pertrusionpercore)
                            {
                                radial_dist_POPO[index_tr*14+sum]++;
                                balance_check = 3;
                            }
                        }
                        
                    }
                }

                if (poly_ion_mov == 1)
                {
                    int sum = local_en_dist[i];
                    int holder = latticepoint_S[checked_index];
                    if (holder == 1 || holder >= 5)
                    {
                        if( holder==5+numbercores*pertrusionpercore)
                        {
                            radial_dist_EuPO[index_tr*14+sum]++;
                            balance_check = 2;
                            
                            if(sum<=10) no_func_int++;
                        }
                    }

                    if (holder == 2)
                    {
                        radial_dist_EuO[index_tr*14+sum]++;
                        balance_check = 2;
                        
                        if(sum<=10) no_solvent_int++;
                        
                    }
                    
                   
                    
                }
            //}
             
        }

        if(poly_ion_mov==1)
        {
            
            if(poly_ion_mov==1 &&  no_func_int==3)
            {
                part_seq[index_tr]++;
                func_per_seq[index_tr]+=no_func_int;
            }
            else{
                if(poly_ion_mov==1 &&  no_func_int>=2 &&  no_solvent_int>1 && 5>=no_func_int+no_solvent_int)
                {
                    perf_seq[index_tr]++;
                    func_per_seq[index_tr]+=no_func_int;
                }
            }
        }
        
        if(balance_check+balance_check2>0)
        {
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
}

int propose_update()    // returns energy change if a spin would be updated
{
    //int ep = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc);    // compares the energy before and after a movement has been preformed

    loc_en_before=0;
    
    //int ep = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);    // compares the energy before and after a movement has been preformed

    no_OH_electrostatic_inter[0]=0;
    no_OEu_electrostatic_inter[0]=0;
    no_OO_electrostatic_inter[0]=0;
    no_HH_electrostatic_inter[0]=0;
    no_EuEu_electrostatic_inter[0]=0;
    no_HEu_electrostatic_inter[0]=0;
    
    poly_mov();

//        uni_counter++;
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "%i %i\n",uni_counter,latticepoint_S[80]);
//        fclose(stdoutlog);
    
    //        stdoutlog = fopen(stdoutlogname, "a");
    //        fprintf(stdoutlog, "Left polymov()\n");
    //        fclose(stdoutlog);

    int loc_energy_before = 0;
    int loc_energy_after = 0;


    
//        uni_counter++;
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "f%i %i\n",uni_counter,latticepoint_S[80]);
//        fclose(stdoutlog);
    
    if (MoveProposal == 1)
    {
        loc_energy_before = loc_en_before;

        if (poly_solvent_mov == 0 && poly_ion_mov == 0)
        {
            lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        }

        if (poly_solvent_mov == 1)
        {
            solvent_reset2(solvent_loc, ph_solvent_loc, 2);
        }

        if (poly_ion_mov == 1)
        {
            solvent_reset(ion_loc, ph_ion_loc, 3);
        }

        no_OH_electrostatic_inter[1]=0;
        no_OEu_electrostatic_inter[1]=0;
        no_OO_electrostatic_inter[1]=0;
        no_HH_electrostatic_inter[1]=0;
        no_EuEu_electrostatic_inter[1]=0;
        no_HEu_electrostatic_inter[1]=0;
        
        loc_energy_after = local_energy(ph_poly_lattice_coordinates, ph_poly_lattice_indexes, ph_solvent_loc, ph_ion_loc,1);

    }

    return (loc_energy_after - loc_energy_before);
}

void pseudowl()    // A Fake WL Function used to explore the energy landscape before the main WL operation to better enesure all lnge values of the system have been explored
{
    int maxe = 0;    //holds the index of the maximum energy found by the system used to reset the for loop counter if a higher energy is found

    int qhold = q;    // Holds the value of q from before (could be better coded but as far as I know would require significant retooling of some functions)
    q = 2;    // minimum q needed for system to explore all of the possible values of the energy landscape (shape of lng(e) is not important in this test

    double lnf = 1.0;

    double wsk, dice;    // wsk: probability
    int wiggle;
    int wiggletwo;

    // terminal modification factor
    double check_flatness_every = 100000;    // in number of sweeps changed to a large number to encourage visitation of all levels (reduces time in the long run)
    //double check_flatness_every = 10;      // in number of sweeps changed to a large number to encourage visitation of all levels (reduces time in the long run)
    int backup;
    int backuptwo;

    int eold, energie;
    
    energie = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
    //energie -= Eglobalmin;
    energie -=Eglobalmin;
    eold = energie;

    int swtch;
    int found = 0;

    iqq = 0;

    for (int i = 0; i < hist_size; i++) HE[i] = 0;    //init H(E)

    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
        ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
        wl_pseudo_chain[i] = poly_lattice_indexes[i];
        wl_pseudo_chain_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        wl_pseudo_chain_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        wl_pseudo_chain_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
    }

    for (int i = 0; i < 3*no_solvent_sites; i++)
    {
        wl_solvent_loc[i] = solvent_loc[i];
        ph_solvent_loc[i] = solvent_loc[i];
    }

    for (int i = 0; i < no_ion_sites; i++)
    {
        ph_ion_loc[i] = ion_loc[i];
        wl_ion_loc[i] = ion_loc[i];
    }

    int within = 0;
    
    if(((energie < Emaxindex) && (energie > Eminindex)))
    {
        within=1;
    } 
    
    stdoutlog = fopen(stdoutlogname, "a");
       fprintf(stdoutlog, "iqq: %lli \n within check: %i  %i %i  %i \n",iqq,within,energie,Eminindex,Emaxindex);
       fclose(stdoutlog);
    
    for (int k = 0; k < check_flatness_every; k++)
    {
        for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + (backbone_length - numbercores); i++)    // this does 1 MC sweep
        {
            energie = eold + propose_update();    // calculate energy of updated configuration

            // reject because you got out of bounds of energy window

            if (MoveProposal == 0 || (within==1 && ((energie > Emaxindex) || (energie < Eminindex))) || (within==0 && myid<min_max_id && (energie > Emaxindex)) || (within==0 && myid>=min_max_id && (energie < Eminindex)))   // boundary check
            {
                //                stdoutlog = fopen(stdoutlogname, "a");
                //                fprintf(stdoutlog, "\nbc Move Proposal %i\n",MoveProposal);
                //                fclose(stdoutlog);

                if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                {
                    //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                    lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                    poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
                    //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

//                    if (reset_indexes[4] == 2)
//                    {
//                        detachment(random_solvent_site);
//                        reset_indexes[4] = 0;
//                        reset_indexes[5] = 0;
//                        reset_indexes[6] = 0;
//                    }
                }

                if (poly_solvent_mov == 1)
                {
                    solvent_reset2(solvent_loc, wl_solvent_loc, 2);
                    solvent_reset2(ph_solvent_loc, solvent_loc, 2);
                }

                if (poly_ion_mov == 1)
                {
                    solvent_reset(ion_loc, wl_ion_loc, 3);
                    solvent_reset(ph_ion_loc, ion_loc, 3);

//                    if (reset_indexes[4] == 1)
//                    {
//                        reattachment(reset_indexes[5], reset_indexes[6]);
//                        reset_indexes[4] = 0;
//                        reset_indexes[5] = 0;
//                        reset_indexes[6] = 0;
//                    }
//
//                    if (reset_indexes[4] == 2)
//                    {
//                        detachment(reset_indexes[5]);
//                        reset_indexes[4] = 0;
//                        reset_indexes[5] = 0;
//                        reset_indexes[6] = 0;
//                    }

                    //reset_indexes[4]=0;
                }

                energie = eold;
                lngE[energie-energie%500] = 1;    // the real lngE of the system is set to a value of 1 if visitation occured at a specific energy (can be subtracted later) to force visitation in the real wl process
            }
            else    // calculate acceptance propbability
            {
                if(within)
                {
                lngE[energie-energie%500] = 1;
                }
                //                lngE[energie] = 1;
                //                if (energie > maxe)    // checks to see if new energy is higher than max e
                //                {
                //                    maxe = energie;
                //                    k = 0;    // reset for loop iterative
                //                }

                //
                //                   // the regular wl process is kept using a place holder pseudolngE to force the system to not move to move to higher energies
                if(within)
                {
                dice = (1.0* rand() / (RAND_MAX + 1.0));    // roll the dice
                wsk = exp(pseudolngE[eold-eold%500] - pseudolngE[energie-energie%500]);    // calculate acceptance probability
                    }
                                else{
                                    wsk=1;
                                    dice=0;
                                }
                                
                    
                if (dice > wsk)    // reject
                {
                    if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                    {
                        //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
                        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

//                        if (reset_indexes[4] == 2)
//                        {
//                            detachment(random_solvent_site);
//                            reset_indexes[4] = 0;
//                            reset_indexes[5] = 0;
//                            reset_indexes[6] = 0;
//                        }
                    }

                    if (poly_solvent_mov == 1)
                    {
                        solvent_reset2(solvent_loc, wl_solvent_loc, 2);
                        solvent_reset2(ph_solvent_loc, solvent_loc, 2);
                    }

                    if (poly_ion_mov == 1)
                    {
                        solvent_reset(ion_loc, wl_ion_loc, 3);
                        solvent_reset(ph_ion_loc, ion_loc, 3);

//                        if (reset_indexes[4] == 1)
//                        {
//                            reattachment(reset_indexes[5], reset_indexes[6]);
//                            reset_indexes[4] = 0;
//                            reset_indexes[5] = 0;
//                            reset_indexes[6] = 0;
//                        }
//
//                        if (reset_indexes[4] == 2)
//                        {
//                            detachment(reset_indexes[5]);
//                            reset_indexes[4] = 0;
//                            reset_indexes[5] = 0;
//                            reset_indexes[6] = 0;
//                        }
//
//                        reset_indexes[4] = 0;
                    }

                    energie = eold;

                }
                else
                {
                    eold = energie;    // accept
                    if(((energie < Emaxindex) && (energie > Eminindex)))
                                       {
                                           within=1;
                                       }
//                    if(reset_indexes[4] == 1)
//                    {
//                        current_attach--;
//                    }
//                    if(reset_indexes[4] == 2)
//                    {
//                        current_attach++;
//                    }
                    
                    if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                    {
                        //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

                        

                        reset_indexes[4] = 0;
                    }

                    if (poly_solvent_mov == 1)
                    {
                        solvent_reset2(solvent_loc, ph_solvent_loc, 2);
                        solvent_reset2(wl_solvent_loc, solvent_loc, 2);

                    }

                    if (poly_ion_mov == 1)
                    {
                        solvent_reset(ion_loc, ph_ion_loc, 3);
                        solvent_reset(wl_ion_loc, ion_loc, 3);
                        reset_indexes[4] = 0;
                        reset_indexes[5] = 0;
                        reset_indexes[6] = 0;
                    }

                    //
                    //                    if(poly_solvent_mov==0)
                    //                    {
                    //                           //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    //                        lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    //                        poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                    //
                    //                        lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    //                        poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                    //
                    //
                    //                        if(core_movement == 1)
                    //                        {
                    //                            lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[2], reset_indexes[3]);
                    //                            poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[2], reset_indexes[3]);
                    //
                    //                            lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[2], reset_indexes[3]);
                    //                            poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[2], reset_indexes[3]);
                    //
                    //                        }

                    //                        reset_indexes[4]=0;
                    //                    }

                    //                    if(poly_solvent_mov==1)
                    //                    {
                    //                        solvent_reset(solvent_loc,ph_solvent_loc,2);
                    //                        solvent_reset(wl_solvent_loc,solvent_loc,2);
                    //                    }

                    //                    if(poly_ion_mov==1)
                    //                    {
                    //                        solvent_reset(ion_loc,ph_ion_loc,3);
                    //                        solvent_reset(wl_ion_loc,ion_loc,3);
                    //                        reset_indexes[4]=0;
                    //                    }

                    ////                    lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    ////                    lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    ////                    poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                    ////                    solvent_reset(solvent_loc,ph_solvent_loc);
                    ////
                    ////
                    ////                    lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    ////                    lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                    ////                    poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                    ////                    solvent_reset(wl_solvent_loc,wl_solvent_loc);
                    //                    if(myid==0 && energie == minimum_en_config && printed == 0)
                    //                    {
                    //                        printed = 1;
                    //
                    //                        sprintf(filename, "Minimum_Config_Poly.txt");
                    //                        if ((file = fopen(filename, "w")) == NULL)
                    //                        {
                    //                            stdoutlog = fopen(stdoutlogname, "a");
                    //                            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                    //                            fclose(stdoutlog);
                    //                            MPI_Abort(MPI_COMM_WORLD, 1);
                    //                        }

                    //                        else
                    //                        {
                    //                            for (int i = 0; i <= lengthpercore*pertrusionpercore*lengthpercore+ backbone_length; i++)
                    //                            {
                    //                                fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
                    //                            }

                    //                        }

                    //                        fclose(file);
                    //
                    //                        sprintf(filename, "Minimum_Config_Solvent.txt");
                    //                        if ((file = fopen(filename, "w")) == NULL)
                    //                        {
                    //                            stdoutlog = fopen(stdoutlogname, "a");
                    //                            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                    //                            fclose(stdoutlog);
                    //                            MPI_Abort(MPI_COMM_WORLD, 1);
                    //                        }

                    //                        else
                    //                        {
                    //                            for (int i = 0; i <= no_solvent_sites; i++)
                    //                            {
                    //                                fprintf(file, "%i\t%i\n", i, solvent_loc[i]);
                    //                            }

                    //                        }

                    //                        fclose(file);
                    //                    }

                    //sysparam(energie);
                    //                }

                    //
                    //                   //if ((energie > eold))
                    //                   //{
                    //                   //    latticepoint[wiggle] = backup;
                    //                   //    latticepoint[wiggletwo] = backuptwo;
                    //                   //    energie = eold;
                    //                   //}

                    //                   //if ((energie < eold))
                    //                   //{
                    //                   //    k = 0;
                    //                   //    eold = energie;
                }
            }
//            current_attach=0;
//            if(vicinity_check)
//            {
//                for(int i=0;i<numbercores*pertrusionpercore;i++)
//                {
//                    attach_track[energie]+=current_attachment[i*lengthpercore+lengthpercore-1];
//                }
//
//            //attach_track[ep] += current_attach;
//            attach_track_no[energie] ++;
//            }

            
//            if (iqq ==22981)
//            {
//                int ep = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
//                if(ep-Eglobalmin!=energie)
//            {
//
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "Difference %i   %i    %i   %i  %lli   %i\n",MoveProposal,ep,energie,chosen_mov,iqq,i);
//                fclose(stdoutlog);
//                MPI_Abort(MPI_COMM_WORLD, 1);
//                }
//            }
            if(within)
                                      {
            pseudolngE[energie-energie%500] += lnf;
            lngE[energie-energie%500] = 1;
        }
            if (energie == 0)
            {
                stdoutlog = fopen(stdoutlogname, "a");
                for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
                {
                    fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2], ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
                }

                fprintf(stdoutlog, "\n");
                fclose(stdoutlog);
            }
        }

        if (((iqq % 200000 == 0) ))
        {
            sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim_S, q, myid, -1);
            if ((file = fopen(filename, "w")) == NULL)
            {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                for (int i = Eminindex; i <= Emaxindex; i++)
                {
                    if (pseudolngE[i] > .5)
                    {
                        fprintf(file, "%i\t%i\t%e\t%e\t%e\t%f\t%f\n", i, (i + (Eglobalmin)), 0.0, lngE[i], pseudolngE[i], rog[i] / visits[i], tortuosity[i] / visits[i]);
                    }
                }
            }

            fclose(file);

             stdoutlog = fopen(stdoutlogname, "a");
                       fprintf(stdoutlog, "iqq: %lli \n within check: %i  %i %i  %i \n",iqq,within,energie,Eminindex,Emaxindex);
            for (int i = 0; i < numberspins_S; i++)
            {
                fprintf(stdoutlog, "%4i", latticepoint_S[i]);
                if ((i + 1) % L1dim_S_xy == 0)
                    fprintf(stdoutlog, "\n");
                if ((i + 1) % (L1dim_S_xy *L1dim_S_xy) == 0)
                    fprintf(stdoutlog, "\n");

            }

            fprintf(stdoutlog, "\n");
            fprintf(stdoutlog, "\n");
            fclose(stdoutlog);

            stdoutlog = fopen(stdoutlogname, "a");
            for (int i = 0; i < numberspins_L; i++)
            {
                fprintf(stdoutlog, "%4i", latticepoint_L[i]);
                if ((i + 1) % L1dim_L_xy == 0)
                    fprintf(stdoutlog, "\n");
                if ((i + 1) % (L1dim_L_xy *L1dim_L_xy) == 0)
                    fprintf(stdoutlog, "\n");

            }

            fclose(stdoutlog);

            stdoutlog = fopen(stdoutlogname, "a");
            for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
            {
                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2], ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
            }

            fprintf(stdoutlog, "\n");
            fclose(stdoutlog);
            
            int ep = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
            if (energie != ep-Eglobalmin)
            {
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Difference %i   %i    %i   %i\n",MoveProposal,ep,energie,chosen_mov);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        iqq++;
    }

    q = qhold;    // when the pseudowl process is done q is reset to the system value

    // The occucpied placeholder lnge and the values of the real lnge are out putted into a prelim file
    sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim_S, q, myid, 0);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i <= Emaxindex; i++)
        {
            fprintf(file, "%i\t%i\t%e\t%e\t%e\t%f\t%f\n", i, (i + (Eglobalmin)), 0.0, lngE[i], pseudolngE[i], rog[i] / visits[i], tortuosity[i] / visits[i]);
        }

        fclose(file);
    }
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "OK  pseudowl  %i  %i\n",eold,energie);
    fclose(stdoutlog);
}

int histflat(int imin, int imax, double ratio)
{
    // checks flatness of histograms for _all_ walkers in energy window
    int myflat, otherflat;
    int flatproc, ioffset;
    int flatprocr, ioffsetr;
    int flatness_crit = 1;

    int multimyflat = 0;    // bool value of wheter in of the concurrent walkers in a window are flat

    int merge_crit = 1;    // merge_hist criteria

    // check own flatness first
    if (flatness_crit == 0)    // Zhou criterion
    {
        // NOT AVAILABLE YET!
        // TODO: CALCULATE histmin FOR ZHOU CRITERION
        //    myflat=1;
        //    for (int x=imin;x<=imax;x++)
        //    if (HE[x] < histmin) myflat=0;
    }
    else if (flatness_crit == 1)    // "Original" percentage criterion
    {
        myflat = 1;
        double minval;
        minval = HE[imin];    // take GS hits as reference

        /*if (myid == 0)
        {
            minval = HE[10];
        }*/

        minval = DBL_MAX;    // Ridiculously large arbitrary value
        for (int x = imin; x <= imax; x++)
        {
            if ((lngE[x] > 0 || HE[x] > 0) && HE[x] <= minval)
            {
                minval = HE[x];
            }
        }

        double average = 0.0;
        double count = 0.0;

        for (int x = imin; x <= imax; x++)
        {
            /*if (((x > 5) || (x == 4)) && (HE[x] < minval))
                minval = HE[x]; */
            //(I am not sure right now why I included the first condition ...)

            if (lngE[x] > 0 || HE[x] > 0)
            {
                average += HE[x];
                count++;
            }
        }

        average /= count;

        flatratio = (ratio *average);
        flatmin = (minval);

        if (minval < (ratio *average))
            myflat = 0;
    }

    // now talk to all the other walkers in the energy window
    // (! this whole thing can be reduced to an MPI_Allreduce once there
    // are separate communicators for energy windows !)
//    if (local_dup_group == localgroup)
//    {
//        if (merge_crit == 1 && merge_hists == 1)    // check flatness of other walkers in window
//        {
//            if (myid == dup_headproc)    // 'root' in energy window, receive individual flatnesses
//            {
//                if (myflat == 1)    //Main node multimyflat (checks for the flat process in an energy window
//                {
//                    flatproc = myid;    // id of flat process
//                    ioffset = dup_headproc;    // offset from the main energy window node
//                    multimyflat = 1;
//                }
//
//                for (int i = (dup_headproc + 1); i < numprocs; i++)
//                {
//                    MPI_Recv(&otherflat, 1, MPI_INT, i, 66, MPI_COMM_WORLD, &status);
//
//                    if (otherflat == 1)    //sets the value of the two variable based on information recieved from the process each one is communicating with
//                    {
//                        flatproc = i;
//                        ioffset = myid - dup_headproc;
//                        multimyflat = 1;
//                    }
//
//                    myflat *= otherflat;    // all have to be '1' to give '1' in the end (manual '&&')
//                }
//
//                for (int i = (dup_headproc + 1); i < numprocs; i++)    // and let everybody know
//                {
//                    MPI_Send(&myflat, 1, MPI_INT, i, 88, MPI_COMM_WORLD);
//                    MPI_Send(&multimyflat, 1, MPI_INT, i, 86, MPI_COMM_WORLD);
//                    if (multimyflat == 1)    // if multimyflat is found to be equal to one flatproc and ioffset
//                    {
//                        MPI_Send(&flatproc, 1, MPI_INT, i, 90, MPI_COMM_WORLD);
//                        MPI_Send(&ioffset, 1, MPI_INT, i, 92, MPI_COMM_WORLD);
//                    }
//                }
//            }
//            else    // send individual flatness and receive merged status variables
//            {
//                MPI_Send(&myflat, 1, MPI_INT, dup_headproc, 66, MPI_COMM_WORLD);
//                MPI_Recv(&otherflat, 1, MPI_INT, dup_headproc, 88, MPI_COMM_WORLD, &status);
//                MPI_Recv(&multimyflat, 1, MPI_INT, dup_headproc, 86, MPI_COMM_WORLD, &status);
//                if (multimyflat == 1)
//                {
//                    MPI_Recv(&flatprocr, 1, MPI_INT, dup_headproc, 90, MPI_COMM_WORLD, &status);
//                    flatproc = flatprocr;
//                    MPI_Recv(&ioffsetr, 1, MPI_INT, dup_headproc, 92, MPI_COMM_WORLD, &status);
//                    ioffset = ioffsetr;
//                }
//
//                myflat = otherflat;    // replace individual flatness by merged
//            }
//
//            if (multimyflat == 1)
//            {
//                if (myid != flatproc)    // non flat process recieve merged flat process density of states and the myflat status flagged is called to perform the lnge merging procedures in the wl routine
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "Proc %3i: recievied flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\t flatproc: %3i\n", myid, flatproc, ioffset, multimyflat, flatproc);
//                    fclose(stdoutlog);
//                    MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, flatproc, 17, MPI_COMM_WORLD, &status);
//                    for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];    // overrides density of states of non flat processes
//                    myflat = 1;
//                }
//                else    // flat process sends out its density of states values
//                {
//                    for (int i = (dup_headproc); i < numprocs; i++)
//                    {
//                        if (myid == flatproc && i != flatproc)
//                        {
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "Proc %3i: sent flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\t flatproc: %3i\t %i\n", myid, myid - (myid % multiple) + i, ioffset, multimyflat, flatproc, i);
//                            fclose(stdoutlog);
//                            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, i, 17, MPI_COMM_WORLD);
//                            myflat = 1;
//                        }
//                    }
//                }
//            }
//        }
//    }

    if (multiple > 1)
    {
        if (merge_crit == 0 && merge_hists == 1)    // check flatness of other walkers in window
        {
            if (myid % multiple == 0)    // 'root' in energy window, receive individual flatnesses
            {
                for (int i = 1; i < multiple; i++)
                {
                    MPI_Recv(&otherflat, 1, MPI_INT, myid + i, 66, MPI_COMM_WORLD, &status);
                    myflat *= otherflat;    // all have to be '1' to give '1' in the end (manual '&&')
                }

                for (int i = 1; i < multiple; i++)    // and let everybody know
                {
                    MPI_Send(&myflat, 1, MPI_INT, myid + i, 88, MPI_COMM_WORLD);
                }
            }
            else    // send individual flatness and receive 'merged' flatness
            {
                MPI_Send(&myflat, 1, MPI_INT, myid - (myid % multiple), 66, MPI_COMM_WORLD);
                MPI_Recv(&otherflat, 1, MPI_INT, myid - (myid % multiple), 88, MPI_COMM_WORLD, &status);
                myflat = otherflat;    // replace individual flatness by merged
            }
        }

        if (merge_crit == 1 && merge_hists == 1)    // check flatness of other walkers in window
        {
            if (myid % multiple == 0)    // 'root' in energy window, receive individual flatnesses
            {
                if (myflat == 1)    //Main node multimyflat (checks for the flat process in an energy window
                {
                    flatproc = myid;    // id of flat process
                    ioffset = (myid - (myid % multiple)) - myid;    // offset from the main energy window node
                    multimyflat = 1;
                }

                for (int i = 1; i < multiple; i++)
                {
                    MPI_Recv(&otherflat, 1, MPI_INT, myid + i, 66, MPI_COMM_WORLD, &status);

                    if (otherflat == 1)    //sets the value of the two variable based on information recieved from the process each one is communicating with
                    {
                        flatproc = (myid + i);
                        ioffset = (myid - (myid % multiple)) - myid;
                        multimyflat = 1;
                    }

                    myflat *= otherflat;    // all have to be '1' to give '1' in the end (manual '&&')
                }

                for (int i = 1; i < multiple; i++)    // and let everybody know
                {
                    MPI_Send(&myflat, 1, MPI_INT, myid + i, 88, MPI_COMM_WORLD);
                    MPI_Send(&multimyflat, 1, MPI_INT, myid + i, 86, MPI_COMM_WORLD);
                    if (multimyflat == 1)    // if multimyflat is found to be equal to one flatproc and ioffset
                    {
                        MPI_Send(&flatproc, 1, MPI_INT, myid + i, 90, MPI_COMM_WORLD);
                        MPI_Send(&ioffset, 1, MPI_INT, myid + i, 92, MPI_COMM_WORLD);
                    }
                }
            }
            else    // send individual flatness and receive merged status variables
            {
                MPI_Send(&myflat, 1, MPI_INT, myid - (myid % multiple), 66, MPI_COMM_WORLD);
                MPI_Recv(&otherflat, 1, MPI_INT, myid - (myid % multiple), 88, MPI_COMM_WORLD, &status);
                MPI_Recv(&multimyflat, 1, MPI_INT, myid - (myid % multiple), 86, MPI_COMM_WORLD, &status);
                if (multimyflat == 1)
                {
                    MPI_Recv(&flatprocr, 1, MPI_INT, myid - (myid % multiple), 90, MPI_COMM_WORLD, &status);
                    flatproc = flatprocr;
                    MPI_Recv(&ioffsetr, 1, MPI_INT, myid - (myid % multiple), 92, MPI_COMM_WORLD, &status);
                    ioffset = ioffsetr;
                }

                myflat = otherflat;    // replace individual flatness by merged
            }

            if (multimyflat == 1)
            {
                if (myid != flatproc)    // non flat process recieve merged flat process density of states and the myflat status flagged is called to perform the lnge merging procedures in the wl routine
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "Proc %3i: recievied flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\n", myid, flatproc, ioffset, multimyflat);
                    fclose(stdoutlog);
                    MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, flatproc, 77, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];    // overrides density of states of non flat processes
                    myflat = 1;
                }
                else    // flat process sends out its density of states values
                {
                    for (int i = 0; i < multiple; i++)
                    {
                        if (myid == flatproc && myid - (myid % multiple) + i != flatproc)
                        {
                            stdoutlog = fopen(stdoutlogname, "a");
                            fprintf(stdoutlog, "Proc %3i: sent flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\n", myid, myid - (myid % multiple) + i, ioffset, multimyflat);
                            fclose(stdoutlog);
                            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple) + i, 77, MPI_COMM_WORLD);
                            myflat = 1;
                        }
                    }
                }
            }
        }
    }

    return (myflat);
    // note: by now, myflat refers to the 'collective' flatness in the energy window,
    // not the flatness of an individual walker
}

void eye()
{
    
//    for(int i=0;i<numberspins_S;i++)
//    {
//        if(latticepoint_S[i]==2)
//        {
//            int yes=0;
//            for(int ii=0;ii<no_solvent_sites;ii++)
//            {
//                if(i==solvent_loc[ii])
//                {yes=1;}
//
//
//            }
//            if(yes==0)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "solvent no\n");
//                fclose(stdoutlog);
//
//                MPI_Abort(MPI_COMM_WORLD,1);
//            }
//        }
//
//        if(latticepoint_S[i]==1)
//        {
//            int yes=0;
//            for(int ii=0;ii<lengthpercore*pertrusionpercore*numbercores;ii++)
//            {
//                if(i==poly_lattice_indexes[ii])
//                {yes=1;}
//            }
//            if(yes==0)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "poly no\n");
//                fclose(stdoutlog);
//
//                MPI_Abort(MPI_COMM_WORLD,1);
//            }
//        }
//    }
    
    for(int i=0;i<lengthpercore*pertrusionpercore*numbercores+backbone_length;i++)
    {

        if(i<lengthpercore*pertrusionpercore*numbercores && i%lengthpercore!=0)
        {
        if (latticepoint_S[poly_lattice_indexes[i]] != poly_spin_reference[i] )
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "veye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i poly_spin_ref %i   %i",i,latticepoint_S[poly_lattice_indexes[i]],poly_spin_reference[i],poly_lattice_indexes[i]);
            fclose(stdoutlog);
            
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            
            if (latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]] != poly_spin_reference[i])
            {
                int x= latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]];
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "neye iii:%i  i:%i\t %i\t %i\t %i\n",iii, i,neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii],x,latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]]);
                fclose(stdoutlog);

                MPI_Abort(MPI_COMM_WORLD,1);
            }
        }
        }
        for(int j=0;j<lengthpercore*pertrusionpercore*numbercores+backbone_length;j++)
        {
            if(j<lengthpercore*pertrusionpercore*numbercores && j%lengthpercore!=0)

            {
                if((int_point_distance(poly_lattice_indexes[i],poly_lattice_indexes[j])<4 || poly_lattice_indexes[i]==poly_lattice_indexes[j]) && i!=j )
                {
                    float x = point_distance(poly_lattice_indexes[i],poly_lattice_indexes[j]);
                    int y = int_point_distance(poly_lattice_indexes[i],poly_lattice_indexes[j]);
                    
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],poly_lattice_indexes[j]),poly_lattice_indexes[i],poly_lattice_indexes[j],x,y);
                    fclose(stdoutlog);
                    
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],poly_lattice_indexes[j]),poly_lattice_indexes[i],poly_lattice_indexes[j],x,y);
                    fclose(stdoutlog);
                    MPI_Abort(MPI_COMM_WORLD,1);
                }
            
            
            }
        }
        
        for (int j = 0; j < 3*no_solvent_sites; j++) // calculatw=es energy for solvent interactions
        {
            if(( int_point_distance( poly_lattice_indexes[i],solvent_loc[j])<4 || poly_lattice_indexes[i]==solvent_loc[j]))
            {
                float x = point_distance(poly_lattice_indexes[i],solvent_loc[j]);
                int y = int_point_distance(poly_lattice_indexes[i],solvent_loc[j]);
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "sol i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],solvent_loc[j]),poly_lattice_indexes[i],solvent_loc[j],x,y);
                fclose(stdoutlog);
      
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            if (latticepoint_S[solvent_loc[j]] != 2 && j%3 == 0)
            {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "sveye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i",i,latticepoint_S[solvent_loc[i]]);
                fclose(stdoutlog);
                
                MPI_Abort(MPI_COMM_WORLD,1);
                
            }
            if (latticepoint_S[solvent_loc[j]] != -2 && j%3 !=0)
            {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "s2veye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i",i,latticepoint_S[solvent_loc[i]]);
                fclose(stdoutlog);
                
                MPI_Abort(MPI_COMM_WORLD,1);
                
            }
//            for (int iii = 0; iii < numberneighbors_L; iii++)
//            {
//
//                if (latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[j] + iii]] != 2)
//                {
//                    int x= latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[j] + iii]];
//
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "sneye iii:%i  i:%i\t %i\t %i\t %i  %i\n",iii, i,neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii],x,latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[i] + iii]],latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[j] + iii]]);
//                    fclose(stdoutlog);
//
//                    MPI_Abort(MPI_COMM_WORLD,1);
//                }
//            }
            
        }
        for (int j = 0; j < no_ion_sites; j++) // calculatw=es energy for solvent interactions
          {
              if(( int_point_distance( poly_lattice_indexes[i],ion_loc[j])<4 || poly_lattice_indexes[i]==ion_loc[j]))
              {
                  float x = point_distance(poly_lattice_indexes[i],ion_loc[j]);
                  int y = int_point_distance(poly_lattice_indexes[i],ion_loc[j]);
                  
                  stdoutlog = fopen(stdoutlogname, "a");
                  fprintf(stdoutlog, "iol i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],ion_loc[j]),poly_lattice_indexes[i],ion_loc[j],x,y);
                  fclose(stdoutlog);
        
                  MPI_Abort(MPI_COMM_WORLD,1);
              }
              if (latticepoint_S[ion_loc[j]] != 3)
              {
                  stdoutlog = fopen(stdoutlogname, "a");
                  fprintf(stdoutlog, "iveye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i",j,latticepoint_S[ion_loc[j]]);
                  fclose(stdoutlog);
                  
                  MPI_Abort(MPI_COMM_WORLD,1);
              }
              for (int iii = 0; iii < numberneighbors_L; iii++)
              {
                  
                  if (latticepoint_L[neighbor_L[numberneighbors_L * ion_loc[j] + iii]] != 3)
                  {
                      int x= latticepoint_L[neighbor_L[numberneighbors_L * ion_loc[j] + iii]];
                      
                      stdoutlog = fopen(stdoutlogname, "a");
                      fprintf(stdoutlog, "ineye iii:%i  i:%i\t %i\t %i\t %i\n",iii, i,neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii],x,latticepoint_L[neighbor_L[numberneighbors_L * ion_loc[i] + iii]]);
                      fclose(stdoutlog);

                      MPI_Abort(MPI_COMM_WORLD,1);
                  }
              }
              
          }
    }
    
}

void init_lattice(double emin, double emax)    // Changes made to bctype check and made latticepoint[numbersign equal to 0
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to init lattice()\n");
    fclose(stdoutlog);
    int e, r;

    latticepoint_S = (int*) malloc((numberspins_S + 2 + 1) *sizeof(int));
    latticepoint_L = (int*) malloc((numberspins_L + 2 + 1) *sizeof(int));
    lattice_polymer = (int*) malloc(((numberspins_S *3) + 2 + 1) *sizeof(int));
    // Note: we reserve space for 3 extra 'spins':
    // 2 extra to store fix values used for certain boundary conditions
    // 1 extra to carry an replica-id flag

    // find a fast way to create valid initial configurations
    // start at either Eglobalmin or Eglobalmax and change spins randomly
    // until energy is in the middle third of local energy range
    // global maximum of g(E) id at E/N=-0.2

    //initialize lattice
    for (int i = 0; i < numberspins_S; i++)
    {
        latticepoint_S[i] = 0;
        lattice_polymer[i *3] = i;
        lattice_polymer[i *3 + 1] = -1;
        lattice_polymer[i *3 + 2] = -1;
    }

    for (int i = 0; i < numberspins_L; i++)
    {
        latticepoint_L[i] = 0;
    }

    init_poly_cores();
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "1st\n");
    fclose(stdoutlog);
    
    stdoutlog = fopen(stdoutlogname, "a");
    for (int i = 0; i < numberspins_S; i++)
    {
        fprintf(stdoutlog, "%4i",latticepoint_S[i]);
        if ((i + 1) % L1dim_S_xy == 0)
            fprintf(stdoutlog, "\n");
        if ((i + 1) % (L1dim_S_xy *L1dim_S_xy) == 0)
            fprintf(stdoutlog, "\n");

    }

    fprintf(stdoutlog, "\n");
    fprintf(stdoutlog, "\n");
    fclose(stdoutlog);

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "2nd\n");
    fclose(stdoutlog);
    
    stdoutlog = fopen(stdoutlogname, "a");
    for (int i = 0; i < numberspins_L; i++)
    {
        fprintf(stdoutlog, "%4i", latticepoint_L[i]);
        if ((i + 1) % L1dim_L_xy == 0)
            fprintf(stdoutlog, "\n");
        if ((i + 1) % (L1dim_L_xy *L1dim_L_xy) == 0)
            fprintf(stdoutlog, "\n");

    }
    fclose(stdoutlog);
    
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Left init polycore()\n");
    fclose(stdoutlog);

    int array[3] = { 0, 0, 0 };

    for (int i = 0; i < numberspins_S; i++)
    {
        general_coord(array, i);

        //stdoutlog = fopen(stdoutlogname, "a");
        //fprintf(stdoutlog, "array: %i \t%i\t%i \n",array[0],array[1],array[2]);
        //fclose(stdoutlog);

        
    }

    e = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc. %i: Initialized lattice with energy e=%i, create setup with %lf<e<%lf\n", myid, e, (emin + (emax - emin) / 3), (emin + 2 *(emax - emin) / 3));
    fclose(stdoutlog);

    //for (int i = 0; i < numberspins; i++)
    //{
    //    latticepoint[i] = (i + count) % (q - 1);
    //    if ((i + 1) % L1dim == 0)
    //    {
    //        int rhld = rand() % (q - 1);
    //        while (rhld == count)
    //        {
    //            rhld = rand() % (q - 1);
    //            count = rhld;
    //        }

    //    }

    //}

    stdoutlog = fopen(stdoutlogname, "a");
    for (int i = 0; i < numberspins_S; i++)
    {
        fprintf(stdoutlog, "%i: %i ",i,latticepoint_S[i]);
        if ((i + 1) % L1dim_S_xy == 0)
            fprintf(stdoutlog, "\n");
        if ((i + 1) % (L1dim_S_xy *L1dim_S_xy) == 0)
            fprintf(stdoutlog, "\n");

    }

    fprintf(stdoutlog, "\n");
    fprintf(stdoutlog, "\n");
    fclose(stdoutlog);

    stdoutlog = fopen(stdoutlogname, "a");
    for (int i = 0; i < numberspins_L; i++)
    {
        fprintf(stdoutlog, "%4i", latticepoint_L[i]);
        if ((i + 1) % L1dim_L_xy == 0)
            fprintf(stdoutlog, "\n");
        if ((i + 1) % (L1dim_L_xy *L1dim_L_xy) == 0)
            fprintf(stdoutlog, "\n");

    }

    int array_inc = 0;
    int sign = -1;
    int column_inc = -1;
    int core_backbone_array[5] = { 1, 3, 4, 6 , 7 };
    array_inc = 0;
    for (int i = 0; i < backbone_length; i++)
    {
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "%i  %i\n", poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i],latticepoint_S[poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i]]);
        fclose(stdoutlog);

        if (core_backbone_array[array_inc] == i)
        {
            for (int h = 0; h < lengthpercore; h++)
            {
                

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "\t%i   %i\n", poly_lattice_indexes[array_inc *lengthpercore+h],latticepoint_S[poly_lattice_indexes[array_inc *lengthpercore+h]]);
                fclose(stdoutlog);

            }

            array_inc++;

        }
    }
    
    
    fclose(stdoutlog);
    // run to a valid energy for the energy window

    // as long as energy is outside the middle third of local energy range
    long long increm = 0;
    int ep = 0;
    int loc_energy_before = 0;
    int loc_energy_after = 0;
    e = -100;

    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)    //defines wl pseudo_chain for use in the below local energy calculation
    {
        wl_pseudo_chain[i] = poly_lattice_indexes[i];
        wl_pseudo_chain_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        wl_pseudo_chain_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        wl_pseudo_chain_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
    }
    
    int min_e = 0;
    int max_e = 0;

    if (1)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\n");
        fprintf(stdoutlog, "\n");
        fprintf(stdoutlog, "increment: %lli", increm);

        //stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\n");
        //            for (int i = 0; i < numberspins; i++)
        //            {
        //                fprintf(stdoutlog, "%4i", latticepoint[i]);
        //                if ((i + 1) % L1dim == 0)
        //                    fprintf(stdoutlog, "\n");
        //
        //            }

        fprintf(stdoutlog, "\n");
        fprintf(stdoutlog, "%i\n", e);
        fprintf(stdoutlog, "\n");
        for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
        {
            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t \t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2], ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
        }

        //
        fprintf(stdoutlog, "\n");
        for (int i = 0; i < 3*no_solvent_sites; i++)
        {
            fprintf(stdoutlog, "%4i\t\n", solvent_loc[i]);
        }

        fprintf(stdoutlog, "\n");
        fprintf(stdoutlog, "\n");
        for (int i = 0; i < no_ion_sites; i++)
        {
            fprintf(stdoutlog, "%4i\t%4i\n", ion_loc[i], ph_ion_loc[i]);
        }

        fclose(stdoutlog);

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\n");
        for (int i = 0; i < numberspins_S; i++)
        {
            fprintf(stdoutlog, "%4i", latticepoint_S[i]);
            if ((i + 1) % L1dim_S_xy == 0)
                fprintf(stdoutlog, "\n");
            if ((i + 1) % (L1dim_S_xy *L1dim_S_xy) == 0)
                fprintf(stdoutlog, "\n");

        }

        fprintf(stdoutlog, "\n");
        fprintf(stdoutlog, "\n");
        fclose(stdoutlog);
        stdoutlog = fopen(stdoutlogname, "a");
        for (int i = 0; i < numberspins_L; i++)
        {
            fprintf(stdoutlog, "%4i", latticepoint_L[i]);
            if ((i + 1) % L1dim_L_xy == 0)
                fprintf(stdoutlog, "\n");
            if ((i + 1) % (L1dim_L_xy *L1dim_L_xy) == 0)
                fprintf(stdoutlog, "\n");

        }

        fclose(stdoutlog);
    }

    //MPI_Abort(MPI_COMM_WORLD, 1);

    ep = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
    e = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "ajkcuc,\n\n\n\n\n");
    fclose(stdoutlog);
//
//    current_attach=0;
//    for(int i=0;i<numbercores*pertrusionpercore;i++)
//    {
//        for(int j=0;j<no_ion_sites;j++)
//        {
//            if(int_point_distance(poly_lattice_indexes[i*lengthpercore+lengthpercore-1],ion_loc[j])<=9)
//            {
//                current_attach++;
//                current_attachment[i*lengthpercore+lengthpercore-1]++;
//            }
//        }
//    }
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "bbbbaoubobuobuc,\n\n\n\n\n");
    fclose(stdoutlog);
//    for(int i=0;i<numbercores*pertrusionpercore;i++)
//    {
//        attach_track[ep]+=current_attachment[i*lengthpercore+lengthpercore-1];
//    }
//    //attach_track[ep] += current_attach;
//    attach_track_no[ep] ++;
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "aoubobuobuc,\n\n\n\n\n");
    fclose(stdoutlog);
    
//    while((e < (emin)) || (e > (emin +(emax - emin))))
//    {
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "sveye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i  %i",0,latticepoint_S[solvent_loc[0]],solvent_loc[0]);
////        fclose(stdoutlog);
//        //eye();
//        uni_counter=0;
//        loc_en_before=0;
//
//        confused2= 0;
//        //ep = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);    // compares the energy before and after a movement has been preformed
//
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "left q2\n");
////        fclose(stdoutlog);
//
//        //no_neg_electrostatic_inter[0]=0;
//        //no_pos_electrostatic_inter[0]=0;
//
//        no_OH_electrostatic_inter[0]=0;
//        no_OEu_electrostatic_inter[0]=0;
//        no_OO_electrostatic_inter[0]=0;
//        no_HH_electrostatic_inter[0]=0;
//        no_EuEu_electrostatic_inter[0]=0;
//        no_HEu_electrostatic_inter[0]=0;
//
//        poly_mov();
//
////        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "%i %i\n",uni_counter,latticepoint_S[80]);
////        fclose(stdoutlog);
//
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "Left polymov()  %i\n",latticepoint_L[0]);
////                fclose(stdoutlog);
//
//        loc_energy_before = 0;
//        loc_energy_after = 0;
//
//
//
//        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "f%i %i\n",uni_counter,latticepoint_L[0]);
////        fclose(stdoutlog);
//
//        confused2=0;
//
//        if (MoveProposal == 1)
//        {
//            loc_energy_before = loc_en_before;
////
////            if (reset_indexes[4] > 0)
////            {
////                if (reset_indexes[4] == 1)
////                {
////                    //                    stdoutlog = fopen(stdoutlogname, "a");
////                    //                    fprintf(stdoutlog, "dettach\n");
////                    //                    fclose(stdoutlog);
////                    reset_indexes[5] = random_solvent_site;
////                    reset_indexes[6] = attachment_ion_ref[2 *random_solvent_site + 1];
////
////                    detachment(random_solvent_site);
////
////                    reset_indexes[4] = 1;
////
////                    current_attach--;
////
////                }
////
////                if (reset_indexes[4] == 2)
////                {
////                    //                    stdoutlog = fopen(stdoutlogname, "a");
////                    //                    fprintf(stdoutlog, "attach tried\n");
////                    //                    fclose(stdoutlog);
////                    int success = 0;
////                    if (poly_ion_mov == 1)
////                    {
////                        success = attachment(random_solvent_site, random_index);
////                    }
////
////                    if (poly_ion_mov == 0 && poly_solvent_mov == 0)
////                    {
////                        reattachment(random_solvent_site, reset_indexes[0] + 1);
////                        success = 1;
////                    }
////
////                    //                    stdoutlog = fopen(stdoutlogname, "a");
////                    //                    fprintf(stdoutlog, "\nattach  %i  %i  %i  %i\t %i  %i f:%i \n",success,random_solvent_site,ph_ion_loc[random_solvent_site],ion_loc[random_solvent_site],ph_poly_lattice_indexes[attachment_ion_ref[2*random_solvent_site+1]],attachment_ion_ref[2*random_solvent_site+1],attachment_ref[2*attachment_ion_ref[2*random_solvent_site+1]+1]);
////                    //                    fclose(stdoutlog);
////
////                    if (success)
////                    {
////                        reset_indexes[4] = 2;
////                        reset_indexes[5] = random_solvent_site;
////                        reset_indexes[6] = attachment_ion_ref[2 *random_solvent_site + 1];
////
////                        current_attach++;
////                    }
////                    else
////                    {
////                        //                        stdoutlog = fopen(stdoutlogname, "a");
////                        //                        fprintf(stdoutlog, "attach failed\n");
////                        //                        fclose(stdoutlog);
////                        reset_indexes[4] = 0;
////                        MoveProposal = 0;
////                        if (poly_ion_mov == 1)
////                        {
////                            solvent_reset(ph_ion_loc, ion_loc, 3);
////                        }
////
////                        if (poly_ion_mov == 0 && poly_solvent_mov == 0)
////                        {
////                            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
////                        }
////                    }
////                }
////            }
//
//            if (poly_solvent_mov == 0 && poly_ion_mov == 0)
//            {
//                //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
////                uni_counter++;
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "jj%i %i  %i\n",uni_counter,latticepoint_S[80],poly_spin_reference[reset_indexes[0]+1]);
////                fclose(stdoutlog);
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "33f%i %i\n",uni_counter,latticepoint_L[0]);
////                fclose(stdoutlog);
//                lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
////                uni_counter++;
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "j%i %i\n",uni_counter,latticepoint_S[80]);
////                fclose(stdoutlog);
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "34f%i %i\n",uni_counter,latticepoint_L[0]);
////                fclose(stdoutlog);
//                poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "35f%i %i\n",uni_counter,latticepoint_L[0]);
////                fclose(stdoutlog);
//                uni_counter++;
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "j%i %i\n",uni_counter,latticepoint_S[80]);
////                fclose(stdoutlog);
//
//                reset_indexes[4] = 0;
//            }
//
//            if (poly_solvent_mov == 1)
//            {
//                solvent_reset2(solvent_loc, ph_solvent_loc, 2);
//
//                //                stdoutlog = fopen(stdoutlogname, "a");
//                //                fprintf(stdoutlog, "solvent_loc[i]: %i ph_solvent_loc[i]: %i\n",solvent_loc[reset_indexes[0]+1],ph_solvent_loc[reset_indexes[0]+1]);
//                //                fclose(stdoutlog);
//            }
//
//            if (poly_ion_mov == 1)
//            {
//                solvent_reset(ion_loc, ph_ion_loc, 3);
//                reset_indexes[4] = 0;
//                //                stdoutlog = fopen(stdoutlogname, "a");
//                //                fprintf(stdoutlog, "solvent_loc[i]: %i ph_solvent_loc[i]: %i\n",solvent_loc[reset_indexes[0]+1],ph_solvent_loc[reset_indexes[0]+1]);
//                //                fclose(stdoutlog);
//            }
//
//            vicinity_check = 0;
//            //no_neg_electrostatic_inter[1]=0;
//            //no_pos_electrostatic_inter[1]=0;
//
//            no_OH_electrostatic_inter[1]=0;
//            no_OEu_electrostatic_inter[1]=0;
//            no_OO_electrostatic_inter[1]=0;
//            no_HH_electrostatic_inter[1]=0;
//            no_EuEu_electrostatic_inter[1]=0;
//            no_HEu_electrostatic_inter[1]=0;
//
//            confused= 0;
//            loc_energy_after = local_energy(ph_poly_lattice_coordinates, ph_poly_lattice_indexes, ph_solvent_loc, ph_ion_loc,1);
//
//        }
//        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "f%i %i\n",uni_counter,latticepoint_L[0]);
////        fclose(stdoutlog);
//        increm = increm + 1;
//
//        e += loc_energy_after-loc_energy_before;
//       // uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "2f%i %i\n",uni_counter,latticepoint_L[0]);
////        fclose(stdoutlog);
////                stdoutlog = fopen(stdoutlogname, "a");
////                fprintf(stdoutlog, "Left totalenergyx2 increm: %lli\n   ep: %i   e:%i   ep-e: %i     %i\n",increm%1,ep,e,ep-e, MoveProposal);
////                fprintf(stdoutlog, "local_energy_before: %i   local_energy_after: %i   loc_en_b-loc_en_a: %i\n",loc_energy_before,loc_energy_after,loc_energy_before-loc_energy_after);
////                fclose(stdoutlog);
////
//
//
////        if (loc_energy_before - loc_energy_after != ep - e)
////        {
////
////            stdoutlog = fopen(stdoutlogname, "a");
////            fprintf(stdoutlog, "Difference %i\n",MoveProposal);
////            fclose(stdoutlog);
////            MPI_Abort(MPI_COMM_WORLD, 1);
////        }
//
//        //        if(min_e>ep)
//        //        {
//        //            min_e=ep;
//        //            stdoutlog = fopen(stdoutlogname, "a");
//        //            fprintf(stdoutlog, "min_e: %i\n",ep);
//        //            fclose(stdoutlog)
//        //        }
//
//        //
//        //        if(max_e < ep)
//        //        {
//        //            max_e=ep;
//        //            stdoutlog = fopen(stdoutlogname, "a");
//        //            fprintf(stdoutlog, "max_e: %i\n",ep);
//        //            fclose(stdoutlog);
//        //        }
//
////        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "a%i %i\n",uni_counter,latticepoint_S[80]);
////        fclose(stdoutlog);
//
////        if(vicinity_check)
////        {
////            for(int i=0;i<numbercores*pertrusionpercore;i++)
////            {
////                attach_track[e]+=current_attachment[i*lengthpercore+lengthpercore-1];
////            }
////
////            //attach_track[ep] += current_attach;
////            attach_track_no[e] ++;
////        }
////        current_attach=0;
//
////        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "j%i %i\n",uni_counter,latticepoint_S[80]);
////        fclose(stdoutlog);
////
////        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "3f%i %i\n",uni_counter,latticepoint_L[0]);
////        fclose(stdoutlog);
//        if ((increm % 10000) == 0)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "\n");
//            fprintf(stdoutlog, "\n");
//            fprintf(stdoutlog, "increment: %lli", increm);
//
//            //stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "\n");
//            //            for (int i = 0; i < numberspins; i++)
//            //            {
//            //                fprintf(stdoutlog, "%4i", latticepoint[i]);
//            //                if ((i + 1) % L1dim == 0)
//            //                    fprintf(stdoutlog, "\n");
//            //
//            //            }
//
//            fprintf(stdoutlog, "\n");
//            fprintf(stdoutlog, "%i\n", e);
//            fprintf(stdoutlog, "\n");
////            for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
////            {
////                fprintf(stdoutlog, "%i: %4i\t\t%4i\t%4i\t%4i\t\t%4i\t\t%4i\t%4i\t%4i\t\t%i\t%i\t%i\n",i, poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2], ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2], poly_spin_reference[i], attachment_ref[2 *i], attachment_ref[2 *i + 1]);
////            }
////
////            //
////            fprintf(stdoutlog, "\n");
////            for (int i = 0; i < 2*no_solvent_sites; i++)
////            {
////                fprintf(stdoutlog, "%4i\t%i\n", solvent_loc[i], ph_solvent_loc[i]);
////            }
////
////            fprintf(stdoutlog, "\n");
////            fprintf(stdoutlog, "\n");
////            for (int i = 0; i < no_ion_sites; i++)
////            {
////                fprintf(stdoutlog, "%4i\t%4i\t %i\t%i\n", ion_loc[i], ph_ion_loc[i], attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
////            }
////
////            fclose(stdoutlog);
////
////            stdoutlog = fopen(stdoutlogname, "a");
////            fprintf(stdoutlog, "\n");
////            for (int i = 0; i < numberspins_S; i++)
////            {
////                fprintf(stdoutlog, "%4i", latticepoint_S[i]);
////                if ((i + 1) % L1dim_S_xy == 0)
////                    fprintf(stdoutlog, "\n");
////                if ((i + 1) % (L1dim_S_xy *L1dim_S_xy) == 0)
////                    fprintf(stdoutlog, "\n");
////
////            }
////
////            fprintf(stdoutlog, "\n");
////            fprintf(stdoutlog, "\n");
////            fclose(stdoutlog);
////            stdoutlog = fopen(stdoutlogname, "a");
////            for (int i = 0; i < numberspins_L; i++)
////            {
////                fprintf(stdoutlog, "%4i", latticepoint_L[i]);
////                if ((i + 1) % L1dim_L_xy == 0)
////                    fprintf(stdoutlog, "\n");
////                if ((i + 1) % (L1dim_L_xy *L1dim_L_xy) == 0)
////                    fprintf(stdoutlog, "\n");
////
////            }
//
//            fclose(stdoutlog);
//        }
//
//
////        uni_counter++;
////        stdoutlog = fopen(stdoutlogname, "a");
////        fprintf(stdoutlog, "5f%i %i\n",uni_counter,latticepoint_L[0]);
////        fclose(stdoutlog);
//        //        ep = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc,ion_loc);    // compares the energy before and after a movement has been preformed
//        //        e = total_energy(ph_poly_lattice_coordinates,ph_poly_lattice_indexes,ph_solvent_loc,ph_ion_loc);
//        //        stdoutlog = fopen(stdoutlogname, "a");
//        //        fprintf(stdoutlog, "\n\n\nLeft totalenergyx2 increm: %lli\n   ep: %i   e:%i   ep-e: %i     %i\n",increm%1,ep,e,ep-e, MoveProposal);
//        //        fclose(stdoutlog);
//    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "passed range finder\n");
    fclose(stdoutlog);

    // use last spin to store replica-id; this is just a marker to
    // allows us to keep track of replicas as they are exchanged
    latticepoint_S[numberspins_S + 2] = myid;
    latticepoint_L[numberspins_L + 2] = myid;

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned lattice points\n");
    fclose(stdoutlog);

    if (fopen("Access_Levels.txt", "r") == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
        fclose(stdoutlog);
    }
    else    // rewrite of prior file reading system to be more portable
    {
        FILE * fptr;
        int value;
        int valuee;

        fptr = fopen("Access_Levels.txt", "r");

        while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
        {
            if ((e >= (emin) && (e <= (emin + (emax - emin)))))
            {
                lngE[value-Eglobalmin] = valuee;
                pseudolngE[value-Eglobalmin] = valuee;
            }
        }

        fclose(fptr);
    }

    // Pseudo WL Process is started for all processes to explore the energy levels
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "\nPseudo WL Proces has started for process: %i\n", myid);
    fclose(stdoutlog);
    pseudowl();
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "\nLeft Pseudo WL Proces has started for process: %i\n", myid);
    fclose(stdoutlog);
    

    
    stdoutlog = fopen(stdoutlogname, "a");

    fprintf(stdoutlog, "local dup group: %i  localgroup: %i\n", local_dup_group, localgroup);

//    if (local_dup_group == localgroup && merge_hists == 1)    // merge g(E) estimators from multiple walkers in the same energy window
//    {
//        fprintf(stdoutlog, "ugig");
//        if (myid == localgroup)    // 'root' in energy window, receive individual g(E) and send merged g(E)
//        {
//            fprintf(stdoutlog, "upig");
//            for (int i = (dup_headproc + 1); i < numprocs; i++)
//            {
//                MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, i, 76, MPI_COMM_WORLD, &status);    // get other dens. of states
//                fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, i);
//                for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j];    // sum up for average
//            }
//
//            for (int i = (dup_headproc + 1); i < numprocs; i++)
//            {
//                MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, i, 98, MPI_COMM_WORLD);
//                fprintf(stdoutlog, "Proc %i: Sent combined lngE to Proc. %i\n", myid, i);
//            }
//        }
//        else    // send individual g(E) and receive merged g(E)
//        {
//            fprintf(stdoutlog, "ufig");
//            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, dup_headproc, 76, MPI_COMM_WORLD);
//            fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, dup_headproc);
//            MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, dup_headproc, 98, MPI_COMM_WORLD, &status);
//            fprintf(stdoutlog, "Proc %i: Received combined lngE from Proc. %i\n", myid, dup_headproc);
//            for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];    // replace individual lngE (could be done directly, yes)
//        }
//    }
//
//    fclose(stdoutlog);

    if (multiple > 1 && merge_hists == 1)    // merge g(E) estimators from multiple walkers in the same energy window
    {
        stdoutlog = fopen(stdoutlogname, "a");
        if (myid % multiple == 0)    // 'root' in energy window, receive individual g(E) and send merged g(E)
        {
            for (int i = 1; i < multiple; i++)
            {
                MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid + i, 76, MPI_COMM_WORLD, &status);    // get other dens. of states
                fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, myid + i);
                for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j];    // sum up for average
            }

            for (int i = 1; i < multiple; i++)
            {
                MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid + i, 98, MPI_COMM_WORLD);
                fprintf(stdoutlog, "Proc %i: Sent combined lngE to Proc. %i\n", myid, myid + i);
            }
        }
        else    // send individual g(E) and receive merged g(E)
        {
            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 76, MPI_COMM_WORLD);
            fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, myid - (myid % multiple));
            MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 98, MPI_COMM_WORLD, &status);
            fprintf(stdoutlog, "Proc %i: Received combined lngE from Proc. %i\n", myid, myid - (myid % multiple));
            for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];    // replace individual lngE (could be done directly, yes)
        }

        fclose(stdoutlog);
    }

    accessiblelevels();

    //Print lattice to logfile
    stdoutlog = fopen(stdoutlogname, "a");
    for (int i = 0; i < numberspins_S; i++)
    {
        fprintf(stdoutlog, "%4i", latticepoint_S[i]);
        if ((i + 1) % L1dim_S_xy == 0)
            fprintf(stdoutlog, "\n");
    }

    fprintf(stdoutlog, "\n");
    //fprintf(stdoutlog, "%i", total_energy());
    fprintf(stdoutlog, " done \n");
    fclose(stdoutlog);

    //MPI_Abort(MPI_COMM_WORLD, 1);
}

void init_hists()    // initialize histograms
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to init hist()\n");
    fclose(stdoutlog);

    lngE = (double*) malloc(hist_size* sizeof(double));
    lngE_buf = (double*) malloc(hist_size* sizeof(double));
    real_lngE = (double*) malloc(hist_size* sizeof(double));
    real_lngE_buf = (double*) malloc(hist_size* sizeof(double));
    microT = (double*) malloc(hist_size* sizeof(double));
    microT_buf = (double*) malloc(hist_size* sizeof(double));
    pseudolngE = (double*) malloc(hist_size* sizeof(double));    // added for pseudo wl process
    HE = (double*) malloc(hist_size* sizeof(double));

    
    for (int i = 0; i < hist_size; i++)
    {
        lngE[i] = 0.0;
        lngE_buf[i] = 0.0;
        real_lngE[i] = 0.0;
        real_lngE_buf[i] = 0.0;
        microT[i] = 0.0;
        pseudolngE[i] = 0.0;
        microT_buf[i] = 0.0;
    }
}

void partial_init_hists()    // initialize histograms
{
    for (int i = 0; i < hist_size; i++)
    {
        lngE_buf[i] = 0.0;
        real_lngE[i] = 0.0;
        real_lngE_buf[i] = 0.0;
        microT[i] = 0.0;
        microT_buf[i] = 0.0;
    }
}

// set boundaries for energy windows:
// either read from file (if precomputed for balanced REWL)
// or just split energy range evenly according to avaliable waklers
int find_local_energy_range(double Eglobmin, double Eglobmax, double overlap, int N)
{
    stdoutlog = fopen(stdoutlogname, "a");

    // consistency check
    if (N % multiple != 0)
    {
        fprintf(stdoutlog, "Total number of processes (%i) must be a multiple of the number of walkers per energy range (%i)!\n", N, multiple);
        fclose(stdoutlog);
        return (1);
    }

    FILE * window_init;

    // check if there is a file containing local energy window boundaries
    // Must be named Ewindows.dat !

    //fprintf(stdoutlog, "\nProc %i: Can't find file Ewindows.dat. Will calculate equal-size windows with overlap %lf\n", myid, overlap);
    //double Ewidth = (Eglobmax - Eglobmin) / (1.0 + ((double)(N / multiple) - 1.0) *(1.0 - overlap));

    //Emin = Eglobmin + (double)(myid / multiple) *(1.0 - overlap) *Ewidth;
    //Emax = Emin + Ewidth;

    //Eminindex = floor(Emin + -Eglobalmin);
    //Emaxindex = ceil(Emax + -Eglobalmin);

    if (fopen("Ewindows_991.dat", "r") == NULL)
    {
        fprintf(stdoutlog, "\nProc %i: Can't find file Ewindows.dat. Will calculate equal-size windows with overlap %lf\n", myid, overlap);
        double Ewidth = (Eglobmax - Eglobmin) / (1.0 + ((double)(N / multiple) - 1.0) *(1.0 - overlap));

        Emin = Eglobmin + (double)(myid / multiple) *(1.0 - overlap) *Ewidth;
        Emax = Emin + Ewidth;

        Eminindex = floor(Emin + -Eglobalmin);
        Emaxindex = ceil(Emax + -Eglobalmin);

        time(&timenow);
        fprintf(stdoutlog, "Proc %3i: Parameter: Eglobmin -- Eglobmax: %lf -- %lf; overlap=%i percent, Ewindowwidth=%lf, %s", myid, Eglobmin, Eglobmax, (int)(overlap *100.0), Ewidth, ctime(&timenow));
        //fprintf(stdoutlog, "Proc %3i: Parameter: Eglobmin -- Eglobmax: %lf -- %lf; %s", myid, Emin, Emax, ctime(&timenow));
    }
    else    // rewrite of prior file reading system to be more portable
    {
        FILE * fptr;
        int value;
        int valuee;
        int valueee;
        int valueeee;

        fptr = fopen("Ewindows_991.dat", "r");

        while (fscanf(fptr, "%i,%i,%i,%i\n", &value, &valuee, &valueee, &valueeee) > 0)
        {
            if (value == myid / multiple)
            {
                Eminindex = valuee-Eglobalmin;
                Emaxindex = valueee-Eglobalmin;
                localgroup = valueeee;
                fprintf(stdoutlog, "%i %i %i %i\n", value, Eminindex, Emaxindex, localgroup);
            }

            Emin =Eglobalmin + Eminindex;
            Emax =Eglobalmin + Emaxindex;
        }

        fclose(fptr);
    }

    fclose(stdoutlog);

    return (0);
}

// THIS IS THE MASTER RE / SWAP FUNCTION
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int replica_exchange(int swap_direction, int index_akt)
{
    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "1\n");
    //    fclose(stdoutlog);

    int i_new;    // histogram index of my configuration
    int Ecur;    // current energy

    // frac refers to local exchange probability
    // wk is combined exchange probability
    double myfrac, otherfrac, randx, wk;

    int change = 0;    // boolean: 0 for not exchanging, 1 for exchanging
    int swap_partner = -1;    // id of swap partner (receive_buffer)

    swap_direction = swap_direction % 2;    // comes actually as number of swap attempt

    // everyone has to find its swap-partner

    int *pairs;    // array containing the partners of each process (send_buffer)
    pairs = (int*) malloc(2 *multiple* sizeof(int));

    if (mylocalid[swap_direction] == 0)    // 'head-node' in the energy window determines pairs of flippartners
    {
        //stdoutlog = fopen(stdoutlogname, "a");
        //fprintf(stdoutlog, "2\n");
        //fclose(stdoutlog);

        int choose_from = multiple;    // number of free partners in higher window of communicator
        int select;    // storage for random number

        int *libre;    // list of free partners from higher window in communicator
        libre = (int*) malloc(multiple* sizeof(int));

        for (int i = 0; i < multiple; i++) libre[i] = multiple + i;    // initialise

        // idea: processes from the lower window choose someone from the higher window at random
        // of course, the chosen walker can't have an exchange partner yet
        for (int i = 0; i < multiple; i++)    // loop over processes in the lower window
        {
            select = rand() % choose_from;
            pairs[i] = libre[select];
            pairs[libre[select]] = i;    // the 'vice-versa pair'
            // update list
            choose_from--;
            for (int j = select; j < choose_from; j++)
                libre[j] = libre[j + 1];
        }

        //       stdoutlog=fopen(stdoutlogname,"a");
        //       fprintf(stdoutlog,"Proc %3i: Drew the following swap partners:\n",myid);
        //       for (int i=0;i<2*multiple;i++)
        //     fprintf(stdoutlog,"Proc %3i: %i -- %i (local ids in communicator)\n",myid,i,pairs[i]);
        //       fclose(stdoutlog);

        free(libre);
    }

    // at this point, every walker has a swap partner assigned, now they must be communicated
    if ((swap_direction == 0) && (myid < (numprocs - multiple)))    // the walkers from the last node should not swap
    {
        comm_id = 2 *(myid / (2 *multiple));    // ! all integer, the '/' is a (div) ! Not the same as myid/multiple !
        MPI_Scatter(pairs, 1, MPI_INT, &swap_partner, 1, MPI_INT, 0, mpi_local_comm[comm_id]);
    }

    if ((swap_direction == 1) && (myid >= multiple))    // the walkers from the zero-node should not swap
    {
        comm_id = ((myid - multiple) / (2 *multiple)) *2 + 1;    // ! all integer, the '/' is a (div) ! See above
        MPI_Scatter(pairs, 1, MPI_INT, &swap_partner, 1, MPI_INT, 0, mpi_local_comm[comm_id]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(pairs);

    if (swap_partner != -1)    // i.e. if there is a swap-partner for me (if I am at a boundary, I might not have a swap partner this time)
    {
        //stdoutlog = fopen(stdoutlogname, "a");
        //fprintf(stdoutlog, "3\n");
        //fclose(stdoutlog);

        // statistics
        if (swap_partner > mylocalid[swap_direction]) tryright++;
        else tryleft++;

        // safety cross check
        Ecur = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
        if (Ecur + (-Eglobalmin) != index_akt)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Proc %3i, replica_exchange(): Something is wrong here! Received index=%i, calculated index=%i totalenergy=%i. Abort.\n", myid, index_akt, Ecur + (-Eglobalmin), total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc));
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        i_new = index_akt;

        // get histogram index from my swap partner
        MPI_Sendrecv_replace(&i_new, 1, MPI_INT, swap_partner, 1, swap_partner, 1, mpi_local_comm[comm_id], &status);

        if ((i_new > Emaxindex) || (i_new < Eminindex))    // energyranges must overlap!
        {
            myfrac = -1.0;
        }
        else
        {
            // calculate my part of the exchange probability
            myfrac = exp(lngE[index_akt] - lngE[i_new]);    // g(myE)/g(otherE)
        }

        if (mylocalid[swap_direction] < multiple)    // I am receiver and calculator
        {
            // get my partners part of the exchange probability
            MPI_Recv(&otherfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_local_comm[comm_id], &status);

            // calculate combined exchange probability and roll the dice
            if ((myfrac > 0.0) && (otherfrac > 0.0))
            {
                randx = (1.0* rand() / (RAND_MAX + 1.0));
                wk = myfrac * otherfrac;
                if (randx < wk) change = 1;
            }

            // tell my swap partner whether to exchange or not
            MPI_Send(&change, 1, MPI_INT, swap_partner, 3, mpi_local_comm[comm_id]);
        }
        else    // I just send my part of exchange probability and await decision
        {
            MPI_Send(&myfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_local_comm[comm_id]);
            MPI_Recv(&change, 1, MPI_INT, swap_partner, 3, mpi_local_comm[comm_id], &status);
        }

        // if decision was made to exchange configurations
        if (change == 1)
        {
            //stdoutlog = fopen(stdoutlogname, "a");
            //fprintf(stdoutlog, "4\n");
            //fclose(stdoutlog);
            // exchange spin conformations (incl. the 3 'special' lattice sites)
            MPI_Sendrecv_replace(&poly_lattice_indexes[0], numbercores *lengthpercore *pertrusionpercore + backbone_length, MPI_INT, swap_partner, 20, swap_partner, 20, mpi_local_comm[comm_id], &status);
            MPI_Sendrecv_replace(&poly_lattice_coordinates[0], (numbercores *lengthpercore *pertrusionpercore + backbone_length) *3, MPI_INT, swap_partner, 10, swap_partner, 10, mpi_local_comm[comm_id], &status);
            MPI_Sendrecv_replace(&solvent_loc[0], 3*no_solvent_sites, MPI_INT, swap_partner, 30, swap_partner, 30, mpi_local_comm[comm_id], &status);

            MPI_Sendrecv_replace(&ion_loc[0], no_ion_sites, MPI_INT, swap_partner, 40, swap_partner, 40, mpi_local_comm[comm_id], &status);

            //            stdoutlog = fopen(stdoutlogname, "a");
            //            fprintf(stdoutlog, "5\n");
            //            fclose(stdoutlog);

            switched = i_new;

            // statistics
            if (swap_partner > mylocalid[swap_direction]) exchangeright++;
            else exchangeleft++;
        }
    }

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "6\n");
    //fclose(stdoutlog);
    return (change);
    // returns whether or not configs were actually exchanged
}

void sysparam(int index)
{
    
  radial_counter(poly_lattice_indexes, solvent_loc, ion_loc,index);
    
    double rsum[numbercores *pertrusionpercore] = { 0 };
    double radius = 0;

    
    
    int dx[numbercores *pertrusionpercore *lengthpercore] = { 0 };
    int dy[numbercores *pertrusionpercore *lengthpercore] = { 0 };
    int dz[numbercores *pertrusionpercore *lengthpercore] = { 0 };

    int s_dx[numbercores *pertrusionpercore *lengthpercore] = { 0 };
    int s_dy[numbercores *pertrusionpercore *lengthpercore] = { 0 };
    int s_dz[numbercores *pertrusionpercore *lengthpercore] = { 0 };

    int dxbb[backbone_length+numbercores] = { 0 };
    int dybb[backbone_length+numbercores] = { 0 };
    int dzbb[backbone_length+numbercores] = { 0 };

    int s_dxbb[backbone_length+numbercores] = { 0 };
    int s_dybb[backbone_length+numbercores] = { 0 };
    int s_dzbb[backbone_length+numbercores] = { 0 };

    for (int i = (0); i < numbercores * pertrusionpercore; i++)    // performs the rotations
    {
        double centerofmass[3] = { 0 };

        s_dx[i *lengthpercore] = poly_lattice_coordinates[3 *(i *lengthpercore)];
        s_dy[i *lengthpercore] = poly_lattice_coordinates[3 *(i *lengthpercore) + 1];
        s_dz[i *lengthpercore] = poly_lattice_coordinates[3 *(i *lengthpercore) + 2];

        centerofmass[0] += s_dx[i *lengthpercore];
        centerofmass[1] += s_dy[i *lengthpercore];
        centerofmass[2] += s_dz[i *lengthpercore];

        for (int j = 1; j < lengthpercore; j++)    // performs the rotations
        {
            dx[i *lengthpercore + j] = poly_lattice_coordinates[3 *(i *lengthpercore + j)] - poly_lattice_coordinates[3 *(i *lengthpercore + j - 1)];

            dy[i *lengthpercore + j] = poly_lattice_coordinates[3 *(i *lengthpercore + j) + 1] - poly_lattice_coordinates[3 *(i *lengthpercore + j - 1) + 1];

            dz[i *lengthpercore + j] = poly_lattice_coordinates[3 *(i *lengthpercore + j) + 2] - poly_lattice_coordinates[3 *(i *lengthpercore + j - 1) + 2];

            //            stdoutlog = fopen(stdoutlogname, "a");
            //            fprintf(stdoutlog, "d={%i,%i,%i}    d-1={%i,%i,%i}\n\tp_l_c_diff={%i,%i,%i}\n",dx[i+j],dy[i+j],dz[i+j],dx[i+j-1],dy[i+j-1],dz[i+j-1],poly_lattice_coordinates[3 *(i+j)]-poly_lattice_coordinates[3 *(i+j-1)] ,poly_lattice_coordinates[3 *(i+j)+1] -poly_lattice_coordinates[3 *(i+j-1) + 1],poly_lattice_coordinates[3 *(i+j)+2]-poly_lattice_coordinates[3 *(i+j-1)+2]);
            //            fclose(stdoutlog);

            dx[i *lengthpercore + j] = dx[i *lengthpercore + j - 1] - (dx[i *lengthpercore + j] - anint(dx[i *lengthpercore + j], L1dim_S_xy) *L1dim_S_xy);
            dy[i *lengthpercore + j] = dy[i *lengthpercore + j - 1] - (dy[i *lengthpercore + j] - anint(dy[i *lengthpercore + j], L1dim_S_xy) *L1dim_S_xy);
            dz[i *lengthpercore + j] = dz[i *lengthpercore + j - 1] - (dz[i *lengthpercore + j] - anint(dz[i *lengthpercore + j], L1dim_S_z) *L1dim_S_z);

            s_dx[i *lengthpercore + j] += s_dx[i *lengthpercore] + dx[i *lengthpercore + j];
            s_dy[i *lengthpercore + j] += s_dy[i *lengthpercore] + dy[i *lengthpercore + j];
            s_dz[i *lengthpercore + j] += s_dz[i *lengthpercore] + dz[i *lengthpercore + j];

            centerofmass[0] += s_dx[i *lengthpercore + j];
            centerofmass[1] += s_dy[i *lengthpercore + j];
            centerofmass[2] += s_dz[i *lengthpercore + j];

            //            stdoutlog = fopen(stdoutlogname, "a");
            //            fprintf(stdoutlog, "d={%i,%i,%i}  cm={%f,%f,%f}\n\tp_l_c[i+j]={%i,%i,%i}     p_l_c[i+j-1]={%i,%i,%i}\n",dx[i+j],dy[i+j],dz[i+j],centerofmass[0],centerofmass[1],centerofmass[2],poly_lattice_coordinates[3 *(i+j)] ,poly_lattice_coordinates[3 *(i+j)+1] ,poly_lattice_coordinates[3 *(i+j)+2],poly_lattice_coordinates[3 *(i+j-1)],poly_lattice_coordinates[3 *(i+j-1) + 1],poly_lattice_coordinates[3 *(i+j-1)+2]);
            //            fclose(stdoutlog);
        }

        centerofmass[0] /= lengthpercore;
        centerofmass[1] /= lengthpercore;
        centerofmass[2] /= lengthpercore;

        for (int j = 0; j < lengthpercore; j++)    // performs the rotations
        {
            rsum[i] += ((double) s_dx[i *lengthpercore + j] - centerofmass[0]) *((double) s_dx[i *lengthpercore + j] - centerofmass[0]) + ((double) s_dy[i *lengthpercore + j] - centerofmass[1]) *((double) dy[i *lengthpercore + j] - centerofmass[1]) + ((double) s_dz[i *lengthpercore + j] - centerofmass[2]) *((double) s_dz[i *lengthpercore + j] - centerofmass[2]);

            //            stdoutlog = fopen(stdoutlogname, "a");
            //            fprintf(stdoutlog, "b%i rsum: %f  %f    %f    %f\t",j,rsum[i],((double)s_dx[i*lengthpercore+j] - centerofmass[0]) *((double)s_dx[i*lengthpercore+j] - centerofmass[0]),((double)s_dy[i*lengthpercore+j] - centerofmass[1]) *((double)dy[i*lengthpercore+j] - centerofmass[1]),((double)s_dz[i*lengthpercore+j] - centerofmass[2]) *((double)s_dz[i*lengthpercore+j] - centerofmass[2]));
            //            fclose(stdoutlog);
        }
    }

    for (int i = 0; i < numbercores * pertrusionpercore; i++)
    {
        //        stdoutlog = fopen(stdoutlogname, "a");
        //        fprintf(stdoutlog, "\nrsum: %f   %i\t\n",rsum[i],lengthpercore);
        //        fclose(stdoutlog);
        radius += sqrt(rsum[i] / lengthpercore);
    }

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "radius: %f\n",radius);
    //    fclose(stdoutlog);
    //
    rog[index] += radius / (numbercores *pertrusionpercore);

    //    stdoutlog = fopen(stdoutlogname, "a");
    //    fprintf(stdoutlog, "rog: %f\n",rog[index]);
    //    fclose(stdoutlog);
    double radiusbb = 0;

    int func_counter = 0;
    int func_counter2=0;
    
    
    for (int jj = 0; jj < numberofchains; jj++)    // performs the rotations
    {
        int func_counter3=0;
        double rsumbb = 0;
            
            double centerofmass[3] = { 0 };
        
            s_dxbb[jj*lengthperchain] = poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain)];
            s_dybb[jj*lengthperchain] = poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain) + 1];
            s_dzbb[jj*lengthperchain] = poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain) + 2];
        
            centerofmass[0] += s_dxbb[jj*lengthperchain];
            centerofmass[1] += s_dybb[jj*lengthperchain];
            centerofmass[2] += s_dzbb[jj*lengthperchain];
        
        if(poly_spin_reference[(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain)]!=1)
        {
            int ref=poly_spin_reference[(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain)]-5;
            
            dxbb[backbone_length+func_counter] = poly_lattice_coordinates[3 *(ref*lengthpercore+1)] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain)];

            dybb[backbone_length+func_counter] = poly_lattice_coordinates[3 *(ref*lengthpercore+1 )+ 1] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore + jj*lengthperchain) + 1];

            dzbb[backbone_length+func_counter] = poly_lattice_coordinates[3 * (ref*lengthpercore+1)+ 2] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore + jj*lengthperchain) + 2];
            
            s_dxbb[backbone_length+func_counter] += s_dxbb[jj*lengthperchain] + dxbb[backbone_length+func_counter];
            s_dybb[backbone_length+func_counter] += s_dybb[jj*lengthperchain] + dybb[backbone_length+func_counter];
            s_dzbb[backbone_length+func_counter] += s_dzbb[jj*lengthperchain] + dzbb[backbone_length+func_counter];
            
            centerofmass[0] += s_dxbb[backbone_length+func_counter];
            centerofmass[1] += s_dybb[backbone_length+func_counter];
            centerofmass[2] += s_dzbb[backbone_length+func_counter];
            
//            int j=0;
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "1d={%i,%i,%i}  cm={%f,%f,%f}   s_d={%i,%i,%i}\n\tp_l_c[i+j]={%i,%i,%i}     p_l_c[i+j-1]={%i,%i,%i}  %i   %i\n\n",dxbb[backbone_length+func_counter],dybb[backbone_length+func_counter],dzbb[backbone_length+func_counter],centerofmass[0],centerofmass[1],centerofmass[2],s_dxbb[backbone_length+func_counter],s_dybb[backbone_length+func_counter],s_dzbb[backbone_length+func_counter],poly_lattice_coordinates[3 *(ref*lengthpercore+1)] ,poly_lattice_coordinates[3 *(ref*lengthpercore+1)+1] ,poly_lattice_coordinates[3 *(ref*lengthpercore+1)+2],poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain)],poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain) + 1],poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain)+2],ref,func_counter);
//            fclose(stdoutlog);
            
             func_counter++;
        }
        
            for (int j = 1; j < lengthperchain; j++)    // performs the rotations
            {
                dxbb[jj*lengthperchain+j] = poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j)] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j - 1)];

                dybb[jj*lengthperchain+j] = poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j) + 1] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore + jj*lengthperchain + j - 1) + 1];

                dzbb[jj*lengthperchain+j] = poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j) + 2] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore + jj*lengthperchain + j - 1) + 2];

                //            stdoutlog = fopen(stdoutlogname, "a");
                //            fprintf(stdoutlog, "d={%i,%i,%i}    d-1={%i,%i,%i}\n\tp_l_c_diff={%i,%i,%i}\n",dx[i+j],dy[i+j],dz[i+j],dx[i+j-1],dy[i+j-1],dz[i+j-1],poly_lattice_coordinates[3 *(i+j)]-poly_lattice_coordinates[3 *(i+j-1)] ,poly_lattice_coordinates[3 *(i+j)+1] -poly_lattice_coordinates[3 *(i+j-1) + 1],poly_lattice_coordinates[3 *(i+j)+2]-poly_lattice_coordinates[3 *(i+j-1)+2]);
                //            fclose(stdoutlog);

                dxbb[jj*lengthperchain+j] = dxbb[jj*lengthperchain+j - 1] - (dxbb[jj*lengthperchain+j] - anint(dxbb[jj*lengthperchain+j], L1dim_S_xy) *L1dim_S_xy);
                dybb[jj*lengthperchain+j] = dybb[jj*lengthperchain+j - 1] - (dybb[jj*lengthperchain+j] - anint(dybb[jj*lengthperchain+j], L1dim_S_xy) *L1dim_S_xy);
                dzbb[jj*lengthperchain+j] = dzbb[jj*lengthperchain+j - 1] - (dzbb[jj*lengthperchain+j] - anint(dzbb[jj*lengthperchain+j], L1dim_S_z) *L1dim_S_z);

                s_dxbb[jj*lengthperchain+j] += s_dxbb[jj*lengthperchain] + dxbb[jj*lengthperchain+j];
                s_dybb[jj*lengthperchain+j] += s_dybb[jj*lengthperchain] + dybb[jj*lengthperchain+j];
                s_dzbb[jj*lengthperchain+j] += s_dzbb[jj*lengthperchain] + dzbb[jj*lengthperchain+j];

                centerofmass[0] += s_dxbb[jj*lengthperchain+j];
                centerofmass[1] += s_dybb[jj*lengthperchain+j];
                centerofmass[2] += s_dzbb[jj*lengthperchain+j];

                if(poly_spin_reference[(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain+j)]!=1)
                {
                    int ref=poly_spin_reference[(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain+j)]-5;
                    
                    dxbb[backbone_length+func_counter] = poly_lattice_coordinates[3 *(ref*lengthpercore+1)] -  poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j)];
                    
                    dybb[backbone_length+func_counter] = poly_lattice_coordinates[3 *(ref*lengthpercore+1) + 1] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j)+1];

                    dzbb[backbone_length+func_counter] = poly_lattice_coordinates[3 * (ref*lengthpercore+1)+ 2] - poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j)+2];
                    
                    s_dxbb[backbone_length+func_counter] += s_dxbb[jj*lengthperchain+j] + dxbb[backbone_length+func_counter];
                    s_dybb[backbone_length+func_counter] += s_dybb[jj*lengthperchain+j] + dybb[backbone_length+func_counter];
                    s_dzbb[backbone_length+func_counter] += s_dzbb[jj*lengthperchain+j] + dzbb[backbone_length+func_counter];
                    
                    centerofmass[0] += s_dxbb[backbone_length+func_counter];
                    centerofmass[1] += s_dybb[backbone_length+func_counter];
                    centerofmass[2] += s_dzbb[backbone_length+func_counter];
                    
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "d={%i,%i,%i}  cm={%f,%f,%f}   s_d={%i,%i,%i}\n\tp_l_c[i+j]={%i,%i,%i}     p_l_c[i+j-1]={%i,%i,%i}  %i   %i\n\n",dxbb[backbone_length+func_counter],dybb[backbone_length+func_counter],dzbb[backbone_length+func_counter],centerofmass[0],centerofmass[1],centerofmass[2],s_dxbb[backbone_length+func_counter],s_dybb[backbone_length+func_counter],s_dzbb[backbone_length+func_counter],poly_lattice_coordinates[3 *(ref*lengthpercore+1)] ,poly_lattice_coordinates[3 *(ref*lengthpercore+1)+1] ,poly_lattice_coordinates[3 *(ref*lengthpercore+1)+2],poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j)],poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j) + 1],poly_lattice_coordinates[3 *(numbercores *pertrusionpercore *lengthpercore +jj*lengthperchain + j)+2],ref,func_counter);
//                    fclose(stdoutlog);
                    
                    func_counter++;
                }
                
                            
            }

            centerofmass[0] /= (lengthperchain+func_counter);
            centerofmass[1] /= (lengthperchain+func_counter);
            centerofmass[2] /=(lengthperchain+func_counter);

//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "cm={%f,%f,%f}\n",centerofmass[0],centerofmass[1],centerofmass[2]);
//            fclose(stdoutlog);
        
            for (int j = 0; j < lengthperchain; j++)    // performs the rotations
            {
                rsumbb += (s_dxbb[jj*lengthperchain+j] - centerofmass[0]) *(s_dxbb[jj*lengthperchain+j] - centerofmass[0]) + (s_dybb[jj*lengthperchain+j] - centerofmass[1]) *(s_dybb[jj*lengthperchain+j] - centerofmass[1]) + (s_dzbb[jj*lengthperchain+j] - centerofmass[2]) *(s_dzbb[jj*lengthperchain+j] - centerofmass[2]);

                if(poly_spin_reference[(numbercores *pertrusionpercore *lengthpercore + jj*lengthperchain + j)]!=1)
                {
                    rsumbb += (s_dxbb[backbone_length+func_counter2] - centerofmass[0]) * (s_dxbb[backbone_length+func_counter2] - centerofmass[0]) + (s_dybb[backbone_length+func_counter2] - centerofmass[1]) * (s_dybb[backbone_length+func_counter2] - centerofmass[1]) + (s_dzbb[backbone_length+func_counter2] - centerofmass[2]) * (s_dzbb[backbone_length+func_counter2] - centerofmass[2]);
                    
                    func_counter2++;
                    func_counter3++;
                }
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "b%i rsum: %f  %f   %f   %f\t",j,rsumbb,((double)s_dxbb[jj*lengthperchain+j] - centerofmass[0]) *((double)s_dxbb[jj*lengthperchain+j] - centerofmass[0]),((double)s_dybb[jj*lengthperchain+j] - centerofmass[1]) *((double)s_dybb[jj*lengthperchain+j] - centerofmass[1]),((double)s_dzbb[jj*lengthperchain+j] - centerofmass[2]) *((double)s_dzbb[jj*lengthperchain+j] - centerofmass[2]));
//                fclose(stdoutlog);
            }
        
            radiusbb += sqrt((rsumbb / (lengthperchain+func_counter3)));

        }
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "radius: %f\n",radiusbb);
//                fclose(stdoutlog);
    
            rogbb[index] += radiusbb/numberofchains;
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "rog[index]: %f",rogbb[index]/visits[index]);
//                fclose(stdoutlog);

            double ssum[numbercores *pertrusionpercore] = { 0 };
            double tort = 0;

            int x1, y1, z1;
            int x2, y2, z2;

            for (int i = (0); i < numbercores * pertrusionpercore; i++)    // performs the rotations
            {
                int s_vector[3 *lengthpercore] = { 0 };
                double s_avg[3] = { 0 };

                x1 = s_dx[i *lengthpercore + 0] - s_dx[i *lengthpercore + 0 + 1];
                x2 = s_dx[i *lengthpercore + 0] - s_dx[i *lengthpercore + 0 + 2];

                y1 = s_dy[i *lengthpercore + 0] - s_dy[i *lengthpercore + 0 + 1];
                y2 = s_dy[i *lengthpercore + 0] - s_dy[i *lengthpercore + 0 + 2];

                z1 = s_dz[i *lengthpercore + 0] - s_dz[i *lengthpercore + 0 + 1];
                z2 = s_dz[i *lengthpercore + 0] - s_dz[i *lengthpercore + 0 + 2];

                s_vector[3 * 0 + 0] = y1 *z2 - z1 * y2;
                s_vector[3 * 0 + 1] = z1 *x2 - x1 * z2;
                s_vector[3 * 0 + 2] = x1 *y2 - y1 * x2;

                for (int j = 1; j < lengthpercore - 2; j++)    // performs the rotations
                {
        //            stdoutlog = fopen(stdoutlogname, "a");
        //            fprintf(stdoutlog, "\n{%f,%f,%f}\n", s_vector[3 *j + 0], s_vector[3 *j + 1], s_vector[3 *j + 2]);
        //            fclose(stdoutlog);
                    x1 = s_dx[i *lengthpercore + j] - s_dx[i *lengthpercore + j + 1];
                    x2 = s_dx[i *lengthpercore + j] - s_dx[i *lengthpercore + j + 2];

                    y1 = s_dy[i *lengthpercore + j] - s_dy[i *lengthpercore + j + 1];
                    y2 = s_dy[i *lengthpercore + j] - s_dy[i *lengthpercore + j + 2];

                    z1 = s_dz[i *lengthpercore + j] - s_dz[i *lengthpercore + j + 1];
                    z2 = s_dz[i *lengthpercore + j] - s_dz[i *lengthpercore + j + 2];

                    s_vector[3 *j + 0] = y1 *z2 - z1 *y2 + s_vector[3 *(j - 1) + 0];
                    s_vector[3 *j + 1] = z1 *x2 - x1 *z2 + s_vector[3 *(j - 1) + 1];
                    s_vector[3 *j + 2] = x1 *y2 - y1 *x2 + s_vector[3 *(j - 1) + 2];
        //            stdoutlog = fopen(stdoutlogname, "a");
        //            fprintf(stdoutlog, "{%f,%f,%f}\n", s_vector[3 *j + 0], s_vector[3 *j + 1], s_vector[3 *j + 2]);
        //            fclose(stdoutlog);
                }

                for (int j = 0; j < lengthpercore - 2; j++)    // performs the rotations
                {
                    s_avg[0] += s_vector[3 *j + 0];
                    s_avg[1] += s_vector[3 *j + 1];
                    s_avg[2] += s_vector[3 *j + 2];
                }

                s_avg[0] /= ((double)((lengthpercore - 2)));
                s_avg[1] /= ((double)((lengthpercore - 2)));
                s_avg[2] /= ((double)((lengthpercore - 2)));

                for (int j = 0; j < lengthpercore - 2; j++)    // performs the rotations
                {
                    ssum[i] = (((double)(s_vector[3 *(j)])) - s_avg[0]) *(((double)(s_vector[3 *(j)])) - s_avg[0]) + (((double)(s_vector[3 *(j) + 1])) - s_avg[1]) *(((double)(s_vector[3 *(j) + 1])) - s_avg[1]) + (((double)(s_vector[3 *(j) + 2])) - s_avg[2]) *(((double)(s_vector[3 *(j) + 2])) - s_avg[2]);
                }
            }

            for (int i = 0; i < numbercores * pertrusionpercore; i++)
            {
                tort += sqrt(ssum[i] / (lengthpercore - 2));
            }

            tortuosity[index] += tort / (numbercores *pertrusionpercore);

            double tortbb = 0;
            for (int jj = 0; jj < numberofchains; jj++)    // performs the rotations
            {
            double ssumbb = 0;
            
            int s_vector[3 *lengthperchain*numberofchains] = { 0 };
            double s_avg[3] = { 0 };

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "s_v0 = {%i,%i,%i}  \n",s_vector[3 *(jj*lengthperchain+0) + 0], s_vector[3 *(jj*lengthperchain+0) + 1], s_vector[3 *(jj*lengthperchain+0) + 2]);
//                fclose(stdoutlog);
            x1 = s_dxbb[jj*lengthperchain] - s_dxbb[jj*lengthperchain + 1];
            x2 = s_dxbb[jj*lengthperchain] - s_dxbb[jj*lengthperchain + 2];

            y1 = s_dybb[jj*lengthperchain] - s_dybb[jj*lengthperchain + 1];
            y2 = s_dybb[jj*lengthperchain] - s_dybb[jj*lengthperchain + 2];

            z1 = s_dzbb[jj*lengthperchain] - s_dzbb[(jj*lengthperchain) + 1];
            z2 = s_dzbb[jj*lengthperchain] - s_dzbb[(jj*lengthperchain) + 2];

            s_vector[3 * 0 + 0] = y1 *z2 - z1 * y2;
            s_vector[3 * 0 + 1] = z1 *x2 - x1 * z2;
            s_vector[3 * 0 + 2] = x1 *y2 - y1 * x2;
                //int j=0;
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "j%i jj%i x1 = %d x2 = %d\ty1 = %d y2 = %d\tz1 = %d z2 = %d \n s_v = {%i,%i,%i}  \n", j,jj,x1,x2,y1,y2,z1,z2,s_vector[3 *(jj*lengthperchain+j) + 0], s_vector[3 *(jj*lengthperchain+j) + 1], s_vector[3 *(jj*lengthperchain+j) + 2]);
//                fclose(stdoutlog);
                
            for (int j = 1; j < lengthperchain - 2; j++)    // performs the rotations
            {
                x1 = s_dxbb[jj*lengthperchain+j] - s_dxbb[jj*lengthperchain+j + 1];
                x2 = s_dxbb[jj*lengthperchain+j] - s_dxbb[jj*lengthperchain+j + 2];

                y1 = s_dybb[jj*lengthperchain+j] - s_dybb[jj*lengthperchain+j + 1];
                y2 = s_dybb[jj*lengthperchain+j] - s_dybb[jj*lengthperchain+j + 2];

                z1 = s_dzbb[jj*lengthperchain+j] - s_dzbb[jj*lengthperchain+j + 1];
                z2 = s_dzbb[jj*lengthperchain+j] - s_dzbb[jj*lengthperchain+j + 2];

                s_vector[3 *(jj*lengthperchain+j) + 0] = y1 *z2 - z1 *y2 + s_vector[3 *(jj*lengthperchain+j - 1) + 0];
                s_vector[3 *(jj*lengthperchain+j) + 1] = z1 *x2 - x1 *z2 + s_vector[3 *(jj*lengthperchain+j - 1) + 1];
                s_vector[3 *(jj*lengthperchain+j) + 2] = x1 *y2 - y1 *x2 + s_vector[3 *(jj*lengthperchain+j - 1) + 2];
                
//                            stdoutlog = fopen(stdoutlogname, "a");
//                             fprintf(stdoutlog, "j%i jj%i x1 = %d x2 = %d\ty1 = %d y2 = %d\tz1 = %d z2 = %d \n s_v = {%i,%i,%i}  s_v-1 = {%i,%i,%i}\n", j,jj,x1,x2,y1,y2,z1,z2,s_vector[3 *(jj*lengthperchain+j) + 0], s_vector[3 *(jj*lengthperchain+j) + 1], s_vector[3 *(jj*lengthperchain+j) + 2],s_vector[3 *(jj*lengthperchain+j-1) + 0], s_vector[3 *(jj*lengthperchain+j-1) + 1], s_vector[3 *(jj*lengthperchain+j-1) + 2]);
//                             fclose(stdoutlog);
            }

            for (int j = 0; j < lengthperchain - 2; j++)    // performs the rotations
            {
                s_avg[0] += s_vector[3 *(jj*lengthperchain+j) + 0];
                s_avg[1] += s_vector[3 *(jj*lengthperchain+j) + 1];
                s_avg[2] += s_vector[3 *(jj*lengthperchain+j) + 2];
            }

            s_avg[0] /= ((double)((lengthperchain - 2)));
            s_avg[1] /= ((double)((lengthperchain - 2)));
            s_avg[2] /= ((double)((lengthperchain - 2)));

//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "\ns_avg = {%f,%f,%f} \n\n", s_avg[0],s_avg[1],s_avg[2]);
//                fclose(stdoutlog);
                
            for (int j = 0; j < lengthperchain - 2; j++)    // performs the rotations
            {
                ssumbb += (((double)(s_vector[3 *(jj*lengthperchain+j)])) - s_avg[0]) *(((double)(s_vector[3 *(jj*lengthperchain+j)])) - s_avg[0]) +
                (((double)(s_vector[3 *(jj*lengthperchain+j) + 1])) - s_avg[1]) *(((double)(s_vector[3 *(jj*lengthperchain+j) + 1])) - s_avg[1]) +
                (((double)(s_vector[3 *(jj*lengthperchain+j) + 2])) - s_avg[2]) *(((double)(s_vector[3 *(jj*lengthperchain+j) + 2])) - s_avg[2]);
                
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "j%i ssumbb = %f \n", j,ssumbb);
//                fclose(stdoutlog);
            }

            tortbb += sqrt(ssumbb / (backbone_length - 2));
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "jj%i tortbb = %f \n", jj,tortbb);
//                fclose(stdoutlog);

    }
    tortuositybb[index] += tortbb/numberofchains;
    


//    FILE* fptr;
//    sprintf(filename2, "Functional_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
//    fptr = fopen(filename2, "w");
//    fclose(fptr);
//        int summed_diff[3] = {0,0,0};
//                 int func_diff[3] = {0,0,0};
//            int arm = 0;
//            sprintf(filename, "Chain_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
//            if ((file = fopen(filename, "w")) == NULL)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
//                fclose(stdoutlog);
//                MPI_Abort(MPI_COMM_WORLD, 1);
//            }
//            else
//            {
//                for (int i = numbercores*pertrusionpercore*lengthpercore; i < numbercores*pertrusionpercore*lengthpercore+backbone_length; i++)
//                {
//                    summed_diff[0]=poly_lattice_coordinates[3 * (arm*lengthperchain)];
//                    summed_diff[1]=poly_lattice_coordinates[3 * (arm*lengthperchain) + 1];
//                    summed_diff[2]=poly_lattice_coordinates[3 * (arm*lengthperchain) + 2];
//
//                    summed_diff[0]+=dxbb[i-numbercores*pertrusionpercore*lengthpercore];
//                    summed_diff[1]+=dybb[i-numbercores*pertrusionpercore*lengthpercore];
//                    summed_diff[2]+=dzbb[i-numbercores*pertrusionpercore*lengthpercore];
//
//                    fprintf(file, "%i\t%i\t%i\n",  summed_diff[0], summed_diff[1], summed_diff[2]);
//
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "%i\t%i\t%i   \t%i   %i\n",  summed_diff[0], summed_diff[1], summed_diff[2],i,poly_spin_reference[i]);
//                    fclose(stdoutlog);
//
//                    if(poly_spin_reference[i]>=5)
//                    {
//                        sprintf(filename2, "Functional_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
//                        if ((fptr = fopen(filename2, "a")) == NULL)
//                        {
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
//                            fclose(stdoutlog);
//                            MPI_Abort(MPI_COMM_WORLD, 1);
//                        }
//                        else{
//                            for (int j = 0; j < lengthpercore; j++)
//                            {
//                                 func_diff[0]=summed_diff[0];
//                                 func_diff[1]=summed_diff[1];
//                                func_diff[2]=summed_diff[2];
//
//                            func_diff[0]+=dx[lengthpercore*pertrusionpercore*(poly_spin_reference[i]-5)+j];
//                                func_diff[1]+=dy[lengthpercore*pertrusionpercore*(poly_spin_reference[i]-5)+j];
//                                func_diff[2]+=dz[lengthpercore*pertrusionpercore*(poly_spin_reference[i]-5)+j];
//
//                            fprintf(fptr, "%i\t%i\t%i\n",  func_diff[0], func_diff[1], func_diff[2]);
//
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog,  "\t%i\t%i\t%i\n",  func_diff[0], func_diff[1], func_diff[2]);
//                                fclose(stdoutlog);
//                            }
//                        }
//
//                        fclose(fptr);
//                    }
//
//                    if(((i-numbercores*pertrusionpercore*lengthpercore)+1)%lengthperchain==0)
//                    {
//                        arm++;
//                    }
//                }
//            }
//            fclose(file);
//
//    int s_array[3]={0};
//    arm =0;
//    int s_dist_array[3]={0};
//                sprintf(filename, "Solvent_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
//                if ((file = fopen(filename, "w")) == NULL)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
//                    fclose(stdoutlog);
//                    MPI_Abort(MPI_COMM_WORLD, 1);
//                }
//                else
//                {
//                    for (int i = 0; i < 3*no_solvent_sites; i++)
//                    {
//                        general_coord(s_array,solvent_loc[3*arm]);
//
//                        point_distance(solvent_loc[3*arm],solvent_loc[i],s_dist_array);
//
//                        s_array[0]+=s_dist_array[0];
//                        s_array[1]+=s_dist_array[1];
//                        s_array[2]+=s_dist_array[2];
//
//                        fprintf(file, "%i\t%i\t%i\n",  s_array[0], s_array[1], s_array[2]);
//
//                        if((i+1)%3==0)
//                        {
//                            arm++;
//                        }
//                    }
//                }
//                fclose(file);
//
//    sprintf(filename, "Ion_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
//                if ((file = fopen(filename, "w")) == NULL)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
//                    fclose(stdoutlog);
//                    MPI_Abort(MPI_COMM_WORLD, 1);
//                }
//                else
//                {
//                    for (int i = 0; i < no_ion_sites; i++)
//                    {
//                            general_coord(s_array,ion_loc[i]);
//
//                        fprintf(file, "%i\t%i\t%i\n",  s_array[0], s_array[1], s_array[2]);
//
//                    }
//                }
//                fclose(file);
//
    
    visits[index]++;
//    if(myid==6)MPI_Abort(MPI_COMM_WORLD,1);
}

void accessiblelevels()
{
    //stdoutlog = fopen(stdoutlogname, "a");

    if (myid == 0)    // 'root' in energy window, receive individual g(E) and send merged g(E)
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.005)
            {
                real_lngE[i] = 1;
            }
        }

        MPI_Recv(&maximum_en_config, 1, MPI_INT, numprocs - 1, 71, MPI_COMM_WORLD, &status);
        MPI_Recv(&max_printed, 1, MPI_INT, numprocs - 1, 72, MPI_COMM_WORLD, &status);

        for (int ii = 1; ii < numprocs; ii++)
        {
            MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, ii, 76, MPI_COMM_WORLD, &status);    // get other dens. of states
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\nProc %i: Received lngE from Proc. %i\n", myid, ii);
            fclose(stdoutlog);
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (lngE_buf[i] > 0.005)
                {
                    real_lngE[i] = 1;
                }
            }

            for (int i = 0; i < hist_size; i++)
            {
                lngE_buf[i] = 0.0;
            }
        }

        fclose(file);

        int w_l = 0;

        sprintf(filename, "Access_Levels.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i < Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.005)
                {
                    if (w_l == 0)
                    {
                        if (i != minimum_en_config)
                        {
                            printed = 0;
                        }

                        minimum_en_config = i;
                        w_l = 1;
                    }

                    if (w_l == 1)
                    {
                        if (i > maximum_en_config)
                        {
                            max_printed = 0;
                        }

                        maximum_en_config = i;
                    }

                    fprintf(file, "%i\t%f\n", (i + (Eglobalmin)), real_lngE[i]);
                }
            }
        }

        fclose(file);

        MPI_Send(&maximum_en_config, 1, MPI_INT, numprocs - 1, 73, MPI_COMM_WORLD);
        MPI_Send(&max_printed, 1, MPI_INT, numprocs - 1, 74, MPI_COMM_WORLD);
    }
    else    // send individual g(E) and receive merged g(E)
    {
        if (myid == numprocs - 1)
        {
            MPI_Send(&maximum_en_config, 1, MPI_INT, 0, 71, MPI_COMM_WORLD);
            MPI_Send(&max_printed, 1, MPI_INT, 0, 72, MPI_COMM_WORLD);
        }

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "ufig");
        fclose(stdoutlog);
        MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, 0, 76, MPI_COMM_WORLD);
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, 0);
        fclose(stdoutlog);

        if (myid == numprocs - 1)
        {
            MPI_Recv(&maximum_en_config, 1, MPI_INT, 0, 73, MPI_COMM_WORLD, &status);
            MPI_Recv(&max_printed, 1, MPI_INT, 0, 74, MPI_COMM_WORLD, &status);
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Recieved max details %i %i\n", max_printed, maximum_en_config);
            fclose(stdoutlog);
        }
    }
}

void recombine(double countd)
{
   // stdoutlog = fopen(stdoutlogname, "a");

    sprintf(filename, "SysParam_%i.txt",myid);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5) fprintf(file, "%i\t%1.13f\t%1.13f\t%1.13f\t%1.13f\t%lli\n", i, (float) rog[i], (float) tortuosity[i] , (float) rogbb[i] , (float) tortuositybb[i] , visits[i]);
        }
        
    }
    fclose(file);
    
    sprintf(filename, "POPO_%i.txt",myid);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5){ fprintf(file, "%i\t", i);
            fprintf(file, "%f\t", radial_visits_poly[i]);
            //for(int j =0;j<no_ion_sites;j++)
            //{
                for(int k =0;k<14;k++)
                {
                    fprintf(file, "%f\t", ((double)radial_dist_POPO[i*14+k]));
                }
            //}
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }

    sprintf(filename, "PMPO_%i.txt",myid);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5){ fprintf(file, "%i\t", i);
            fprintf(file, "%f\t", radial_visits_poly[i]);
            //for(int j =0;j<no_ion_sites;j++)
            //{
                for(int k =0;k<14;k++)
                {
                    fprintf(file, "%f\t", ((double)radial_dist_PMPO[i*14+k]));
                }
            //}
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
    
    sprintf(filename, "EuO_%i.txt",myid);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i < Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5){ fprintf(file, "%i\t", i);
                fprintf(file, "%f\t", radial_visits_ion[i]);
            
            //for(int j =0;j<no_ion_sites;j++)
            //{
                for(int k =0;k<14;k++)
                {
                    fprintf(file, "%f\t", ((double)radial_dist_EuO[i*14+k]));
                }
            //}
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
    
    sprintf(filename, "EuPO_%i.txt",myid);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5){ fprintf(file, "%i\t", i);
            fprintf(file, "%f\t", radial_visits_ion[i]);
            //for(int j =0;j<no_ion_sites;j++)
            //{
                for(int k =0;k<14;k++)
                {
                    fprintf(file, "%f\t",((double)radial_dist_EuPO[i*14+k]));
                }
            //}
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
    
    sprintf(filename, "EuEu_%i.txt",myid);
    if ((file = fopen(filename, "w")) == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5){ fprintf(file, "%i\t", i);
            fprintf(file, "%f\t", radial_visits_ion[i]);
            //for(int j =0;j<no_ion_sites;j++)
            //{
                for(int k =0;k<14;k++)
                {
                    fprintf(file, "%f\t",((double)radial_dist_EuEu[i*14+k]));
                }
            //}
                fprintf(file, "\n");
            }
        }
        fclose(file);
    }
	
    double init_dos = 19962000000;

    int init_check = 0;
    int init_check_2 = 0;
    double rec_lnge;
    double rec_real_lnge;
    int rec_en;
    int rec_en_other;
    int rec_en_2;
    int rec_en_other_2;
    double microT_compare;

    if (myid == 0)    // 'root' in energy window, receive individual g(E) and send merged g(E)
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.5 && init_check == 1)
            {
                real_lngE[i] = rec_real_lnge + (lngE[i] - rec_lnge);
                microT[i] = rec_real_lnge + ((rec_lnge - lngE[i]) / (rec_en - i));

                rec_real_lnge = real_lngE[i];
                rec_lnge = lngE[i];
                rec_en = i;
            }

            if (lngE[i] > 0.5 && init_check == 0)
            {
                real_lngE[i] = init_dos;
                init_check = 1;
                rec_real_lnge = real_lngE[i];
                rec_lnge = lngE[i];
                rec_en = i;
            }
        }

        //test output 1 for 0
        sprintf(filename, "TestL%iq%i.HE.proc%04i.iter0", L1dim_S_xy, q, myid);
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (lngE[i] > 0.5) fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + (Eglobalmin)), lngE[i], real_lngE[i], microT[i]);
            }

            fclose(file);
        }

        for (int ii = 1; ii < numprocs; ii++)
        {
            MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, ii, 76, MPI_COMM_WORLD, &status);    // get other dens. of states
            stdoutlog = fopen(stdoutlogname, "a");
            
            fprintf(stdoutlog, "\nProc %i: Received lngE from Proc. %i\n", myid, ii);
            fclose(stdoutlog);
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5 && lngE_buf[i] > 0.5 && init_check_2 == 2)
                {
                    microT_buf[i] = real_lngE[rec_en_2] + ((rec_lnge - lngE_buf[i]) / (rec_en_2 - i));

                    if (abs(microT[i] - microT_buf[i]) < microT_compare)
                    {
                        microT_compare = abs(microT[i] - microT_buf[i]);
                        rec_en_other_2 = i;
                    }

                    rec_lnge = lngE_buf[i];
                    rec_en_2 = i;

                    /*fprintf(stdoutlog,"\n%i\n",i);
                    fprintf(stdoutlog,"\n%i\n",rec_en_other_2);
                    fprintf(stdoutlog,"\n%e\n",lngE_buf[i]);
                    fprintf(stdoutlog,"\n%f\n",microT[i]);
                    fprintf(stdoutlog,"\n%f\n",microT_buf[i]);
                    fprintf(stdoutlog,"\n%f\n", microT_compare);
                    */
                }

                if (real_lngE[i] > 0.5 && lngE_buf[i] > 0.5 && init_check_2 == 1)
                {
                    init_check_2 = 2;

                    microT_buf[i] = real_lngE[rec_en_2] + ((rec_lnge - lngE_buf[i]) / (rec_en_2 - i));
                    microT_compare = abs(microT[i] - microT_buf[i]);

                    rec_lnge = lngE_buf[i];

                    /*
                    fprintf(stdoutlog,"\n%i\n",i);
                    fprintf(stdoutlog,"\n%i\n",rec_en_2);
                    fprintf(stdoutlog,"\n%e\n",lngE_buf[i]);
                    fprintf(stdoutlog,"\n%f\n",microT[i]);
                    fprintf(stdoutlog,"\n%f\n",microT_buf[i]);
                    fprintf(stdoutlog,"\n%f\n", microT_compare);
                     */
                    rec_en_other_2 = i;
                    rec_en_2 = i;
                }

                if (real_lngE[i] > 0.5 && lngE_buf[i] > 0.5 && init_check_2 == 0)
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    
                    
                    fprintf(stdoutlog, "\n%i\n", i);
                    fclose(stdoutlog);
                    rec_en_2 = i;
                    init_check_2 = 1;
                    rec_lnge = lngE_buf[i];
                }
            }

            init_check_2 = 0;

            for (int i = rec_en_other_2; i <= Eglobalwidth; i++)
            {
                if (lngE_buf[i] > 0.5 && init_check_2 == 1)
                {
                    real_lngE[i] = rec_real_lnge + (lngE_buf[i] - rec_lnge);
                    microT[i] = rec_real_lnge + ((rec_lnge - lngE_buf[i]) / (rec_en_2 - i));

                    /*
                    fprintf(stdoutlog,"\ng %i\n",i);
                    fprintf(stdoutlog,"\n%i\n",rec_en_2);
                    fprintf(stdoutlog,"\n%e\n",rec_real_lnge);
                    fprintf(stdoutlog,"\n%e\n",lngE_buf[i]);
                    fprintf(stdoutlog,"\n%e\n",microT[i]);
                    fprintf(stdoutlog,"\n%f\n",rec_real_lnge+((rec_lnge-lngE_buf[i])/(rec_en_2-i)));
                    */

                    rec_real_lnge = real_lngE[i];
                    rec_lnge = lngE_buf[i];
                    rec_en_2 = i;
                }

                if (lngE_buf[i] > 0.5 && init_check_2 == 0)
                {
                    rec_real_lnge = real_lngE[i];
                    microT[i] = microT_buf[i];

                    rec_en_2 = i;
                    init_check_2 = 1;
                    rec_lnge = lngE_buf[i];
                }
            }

            for (int i = 0; i < hist_size; i++)
            {
                microT_buf[i] = 0.0;
                lngE_buf[i] = 0.0;
            }
        }

        sprintf(filename, "Test2L%iq%i.HE.proc%04i.iter0", L1dim_S_xy, q, myid);
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5) fprintf(file, "%i\t%i\t%e\t%f\t%e\n", i, (i + (Eglobalmin)), lngE[i], real_lngE[i], microT[i]);
            }

            fclose(file);
        }

        sprintf(filename, "Recombined_Output_%e.txt", countd);
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                
                if (real_lngE[i] > 0.5) fprintf(file, "%i\t%f\n", (i + (Eglobalmin)), real_lngE[i]);
            }

           fclose(file);
        }

        double* rd_PMPO = (double*) malloc(hist_size*14* sizeof(double));
        double* rd_POPO = (double*) malloc(hist_size*14* sizeof(double));
        double* rd_EuPO = (double*) malloc(hist_size*14* sizeof(double));
        double* rd_EuO = (double*) malloc(hist_size*14* sizeof(double));
		double* rd_EuEu = (double*) malloc(hist_size*14* sizeof(double));
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (real_lngE[i] > 0.5)
            {
                
                    for(int k =0;k<14;k++)
                    {
                        rd_EuPO[i*14+k]=0;
						rd_EuEu[i*14+k]=0;
						rd_POPO[i*14+k]=0;
                        rd_EuO[i*14+k]=0;
                        rd_PMPO[i*14+k]=0;
                    }
            }
        }
        for (int ii = 1; ii < numprocs; ii++)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            
            
            fprintf(stdoutlog, "\nProcess Paramw %i\n", ii);
            fclose(stdoutlog);
            
            MPI_Recv(&tortuosity_buf[0], hist_size, MPI_DOUBLE, ii, 83, MPI_COMM_WORLD, &status);
            MPI_Recv(&rog_buf[0], hist_size, MPI_DOUBLE, ii, 84, MPI_COMM_WORLD, &status);
            MPI_Recv(&visits_buf[0], hist_size, MPI_LONG_LONG, ii, 85, MPI_COMM_WORLD, &status);

            MPI_Recv(&tortuositybb_buf[0], hist_size, MPI_DOUBLE, ii, 86, MPI_COMM_WORLD, &status);
            MPI_Recv(&rogbb_buf[0], hist_size, MPI_DOUBLE, ii, 87, MPI_COMM_WORLD, &status);
           
            MPI_Recv(&radial_dist_EuPO_buf[0], hist_size*14, MPI_FLOAT, ii, 390, MPI_COMM_WORLD, &status);
            MPI_Recv(&radial_dist_EuO_buf[0], hist_size*14, MPI_FLOAT, ii, 391, MPI_COMM_WORLD, &status);
            MPI_Recv(&radial_dist_PMPO_buf[0], hist_size*14, MPI_FLOAT, ii, 392, MPI_COMM_WORLD, &status);
             MPI_Recv(&radial_dist_POPO_buf[0], hist_size*14, MPI_FLOAT, ii, 3923, MPI_COMM_WORLD, &status);
			 MPI_Recv(&radial_dist_EuEu_buf[0], hist_size*14, MPI_FLOAT, ii, 3924, MPI_COMM_WORLD, &status);
            MPI_Recv(&radial_visits_ion_buf[0], hist_size, MPI_FLOAT, ii, 393, MPI_COMM_WORLD, &status);
            MPI_Recv(&radial_visits_poly_buf[0], hist_size, MPI_FLOAT, ii, 394, MPI_COMM_WORLD, &status);
            
             MPI_Recv(&perf_seq_buf[0], hist_size, MPI_INT, ii, 494, MPI_COMM_WORLD, &status);
             MPI_Recv(&func_per_seq_buf[0], hist_size, MPI_FLOAT, ii, 495, MPI_COMM_WORLD, &status);
             MPI_Recv(&part_seq_buf[0], hist_size, MPI_INT, ii, 496, MPI_COMM_WORLD, &status);
            
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5)
                {
                    tortuosity[i] += tortuosity_buf[i];
                    rog[i] += rog_buf[i];
                    tortuositybb[i] += tortuositybb_buf[i];
                    rogbb[i] += rogbb_buf[i];
                    visits[i] += visits_buf[i];
                    
                    
                    
                    radial_visits_ion[i]+=radial_visits_ion_buf[i];
                        
                    for(int k =0;k<14;k++)
                    {
                        radial_dist_EuPO[i*14+k]+=radial_dist_EuPO_buf[i*14+k];
                        radial_dist_EuO[i*14+k]+=radial_dist_EuO_buf[i*14+k];
						radial_dist_EuEu[i*14+k]+=radial_dist_EuEu_buf[i*14+k];
                    }
                        
                   
                    for(int k =0;k<14;k++)
                    {
                        radial_dist_PMPO[i*14+k]+=radial_dist_PMPO_buf[i*14+k];
                        radial_dist_POPO[i*14+k]+=radial_dist_POPO_buf[i*14+k];
                    }
                    radial_visits_poly[i]+=radial_visits_poly_buf[i];
                    
                    perf_seq[i] += perf_seq_buf[i];
                    func_per_seq[i] += func_per_seq_buf[i];
                    part_seq[i] += part_seq_buf[i];
                    
                    
                }
            }
        }
        
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if(radial_visits_ion[i]>2)
            {
                radial_visits_ion[i]-=2;
            }
            
            if(radial_visits_poly[i]>2)
            {
                radial_visits_poly[i]-=2;
            }
            
            if (real_lngE[i] > 0.5)
            {
                for(int j =0;j<no_ion_sites;j++)
                {
                    for(int k =0;k<14;k++)
                    {
                    rd_EuPO[i*14+k]+=((double)radial_dist_EuPO[i*14+k])/((double)radial_visits_ion[i]);
                    rd_EuO[i*14+k]+=((double)radial_dist_EuO[i*14+k])/((double)radial_visits_ion[i]);
					rd_EuEu[i*14+k]+=((double)radial_dist_EuEu[i*14+k])/((double)radial_visits_ion[i]);
                    }
                }
                
                    for(int k =0;k<14;k++)
                    { rd_PMPO[i*14+k]+=((double)radial_dist_PMPO[i*14+k])/((double)radial_visits_poly[i]);
                        
                        rd_POPO[i*14+k]+=((double)radial_dist_POPO[i*14+k])/((double)radial_visits_poly[i]);
                    }
                
            }
        }

        sprintf(filename, "SysParam_%e.txt", countd);
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5) fprintf(file, "%i\t%1.13f\t%1.13f\t%1.13f\t%1.13f\t%lli\t%1.13f\t%1.13f\t%1.13f\n", i, (float) rog[i] / ((double) visits[i]), (float) tortuosity[i] / ((double) visits[i]), (float) rogbb[i] / ((double) visits[i]), (float) tortuositybb[i] / ((double) visits[i]), visits[i],(float) part_seq[i]/ ((double) visits[i]),(float)perf_seq[i]/ ((double) visits[i]),(func_per_seq[i]/(((float) part_seq[i])+((float)perf_seq[i]))));
            }
            fclose(file);
        }
        
        sprintf(filename, "POPO.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5){ fprintf(file, "%i\t", i);
                fprintf(file, "%f\t", radial_visits_poly[i]);
                //for(int j =0;j<no_ion_sites;j++)
                //{
                    for(int k =0;k<14;k++)
                    {
                        fprintf(file, "%f\t", rd_POPO[i*14+k]/numbercores);
                    }
                //}
                    fprintf(file, "\n");
                }
            }
            fclose(file);
        }

        sprintf(filename, "PMPO.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5){ fprintf(file, "%i\t", i);
                fprintf(file, "%f\t", radial_visits_poly[i]);
                //for(int j =0;j<no_ion_sites;j++)
                //{
                    for(int k =0;k<14;k++)
                    {
                        fprintf(file, "%f\t", rd_PMPO[i*14+k]/numbercores);
                    }
                //}
                    fprintf(file, "\n");
                }
            }
            fclose(file);
        }
        
        sprintf(filename, "EuO.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i < Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5){ fprintf(file, "%i\t", i);
                    fprintf(file, "%f\t", radial_visits_ion[i]);
                
                //for(int j =0;j<no_ion_sites;j++)
                //{
                    for(int k =0;k<14;k++)
                    {
                        fprintf(file, "%f\t", rd_EuO[i*14+k]/no_ion_sites);
                    }
                //}
                    fprintf(file, "\n");
                }
            }
            fclose(file);
        }
        
        sprintf(filename, "EuPO.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5){ fprintf(file, "%i\t", i);
                fprintf(file, "%f\t", radial_visits_ion[i]);
                //for(int j =0;j<no_ion_sites;j++)
                //{
                    for(int k =0;k<14;k++)
                    {
                        fprintf(file, "%f\t", rd_EuPO[i*14+k]/no_ion_sites);
                    }
                //}
                    fprintf(file, "\n");
                }
            }
            fclose(file);
        }
        
		sprintf(filename, "EuEu.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5){ fprintf(file, "%i\t", i);
                fprintf(file, "%f\t", radial_visits_ion[i]);
                //for(int j =0;j<no_ion_sites;j++)
                //{
                    for(int k =0;k<14;k++)
                    {
                        fprintf(file, "%f\t", rd_EuEu[i*14+k]/no_ion_sites);
                    }
                //}
                    fprintf(file, "\n");
                }
            }
            fclose(file);
        }
		
        for (int ii = 1; ii < numprocs; ii++)
        {
            MPI_Recv(&tortuosity_buf[0], hist_size, MPI_DOUBLE, ii, 883, MPI_COMM_WORLD, &status);
            MPI_Recv(&rog_buf[0], hist_size, MPI_DOUBLE, ii, 884, MPI_COMM_WORLD, &status);
            MPI_Recv(&visits_buf[0], hist_size, MPI_LONG_LONG, ii, 885, MPI_COMM_WORLD, &status);

            MPI_Recv(&tortuositybb_buf[0], hist_size, MPI_DOUBLE, ii, 886, MPI_COMM_WORLD, &status);
            MPI_Recv(&rogbb_buf[0], hist_size, MPI_DOUBLE, ii, 887, MPI_COMM_WORLD, &status);
           
            MPI_Recv(&radial_dist_EuPO_buf[0], hist_size*14, MPI_FLOAT, ii, 3390, MPI_COMM_WORLD, &status);
                       MPI_Recv(&radial_dist_EuO_buf[0], hist_size*14, MPI_FLOAT, ii, 3391, MPI_COMM_WORLD, &status);
                       MPI_Recv(&radial_dist_PMPO_buf[0], hist_size*14, MPI_FLOAT, ii, 3392, MPI_COMM_WORLD, &status);
             MPI_Recv(&radial_dist_POPO_buf[0], hist_size*14, MPI_FLOAT, ii, 33923, MPI_COMM_WORLD, &status);
			 MPI_Recv(&radial_dist_EuEu_buf[0], hist_size*14, MPI_FLOAT, ii, 33924, MPI_COMM_WORLD, &status);
                       MPI_Recv(&radial_visits_ion_buf[0], hist_size, MPI_FLOAT, ii, 3393, MPI_COMM_WORLD, &status);
                       MPI_Recv(&radial_visits_poly_buf[0], hist_size, MPI_FLOAT, ii, 3394, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&perf_seq_buf[0], hist_size, MPI_INT, ii, 4494, MPI_COMM_WORLD, &status);
                        MPI_Recv(&func_per_seq_buf[0], hist_size, MPI_FLOAT, ii, 4495, MPI_COMM_WORLD, &status);
                        MPI_Recv(&part_seq_buf[0], hist_size, MPI_INT, ii, 4496, MPI_COMM_WORLD, &status);
            
            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5)
                {
                    tortuosity[i] -= tortuosity_buf[i];
                    rog[i] -= rog_buf[i];
                    tortuositybb[i] -= tortuositybb_buf[i];
                    rogbb[i] -= rogbb_buf[i];
                    visits[i] -= visits_buf[i];
                 
                        for(int k =0;k<14;k++)
                        {
                            radial_dist_EuPO[i*14+k]-=radial_dist_EuPO_buf[i*14+k];
                            radial_dist_EuO[i*14+k]-=radial_dist_EuO_buf[i*14+k];
							 radial_dist_EuEu[i*14+k]-=radial_dist_EuEu_buf[i*14+k];
                        }
                        
                        radial_visits_ion[i]-=radial_visits_ion_buf[i];
                    
                    
                        for(int k =0;k<14;k++)
                        {
                            radial_dist_PMPO[i*14+k]-=radial_dist_PMPO_buf[i*14+k];
                             radial_dist_POPO[i*14+k]-=radial_dist_POPO_buf[i*14+k];
                        }
                    
                        radial_visits_poly[i]-=radial_visits_poly_buf[i];
                    
                    if(radial_visits_ion[i]>2)
                    {
                        radial_visits_ion[i]+=2;
                    }
                    
                    if(radial_visits_poly[i]>2)
                    {
                        radial_visits_poly[i]+=2;
                    }
                    
                    perf_seq[i] -= perf_seq_buf[i];
                    func_per_seq[i] -= func_per_seq_buf[i];
                    part_seq[i] -= part_seq_buf[i];
                }
            }
        }
        
        free(rd_EuPO);
		free(rd_EuEu);
        free(rd_EuO);
        free(rd_PMPO);
        free(rd_POPO);
        
    }
    else    // send individual g(E) and receive merged g(E)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "reufig");
        fclose(stdoutlog);
        MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, 0, 76, MPI_COMM_WORLD);
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, 0);

        fclose(stdoutlog);
        MPI_Send(&tortuosity[0], hist_size, MPI_DOUBLE, 0, 83, MPI_COMM_WORLD);
        MPI_Send(&rog[0], hist_size, MPI_DOUBLE, 0, 84, MPI_COMM_WORLD);
        MPI_Send(&visits[0], hist_size, MPI_LONG_LONG, 0, 85, MPI_COMM_WORLD);

        MPI_Send(&tortuositybb[0], hist_size, MPI_DOUBLE, 0, 86, MPI_COMM_WORLD);
        MPI_Send(&rogbb[0], hist_size, MPI_DOUBLE, 0, 87, MPI_COMM_WORLD);
       
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "first old batch sent");
        fclose(stdoutlog);
        
        MPI_Send(&radial_dist_EuPO[0], hist_size*14, MPI_FLOAT, 0, 390, MPI_COMM_WORLD);
        MPI_Send(&radial_dist_EuO[0], hist_size*14, MPI_FLOAT, 0, 391, MPI_COMM_WORLD);
        MPI_Send(&radial_dist_PMPO[0], hist_size*14, MPI_FLOAT, 0, 392, MPI_COMM_WORLD);
        MPI_Send(&radial_dist_POPO[0], hist_size*14, MPI_FLOAT, 0, 3923, MPI_COMM_WORLD);
		MPI_Send(&radial_dist_EuEu[0], hist_size*14, MPI_FLOAT, 0, 3924, MPI_COMM_WORLD);
        MPI_Send(&radial_visits_ion[0], hist_size, MPI_FLOAT, 0, 393, MPI_COMM_WORLD);
        MPI_Send(&radial_visits_poly[0], hist_size, MPI_FLOAT, 0, 394, MPI_COMM_WORLD);
        
        MPI_Send(&perf_seq[0], hist_size, MPI_INT, 0, 494, MPI_COMM_WORLD);
        MPI_Send(&func_per_seq[0], hist_size, MPI_FLOAT, 0, 495, MPI_COMM_WORLD);
        MPI_Send(&part_seq[0], hist_size, MPI_INT, 0, 496, MPI_COMM_WORLD);
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "first new batch sent");
        fclose(stdoutlog);
        
        MPI_Send(&tortuosity[0], hist_size, MPI_DOUBLE, 0, 883, MPI_COMM_WORLD);
        MPI_Send(&rog[0], hist_size, MPI_DOUBLE, 0, 884, MPI_COMM_WORLD);
        MPI_Send(&visits[0], hist_size, MPI_LONG_LONG, 0, 885, MPI_COMM_WORLD);

        MPI_Send(&tortuositybb[0], hist_size, MPI_DOUBLE, 0, 886, MPI_COMM_WORLD);
        MPI_Send(&rogbb[0], hist_size, MPI_DOUBLE, 0, 887, MPI_COMM_WORLD);
        
     
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "second old batch sent");
        fclose(stdoutlog);
        
        MPI_Send(&radial_dist_EuPO[0], hist_size*14, MPI_FLOAT, 0, 3390, MPI_COMM_WORLD);
               MPI_Send(&radial_dist_EuO[0], hist_size*14, MPI_FLOAT, 0, 3391, MPI_COMM_WORLD);
               MPI_Send(&radial_dist_PMPO[0], hist_size*14, MPI_FLOAT, 0, 3392, MPI_COMM_WORLD);
        MPI_Send(&radial_dist_POPO[0], hist_size*14, MPI_FLOAT, 0, 33923, MPI_COMM_WORLD);
         MPI_Send(&radial_dist_EuEu[0], hist_size*14, MPI_FLOAT, 0, 33924, MPI_COMM_WORLD);
               MPI_Send(&radial_visits_ion[0], hist_size, MPI_FLOAT, 0, 3393, MPI_COMM_WORLD);
               MPI_Send(&radial_visits_poly[0], hist_size, MPI_FLOAT, 0, 3394, MPI_COMM_WORLD);
        
        MPI_Send(&perf_seq[0], hist_size, MPI_INT, 0, 4494, MPI_COMM_WORLD);
               MPI_Send(&func_per_seq[0], hist_size, MPI_FLOAT, 0, 4495, MPI_COMM_WORLD);
               MPI_Send(&part_seq[0], hist_size, MPI_INT, 0, 4496, MPI_COMM_WORLD);
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "second new batch sent");
        fclose(stdoutlog);
    }

    
}

void restart(int(&energie),int(&eold),int(&iteration),double(&lnf),int(&swap_every),double(&sweep))
{
    int stop =0;
            sprintf(resetbuffer, "SysStatus.L%i.proc%04i.copy5", L1dim_S_xy, myid);    // Preliminary output file that allows for process inspection
               FILE* rfptr;
               long long phldr;
    long long counter;
    long long counterp=0;
               int Ewidth = Emaxindex - Eminindex;
               int fakeEmin1 = Eminindex;
               int fakeEmin2 = Eminindex;
       
           int fakeEmin3 = Eminindex;
           int fakeEmin4 = Eminindex;
           int fakeEmin5 = Eminindex;
           int fakeEmin6 = Eminindex;
           int fakeEmin7 = Eminindex;
           int fakeEmin8 = 0;
           int fakeEmin9 = 0;
           int fakeEmin10 = 0;
           int fakeEmin11 = ((Eminindex - 1) * 14) + 14;
           int fakeEmin12 = ((Eminindex - 1) * 14) + 14;
           int fakeEmin13 = ((Eminindex - 1) * 14) + 14;
           int fakeEmin14 = Eminindex;
           int fakeEmin15 = Eminindex;
    int fakeEmin17 = Eminindex;
    int fakeEmin18 = Eminindex;
    int fakeEmin19 = Eminindex;
    int fakeEmin20 = Eminindex;
    int fakeEmin16 = ((Eminindex - 1) * 14) + 14;
    int chain_increment=0;
    
           int indexer =1;
       
               
               rfptr = fopen(resetbuffer, "r");
       
               int resetiter = 1;
       
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "starting\n");
    fclose(stdoutlog);
               while (fscanf(rfptr, "%lld,%lld\n",&counter, &phldr) > 0)
               {
                   
                   if (resetiter == 1)
                   {
                          //ratiovalue =0.9956;
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "1: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter == 2)
                   {
                       energie = phldr;
       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "2: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter == 3)
                   {
                       eold = phldr;
       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "3: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter == 4)
                   {
                       sweep = phldr;
       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "4: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter == 5)
                   {
                       swap_every = phldr;
       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "5: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter == 6)
                   {
                       lnf = ((double)phldr) / pow(10, 15);
       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "6: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter == 7)
                   {
                       iteration = phldr;
       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "7: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8) && resetiter < (8 + Ewidth+1))
                   {
                       HE[fakeEmin1] = ((int)phldr);
       
                       fakeEmin1++;
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "8: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

       
                   if (resetiter >= (8 + (Ewidth+1)) && resetiter < (8 + (Ewidth+1)*2))
                   {
                       lngE[fakeEmin2] = ((double)phldr) / pow(10, 9);
       
                       fakeEmin2++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "9: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*2) && resetiter < (8 + (Ewidth+1)*3))
                   {
                       rog[fakeEmin3] = ((double)phldr) / pow(10, 6);
       
                       fakeEmin3++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "10: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*3) && resetiter < (8 + (Ewidth+1)*4))
                   {
                       rogbb[fakeEmin4] = ((double)phldr) / pow(10, 6);
       
                       fakeEmin4++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "11: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*4) && resetiter < (8 + (Ewidth+1)*5))
                   {
                       tortuosity[fakeEmin5] = ((double)phldr) / pow(10, 6);
                       fakeEmin5++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "12: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*5) && resetiter < (8 + (Ewidth+1)*6))
                   {
                       tortuositybb[fakeEmin6] = ((double)phldr) / pow(10, 6);
                       fakeEmin6++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "13: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*6) && resetiter < (8 + (Ewidth+1)*7))
                   {
                       visits[fakeEmin7] = ((long long)phldr);
                       fakeEmin7++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "14: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*7) && resetiter < (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1))
                   {
                       poly_lattice_indexes[chain_increment] = ((int)phldr);
                       
                       
                       
                       stdoutlog = fopen(stdoutlogname, "a");
                       fprintf(stdoutlog, "pl: %i,%i,%i\n", chain_increment,poly_lattice_indexes[chain_increment],((int)phldr));
                       fclose(stdoutlog);
                       chain_increment++;
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "15: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1) && resetiter < (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites))
                   {
                       solvent_loc[fakeEmin9] = ((int)phldr);
                       
                       
                       stdoutlog = fopen(stdoutlogname, "a");
                       fprintf(stdoutlog, "sl: %i,%i,%i\n", fakeEmin9,solvent_loc[fakeEmin9],((int)phldr));
                       fclose(stdoutlog);
                       fakeEmin9++;
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "16: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }

                   if (resetiter >= (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites) && resetiter < (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites))
                   {
                       ion_loc[fakeEmin10] = ((int)phldr);
                       
                       
                       stdoutlog = fopen(stdoutlogname, "a");
                       fprintf(stdoutlog, "il: %i,%i,%i\n", fakeEmin10,ion_loc[fakeEmin10],((int)phldr));
                       fclose(stdoutlog);
                       fakeEmin10++;
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "17:%lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }
                   
                   if (resetiter >= (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites)+no_ion_sites && resetiter < (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14)))
                   {
                       radial_dist_EuPO[fakeEmin11] = ((int)phldr);
                       fakeEmin11++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "18: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }
                   
                   if (resetiter >= (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites)+no_ion_sites+((Ewidth + 1) * (14))*1 && resetiter < (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*2)
                   {
                       radial_dist_EuO[fakeEmin12] = ((int)phldr);
                       fakeEmin12++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "19: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                   }
                   
                   if (resetiter >= (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites)+no_ion_sites+((Ewidth + 1) * (14))*2 && resetiter < (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*3)
                   {
                                         radial_dist_PMPO[fakeEmin13] = ((int)phldr);
                                         fakeEmin13++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "20: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                    }
                   
                   if (resetiter >= (8 + (Ewidth+1)*7+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites)+no_ion_sites+((Ewidth + 1) * (14))*3 && resetiter < (8 + (Ewidth+1)*8+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*3)
                   {
                                         radial_visits_ion[fakeEmin14] = ((int)phldr);
                                         fakeEmin14++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "21: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                    }
                   if (resetiter >= (8 + (Ewidth+1)*8+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites)+no_ion_sites+((Ewidth + 1) * (14))*3 && resetiter < (8 + (Ewidth+1)*9+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*3)
                   {
                       
                                         radial_visits_poly[fakeEmin15] = ((int)phldr);
                                         fakeEmin15++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "224: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                    }
                   
                   if (resetiter >= (8 + (Ewidth+1)*9+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*3 && resetiter < (8 + (Ewidth+1)*9+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4)
                   {
                       
                                         radial_dist_POPO[fakeEmin16] = ((int)phldr);
                                         fakeEmin16++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "225: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                    }
                   if (resetiter >= (8 + (Ewidth+1)*9+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4 && resetiter < (8 + (Ewidth+1)*10+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4)
                   {
                       
                                         part_seq[fakeEmin17] = ((int)phldr);
                                         fakeEmin17++;
                       
                       if(counterp!=counter)
                       {
                           stdoutlog = fopen(stdoutlogname, "a");
                           fprintf(stdoutlog, "225: %lld,%lld\t%lld\n",counter,phldr,counterp);
                           fclose(stdoutlog);
                           MPI_Abort(MPI_COMM_WORLD,1);
                       }
                       counterp++;
                    }
                   if (resetiter >= (8 + (Ewidth+1)*10+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4 && resetiter < (8 + (Ewidth+1)*11+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4)
                                     {
                                         
                                                func_per_seq[fakeEmin18] = ((float)phldr);
                                                           fakeEmin18++;
                                         
                                         if(counterp!=counter)
                                         {
                                             stdoutlog = fopen(stdoutlogname, "a");
                                             fprintf(stdoutlog, "225: %lld,%lld\t%lld\n",counter,phldr,counterp);
                                             fclose(stdoutlog);
                                             MPI_Abort(MPI_COMM_WORLD,1);
                                         }
                                         counterp++;
                                      }
                   if (resetiter >= (8 + (Ewidth+1)*11+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4 && resetiter < (8 + (Ewidth+1)*12+(numbercores *lengthpercore *pertrusionpercore+backbone_length)*1+3*no_solvent_sites+no_ion_sites)+((Ewidth + 1) * (14))*4)
                                     {
                                         
                                                           perf_seq[fakeEmin19] = ((int)phldr);
                                                           fakeEmin19++;
                                         
                                         if(counterp!=counter)
                                         {
                                             stdoutlog = fopen(stdoutlogname, "a");
                                             fprintf(stdoutlog, "225: %lld,%lld\t%lld\n",counter,phldr,counterp);
                                             fclose(stdoutlog);
                                             MPI_Abort(MPI_COMM_WORLD,1);
                                         }
                                         counterp++;
                                      }
                       
                   

                   

       
                   resetiter++;
               }

               fclose(rfptr);
       
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "starting2\n");
    fclose(stdoutlog);
    
    int array_inc=0;
     int core_backbone_array[backbone_length]  ={1,1,1,0,1,1,1,0,0,1,1,1,1,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,0,1,1,0};
    
    for (int i = 0; i < numberspins_S; i++)
    {
        latticepoint_S[i] = 0;
        lattice_polymer[i *3] = i;
        lattice_polymer[i *3 + 1] = -1;
        lattice_polymer[i *3 + 2] = -1;
    }

    for (int i = 0; i < numberspins_L; i++)
    {
        latticepoint_L[i] = 0;
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Assigned index init poly_core()%i\n",poly_lattice_indexes[36]);
    fclose(stdoutlog);

    // writes the chain occupation to lattice point
    int poly_index = -1;
    int spin_inc=0;
    for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
    {
        int value =1;
        
        if((i+1)%lengthpercore==0) value = 5+numbercores*pertrusionpercore;
        if((i)%lengthpercore==0){ value = 5+spin_inc;spin_inc++;}
        
        latticepoint_S[poly_lattice_indexes[i]] = value;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = value;
        }

        poly_spin_reference[i] = value;
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "psr:%i  %i\n",i,value);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    int spin_increment = 0;
    for (int i = numbercores *lengthpercore * pertrusionpercore; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int spin = 1;

        for (int j = 0; j < numbercores *lengthpercore * pertrusionpercore; j++)
        {
            if (poly_lattice_indexes[i] == poly_lattice_indexes[j])
            {
                spin = 5 + spin_increment;
                spin_increment++;
            }
        }

        latticepoint_S[poly_lattice_indexes[i]] = spin;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = spin;
        }

        poly_spin_reference[i] = spin;

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "psr2:%i  %i %i\n",i,spin,spin_increment);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "\n");
    fclose(stdoutlog);

    
    array_inc = 0;
    for (int i = 0; i < backbone_length; i++)
    {
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "%i  %i\n", poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i],latticepoint_S[poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i]]);
        fclose(stdoutlog);

        if (core_backbone_array[array_inc] == i)
        {
            for (int h = 0; h < lengthpercore; h++)
            {
                

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "\t%i   %i\n", poly_lattice_indexes[array_inc *lengthpercore+h],latticepoint_S[poly_lattice_indexes[array_inc *lengthpercore+h]]);
                fclose(stdoutlog);

            }

            array_inc++;

        }
    }
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice point init poly_core()%i\n",poly_lattice_indexes[36]);
    fclose(stdoutlog);

    // writes the lattice polymer data
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int chain_position = i % lengthpercore;

        lattice_polymer[3 *poly_lattice_indexes[i]] = poly_lattice_indexes[i];

        if (chain_position == 0)    // core definition needs to be removed or ignored if there are multiple petrusions from a single core
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = -100;
        }
        else if (chain_position == (lengthpercore - 1))    // head of the poly chain
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = -200;
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
        else
        {
            // aqny intermediate monomer
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice polymer init poly_core()\n");
    fclose(stdoutlog);

    // fills poly_lattice_coordinates
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        poly_coord_decompose(poly_lattice_indexes[i], i);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned coord init poly_core()\n");
    fclose(stdoutlog);

    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
        ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
    }
    
    for (int i = 0; i < (pertrusionpercore *lengthpercore *numbercores) + backbone_length; i++)
    {
        attachment_ref[2 *i] = 0;
        attachment_ref[2 *i + 1] = 0;
    }

    for (int i = 0; i < no_ion_sites; i++)
    {
        attachment_ion_ref[2 *i] = 0;
        attachment_ion_ref[2 *i + 1] = 0;
    }

    reset_indexes[0]=-1;
    reset_indexes[1]=3*no_solvent_sites;
    solvent_reset2(ph_solvent_loc, solvent_loc, 2);
    solvent_reset2(wl_solvent_loc, solvent_loc, 2);
    reset_indexes[0]=-1;
    reset_indexes[1]=no_ion_sites;
    solvent_reset(ph_ion_loc, ion_loc, 3);
    solvent_reset(wl_ion_loc, ion_loc, 3);
    if (myid == 0)
    {
        printed = 1;

        sprintf(filename, "Minimum_Config_Poly3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
            }
        }

        fclose(file);
        sprintf(filename, "Minimum_Config_Poly_Attach3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ref[2 *i], attachment_ref[2 *i + 1]);
            }
        }

        fclose(file);
        sprintf(filename, "Minimum_Config_Solvent3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_solvent_sites; i++)
            {
                fprintf(file, "%i\t%i\t%i\t%i\n", i, solvent_loc[3*i],solvent_loc[3*i+1],solvent_loc[3*i+2]);
            }
        }

        fclose(file);

        sprintf(filename, "Minimum_Config_Ion3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\n", i, ion_loc[i]);
            }
        }

        fclose(file);
        
        sprintf(filename, "Minimum_Config_Ion_Attach3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
            }
        }

        fclose(file);
    }

    if (myid == numprocs - 1)
    {
        max_printed = 1;

        sprintf(filename, "Maximum_Config_Poly3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
            }
        }

        fclose(file);

        sprintf(filename, "Maximum_Config_Poly_Attach3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ref[2 *i], attachment_ref[2 *i + 1]);
            }
        }

        fclose(file);

        sprintf(filename, "Maximum_Config_Solvent3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_solvent_sites; i++)
            {
               fprintf(file, "%i\t%i\t%i\t%i\n", i, solvent_loc[3*i],solvent_loc[3*i+1],solvent_loc[3*i+2]);
            }
        }

        fclose(file);
        sprintf(filename, "Maximum_Config_Ion3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\n", i, ion_loc[i]);
            }
        }

        fclose(file);
        sprintf(filename, "Maximum_Config_Ion_Attach3.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = 0; i < no_ion_sites; i++)
            {
                fprintf(file, "%i\t%i\t%i\n", i, attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
            }
        }

        fclose(file);

    }
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "done 3\n");
    fclose(stdoutlog);
    for (int i = 0; i < numberspins_S; i++)
    {
        latticepoint_S[i] = 0;
        lattice_polymer[i *3] = i;
        lattice_polymer[i *3 + 1] = -1;
        lattice_polymer[i *3 + 2] = -1;
    }

    for (int i = 0; i < numberspins_L; i++)
    {
        latticepoint_L[i] = 0;
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Assigned index init poly_core()\n");
    fclose(stdoutlog);

    // writes the chain occupation to lattice point
    poly_index = -1;
    spin_inc=0;
    for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
    {
        int value =1;
        
        if((i+1)%lengthpercore==0) value = 5+numbercores*pertrusionpercore;
        if((i)%lengthpercore==0){ value = 5+spin_inc;spin_inc++;}
        
        latticepoint_S[poly_lattice_indexes[i]] = value;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = value;
        }

        poly_spin_reference[i] = value;
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "2psr:%i  %i\n",i,value);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    spin_increment = 0;
    for (int i = numbercores *lengthpercore * pertrusionpercore; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int spin = 1;

        for (int j = 0; j < numbercores *lengthpercore * pertrusionpercore; j++)
        {
            if (poly_lattice_indexes[i] == poly_lattice_indexes[j])
            {
                spin = 5 + spin_increment;
                spin_increment++;
            }
        }

        latticepoint_S[poly_lattice_indexes[i]] = spin;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *poly_lattice_indexes[i] + iii]] = spin;
        }

        poly_spin_reference[i] = spin;

        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "2psr2:%i  %i %i\n",i,spin,spin_increment);
        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Doublecheck:%i  \n",latticepoint_S[poly_lattice_indexes[5]]);
    fclose(stdoutlog);
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice point init poly_core()\n");
    fclose(stdoutlog);

    // writes the lattice polymer data
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        int chain_position = i % lengthpercore;

        lattice_polymer[3 *poly_lattice_indexes[i]] = poly_lattice_indexes[i];

        if (chain_position == 0)    // core definition needs to be removed or ignored if there are multiple petrusions from a single core
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = -100;
        }
        else if (chain_position == (lengthpercore - 1))    // head of the poly chain
        {
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = -200;
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
        else
        {
            // aqny intermediate monomer
            lattice_polymer[3 *poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 *poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice polymer init poly_core()\n");
    fclose(stdoutlog);

    // fills poly_lattice_coordinates
    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        poly_coord_decompose(poly_lattice_indexes[i], i);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned coord init poly_core()2\n");
    fclose(stdoutlog);

    for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
    {
        ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
        ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
        wl_pseudo_chain[i] = poly_lattice_indexes[i];
        wl_pseudo_chain_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
        wl_pseudo_chain_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
        wl_pseudo_chain_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];

    }

    for (int i = 0; i < no_solvent_sites; i++)
    {
        wl_solvent_loc[3*i] = solvent_loc[3*i];
        wl_solvent_loc[3*i+1] = solvent_loc[3*i+1];
        wl_solvent_loc[3*i+2] = solvent_loc[3*i+2];
        ph_solvent_loc[3*i] = solvent_loc[3*i];
        ph_solvent_loc[3*i+1] = solvent_loc[3*i+1];
        ph_solvent_loc[3*i+2] = solvent_loc[3*i+2];

        latticepoint_S[solvent_loc[3*i]] = 2;
        latticepoint_S[solvent_loc[3*i+1]] = -2;
        latticepoint_S[solvent_loc[3*i+2]] = -2;
        
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L *solvent_loc[3*i] + iii]] = 2;
            latticepoint_L[neighbor_L[numberneighbors_L *solvent_loc[3*i+1] + iii]] = -2;
            latticepoint_L[neighbor_L[numberneighbors_L *solvent_loc[3*i+2] + iii]] = -2;
        }
    }

    for (int i = 0; i < no_ion_sites; i++)
    {
        ph_ion_loc[i] = ion_loc[i];
        wl_ion_loc[i] = ion_loc[i];

        latticepoint_S[ion_loc[i]] = 3;
        init_solvent_sub(ion_loc[i], 3);
    }
    array_inc = 0;
    for (int i = 0; i < backbone_length; i++)
    {
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "%i  %i\n", poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i],latticepoint_S[poly_lattice_indexes[numbercores *lengthpercore *pertrusionpercore + i]]);
        fclose(stdoutlog);

        if (core_backbone_array[array_inc] == i)
        {
            for (int h = 0; h < lengthpercore; h++)
            {
                

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "\t%i   %i\n", poly_lattice_indexes[array_inc *lengthpercore+h],latticepoint_S[poly_lattice_indexes[array_inc *lengthpercore+h]]);
                fclose(stdoutlog);

            }

            array_inc++;

        }
    }
    
    
    eye();
    
    int eenergie = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc,ion_loc);
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: energy at start=%i\n", myid, energie);
    fclose(stdoutlog);
    //energie -= Eglobalmin; // shift to positive values to use it as array index

    eenergie -= Eglobalmin;
    Estartindex = energie;

    stop = 0;
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: Parameter: Eminindex -- Emaxindex: %i -- %i; Estartindex=%i\t", myid, Eminindex, Emaxindex, Estartindex);
    fclose(stdoutlog);

    if (eenergie != energie)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "2nd Estart out of range!!!  energie %i  eeenergie: %i\n", energie, eenergie);
        fclose(stdoutlog);

        stop = 1;
    }
    else
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "2nd OK\n");
        fclose(stdoutlog);

    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (stop) MPI_Abort(MPI_COMM_WORLD, 1);
}


void paramrestart()
{
    
    
    int plhdr_i;
    float plhdr_rog;
    float plhdr_tort;
    float plhdr_rogbb;
    float plhdr_tortbb;
    long long plhdr_visits;
    FILE* rfptr;
    for(int i = 0;i<numprocs;i++)
    {
    
            sprintf(resetbuffer, "SysParam_%i.txt", i);    // Preliminary output file that allows for process inspection
               
               rfptr = fopen(resetbuffer, "r");
       
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "starting\n");
    fclose(stdoutlog);
    
    while (fscanf(rfptr,"%i\t%f\t%f\t%f\t%f\t%lli\n", &plhdr_i, &plhdr_rog, &plhdr_tort ,  &plhdr_rogbb,  &plhdr_tortbb , &plhdr_visits) > 0)
    {
        if (((plhdr_i > Emaxindex) || (plhdr_i < Eminindex)))    // boundary check
        {
        rog[plhdr_i]+=plhdr_rog;
        tortuosity[plhdr_i]+=plhdr_tort;
        rogbb[plhdr_i]+=plhdr_rogbb;
        tortuositybb[plhdr_i]+=plhdr_tortbb;
        visits[plhdr_i]+=plhdr_visits;
        }
    }
    
    fclose(rfptr);
    }
    
    float plhdr_visitss;
    double value1;
    double value2;
    double value3;
    double value4;
    double value5;
    double value6;
    double value7;
    double value8;
    double value9;
    double value10;
    double value11;
    double value12;
    double value13;
    double value14;
    
    
    for(int i = 0;i<numprocs;i++)
    {
            sprintf(resetbuffer, "POPO_%i copy.txt", i);    // Preliminary output file that allows for process inspection
               rfptr = fopen(resetbuffer, "r");
       
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "starting\n");
    fclose(stdoutlog);
    while (fscanf(rfptr,"%i\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &plhdr_i, &plhdr_visitss, &value1 , &value2 , &value3 , &value4 , &value5 , &value6, &value7 , &value8 , &value9 , &value10 , &value11 , &value12 , &value13 , &value14 ) > 0)
    {
        if (((plhdr_i > Emaxindex) || (plhdr_i < Eminindex)))    // boundary check
        {
        radial_visits_poly[plhdr_i]+=plhdr_visitss;
        radial_dist_POPO[plhdr_i*14+0]+=value1;
        radial_dist_POPO[plhdr_i*14+1]+=value2;
        radial_dist_POPO[plhdr_i*14+2]+=value3;
        radial_dist_POPO[plhdr_i*14+3]+=value4;
        radial_dist_POPO[plhdr_i*14+4]+=value5;
        radial_dist_POPO[plhdr_i*14+5]+=value6;
        radial_dist_POPO[plhdr_i*14+6]+=value7;
        radial_dist_POPO[plhdr_i*14+7]+=value8;
        radial_dist_POPO[plhdr_i*14+8]+=value9;
        radial_dist_POPO[plhdr_i*14+9]+=value10;
        radial_dist_POPO[plhdr_i*14+10]+=value11;
        radial_dist_POPO[plhdr_i*14+11]+=value12;
        radial_dist_POPO[plhdr_i*14+12]+=value13;
        radial_dist_POPO[plhdr_i*14+13]+=value14;
        }
        
    }
    fclose(rfptr);
    }
        for(int i = 0;i<numprocs;i++)
        {
            sprintf(resetbuffer, "PMPO_%i copy.txt", i);    // Preliminary output file that allows for process inspection
               rfptr = fopen(resetbuffer, "r");
       
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "starting8\n");
    fclose(stdoutlog);
    while (fscanf(rfptr,"%i\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &plhdr_i, &plhdr_visitss, &value1 , &value2 , &value3 , &value4 , &value5 , &value6, &value7 , &value8 , &value9 , &value10 , &value11 , &value12 , &value13 , &value14 ) > 0)
    {
        if (((plhdr_i > Emaxindex) || (plhdr_i < Eminindex)))    // boundary check
        {
        radial_visits_poly[plhdr_i]+=plhdr_visitss;
        radial_dist_PMPO[plhdr_i*14+0]+=value1;
        radial_dist_PMPO[plhdr_i*14+1]+=value2;
        radial_dist_PMPO[plhdr_i*14+2]+=value3;
        radial_dist_PMPO[plhdr_i*14+3]+=value4;
        radial_dist_PMPO[plhdr_i*14+4]+=value5;
        radial_dist_PMPO[plhdr_i*14+5]+=value6;
        radial_dist_PMPO[plhdr_i*14+6]+=value7;
        radial_dist_PMPO[plhdr_i*14+7]+=value8;
        radial_dist_PMPO[plhdr_i*14+8]+=value9;
        radial_dist_PMPO[plhdr_i*14+9]+=value10;
        radial_dist_PMPO[plhdr_i*14+10]+=value11;
        radial_dist_PMPO[plhdr_i*14+11]+=value12;
        radial_dist_PMPO[plhdr_i*14+12]+=value13;
        radial_dist_PMPO[plhdr_i*14+13]+=value14;
        }
        
    }
    fclose(rfptr);
        }
        
        for(int i = 0;i<numprocs;i++)
        {
            sprintf(resetbuffer, "EuPO_%i copy.txt", i);    // Preliminary output file that allows for process inspection
               rfptr = fopen(resetbuffer, "r");
       
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "starting9\n");
    fclose(stdoutlog);
    while (fscanf(rfptr,"%i\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &plhdr_i, &plhdr_visitss, &value1 , &value2 , &value3 , &value4 , &value5 , &value6, &value7 , &value8 , &value9 , &value10 , &value11 , &value12 , &value13 , &value14 ) > 0)
    {
        if (((plhdr_i > Emaxindex) || (plhdr_i < Eminindex)))    // boundary check
        {
            radial_visits_ion[plhdr_i]+=plhdr_visitss;
            radial_dist_EuPO[plhdr_i*14+0]+=value1;
            radial_dist_EuPO[plhdr_i*14+1]+=value2;
            radial_dist_EuPO[plhdr_i*14+2]+=value3;
            radial_dist_EuPO[plhdr_i*14+3]+=value4;
            radial_dist_EuPO[plhdr_i*14+4]+=value5;
            radial_dist_EuPO[plhdr_i*14+5]+=value6;
            radial_dist_EuPO[plhdr_i*14+6]+=value7;
            radial_dist_EuPO[plhdr_i*14+7]+=value8;
            radial_dist_EuPO[plhdr_i*14+8]+=value9;
            radial_dist_EuPO[plhdr_i*14+9]+=value10;
            radial_dist_EuPO[plhdr_i*14+10]+=value11;
            radial_dist_EuPO[plhdr_i*14+11]+=value12;
            radial_dist_EuPO[plhdr_i*14+12]+=value13;
            radial_dist_EuPO[plhdr_i*14+13]+=value14;
        }
        
    }
    fclose(rfptr);
        }
      for(int i = 0;i<numprocs;i++)
      {
                    sprintf(resetbuffer, "EuO_%i copy.txt", i);    // Preliminary output file that allows for process inspection
                       rfptr = fopen(resetbuffer, "r");
               
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "starting7\n");
            fclose(stdoutlog);
            while (fscanf(rfptr,"%i\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &plhdr_i, &plhdr_visitss, &value1 , &value2 , &value3 , &value4 , &value5 , &value6, &value7 , &value8 , &value9 , &value10 , &value11 , &value12 , &value13 , &value14 ) > 0)
            {
                if (((plhdr_i > Emaxindex) || (plhdr_i < Eminindex)))    // boundary check
                {
                radial_visits_ion[plhdr_i]+=plhdr_visitss;
                radial_dist_EuO[plhdr_i*14+0]+=value1;
                radial_dist_EuO[plhdr_i*14+1]+=value2;
                radial_dist_EuO[plhdr_i*14+2]+=value3;
                radial_dist_EuO[plhdr_i*14+3]+=value4;
                radial_dist_EuO[plhdr_i*14+4]+=value5;
                radial_dist_EuO[plhdr_i*14+5]+=value6;
                radial_dist_EuO[plhdr_i*14+6]+=value7;
                radial_dist_EuO[plhdr_i*14+7]+=value8;
                radial_dist_EuO[plhdr_i*14+8]+=value9;
                radial_dist_EuO[plhdr_i*14+9]+=value10;
                radial_dist_EuO[plhdr_i*14+10]+=value11;
                radial_dist_EuO[plhdr_i*14+11]+=value12;
                radial_dist_EuO[plhdr_i*14+12]+=value13;
                radial_dist_EuO[plhdr_i*14+13]+=value14;
                }
                
            }
            fclose(rfptr);
      }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char *argv[])
{
    /*time(&timenow);
    std::cout << ctime(&timenow) << std::endl; */

    clock_t tStart = clock();

    // set up MPI and MPI_COMM_WORLD
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // check command line arguments
    sprintf(stdoutlogname, "Error%04i.log", myid);
    if (argc < 5)
    {
        if (myid == 0)
        {
            stdoutlog = fopen(stdoutlogname, "w");
            fprintf(stdoutlog, "Unexpected number of command line arguments!\n");
            fprintf(stdoutlog, "Expect 4 arguments, %d were provided.\n", argc - 1);
            fprintf(stdoutlog, "Syntax: ./WLpotts_mpi[arg1][arg2][arg3][arg4] \n\n");
            fprintf(stdoutlog, "Please provide the following command line arguments:\n");
            fprintf(stdoutlog, "1. Overlap between consecutive windows.[double, 0 <= overlap <= 1]\n");
            fprintf(stdoutlog, "2. Number of walkers per energy subwindow.[integer]\n");
            fprintf(stdoutlog, "3. Number of Monte Carlo steps between replica exchange.[integer]\n");
            fprintf(stdoutlog, "4. Random number seed.[integer]\n\n");
            fclose(stdoutlog);

            printf("ERROR: Unexpected number of command line arguments!\n");
            printf("       Expect 4 arguments, %d were provided.\n", argc - 1);
            printf("Syntax: ./WLpotts_mpi[arg1][arg2][arg3][arg4] \n\n");
            printf("Please provide the following command line arguments:\n");
            printf("1. Overlap between consecutive windows.[double, 0 <= overlap <= 1]\n");
            printf("2. Number of walkers per energy subwindow.[integer]\n");
            printf("3. Number of Monte Carlo steps between replica exchange.[integer]\n");
            printf("4. Random number seed.[integer]\n\n");

        }

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // number of walkers per energy window
    multiple = atoi(argv[2]);

    // at the moment, the code works only for an _odd_ number of energy windows
    // (to make the RE in windows at the outside consistent)
    if ((numprocs / multiple) % 2 == 0)
    {
        if (myid == 0)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "ERROR: Even number of energy windows (%d) requested. Please request an odd number of energy windows.\n\n", numprocs / multiple);
            fclose(stdoutlog);

            printf("ERROR: Even number of energy windows (%d) requested. Please request an odd number of energy windows.\n\n", numprocs / multiple);
        }

        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    // initiate random numbers
    rseed = atoi(argv[4]);
    rseed = time(NULL);
    //rseed = 1651894561;
    srand(rseed + myid);

    int restart_check = atoi(argv[7]);
    
    int swap_every = atoi(argv[3]);    // after this number of sweeps try conformations swap
    int swap_every_init = swap_every;

    // init log file
    sprintf(stdoutlogname, "Proc%04i.sim%i.log", myid, rseed);

    // set local energy range
    ret_status = find_local_energy_range(Eglobalmin, Eglobalmax, atof(argv[1]), numprocs);

    local_dup_group = atoi(argv[5]);
    num_dup_proc = atoi(argv[6]);
    dup_headproc = numprocs - num_dup_proc;

    MPI_Barrier(MPI_COMM_WORLD);

    if (ret_status != 0)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc. %3i: find_local_energy_range() returned %i\n", myid, ret_status);
        fclose(stdoutlog);

        MPI_Abort(MPI_COMM_WORLD, 1);    // something went wrong in find_local_energy_range()
    }

    init_hists();    // moved above init_lattice() for calculation considerations
    init_neighbors();
    init_lattice(Emin, Emax);    // 0 - random; 1 - all equal

    // calculate energy for the first time
    int eold, energie;
    energie = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: energy at start=%i\n", myid, energie);
    fclose(stdoutlog);
    energie -= Eglobalmin;    // shift to positive values to use it as array index

    Estartindex = energie;

    int stop = 0;
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: Parameter: Eminindex -- Emaxindex: %i -- %i; Estartindex=%i\t", myid, Eminindex, Emaxindex, Estartindex);
    fclose(stdoutlog);

    if ((Estartindex > Emaxindex) || (Estartindex < Eminindex))
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Estart out of range!!!\n");
        fclose(stdoutlog);

        stop = 1;
    }
    else
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "OK\n");
        fclose(stdoutlog);

    }

    //MPI_Abort(MPI_COMM_WORLD,1);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (stop) MPI_Abort(MPI_COMM_WORLD, 1);

    //Teststop
    //   MPI_Barrier(MPI_COMM_WORLD);
    //   MPI_Abort(MPI_COMM_WORLD,1);

    // create new groups and communicators for each energy range window
    stdoutlog = fopen(stdoutlogname, "a");
    MPI_Comm_group(MPI_COMM_WORLD, &world);    // get the group of processes in MPI_COMM_WORLD (i.e. all)
    int *ranks;
    ranks = (int*) malloc(2 *multiple* sizeof(int));
    mpi_local_group = (MPI_Group*) malloc(((numprocs / multiple) - 1) *sizeof(MPI_Group));
    mpi_local_comm = (MPI_Comm*) malloc(((numprocs / multiple) - 1) *sizeof(MPI_Comm));

    for (int i = 0; i < ((numprocs / multiple) - 1); i++)    // i is counter for the energy range windows
    {
        for (int j = 0; j < 2 * multiple; j++)
        {
            ranks[j] = i *multiple + j;    // contains the ranks of processes in MPI_COMM_WORLD which should get into local group
            if (myid == 0)
            {
                fprintf(stdoutlog, "Proc %3i: %i will be part of communicator/group %i\n", myid, ranks[j], i);
            }
        }

        MPI_Group_incl(world, 2 *multiple, ranks, &mpi_local_group[i]);    // create local group
        MPI_Comm_create(MPI_COMM_WORLD, mpi_local_group[i], &mpi_local_comm[i]);    // create communicator for that group
    }

    fclose(stdoutlog);
    free(ranks);

    // get my local id (in my local communicators)
    if (myid < numprocs - multiple)
    {
        comm_id = 2 *(myid / (2 *multiple));
        MPI_Comm_rank(mpi_local_comm[comm_id], &mylocalid[0]);
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %3i: I am part of communicator/group %i with local_id[0]=%i\n", myid, comm_id, mylocalid[0]);
        fclose(stdoutlog);

    }
    else
    {
        mylocalid[0] = INT_MAX;    // just to give it a value
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %3i: got local_id[0]=%i\n", myid, mylocalid[0]);
        fclose(stdoutlog);
    }

    if (myid >= multiple)
    {
        comm_id = 2 *((myid - multiple) / (2 *multiple)) + 1;
        MPI_Comm_rank(mpi_local_comm[comm_id], &mylocalid[1]);
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %3i: I am part of communicator/group %i with local_id[1]=%i\n", myid, comm_id, mylocalid[1]);
        fclose(stdoutlog);
    }
    else
    {
        mylocalid[1] = INT_MAX;    // just to give it a value
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %3i: got local_id[1]=%i\n", myid, mylocalid[1]);
        fclose(stdoutlog);
    }

    //fclose(stdoutlog);

    // log 'path'data (the same data is written, just differently sorted: with respect to the sample id and with respect to process id)
    sprintf(wanderlogname, "wanderlog_rseed%i_sample_%i.dat", rseed, latticepoint_S[numberspins_S + 2]);
    wanderlog = fopen(wanderlogname, "w");
    time(&timenow);
    fprintf(wanderlog, "#sweep\t\tproc-id\tconf-id\tenergy\t\ttime\n");
    fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", 0.0, myid, latticepoint_S[numberspins_S + 2], energie, ctime(&timenow));
    fclose(wanderlog);

    sprintf(wanderlogname, "wanderlog_rseed%i_walker_%i.dat", rseed, myid);
    wanderlog = fopen(wanderlogname, "w");
    time(&timenow);
    fprintf(wanderlog, "#sweep\t\tproc-id\tconf-id\tenergy\t\ttime\n");
    fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", 0.0, myid, latticepoint_S[numberspins_S + 2], energie, ctime(&timenow));
    fclose(wanderlog);

    //start simulation
    double wsk, dice;    // wsk: probability
    int wiggle;
    int wiggletwo;    // ID of spin to be updated (only single-spin update implemented)

    // the WL parameter should eventually go into the init file, too
    double countdown = 2;
    double lnf = 1.0;    // my modification factor
    double lnf_slowest = lnf;    // modification factor of slowest walker in my window
    long long recombine_counter = 0;
    long long recombine_counter_c = 0;// modification factor of slowest walker in my window
    //double lnfmin=log(1.000000001);
    double lnfmin = log(1.0000000001);    // terminal modification factor
    double sweep = 0;    // counter for MC sweeps
    int flat;    // 0 - histogram not flat; 1 - flat
    int iteration = 1;    // WL iteration counter
    int iter = 1;
    double check_flatness_every = 3*10000;    // in number of sweeps
    int backup;
    int backuptwo;

    eold = energie;

    int swtch;
    int found = 0;

    int prelim_counter = 0;

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: Start WL iteration\n", myid);
    fclose(stdoutlog);

        for (int i = 0; i < hist_size; i++) HE[i] = 0;    //init H(E)

        if (restart_check)
           {
               restart(energie,eold,iteration,lnf,swap_every,sweep);
           }
           
        //paramrestart();
        
        for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
        {
            ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
            ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
            ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
            ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
            wl_pseudo_chain[i] = poly_lattice_indexes[i];
            wl_pseudo_chain_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
            wl_pseudo_chain_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
            wl_pseudo_chain_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
        }

        for (int i = 0; i < 3*no_solvent_sites; i++)
        {
            wl_solvent_loc[i] = solvent_loc[i];
            ph_solvent_loc[i] = solvent_loc[i];
        }

        for (int i = 0; i < no_ion_sites; i++)
        {
            ph_ion_loc[i] = ion_loc[i];
            wl_ion_loc[i] = ion_loc[i];
        }

        long long restartiter = 0;
        int eee;
        // Main Wang-Landau routine
        while (lnf_slowest > lnfmin)
        {
            if ((restartiter) % 30 == 0)
            {
                 eye();
                eee = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc, ion_loc);
                if (eee-Eglobalmin != energie)
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "energy different  eee:%i  energue:%i\n", eee, energie);
                    fclose(stdoutlog);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                sprintf(resetbuffer, "SysStatus.L%i.proc%04i.copy6", L1dim_S_xy, myid);    // Preliminary output file that allows for process inspection

                FILE * filee;

                if ((filee = fopen(resetbuffer, "w")) == NULL)
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", resetbuffer);
                    fclose(stdoutlog);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                else
                {
                    long long counterr=0;
                    fprintf(filee, "%lld,%lld\n", counterr,(long long)(.99* pow(10, 17)));
                    counterr++;
                    fprintf(filee, "%lld,%i\n", counterr,energie);
                    counterr++;
                    fprintf(filee, "%lld,%i\n", counterr,eold);
                    counterr++;
                    fprintf(filee, "%lld,%lld\n", counterr,(long long)(sweep));
                    counterr++;
                    fprintf(filee, "%lld,%i\n",counterr, swap_every);
                    counterr++;
                    fprintf(filee, "%lld,%lld\n",counterr, (long long)(lnf* pow(10, 15)));
                    counterr++;
                    fprintf(filee, "%lld,%i\n",counterr, iteration);
                    counterr++;
                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long) HE[i]);
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(lngE[i] *pow(10, 9)));
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(rog[i] *pow(10, 6)));
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(rogbb[i] *pow(10, 6)));
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(tortuosity[i] *pow(10, 6)));
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(tortuositybb[i] *pow(10, 6)));
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(visits[i]));
                        counterr++;
                    }

                    for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
                    {
                        fprintf(filee, "%lld,%lld\n", counterr,(long long) poly_lattice_indexes[i]);
    //                    stdoutlog = fopen(stdoutlogname, "a");
    //                    fprintf(stdoutlog, "%i,%i,%i\n",i,poly_lattice_indexes[i],latticepoint_S[poly_lattice_indexes[i]]);
    //                    fclose(stdoutlog);
                        
                        counterr++;
                    }

                    for (int i = 0; i < 3*no_solvent_sites; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long) solvent_loc[i]);
                        
    //                    stdoutlog = fopen(stdoutlogname, "a");
    //                    fprintf(stdoutlog, "%i,%i,%i\n",i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
    //                    fclose(stdoutlog);
                        
                        counterr++;
                    }

                    for (int i = 0; i < no_ion_sites; i++)
                    {
                        fprintf(filee, "%lld,%lld\n",counterr, (long long) ion_loc[i]);
                        
    //                    stdoutlog = fopen(stdoutlogname, "a");
    //                    fprintf(stdoutlog, "%i,%i,%i\n",i,ion_loc[i],latticepoint_S[ion_loc[i]]);
    //                    fclose(stdoutlog);
                        counterr++;
                    }

                    for (int i = Eminindex; i <= Emaxindex; i++){
                        for (int j = 0; j < 14; j++){
                            fprintf(filee, "%lld,%lld\n", counterr,(long long)(radial_dist_EuPO[i*14+j]));
                            counterr++;
                        }
                    }
                    for (int i = Eminindex; i <= Emaxindex; i++){
                        for (int j = 0; j < 14; j++){
                            fprintf(filee, "%lld,%lld\n", counterr,(long long)(radial_dist_EuO[i*14+j]));
                            counterr++;
                        }
                    }
                    for (int i = Eminindex; i <= Emaxindex; i++){
                        for (int j = 0; j < 14; j++){
                            fprintf(filee, "%lld,%lld\n", counterr,(long long)(radial_dist_PMPO[i*14+j]));
                            counterr++;
                        }
                    }
                    for (int i = Eminindex; i <= Emaxindex; i++){
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(radial_visits_ion[i]));
                        counterr++;
                        
                    }
                    for (int i = Eminindex; i <= Emaxindex; i++){
                        fprintf(filee, "%lld,%lld\n",counterr, (long long)(radial_visits_poly[i]));
                        counterr++;
                    }
                    for (int i = Eminindex; i <= Emaxindex; i++){
                        for (int j = 0; j < 14; j++){
                            fprintf(filee, "%lld,%lld\n", counterr,(long long)(radial_dist_POPO[i*14+j]));
                            counterr++;
                        }
                    }
                    fclose(filee);

                }

                restartiter = 0;
            }

            restartiter++;
            //MPI_Abort(MPI_COMM_WORLD, 1);
            for (int k = 0; k < check_flatness_every; k++)
            {
                for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length+no_ion_sites+no_solvent_sites; i++)    // this does 1 MC sweep
                {
                    index_track = energie;
                    
                    energie = eold + propose_update();    // calculate energy of updated configuration

                    
                    // reject because you got out of bounds of energy window

                    if (MoveProposal == 0 || ((energie > Emaxindex) || (energie < Eminindex)))    // boundary check
                    {
                        //                stdoutlog = fopen(stdoutlogname, "a");
                                        //                fprintf(stdoutlog, "\nbc Move Proposal %i\n",MoveProposal);
                                        //                fclose(stdoutlog);
                        if(MoveProposal==1)
                        {
                                        if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                                        {
                                            //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                                            lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                                            poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
                                            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

                        //                    if (reset_indexes[4] == 2)
                        //                    {
                        //                        detachment(random_solvent_site);
                        //                        reset_indexes[4] = 0;
                        //                        reset_indexes[5] = 0;
                        //                        reset_indexes[6] = 0;
                        //                    }
                                        }

                                        if (poly_solvent_mov == 1)
                                        {
                                            solvent_reset2(solvent_loc, wl_solvent_loc, 2);
                                            solvent_reset2(ph_solvent_loc, solvent_loc, 2);
                                        }

                                        if (poly_ion_mov == 1)
                                        {
                                            solvent_reset(ion_loc, wl_ion_loc, 3);
                                            solvent_reset(ph_ion_loc, ion_loc, 3);

                        //                    if (reset_indexes[4] == 1)
                        //                    {
                        //                        reattachment(reset_indexes[5], reset_indexes[6]);
                        //                        reset_indexes[4] = 0;
                        //                        reset_indexes[5] = 0;
                        //                        reset_indexes[6] = 0;
                        //                    }
                        //
                        //                    if (reset_indexes[4] == 2)
                        //                    {
                        //                        detachment(reset_indexes[5]);
                        //                        reset_indexes[4] = 0;
                        //                        reset_indexes[5] = 0;
                        //                        reset_indexes[6] = 0;
                        //                    }

                                            //reset_indexes[4]=0;
                                        }
                        }
                                        energie = eold;
                        
                    }
                    else    // calculate acceptance propbability
                    {
                        // the regular wl process is kept using a place holder pseudolngE to force the system to not move to move to higher energies
                        dice = (1.0* rand() / (RAND_MAX + 1.0));    // roll the dice
                        wsk = exp(lngE[eold-eold%500] - lngE[energie-energie%500]);    // calculate acceptance probability
                        if (dice > wsk)    // reject
                        {
                            if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                                                {
                                                    //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                                                    lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                                                    poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
                                                    //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                                    lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                                    poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

                            //                        if (reset_indexes[4] == 2)
                            //                        {
                            //                            detachment(random_solvent_site);
                            //                            reset_indexes[4] = 0;
                            //                            reset_indexes[5] = 0;
                            //                            reset_indexes[6] = 0;
                            //                        }
                                                }

                                                if (poly_solvent_mov == 1)
                                                {
                                                    solvent_reset2(solvent_loc, wl_solvent_loc, 2);
                                                    solvent_reset2(ph_solvent_loc, solvent_loc, 2);
                                                }

                                                if (poly_ion_mov == 1)
                                                {
                                                    solvent_reset(ion_loc, wl_ion_loc, 3);
                                                    solvent_reset(ph_ion_loc, ion_loc, 3);

                            //                        if (reset_indexes[4] == 1)
                            //                        {
                            //                            reattachment(reset_indexes[5], reset_indexes[6]);
                            //                            reset_indexes[4] = 0;
                            //                            reset_indexes[5] = 0;
                            //                            reset_indexes[6] = 0;
                            //                        }
                            //
                            //                        if (reset_indexes[4] == 2)
                            //                        {
                            //                            detachment(reset_indexes[5]);
                            //                            reset_indexes[4] = 0;
                            //                            reset_indexes[5] = 0;
                            //                            reset_indexes[6] = 0;
                            //                        }
                            //
                            //                        reset_indexes[4] = 0;
                                                }

                                                energie = eold;

                        }
                        else
                        {
                           eold = energie;    // accept

                            //                    if(reset_indexes[4] == 1)
                            //                    {
                            //                        current_attach--;
                            //                    }
                            //                    if(reset_indexes[4] == 2)
                            //                    {
                            //                        current_attach++;
                            //                    }
                                                
                                                if (poly_solvent_mov == 0 && poly_ion_mov == 0)
                                                {
                                                    //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
    //                                                lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
    //                                                poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                                                    //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                                    lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                                    poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

                                                    

                                                    reset_indexes[4] = 0;
                                                }

                                                if (poly_solvent_mov == 1)
                                                {
                                                    //solvent_reset2(solvent_loc, ph_solvent_loc, 2);
                                                    solvent_reset2(wl_solvent_loc, solvent_loc, 2);

                                                }

                                                if (poly_ion_mov == 1)
                                                {
                                                    //solvent_reset(ion_loc, ph_ion_loc, 3);
                                                    solvent_reset(wl_ion_loc, ion_loc, 3);
    //                                                reset_indexes[4] = 0;
    //                                                reset_indexes[5] = 0;
    //                                                reset_indexes[6] = 0;
                                                }

                            //
                            //                    if(poly_solvent_mov==0)
                            //                    {
                            //                           //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            //                        lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            //                        poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                            //
                            //                        lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            //                        poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                            //
                            //
                            //                        if(core_movement == 1)
                            //                        {
                            //                            lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[2], reset_indexes[3]);
                            //                            poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[2], reset_indexes[3]);
                            //
                            //                            lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[2], reset_indexes[3]);
                            //                            poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[2], reset_indexes[3]);
                            //
                            //                        }

                            //                        reset_indexes[4]=0;
                            //                    }

                            //                    if(poly_solvent_mov==1)
                            //                    {
                            //                        solvent_reset(solvent_loc,ph_solvent_loc,2);
                            //                        solvent_reset(wl_solvent_loc,solvent_loc,2);
                            //                    }

                            //                    if(poly_ion_mov==1)
                            //                    {
                            //                        solvent_reset(ion_loc,ph_ion_loc,3);
                            //                        solvent_reset(wl_ion_loc,ion_loc,3);
                            //                        reset_indexes[4]=0;
                            //                    }

                            ////                    lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            ////                    lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            ////                    poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                            ////                    solvent_reset(solvent_loc,ph_solvent_loc);
                            ////
                            ////
                            ////                    lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            ////                    lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                            ////                    poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                            ////                    solvent_reset(wl_solvent_loc,wl_solvent_loc);
                            if (energie < minimum_en_config)
                            {
                                printed = 0;
                                minimum_en_config = energie;
                            }

                            if (myid == 0 && energie <= minimum_en_config && printed == 0)
                            {
                                printed = 1;

                                sprintf(filename, "Minimum_Config_Poly.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
                                    {
                                        fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
                                    }
                                }

                                fclose(file);
                                sprintf(filename, "Minimum_Config_Poly_Attach.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
                                    {
                                        fprintf(file, "%i\t%i\t%i\n", i, attachment_ref[2 *i], attachment_ref[2 *i + 1]);
                                    }
                                }

                                fclose(file);
                                sprintf(filename, "Minimum_Config_Solvent.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < no_solvent_sites; i++)
                                    {
                                        fprintf(file, "%i\t%i\t%i\t%i\n", i, solvent_loc[3*i],solvent_loc[3*i+1],solvent_loc[3*i+2]);
                                    }
                                }

                                fclose(file);

                                stdoutlog = fopen(stdoutlogname, "a");
                                fprintf(stdoutlog, "Printed energie %i\n", energie);
                                fclose(stdoutlog);
                                sprintf(filename, "Minimum_Config_Ion.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < no_ion_sites; i++)
                                    {
                                        fprintf(file, "%i\t%i\n", i, ion_loc[i]);
                                    }
                                }

                                fclose(file);
                                sprintf(filename, "Minimum_Config_Ion_Attach.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < no_ion_sites; i++)
                                    {
                                        fprintf(file, "%i\t%i\t%i\n", i, attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
                                    }
                                }

                                fclose(file);
                            }

                            if (energie > maximum_en_config)
                            {
                                max_printed = 0;
                                maximum_en_config = energie;
                            }

                            if (myid == numprocs - 1 && energie >= maximum_en_config && max_printed == 0)
                            {
                                max_printed = 1;

                                sprintf(filename, "Maximum_Config_Poly.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
                                    {
                                        fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
                                    }
                                }

                                fclose(file);

                                sprintf(filename, "Maximum_Config_Poly_Attach.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < numbercores *pertrusionpercore *lengthpercore + backbone_length; i++)
                                    {
                                        fprintf(file, "%i\t%i\t%i\n", i, attachment_ref[2 *i], attachment_ref[2 *i + 1]);
                                    }
                                }

                                fclose(file);

                                sprintf(filename, "Maximum_Config_Solvent.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < no_solvent_sites; i++)
                                    {
                                        fprintf(file, "%i\t%i\t%i\t%i\n", i, solvent_loc[3*i],solvent_loc[3*i+1],solvent_loc[3*i+2]);
                                    }
                                }

                                fclose(file);
                                sprintf(filename, "Maximum_Config_Ion.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < no_ion_sites; i++)
                                    {
                                        fprintf(file, "%i\t%i\n", i, ion_loc[i]);
                                    }
                                }

                                fclose(file);
                                sprintf(filename, "Maximum_Config_Ion_Attach.txt");
                                if ((file = fopen(filename, "w")) == NULL)
                                {
                                    stdoutlog = fopen(stdoutlogname, "a");
                                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                                    fclose(stdoutlog);
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                else
                                {
                                    for (int i = 0; i < no_ion_sites; i++)
                                    {
                                        fprintf(file, "%i\t%i\t%i\n", i, attachment_ion_ref[2 *i], attachment_ion_ref[2 *i + 1]);
                                    }
                                }

                                fclose(file);
                                stdoutlog = fopen(stdoutlogname, "a");
                                fprintf(stdoutlog, "Printed energie %i\n", energie);
                                fclose(stdoutlog);

                            }

                         sysparam(energie-energie%500);

                            //
                            //                   //if ((energie > eold))
                            //                   //{
                            //                   //    latticepoint[wiggle] = backup;
                            //                   //    latticepoint[wiggletwo] = backuptwo;
                            //                   //    energie = eold;
                            //                   //}

                            //                   //if ((energie < eold))
                            //                   //{
                            //                   //    k = 0;
                            //                   //    eold = energie;
                        }

                        //                        sysparam(energie);
                        //                    lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        //                    lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        //                    poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                        //                    solvent_reset(solvent_loc,ph_solvent_loc);
                        //
                        //
                        //                    lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        //                    lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        //                    poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                        //                    solvent_reset(wl_solvent_loc,wl_solvent_loc);

                    }

    //                if(vicinity_check)
    //                {
    //                    for(int i=0;i<numbercores*pertrusionpercore;i++)
    //                    {
    //                        attach_track[energie]+=current_attachment[i*lengthpercore+lengthpercore-1];
    //                    }
    //
    //                //attach_track[ep] += current_attach;
    //                    attach_track_no[energie] ++;
    //                }

                    //attach_track[ep] += current_attach;
                    // update histograms
                    lngE[energie-energie%500] += lnf;
                    HE[energie-energie%500]++;

                    //                int eee = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
                    //                stdoutlog = fopen(stdoutlogname, "a");
                    //                fprintf(stdoutlog, "eee: %i   energie:%i  s_or_p:%i   r0:%i   r1:%i\n",eee,energie,poly_solvent_mov,reset_indexes[0],reset_indexes[1]);
                    //                fclose(stdoutlog);
                    //
                    //                stdoutlog = fopen(stdoutlogname, "a");
                    //                for (int i = 0; i < numbercores *lengthpercore * pertrusionpercore; i++)
                    //                {
                    //                    fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 *i], poly_lattice_coordinates[3 *i + 1], poly_lattice_coordinates[3 *i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 *i], ph_poly_lattice_coordinates[3 *i + 1], ph_poly_lattice_coordinates[3 *i + 2]);
                    //                }

                    //                for (int i = 0; i < no_solvent_sites; i++)
                    //                {
                    //                    fprintf(stdoutlog, "%4i\t\t%4i\t\t%i\n", solvent_loc[i],ph_solvent_loc[i],wl_solvent_loc[i]);
                    //                }

                    //
                    //                fclose(stdoutlog);
                    //
                    //                if(eee!=energie)
                    //                {
                    //                MPI_Abort(MPI_COMM_WORLD, 1);
                    //                }
                }    // end 1 sweep
    radial_counter(poly_lattice_indexes, solvent_loc, ion_loc,eold);
                sweep++;    // sweep counter

                swap_every--;    // RE counter (countdown)
                if (swap_every == 0)    // ignition! (time to try to exchange configurations)
                {
                    swap_every = swap_every_init;    // reset countdown clock
                    if (replica_exchange(sweep / swap_every, eold) == 1)    // if configs were exchanged ('accept'-branch)
                    {
                        energie = switched;
                        eold = energie;

                        reset_indexes[0] = 0 - 1;
                        reset_indexes[1] = 3*no_solvent_sites;
                        poly_latticepoint_clear(wl_solvent_loc, reset_indexes[0], reset_indexes[1]);

                        reset_indexes[0] = 0 - 1;
                        reset_indexes[1] = no_ion_sites;
                        poly_latticepoint_clear(wl_ion_loc, reset_indexes[0], reset_indexes[1]);

                        reset_indexes[0] = 0 - 1;
                        reset_indexes[1] = numbercores *lengthpercore *pertrusionpercore + backbone_length;
                        poly_latticepoint_clear(wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);

                        for (int i = 0; i < numbercores *lengthpercore *pertrusionpercore + backbone_length; i++)
                        {
                            wl_pseudo_chain[i] = poly_lattice_indexes[i];
                            wl_pseudo_chain_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
                            wl_pseudo_chain_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
                            wl_pseudo_chain_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];

                            ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
                            ph_poly_lattice_coordinates[3 *i] = poly_lattice_coordinates[3 *i];
                            ph_poly_lattice_coordinates[3 *i + 1] = poly_lattice_coordinates[3 *i + 1];
                            ph_poly_lattice_coordinates[3 *i + 2] = poly_lattice_coordinates[3 *i + 2];
                        }

                        for (int i = 0; i < 3*no_solvent_sites; i++)
                        {
                            wl_solvent_loc[i] = solvent_loc[i];
                            ph_solvent_loc[i] = solvent_loc[i];
                        }

                        for (int i = 0; i < no_ion_sites; i++)
                        {
                            wl_ion_loc[i] = ion_loc[i];
                            ph_ion_loc[i] = ion_loc[i];
                        }

                        reset_indexes[0] = 0 - 1;
                        reset_indexes[1] = numbercores *lengthpercore *pertrusionpercore + backbone_length;

                        //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

                        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

                        reset_indexes[0] = 0 - 1;
                        reset_indexes[1] = 3*no_solvent_sites;

                        solvent_reset2(wl_solvent_loc, solvent_loc, 2);
                        solvent_reset2(ph_solvent_loc, solvent_loc, 2);
                        reset_indexes[0] = 0 - 1;
                        reset_indexes[1] = no_ion_sites;

                        solvent_reset(wl_ion_loc, ion_loc, 3);
                        solvent_reset(ph_ion_loc, ion_loc, 3);

                       
                            stdoutlog = fopen(stdoutlogname, "a");
                            fprintf(stdoutlog, "Proc %3i: %i iteration, %e sweeps, Replica Exchange Success\n Energy Calculted: %i \t Recievced Index Energy: %i\n", myid, iteration, sweep, energie, switched);
                            fclose(stdoutlog);
                        

    //                    current_attach=0;
    //                    for(int i=0;i<numbercores*pertrusionpercore;i++)
    //                    {
    //                        current_attachment[i]=0;
    //                        for(int j=0;j<no_ion_sites;j++)
    //                        {
    //                            if(int_point_distance(poly_lattice_indexes[i*lengthpercore+lengthpercore-1],ion_loc[j])<=9)
    //                            {
    //                                current_attach++;
    //                                current_attachment[i*lengthpercore+lengthpercore-1]++;
    //                            }
    //                        }
    //                    }
    //                    for(int i=0;i<numbercores*pertrusionpercore;i++)
    //                    {
    //                        attach_track[energie]+=current_attachment[i*lengthpercore+lengthpercore-1];
    //                    }
    //                    //attach_track[ep] += current_attach;
    //                    attach_track_no[energie] ++;
                    }

                    // update histograms (independently of whether RE happened or not)
                    lngE[energie-energie%500] += lnf;
                    HE[energie-energie%500]++;

                }
            }    // end computation block

                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Prelim\n");
                fclose(stdoutlog);

                sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim_S_xy, q, myid, iteration);    // Preliminary output file that allows for process inspection
                if ((file = fopen(filename, "w")) == NULL)
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                    fclose(stdoutlog);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                else
                {
                    for (int i = Eminindex; i <= Emaxindex; i++)
                    {
                        if (lngE[i] > 0.5)
                        {
                            fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + (Eglobalmin)), 0.0, lngE[i], HE[i]);
                        }
                    }

                    fclose(file);
                }
            

            prelim_counter++;
            //prelim_counter++;

            // check flatness
            flat = histflat(Eminindex, Emaxindex, 0.92);
            if (flat == 0)    // histogram not flat, do nothing except maybe some logging
            {
                
                    stdoutlog = fopen(stdoutlogname, "a");    // Prints out effictively flatness progress for inspection
                    fprintf(stdoutlog, "Proc %3i: %i iteration, %e sweeps, Histogram noch nicht flach\nValue of Tried Minimum Histogram Value: %f    Value of Tried Ratio: %f\n", myid, iteration, sweep, flatmin, flatratio);
                    fclose(stdoutlog);
                    phrepv = flatmin;    // holds the minimum value of the histogram and checks change and prints out the current lattice if the value has reoccured frequently
                    //if (repv == phrepv)
                    //{
                    //    iterr++;
                    //    repv = phrepv;
                    //}

                    //else
                    //{
                    //    repv = phrepv;
                    //}

                    //if (iterr == 20)
                    //{
                    //    iterr = 0;

                    /*
                    fprintf(stdoutlog, "\n");

                    for (int i = 0; i < numberspins_S; i++)
                    {
                        fprintf(stdoutlog, "%4i", latticepoint_S[i]);
                        if ((i + 1) % L1dim_S_xy == 0)
                            fprintf(stdoutlog, "\n");

                    }

                    fprintf(stdoutlog, "\n Sytem Energy Index: %i \n", energie);

                    fprintf(stdoutlog, "\n");
                    */
                    //}

                    //            accessiblelevels();
                    //                partial_init_hists();
                    //                stdoutlog = fopen(stdoutlogname, "a");    // Prints out effictively flatness progress for inspection
                    //                fprintf(stdoutlog, "\tLeft Access\n");
                    //                 fclose(stdoutlog);
                
            }
            else    // histograms of all walkers are 'flat'
            {
                if (lnf > lnfmin)    // as long as modification factor is still larger than terminal value
                {
                    sprintf(filename, "L%iq%i.HE.proc%04i.iter%i", L1dim_S_xy, q, myid, iteration);
                    if ((file = fopen(filename, "w")) == NULL)
                    {
                        stdoutlog = fopen(stdoutlogname, "a");
                        fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                        fclose(stdoutlog);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    else
                    {
                        for (int i = Eminindex; i <= Emaxindex; i++)
                        {
                            if (HE[i] > 0.5) fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + (Eglobalmin)), 0.0, lngE[i], HE[i]);
                        }

                        fclose(file);
                    }
                }

                // decrease modification factor
                // canonical method (reduce by sqrt(2)) implemented
                lnf /= 2.0;

                for (int i = 0; i < hist_size; i++) HE[i] = 0;    //reset H(E)
                iteration++;    // iteration counter

                if (merge_hists == 1)    // merge g(E) estimators from multiple walkers in the same energy window
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    if (myid % multiple == 0)    // 'root' in energy window, receive individual g(E) and send merged g(E)
                    {
                        for (int i = 1; i < multiple; i++)
                        {
                            MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid + i, 77, MPI_COMM_WORLD, &status);    // get other dens. of states
                            fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, myid + i);
                            for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j];    // sum up for average
                        }

                        for (int j = 0; j < hist_size; j++) lngE[j] /= (double) multiple;    // normalize
                        for (int i = 1; i < multiple; i++)
                        {
                            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid + i, 99, MPI_COMM_WORLD);
                            fprintf(stdoutlog, "Proc %i: Sent merged lngE to Proc. %i\n", myid, myid + i);
                        }
                    }
                    else    // send individual g(E) and receive merged g(E)
                    {
                        MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 77, MPI_COMM_WORLD);
                        fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, myid - (myid % multiple));
                        MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 99, MPI_COMM_WORLD, &status);
                        fprintf(stdoutlog, "Proc %i: Received merged lngE from Proc. %i\n", myid, myid - (myid % multiple));
                        for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];    // replace individual lngE (could be done directly, yes)
                    }

                    fclose(stdoutlog);
                }
            }

            recombine_counter++;
            // check progress from all other windows
           // MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allreduce(&lnf, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&recombine_counter, &recombine_counter_c, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            
            // check progress from all other windows
            MPI_Allreduce(&lnf, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            if (lnf_slowest <= (((double)(10000000000 / countdown)) *.0000000001)|| recombine_counter_c>50)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "access start\n");
                fclose(stdoutlog);
                accessiblelevels();
                partial_init_hists();
                stdoutlog = fopen(stdoutlogname, "a");    // Prints out effictively flatness progress for inspection
                fprintf(stdoutlog, "\tLeft Access\n");
                fclose(stdoutlog);
                MPI_Barrier(MPI_COMM_WORLD);
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "recombine\n");
                fclose(stdoutlog);
                
                recombine((((double)(10000000000 / countdown)) *.0000000001));
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "left recombine\n");
                fclose(stdoutlog);
                if(lnf_slowest <= (((double)(10000000000 / countdown)) * .0000000001))
                {
                countdown *= 2;
                }
                recombine_counter=0;
                partial_init_hists();
                //MPI_Barrier(MPI_COMM_WORLD);
            }

            // just some logging
            if (flat == 1)
            {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Proc %3i: Start %i iteration, %e sweeps total so far, lnf now %e (lnfmin=%.2e, lnf_slowest=%e)\n", myid, iteration, sweep, lnf, lnfmin, lnf_slowest);
                fprintf(stdoutlog, "Proc %3i: tryleft: %i, exchangeleft %i (Akzeptanzleft:%.2lf) < --> tryright: %i, exchangeright %i (Akzeptanzright:%.2lf)\n", myid, tryleft, exchangeleft, (double) exchangeleft / (double) tryleft, tryright, exchangeright, (double) exchangeright / (double) tryright);
                fclose(stdoutlog);
            }
        }    // end while(lnf_slowest>lnfmin) -> this terminates the simulation

        //   free(neighbor);
        //   free(latticepoint);

        // normalize results
        double norm = lngE[Eminindex] - log(q);
        for (int i = 0; i <= 2 *(-Eglobalmin); i++) lngE[i] -= norm;

        // write final data to file
        sprintf(filename, "L%iq%i.proc%04i.lngE", L1dim_S_xy, q, myid);
        if ((file = fopen(filename, "w")) != NULL)
        {
            for (int i = Eminindex; i <= Emaxindex; i++)
                fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + Eglobalmin), 0.0, lngE[i], HE[i]);
            fclose(file);
        }

        //   free(HE);
        //   free(lngE);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();

        /*time(&timenow);
        std::cout << ctime(&timenow) << std::endl; */

        printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        return (0);    // you did it!
    }
