//mpirun --oversubscribe -np 19 ./EXE 0.95 1 80000 07891 0
// Fully adapted to periodic boundary condtions utilized Dr. Bai's lecture notes in order to get everthing worked out correctly
// The only changes that should be left are verifying current equations
// Adapting the remaining function (should require no change just adaptations)
// Also running a small test case
// work on a local energy function this time
/*
  Replica Exchange Wang Landau demo code for simulating the 2D Potts model
  (c) Thomas Vogel and Ying Wai Li (2013-2019)

  License: Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
           https://creativecommons.org/licenses/by-sa/4.0/legalcode

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
         - T. Vogel et al., J. Phys.: Conf. Ser. 487 (2014) 012001
             - Y.W. Li et al., J. Phys.: Conf. Ser. 510 (2014) 012012

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
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <mpi.h>

#include <stdlib.h>

#include <stdio.h>

#include <math.h>

#include <time.h>

#include <limits.h>

#include <float.h>

#include <random>


#include "Constants.h"

#include "General.h"

#include "Energy.h"

#include "Mov.h"

#include "EXWL.h"

int switched = INT_MAX;
int cccounter = 0;

int iterr; // Apart of an added debug output. Effictively just a counter for the number of repeated histogram values when flatness is checked
int repv = 0; // Holds the value of the (potentially) repeated minvalue of the histogram
int phrepv; // Holds the value of the previous minimum histogram value to compare with the prior variable
int minmaxid = 12;

int * poly_lattice_indexes = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * sizeof(int));
int * ph_poly_lattice_indexes = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * sizeof(int));
int * wl_pseudo_chain = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * sizeof(int));

//may not need forgot why this is here in the first place aa it will be tracked by the lattice polymer
int * poly_lattice_connection = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * sizeof(int));

//hold the physical coordinates of a given monomer
int * poly_lattice_coordinates = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * 3 * sizeof(int));
int * ph_poly_lattice_coordinates = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * 3 * sizeof(int));
int * wl_pseudo_chain_coordinates = (int * ) malloc(numbercores * lengthpercore * pertrusionpercore * 3 * sizeof(int));

int reset_solvent_index = 0;
int * solvent_loc = (int * ) malloc(no_solvent_sites * sizeof(int)); // initializes main function to track solvent
int * ph_solvent_loc = (int * ) malloc(no_solvent_sites * sizeof(int)); // initializes main function to track solvent
int * wl_solvent_loc = (int * ) malloc(no_solvent_sites * sizeof(int)); // initializes wl function to track solvent vector for replacemt if movement fails

int q = 2;

int * neighbor_S; // list containing indices of neighbors for all spins
// energy histogram
double * HE;
// ln g(E) estimator
double * lngE;

double * lngE_buf; // ln g(E) estimator buffer for exchange
double * pseudolngE;
double * real_lngE;
double * real_lngE_buf;
double * microT;
double * microT_buf;

int rseed; // seed for random number generator
int energy, localenergy;

double Emin, Emax; // doubles as calculation of boundaries in Energy is not in 'int'
int Eminindex, Emaxindex, Estartindex; // local boundaries as index

// MPI; set up with local communicators for replica exchange (RE)
// Needed so that two processes can talk w/o depending on / bothering the others
int numprocs, myid, multiple, comm_id;
// each process belongs to two local groups,
// one to communicate to left neighbor and one to communicate to right neighbor
// each process has different loca IDs in different communicatore, in general
int mylocalid[2]; // id in local communicators
MPI_Comm * mpi_local_comm;
MPI_Group * mpi_local_group;
MPI_Group world;
MPI_Status status;
int merge_hists = 1; // flag whether or not to merge histograms between iterations for multiple walkers on the same energy range

// to keep track of exchange statistics
int tryleft, tryright, exchangeleft, exchangeright;

// File handlers for I/O
FILE * file;
FILE * stdoutlog;
FILE * wanderlog;

FILE * DOSEstimate; //holds values for the dos estimate all values for the rog and tort are placed here to later be deemed unique or not
char dosestimate[50];

char filename[50];
char resetbuffer[50];
char stdoutlogname[128];
char wanderlogname[128];
int ret_status;

double flatratio;
double flatmin;

int en_array[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int en_array_pp[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int en_array_ps[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int en_array_ss[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int en_array_bond[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

int int_en_array[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int int_en_array_pp[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int int_en_array_ps[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int int_en_array_ss[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
int int_en_array_bond[11]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1

int en_array_wall[4][4][4]; // plan to hold the available energy values of the hamiltonian will contain 1's for available energies and 0 for others, made to be a size of three so all values can be looked up other wise Id have to reduce by 1
double dist_array[4][4][4]; // plan to hold the available distances in this array so the square root does not need to be called when values are needed
int int_dist_array[4][4][4]; // plan to hold the available distances in this array so the square root does not need to be called when values are neededxxxxxx
// to track execution times
time_t timenow, timestart, timeend;

int * latticepoint_S = (int * ) malloc((numberspins_S + 2 + 1) * sizeof(int));
int * latticepoint_L = (int * ) malloc((numberspins_L + 2 + 1) * sizeof(int));
int * lattice_polymer = (int * ) malloc(((numberspins_S * 3) + 2 + 1) * sizeof(int));

int * neighbor_L = (int * ) malloc(numberspins_S * numberneighbors_L * sizeof(int));

double * tortuosity = (double * ) malloc(hist_size * sizeof(double));
double * tortuosity_buf = (double * ) malloc(hist_size * sizeof(double));

double * density_profile_buf = (double * ) malloc(hist_size * L1dim_S_z * sizeof(double));
double * density_profile = (double * ) malloc(hist_size * L1dim_S_z * sizeof(double));
double * density_profile_solvent_buf = (double * ) malloc(hist_size * L1dim_S_z * sizeof(double));
double * density_profile_solvent = (double * ) malloc(hist_size * L1dim_S_z * sizeof(double));

double * indiv_rog = (double * ) malloc(hist_size * numbercores * sizeof(double));
double * indiv_rog_buf = (double * ) malloc(hist_size * numbercores * sizeof(double));

double * rog = (double * ) malloc(hist_size * sizeof(double));
double * rog_buf = (double * ) malloc(hist_size * sizeof(double));
double * rog_z = (double * ) malloc(hist_size * sizeof(double));
double * rog_z_buf = (double * ) malloc(hist_size * sizeof(double));
double * rog_xy = (double * ) malloc(hist_size * sizeof(double));
double * rog_xy_buf = (double * ) malloc(hist_size * sizeof(double));
long long * visitdp = (long long * ) malloc(hist_size * sizeof(long long));
long long * visits = (long long * ) malloc(hist_size * sizeof(long long));
long long * visits_buf = (long long * ) malloc(hist_size * sizeof(long long));
double * endtoend = (double * ) malloc(hist_size * sizeof(double));
double * endtoend_buf = (double * ) malloc(hist_size * sizeof(double));

int MoveProposal = 0; // tracks if a move has been performed
//checks to calculate local energy
int poly_solvent_mov = 0; // says wether a chosen mov moved solvent or a monomer
int chosen_mov = 0; // tellw which move was last chosen
int poly_mov_start = 0; // tells which monomwrs were involved with a given move the starting position (on thr chain)
int poly_mov_end = 0; // tells which monomer is the end considered loop

int reset_indexes[2] = {
  0,
  0
};

int what_to_exclude[120 * numberneighbors_L] = {
  120
};
int * local_en_indexes_op = (int * ) malloc(numberspins_S * 120 * sizeof(int));
int local_en_dist[120] = {
  0
};
int local_en_indexes_op_x[120] = {
  0
};
int local_en_indexes_op_y[120] = {
  0
};
int local_en_indexes_op_z[120] = {
  0
};

int * lattice_dist_optomize_2 = (int * ) malloc(L1dim_S_xy * L1dim_S_xy * L1dim_S_z * sizeof(int));

int * lattice_dist_optomize_2_lattice = (int * ) malloc((L1dim_S_xy * L1dim_S_xy) * (L1dim_S_xy * L1dim_S_xy * 7) * sizeof(int));

int * index_optomize_x = (int * ) malloc((L1dim_S_xy * 12) * sizeof(int));
int * index_optomize_y = (int * ) malloc((L1dim_S_xy * 12) * sizeof(int));
int * index_optomize_z = (int * ) malloc((L1dim_S_z * 6) * sizeof(int));

int * calc_index_optomize_x = (int * ) malloc((L1dim_S_xy * 12) * sizeof(int));
int * calc_index_optomize_y = (int * ) malloc((L1dim_S_xy * 12) * sizeof(int));
int * calc_index_optomize_z = (int * ) malloc((L1dim_S_z * 6) * sizeof(int));

int * polymer_optomize = (int * ) malloc(L1dim_S_xy * L1dim_S_xy * L1dim_S_z * sizeof(int));

int * xyz_index = (int * ) malloc((L1dim_S_z + 10) * (L1dim_S_xy + 10) * (L1dim_S_xy + 10) * sizeof(int));

int * index_x = (int * ) malloc(numberspins_S * sizeof(int));
int * index_y = (int * ) malloc(numberspins_S * sizeof(int));
int * index_z = (int * ) malloc(numberspins_S * sizeof(int));

int minimum_en_config = 10000000;
int maximum_en_config = -1000000;
int printed = 0;
int max_printed = 0;

int loc_en_before = 0;

//const int solvent_max = numberspins_S-( L1dim_S_xy*L1dim_S_xy*(plane_z_min+1)) -(L1dim_S_xy*L1dim_S_xy*(L1dim_S_z-plane_z_max));
//const int solvent_offset = L1dim_S_xy*L1dim_S_xy*(plane_z_min+1);

double * output_table;
int table_length = 0;

void keypressed() // just for developing / manual debugging
{
  for (;;)
    if (getchar() == 27) break;
}

//general set up

void sysparam(int);

//const int solvent_max = numberspins_S-( L1dim_S_xy*L1dim_S_xy*(plane_z_min+1)) -(L1dim_S_xy*L1dim_S_xy*(L1dim_S_z-plane_z_max));
//const int solvent_offset = L1dim_S_xy*L1dim_S_xy*(plane_z_min+1);

int local_energy(int * p_l_c, int * p_l_i, int * s_l);

void init_output();
void output_lattice(int index, double rog, double tort, int dx[numbercores * pertrusionpercore * lengthpercore], int dy[numbercores * pertrusionpercore * lengthpercore], int dz[numbercores * pertrusionpercore * lengthpercore]);

void init_output() {
  FILE * fptr;
  fptr = fopen("Combined_58mer.txt", "r");
  int out_counter = 0;
  int value;
  double valuee;
  double valueee;
  while (fscanf(fptr, "%i,%lf,%lf\n", & value, & valuee, & valueee) > 0) {
    fprintf(stdoutlog, "init_output1 %i: %i %f %f\n", out_counter, value, valuee, valueee);
    out_counter++;
  }
  fclose(fptr);

  output_table = (double * ) malloc((out_counter + 1) * 4 * sizeof(double));

  fptr = fopen("Combined_58mer.txt", "r");
  out_counter = 0;
  value = 0;
  valuee = 0;
  valueee = 0;
  while (fscanf(fptr, "%i,%lf,%lf\n", & value, & valuee, & valueee) > 0) {
    output_table[out_counter * 4 + 0] = value;
    output_table[out_counter * 4 + 1] = valuee;
    output_table[out_counter * 4 + 2] = valueee;
    fprintf(stdoutlog, "init_output2 %i: %f %f %f\n", out_counter, output_table[out_counter * 4 + 0], output_table[out_counter * 4 + 1], valueee);
  }
  fclose(fptr);
}

void output_lattice(int index, double rog, double tort, int dx[numbercores * pertrusionpercore * lengthpercore], int dy[numbercores * pertrusionpercore * lengthpercore], int dz[numbercores * pertrusionpercore * lengthpercore]) {
  int break_statement = 0;
  int i = 0;
  while (break_statement == 0 && i < table_length) {
    if (((int)(((double) index) * .95)) > output_table[4 * i + 0] && ((int)(((double) index) * 1.05)) < output_table[4 * i + 0]) {
      if (((int) output_table[4 * i + 3]) % 100 == 0) {
        if (((int)(((double) rog) * .9)) > output_table[4 * i + 1] && ((int)(((double) rog) * 1.1)) < output_table[4 * i + 1]) {
          if (((int)(((double) tort) * .9)) > output_table[4 * i + 2] && ((int)(((double) tort) * 1.1)) < output_table[4 * i + 2]) {
            int summed_diff[3] = {
              0,
              0,
              0
            };
            int arm = 0;
            sprintf(filename, "Output_%i_%f_Chain_CCords_Proc%i_%i_%i.txt", i, output_table[4 * i + 0], myid, numbercores, lengthpercore);
            if ((file = fopen(filename, "w")) == NULL) {
              stdoutlog = fopen(stdoutlogname, "a");
              fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
              fclose(stdoutlog);
              MPI_Abort(MPI_COMM_WORLD, 1);
            } else {
              for (int i = 0; i < numbercores * pertrusionpercore * lengthpercore; i++) {
                summed_diff[0] = poly_lattice_coordinates[3 * (arm * lengthpercore)];
                summed_diff[1] = poly_lattice_coordinates[3 * (arm * lengthpercore) + 1];
                summed_diff[2] = poly_lattice_coordinates[3 * (arm * lengthpercore) + 2];

                summed_diff[0] += dx[i];
                summed_diff[1] += dy[i];
                summed_diff[2] += dz[i];

                fprintf(file, "%i\t%i\t%i\n", summed_diff[0], summed_diff[1], summed_diff[2]);

                if ((i + 1) % lengthpercore == 0) {
                  arm++;
                }
              }
            }
            fclose(file);

            int s_array[3] = {
              0
            };
            sprintf(filename, "Output_%i_%f_Solvent_CCords_Proc%i_%i_%i.txt", i, output_table[4 * i + 0], myid, numbercores, lengthpercore);
            if ((file = fopen(filename, "w")) == NULL) {
              stdoutlog = fopen(stdoutlogname, "a");
              fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
              fclose(stdoutlog);
              MPI_Abort(MPI_COMM_WORLD, 1);
            } else {
              for (int i = 0; i < no_solvent_sites; i++) {
                general_coord(s_array, solvent_loc[i]);
                fprintf(file, "%i\t%i\t%i\n", s_array[0], s_array[1], s_array[2]);
              }
            }
            fclose(file);

            break_statement = 1;
          }
        }
      } else {
        output_table[4 * i + 3]++;
      }
    }
    i++;
  }
}

void init_neighbors() // create neighbor list first for smaller lattice
{
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Made it to init neighbors()\n");
  fclose(stdoutlog);
  // neighbor contains the index of the neighboring spin
  // for each spin there are four neighbors in this order: above, right, below, left
  neighbor_S = (int * ) malloc(numberspins_S * numberneighbors_S * sizeof(int));

  for (int i = 0; i < numberspins_S; i++) // in general
  {
    neighbor_S[numberneighbors_S * i] = i - L1dim_S_xy * L1dim_S_xy; // out of plane
    neighbor_S[numberneighbors_S * i + 1] = i - L1dim_S_xy; // above
    neighbor_S[numberneighbors_S * i + 2] = i + 1; // right
    neighbor_S[numberneighbors_S * i + 3] = i + L1dim_S_xy * L1dim_S_xy; // into plane
    neighbor_S[numberneighbors_S * i + 4] = i + L1dim_S_xy; // below
    neighbor_S[numberneighbors_S * i + 5] = i - 1; // left
  }

  if (bctype == 0) // periodic BC
    for (int i = 0; i < numberspins_S; i++) // now treat boundaries separately
  {
    if (i < L1dim_S_xy * L1dim_S_xy) // highest plane
      neighbor_S[numberneighbors_S * i] = i + (L1dim_S_xy * L1dim_S_xy * (L1dim_S_z - 1));
    //neighbor_S[numberneighbors_S * i] = numberspins_S;

    if (i % (L1dim_S_xy * L1dim_S_xy) < L1dim_S_xy) // top row
      neighbor_S[numberneighbors_S * i + 1] = i + (L1dim_S_xy * (L1dim_S_xy - 1));
    //neighbor_S[numberneighbors_S * i] = numberspins_S;

    if ((i + 1) % L1dim_S_xy == 0) // rightmost column
      neighbor_S[numberneighbors_S * i + 2] = i - (L1dim_S_xy - 1);
    //neighbor_S[numberneighbors_S * i] = numberspins_S;

    if (i >= L1dim_S_xy * L1dim_S_xy * (L1dim_S_z - 1)) // deepest row
      neighbor_S[numberneighbors_S * i + 3] = i - (L1dim_S_xy * L1dim_S_xy * (L1dim_S_z - 1));
    //neighbor_S[numberneighbors_S * i] = numberspins_S;

    if ((i % (L1dim_S_xy * L1dim_S_xy)) > (L1dim_S * (L1dim_S - 1) - 1)) // bottom row
      neighbor_S[numberneighbors_S * i + 4] = i - (L1dim_S * (L1dim_S - 1));
    //neighbor_S[numberneighbors_S * i] = numberspins_S;

    if (i % L1dim_S_xy == 0) // leftmost column
      neighbor_S[numberneighbors_S * i + 5] = i + (L1dim_S_xy - 1);
    //neighbor_S[numberneighbors_S * i] = numberspins_S;
  }

  //print neighbor list (for manual debugging)
  //  for (int i=0; i<numberspins*numberneighbors; i++)
  //    {
  //      printf("%4d", neighbor[i]);
  //      if ((i+1)%numberneighbors == 0) printf("\n");
  //    };

  // contains indexes of larger lattice in indexes for the smaller lattice
  int reference = 0;

  for (int i = 0; i < numberspins_L; i++) // in general
  {
    reference = i;

    int next_plane_ref = ((reference + L1dim_L_xy * L1dim_L_xy) % (L1dim_L_xy * L1dim_L_xy * L1dim_L_z));

    neighbor_L[numberneighbors_L * i] = reference; // upper topleft
    neighbor_L[numberneighbors_L * i + 1] = reference + ((reference + 1) % L1dim_L_xy - ((reference) % L1dim_L_xy));; // upper topright
    neighbor_L[numberneighbors_L * i + 2] = (reference + L1dim_L_xy) % (L1dim_L_xy * L1dim_L_xy) + (reference / (L1dim_L_xy * L1dim_L_xy)) * L1dim_L_xy * L1dim_L_xy; // upper bottom left
    neighbor_L[numberneighbors_L * i + 3] = ((reference + L1dim_L_xy) % (L1dim_L_xy * L1dim_L_xy)) + ((reference + 1) % L1dim_L_xy - ((reference) % L1dim_L_xy)) + (reference / (L1dim_L_xy * L1dim_L_xy)) * L1dim_L_xy * L1dim_L_xy; // upper bottom right
    neighbor_L[numberneighbors_L * i + 4] = next_plane_ref; //lower top right
    neighbor_L[numberneighbors_L * i + 5] = next_plane_ref + ((next_plane_ref + 1) % L1dim_L_xy - ((next_plane_ref) % L1dim_L_xy)); // lower topright
    neighbor_L[numberneighbors_L * i + 6] = (next_plane_ref + L1dim_L_xy) % (L1dim_L_xy * L1dim_L_xy) + (next_plane_ref / (L1dim_L_xy * L1dim_L_xy)) * L1dim_L_xy * L1dim_L_xy; // lower bottom left
    neighbor_L[numberneighbors_L * i + 7] = ((next_plane_ref + L1dim_L_xy) % (L1dim_L_xy * L1dim_L_xy)) + ((next_plane_ref + 1) % L1dim_L_xy - ((next_plane_ref) % L1dim_L_xy)) + (next_plane_ref / (L1dim_L_xy * L1dim_L_xy)) * L1dim_L_xy * L1dim_L_xy; // lower bottom right
  }

  //print neighbor list (for manual debugging)
  //    stdoutlog = fopen(stdoutlogname, "a");
  //    int q=0;
  //    fprintf(stdoutlog,"\n%i: ",q);
  //      for (int i=0; i<numberspins_L*numberneighbors_L; i++)
  //        {
  //          fprintf(stdoutlog, " %4d ", neighbor_L[i]);
  //            if ((i+1)%numberneighbors_L == 0) {q++; fprintf(stdoutlog,"\n%i: ",q);};
  //
  //        };
  //
  //    fclose(stdoutlog);
}

void init_lattice(double emin, double emax) // Changes made to bctype check and made latticepoint[numbersign equal to 0
{
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Made it to init lattice()\n");
  fclose(stdoutlog);
  int e, r;

  // Note: we reserve space for 3 extra 'spins':
  // 2 extra to store fix values used for certain boundary conditions
  // 1 extra to carry an replica-id flag

  // find a fast way to create valid initial configurations
  // start at either Eglobalmin or Eglobalmax and change spins randomly
  // until energy is in the middle third of local energy range
  // global maximum of g(E) id at E/N=-0.2

  //initialize lattice
  for (int i = 0; i < numberspins_S; i++) {
    latticepoint_S[i] = 0;
  }
  for (int i = 0; i < numberspins_L; i++) {
    latticepoint_L[i] = 0;
  }
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Made it to polycore()\n");
  fclose(stdoutlog);
  init_poly_cores();
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Left init polycore()\n");
  fclose(stdoutlog);

  int array[3] = {
    0,
    0,
    0
  };

  for (int i = 0; i < numberspins_S; i++) {
    general_coord(array, i);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "array: %i \t%i\t%i \n",array[0],array[1],array[2]);
    //fclose(stdoutlog);

    if (array[2] == plane_z_max && latticepoint_S[i] == 0) {
      latticepoint_S[i] = -1;
    }
  }

  e = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc);

  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Proc. %i: Initialized lattice with energy e=%i, create setup with %lf<e<%lf\n", myid, e, (emin + (emax - emin) / 3), (emin + 2 * (emax - emin) / 3));
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
  for (int i = 0; i < numberspins_S; i++) {
    fprintf(stdoutlog, "%4i", latticepoint_S[i]);
    if ((i + 1) % L1dim_S_xy == 0)
      fprintf(stdoutlog, "\n");
    if ((i + 1) % (L1dim_S_xy * L1dim_S_xy) == 0)
      fprintf(stdoutlog, "\n");

  }
  fprintf(stdoutlog, "\n");
  fprintf(stdoutlog, "\n");
  fclose(stdoutlog);

  stdoutlog = fopen(stdoutlogname, "a");
  for (int i = 0; i < numberspins_L; i++) {
    fprintf(stdoutlog, "%4i", latticepoint_L[i]);
    if ((i + 1) % L1dim_L_xy == 0)
      fprintf(stdoutlog, "\n");
    if ((i + 1) % (L1dim_L_xy * L1dim_L_xy) == 0)
      fprintf(stdoutlog, "\n");

  }

  fclose(stdoutlog);
  // run to a valid energy for the energy window

  // as long as energy is outside the middle third of local energy range
  long long increm = 0;
  int ep = 0;
  int loc_energy_before = 0;
  int loc_energy_after = 0;
  e = -100;

  for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) //defines wl pseudo_chain for use in the below local energy calculation
  {
    wl_pseudo_chain[i] = poly_lattice_indexes[i];
    wl_pseudo_chain_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
    wl_pseudo_chain_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
    wl_pseudo_chain_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];
  }

  int min_e = 0;
  int max_e = 0;

  e = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc);

  //    while ((e < (emin)) || (e > (emin +(emax - emin)))) //while(1)
  //    {
  //
  //        loc_en_before = 0;
  ////        stdoutlog = fopen(stdoutlogname, "a");
  ////        if(increm == 685)
  ////        {
  ////            fprintf(stdoutlog, "686 1 latticepoint 3233 %i\n",latticepoint_L[3233]);
  ////        }
  ////        fclose(stdoutlog);
  //
  //      // ep = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc); // compares the energy before and after a movement has been preformed
  //
  //        poly_mov();
  //
  ////        stdoutlog = fopen(stdoutlogname, "a");
  ////        fprintf(stdoutlog, "Left polymov()\n");
  ////        fclose(stdoutlog);
  //
  //       //e = total_energy(ph_poly_lattice_coordinates,ph_poly_lattice_indexes,ph_solvent_loc);
  //
  //
  //
  //        loc_energy_before = 0;
  //        loc_energy_after = 0;
  //
  //        if(MoveProposal==1)
  //                              {
  //                                  if(poly_solvent_mov==0)
  //                                  {
  //                                      //lattice_polymer_reset(poly_lattice_indexes , ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  //                                      lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  //                                      poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
  //
  //                                  }
  //                                  if(poly_solvent_mov==1)
  //                                  {
  //                                      solvent_reset(solvent_loc,ph_solvent_loc);
  //
  //                      //                stdoutlog = fopen(stdoutlogname, "a");
  //                      //                fprintf(stdoutlog, "solvent_loc[i]: %i ph_solvent_loc[i]: %i\n",solvent_loc[reset_indexes[0]+1],ph_solvent_loc[reset_indexes[0]+1]);
  //                      //                fclose(stdoutlog);
  //                                  }
  //                              }
  //
  //
  //        if(MoveProposal==1)
  //        {
  //            loc_energy_before  = loc_en_before;
  //        loc_energy_after = local_energy(ph_poly_lattice_coordinates,ph_poly_lattice_indexes,ph_solvent_loc);
  //        }
  //        increm = increm + 1;
  //
  //
  ////        stdoutlog = fopen(stdoutlogname, "a");
  ////        fprintf(stdoutlog, "Left totalenergyx2 increm: %lli\n   ep: %i   e:%i   ep-e: %i\n",increm%1,ep,e,ep-e);
  ////        fprintf(stdoutlog, "local_energy_before: %i   local_energy_after: %i   loc_en_b-loc_en_a: %i\n",loc_energy_before,loc_energy_after,loc_energy_before-loc_energy_after);
  ////        fclose(stdoutlog);
  //
  //
  //
  //        e=e+(loc_energy_after - loc_energy_before);
  //
  ////        if(loc_energy_before-loc_energy_after != ep-e)
  ////        {
  ////            stdoutlog = fopen(stdoutlogname, "a");
  ////            fprintf(stdoutlog, "\nMade it to local_energy\n");
  ////            fprintf(stdoutlog, "chosen_mov:%i \t poly_solv:%i\n",chosen_mov,poly_solvent_mov);
  ////            fprintf(stdoutlog, "reset_indexes[0]+1:%i \t reset_indexes[1]:%i\n",reset_indexes[0]+1,reset_indexes[1]);
  ////            fclose(stdoutlog);
  //
  ////            stdoutlog = fopen(stdoutlogname, "a");
  ////            fprintf(stdoutlog, "Difference  %i  %i\n",cccounter,MoveProposal);
  ////            fclose(stdoutlog);
  ////            eye();
  ////                        MPI_Abort(MPI_COMM_WORLD, 1);
  ////       }
  //
  ////        if(min_e>ep)
  ////        {
  ////            min_e=ep;
  ////            stdoutlog = fopen(stdoutlogname, "a");
  ////            fprintf(stdoutlog, "min_e: %i\n",ep);
  ////            fclose(stdoutlog)
  ////        }
  ////
  ////        if(max_e<ep)
  ////        {
  ////            max_e=ep;
  ////            stdoutlog = fopen(stdoutlogname, "a");
  ////            fprintf(stdoutlog, "max_e: %i\n",ep);
  ////            fclose(stdoutlog);
  ////        }
  //
  //
  //
  //        if (increm % 10000 == 0) {
  //
  //
  //            stdoutlog = fopen(stdoutlogname, "a");
  ////            fprintf(stdoutlog, "\n");
  ////            fprintf(stdoutlog, "\n");
  //            fprintf(stdoutlog, "increment: %lli   %i", increm,e);
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
  ////            fprintf(stdoutlog, "\n");
  ////            fprintf(stdoutlog, "%i\n", e);
  ////            fprintf(stdoutlog, "\n");
  ////            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
  ////            {
  ////                fprintf(stdoutlog, "i: %i  %4i\t\t%4i\t%4i\t%4i\t \t%4i\t\t%4i\t%4i\t%4i\n",i, poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
  ////            }
  //////
  ////            fprintf(stdoutlog, "\n");
  ////            for (int i = 0; i < no_solvent_sites; i++)
  ////            {
  ////                fprintf(stdoutlog, "i: %i\t%4i\t\n", i,solvent_loc[i]);
  ////            }
  ////            fclose(stdoutlog);
  ////
  ////            stdoutlog = fopen(stdoutlogname, "a");
  ////
  ////            int counter=0;
  ////             fprintf(stdoutlog, "%i\n",counter);
  ////            for (int i = 0; i < numberspins_S; i++)
  ////            {
  ////                fprintf(stdoutlog, "%4i", latticepoint_S[i]);
  ////                if ((i + 1) % L1dim_S_xy == 0)
  ////                    fprintf(stdoutlog, "\n");
  ////                if ((i + 1) % (L1dim_S_xy*L1dim_S_xy) == 0){
  ////                    fprintf(stdoutlog, "\n%i\n",counter);counter++;
  ////                }
  ////
  ////            }
  ////            fprintf(stdoutlog, "\n");
  ////            fprintf(stdoutlog, "\n");
  ////            fclose(stdoutlog);
  ////            stdoutlog = fopen(stdoutlogname, "a");
  ////            counter=0;
  ////            fprintf(stdoutlog, "%i\n",counter);
  ////            for (int i = 0; i < numberspins_L; i++)
  ////            {
  ////                fprintf(stdoutlog, "%4i", latticepoint_L[i]);
  ////                if ((i + 1) % L1dim_L_xy == 0)
  ////                    fprintf(stdoutlog, "\n");
  ////                if ((i + 1) % (L1dim_L_xy*L1dim_L_xy) == 0){
  ////                    fprintf(stdoutlog, "\n%i\n",counter);counter++;
  ////                }
  ////
  ////            }
  //
  //            fclose(stdoutlog);
  //
  //
  //
  //        }
  //
  //
  //        //eye();
  //        cccounter++;
  //    }

  int ee = total_energy(ph_poly_lattice_coordinates, ph_poly_lattice_indexes, ph_solvent_loc);
  eye();
  if (e != ee) {

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "DDifference  %i\n", cccounter);
    fclose(stdoutlog);
    eye();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

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

  // Pseudo WL Process is started for all processes to explore the energy levels
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "\nPseudo WL Proces has started for process: %i\n", myid);
  fclose(stdoutlog);
  pseudowl();

  if (fopen("Access_Levels.txt", "r") == NULL) {
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
    fclose(stdoutlog);
  } else // rewrite of prior file reading system to be more portable
  {
    FILE * fptr;
    int value;
    int valuee;

    fptr = fopen("Access_Levels.txt", "r");

    while (fscanf(fptr, "%i\t%i\n", & value, & valuee) > 0) {
      if ((e >= (emin) && (e <= (emin + (emax - emin))))) {
        lngE[value - Eglobalmin] = valuee;
        pseudolngE[value - Eglobalmin] = valuee;
      }
    }

    fclose(fptr);
  }

  stdoutlog = fopen(stdoutlogname, "a");

  fclose(stdoutlog);

  if (multiple > 1 && merge_hists == 1) // merge g(E) estimators from multiple walkers in the same energy window
  {
    stdoutlog = fopen(stdoutlogname, "a");
    if (myid % multiple == 0) // 'root' in energy window, receive individual g(E) and send merged g(E)
    {
      for (int i = 1; i < multiple; i++) {
        MPI_Recv( & lngE_buf[0], hist_size, MPI_DOUBLE, myid + i, 76, MPI_COMM_WORLD, & status); // get other dens. of states
        fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, myid + i);
        for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j]; // sum up for average
      }
      for (int i = 1; i < multiple; i++) {
        MPI_Send( & lngE[0], hist_size, MPI_DOUBLE, myid + i, 98, MPI_COMM_WORLD);
        fprintf(stdoutlog, "Proc %i: Sent combined lngE to Proc. %i\n", myid, myid + i);
      }
    } else // send individual g(E) and receive merged g(E)
    {
      MPI_Send( & lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 76, MPI_COMM_WORLD);
      fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, myid - (myid % multiple));
      MPI_Recv( & lngE_buf[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 98, MPI_COMM_WORLD, & status);
      fprintf(stdoutlog, "Proc %i: Received combined lngE from Proc. %i\n", myid, myid - (myid % multiple));
      for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j]; // replace individual lngE (could be done directly, yes)
    }
    fclose(stdoutlog);
  }

  accessiblelevels();

  //Print lattice to logfile
  //    stdoutlog = fopen(stdoutlogname, "a");
  //    for (int i = 0; i < numberspins_S; i++)
  //    {
  //        fprintf(stdoutlog, "%4i", latticepoint_S[i]);
  //        if ((i + 1) % L1dim_S_xy == 0)
  //            fprintf(stdoutlog, "\n");
  //    }
  //    fprintf(stdoutlog, "\n");
  //    //fprintf(stdoutlog, "%i", total_energy());
  //    fprintf(stdoutlog, " done \n");
  //    fclose(stdoutlog);

  //MPI_Abort(MPI_COMM_WORLD, 1);
}

void init_hists() // initialize histograms
{
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Made it to init hist()\n");
  fclose(stdoutlog);

  lngE = (double * ) malloc(hist_size * sizeof(double));
  lngE_buf = (double * ) malloc(hist_size * sizeof(double));
  real_lngE = (double * ) malloc(hist_size * sizeof(double));
  real_lngE_buf = (double * ) malloc(hist_size * sizeof(double));
  microT = (double * ) malloc(hist_size * sizeof(double));
  microT_buf = (double * ) malloc(hist_size * sizeof(double));
  pseudolngE = (double * ) malloc(hist_size * sizeof(double)); // added for pseudo wl process
  HE = (double * ) malloc(hist_size * sizeof(double));

  for (int i = 0; i < hist_size; i++) {
    lngE[i] = 0.0;
    lngE_buf[i] = 0.0;
    real_lngE[i] = 0.0;
    real_lngE_buf[i] = 0.0;
    microT[i] = 0.0;
    pseudolngE[i] = 0.0;
    microT_buf[i] = 0.0;
  }
}

void partial_init_hists() // initialize histograms
{
  for (int i = 0; i < hist_size; i++) {
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
int find_local_energy_range(double Eglobmin, double Eglobmax, double overlap, int N) {
  stdoutlog = fopen(stdoutlogname, "a");

  // consistency check
  if (N % multiple != 0) {
    fprintf(stdoutlog, "Total number of processes (%i) must be a multiple of the number of walkers per energy range (%i)!\n", N, multiple);
    fclose(stdoutlog);
    return (1);
  }

  FILE * window_init;

  // check if there is a file containing local energy window boundaries
  // Must be named Ewindows.dat !

  //fprintf(stdoutlog, "\nProc %i: Can't find file Ewindows.dat. Will calculate equal-size windows with overlap %lf\n", myid, overlap);
  //double Ewidth = (Eglobmax - Eglobmin) / (1.0 + ((double)(N / multiple) - 1.0) * (1.0 - overlap));

  //Emin = Eglobmin + (double)(myid / multiple) * (1.0 - overlap) * Ewidth;
  //Emax = Emin + Ewidth;

  //Eminindex = floor(Emin + -Eglobalmin);
  //Emaxindex = ceil(Emax + -Eglobalmin);

  if (fopen("Ewindows_99.dat", "r") == NULL) {
    fprintf(stdoutlog, "\nProc %i: Can't find file Ewindows.dat. Will calculate equal-size windows with overlap %lf\n", myid, overlap);
    double Ewidth = (Eglobmax - Eglobmin) / (1.0 + ((double)(N / multiple) - 1.0) * (1.0 - overlap));

    Emin = Eglobmin + (double)(myid / multiple) * (1.0 - overlap) * Ewidth;
    Emax = Emin + Ewidth;

    Eminindex = floor(Emin + -Eglobalmin);
    Emaxindex = ceil(Emax + -Eglobalmin);

    time( & timenow);
    fprintf(stdoutlog, "Proc %3i: Parameter: Eglobmin -- Eglobmax: %lf -- %lf; overlap=%i percent, Ewindowwidth=%lf, %s", myid, Eglobmin, Eglobmax, (int)(overlap * 100.0), Ewidth, ctime( & timenow));
    //fprintf(stdoutlog, "Proc %3i: Parameter: Eglobmin -- Eglobmax: %lf -- %lf; %s", myid, Emin, Emax, ctime(&timenow));
  } else // rewrite of prior file reading system to be more portable
  {
    FILE * fptr;
    int value;
    int valuee;
    int valueee;
    int valueeee;

    fptr = fopen("Ewindows_99.dat", "r");

    while (fscanf(fptr, "%i,%i,%i,%i\n", & value, & valuee, & valueee, & valueeee) > 0) {
      if (value == myid / multiple) {
        Eminindex = valuee - Eglobalmin + 0;
        Emaxindex = valueee - Eglobalmin + 0;
        fprintf(stdoutlog, "%i %i %i %i\n", value, Eminindex, Emaxindex, valueeee);
      }
      Emin = Eglobalmin + Eminindex;
      Emax = Eglobalmin + Emaxindex;
    }

    fclose(fptr);
  }

  fclose(stdoutlog);

  return (0);
}

// THIS IS THE MASTER RE / SWAP FUNCTION
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int replica_exchange(int swap_direction, int index_akt) {
  //    stdoutlog = fopen(stdoutlogname, "a");
  //    fprintf(stdoutlog, "1\n");
  //    fclose(stdoutlog);

  int i_new; // histogram index of my configuration
  int Ecur; // current energy

  // frac refers to local exchange probability
  // wk is combined exchange probability
  double myfrac, otherfrac, randx, wk;

  int change = 0; // boolean: 0 for not exchanging, 1 for exchanging
  int swap_partner = -1; // id of swap partner (receive_buffer)

  swap_direction = swap_direction % 2; // comes actually as number of swap attempt

  // everyone has to find its swap-partner

  int * pairs; // array containing the partners of each process (send_buffer)
  pairs = (int * ) malloc(2 * multiple * sizeof(int));

  if (mylocalid[swap_direction] == 0) // 'head-node' in the energy window determines pairs of flippartners
  {
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "2\n");
    //fclose(stdoutlog);

    int choose_from = multiple; // number of free partners in higher window of communicator
    int select; // storage for random number

    int * libre; // list of free partners from higher window in communicator
    libre = (int * ) malloc(multiple * sizeof(int));

    for (int i = 0; i < multiple; i++) libre[i] = multiple + i; // initialise

    // idea: processes from the lower window choose someone from the higher window at random
    // of course, the chosen walker can't have an exchange partner yet
    for (int i = 0; i < multiple; i++) // loop over processes in the lower window
    {
      select = rand() % choose_from;
      pairs[i] = libre[select];
      pairs[libre[select]] = i; // the 'vice-versa pair'
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
  if ((swap_direction == 0) && (myid < (numprocs - multiple))) // the walkers from the last node should not swap
  {
    comm_id = 2 * (myid / (2 * multiple)); // ! all integer, the '/' is a (div) ! Not the same as myid/multiple !
    MPI_Scatter(pairs, 1, MPI_INT, & swap_partner, 1, MPI_INT, 0, mpi_local_comm[comm_id]);
  }

  if ((swap_direction == 1) && (myid >= multiple)) // the walkers from the zero-node should not swap
  {
    comm_id = ((myid - multiple) / (2 * multiple)) * 2 + 1; // ! all integer, the '/' is a (div) ! See above
    MPI_Scatter(pairs, 1, MPI_INT, & swap_partner, 1, MPI_INT, 0, mpi_local_comm[comm_id]);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  free(pairs);

  if (swap_partner != -1) // i.e. if there is a swap-partner for me (if I am at a boundary, I might not have a swap partner this time)
  {
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "3\n");
    //fclose(stdoutlog);

    // statistics
    if (swap_partner > mylocalid[swap_direction]) tryright++;
    else tryleft++;

    // safety cross check
    //        Ecur = total_energy();
    //        if (Ecur + (-Eglobalmin) != index_akt)
    //        {
    //            stdoutlog = fopen(stdoutlogname, "a");
    //            fprintf(stdoutlog, "Proc %3i, replica_exchange(): Something is wrong here! Received index=%i, calculated index=%i totalenergy=%i. Abort.\n", myid, index_akt, Ecur + (-Eglobalmin), total_energy());
    //            fclose(stdoutlog);
    //            MPI_Abort(MPI_COMM_WORLD, 1);
    //        }

    i_new = index_akt;

    // get histogram index from my swap partner
    MPI_Sendrecv_replace( & i_new, 1, MPI_INT, swap_partner, 1, swap_partner, 1, mpi_local_comm[comm_id], & status);

    if ((i_new > Emaxindex) || (i_new < Eminindex)) // energyranges must overlap!
    {
      myfrac = -1.0;
    } else {
      // calculate my part of the exchange probability
      myfrac = exp(lngE[index_akt] - lngE[i_new]); // g(myE)/g(otherE)
    }

    if (mylocalid[swap_direction] < multiple) // I am receiver and calculator
    {
      // get my partners part of the exchange probability
      MPI_Recv( & otherfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_local_comm[comm_id], & status);

      // calculate combined exchange probability and roll the dice
      if ((myfrac > 0.0) && (otherfrac > 0.0)) {
        randx = (1.0 * rand() / (RAND_MAX + 1.0));
        wk = myfrac * otherfrac;
        if (randx < wk) change = 1;
      }

      // tell my swap partner whether to exchange or not
      MPI_Send( & change, 1, MPI_INT, swap_partner, 3, mpi_local_comm[comm_id]);
    } else // I just send my part of exchange probability and await decision
    {
      MPI_Send( & myfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_local_comm[comm_id]);
      MPI_Recv( & change, 1, MPI_INT, swap_partner, 3, mpi_local_comm[comm_id], & status);
    }

    // if decision was made to exchange configurations
    if (change == 1) {
      //stdoutlog = fopen(stdoutlogname, "a");
      //fprintf(stdoutlog, "4\n");
      //fclose(stdoutlog);
      // exchange spin conformations (incl. the 3 'special' lattice sites)
      MPI_Sendrecv_replace( & poly_lattice_indexes[0], numbercores * lengthpercore * pertrusionpercore, MPI_INT, swap_partner, 20, swap_partner, 20, mpi_local_comm[comm_id], & status);
      MPI_Sendrecv_replace( & poly_lattice_coordinates[0], numbercores * lengthpercore * pertrusionpercore * 3, MPI_INT, swap_partner, 10, swap_partner, 10, mpi_local_comm[comm_id], & status);
      MPI_Sendrecv_replace( & solvent_loc[0], no_solvent_sites, MPI_INT, swap_partner, 30, swap_partner, 30, mpi_local_comm[comm_id], & status);

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

void sysparam(int index) {
  visits[index]++;

  visitdp[index]++;
  double rsum[numbercores * pertrusionpercore] = {
    0
  };
  double rsum_z[numbercores * pertrusionpercore] = {
    0
  };
  double rsum_xy[numbercores * pertrusionpercore] = {
    0
  };
  double radius = 0;
  double radius_z = 0;
  double radius_xy = 0;

  //          summed_diff[0]=poly_lattice_coordinates[3 * (arm*lengthpercore)];
  //            summed_diff[1]=poly_lattice_coordinates[3 * (arm*lengthpercore) + 1];
  //            summed_diff[2]=poly_lattice_coordinates[3 * (arm*lengthpercore) + 2];
  //
  //            summed_diff[0]+=dx[i];
  //            summed_diff[1]+=dy[i];
  //            summed_diff[2]+=dz[i];

  int dx[numbercores * pertrusionpercore * lengthpercore] = {
    0
  };
  int dy[numbercores * pertrusionpercore * lengthpercore] = {
    0
  };
  int dz[numbercores * pertrusionpercore * lengthpercore] = {
    0
  };

  int s_dx[numbercores * pertrusionpercore * lengthpercore] = {
    0
  };
  int s_dy[numbercores * pertrusionpercore * lengthpercore] = {
    0
  };
  int s_dz[numbercores * pertrusionpercore * lengthpercore] = {
    0
  };

  for (int i = (0); i < numbercores * pertrusionpercore; i++) // performs the rotations
  {
    double centerofmass[3] = {
      0
    };
    s_dx[i * lengthpercore] = poly_lattice_coordinates[3 * (i * lengthpercore)];
    s_dy[i * lengthpercore] = poly_lattice_coordinates[3 * (i * lengthpercore) + 1];
    s_dz[i * lengthpercore] = poly_lattice_coordinates[3 * (i * lengthpercore) + 2];

    centerofmass[0] += s_dx[i * lengthpercore];
    centerofmass[1] += s_dy[i * lengthpercore];
    centerofmass[2] += s_dz[i * lengthpercore];

    for (int j = 1; j < lengthpercore; j++) // performs the rotations
    {

      dx[i * lengthpercore + j] = poly_lattice_coordinates[3 * (i * lengthpercore + j)] - poly_lattice_coordinates[3 * (i * lengthpercore + j - 1)];

      dy[i * lengthpercore + j] = poly_lattice_coordinates[3 * (i * lengthpercore + j) + 1] - poly_lattice_coordinates[3 * (i * lengthpercore + j - 1) + 1];

      dz[i * lengthpercore + j] = poly_lattice_coordinates[3 * (i * lengthpercore + j) + 2] - poly_lattice_coordinates[3 * (i * lengthpercore + j - 1) + 2];

      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "d={%i,%i,%i}    d-1={%i,%i,%i}\n\tp_l_c_diff={%i,%i,%i}\n",dx[i+j],dy[i+j],dz[i+j],dx[i+j-1],dy[i+j-1],dz[i+j-1],poly_lattice_coordinates[3 * (i+j)]-poly_lattice_coordinates[3 * (i+j-1)] ,poly_lattice_coordinates[3 * (i+j)+1] -poly_lattice_coordinates[3 * (i+j-1) + 1],poly_lattice_coordinates[3 * (i+j)+2]-poly_lattice_coordinates[3 * (i+j-1)+2]);
      //            fclose(stdoutlog);

      dx[i * lengthpercore + j] = dx[i * lengthpercore + j - 1] - (dx[i * lengthpercore + j] - anint(dx[i * lengthpercore + j], L1dim_S_xy) * L1dim_S_xy);
      dy[i * lengthpercore + j] = dy[i * lengthpercore + j - 1] - (dy[i * lengthpercore + j] - anint(dy[i * lengthpercore + j], L1dim_S_xy) * L1dim_S_xy);
      dz[i * lengthpercore + j] = dz[i * lengthpercore + j - 1] - (dz[i * lengthpercore + j] - anint(dz[i * lengthpercore + j], L1dim_S_z) * L1dim_S_z);

      s_dx[i * lengthpercore + j] += s_dx[i * lengthpercore] + dx[i * lengthpercore + j];
      s_dy[i * lengthpercore + j] += s_dy[i * lengthpercore] + dy[i * lengthpercore + j];
      s_dz[i * lengthpercore + j] += s_dz[i * lengthpercore] + dz[i * lengthpercore + j];

      centerofmass[0] += s_dx[i * lengthpercore + j];
      centerofmass[1] += s_dy[i * lengthpercore + j];
      centerofmass[2] += s_dz[i * lengthpercore + j];
      //
      //                stdoutlog = fopen(stdoutlogname, "a");
      //                fprintf(stdoutlog, "\ni:%i  j:%i   d={%i,%i,%i}  cm={%f,%f,%f}\n\tp_l_c[i+j]={%i,%i,%i}     p_l_c[i+j-1]={%i,%i,%i}",i,j,s_dx[i*lengthpercore+j],s_dy[i*lengthpercore+j],s_dz[i*lengthpercore+j],centerofmass[0],centerofmass[1],centerofmass[2],poly_lattice_coordinates[3 * (i+j)] ,poly_lattice_coordinates[3 * (i+j)+1] ,poly_lattice_coordinates[3 * (i+j)+2],poly_lattice_coordinates[3 * (i+j-1)],poly_lattice_coordinates[3 * (i+j-1) + 1],poly_lattice_coordinates[3 * (i+j-1)+2]);
      //                fclose(stdoutlog);
    }
    centerofmass[0] /= lengthpercore;
    centerofmass[1] /= lengthpercore;
    centerofmass[2] /= lengthpercore;

    //            stdoutlog = fopen(stdoutlogname, "a");
    //            fprintf(stdoutlog, "\ncm={%f,%f,%f}\n",centerofmass[0],centerofmass[1],centerofmass[2]);
    //            fclose(stdoutlog);

    for (int j = 0; j < lengthpercore; j++) // performs the rotations
    {
      rsum[i] += (s_dx[i * lengthpercore + j] - centerofmass[0]) * (s_dx[i * lengthpercore + j] - centerofmass[0]) + (s_dy[i * lengthpercore + j] - centerofmass[1]) * (s_dy[i * lengthpercore + j] - centerofmass[1]) + (s_dz[i * lengthpercore + j] - centerofmass[2]) * (s_dz[i * lengthpercore + j] - centerofmass[2]);
      rsum_z[i] += (s_dz[i * lengthpercore + j] - centerofmass[2]) * (s_dz[i * lengthpercore + j] - centerofmass[2]);
      rsum_xy[i] += (s_dx[i * lengthpercore + j] - centerofmass[0]) * (s_dx[i * lengthpercore + j] - centerofmass[0]) + (s_dy[i * lengthpercore + j] - centerofmass[1]) * (s_dy[i * lengthpercore + j] - centerofmass[1]);

      //                stdoutlog = fopen(stdoutlogname, "a");
      //                fprintf(stdoutlog, "b%i rsum: %f\t",j,rsum[i]);
      //                fclose(stdoutlog);
    }
  }

  for (int i = 0; i < numbercores * pertrusionpercore; i++) {
    double rr = sqrt(rsum[i] / lengthpercore);
    double rr_z = sqrt(rsum_z[i] / lengthpercore);
    double rr_xy = sqrt(rsum_xy[i] / lengthpercore);
    radius += rr;
    radius_z += rr_z;
    radius_xy += rr_xy;
    indiv_rog[index * numbercores * pertrusionpercore + i] += (rr) * (rr) / lengthpercore;
  }

  //        stdoutlog = fopen(stdoutlogname, "a");
  //        fprintf(stdoutlog, "\nradius: %f  radius_z: %f  radius_xy:%f\n",radius,radius_z,radius_xy);
  //        fclose(stdoutlog);
  rog[index] += (radius / (numbercores * pertrusionpercore)) * (radius / (numbercores * pertrusionpercore)) / lengthpercore;
  rog_z[index] += (radius_z / (numbercores * pertrusionpercore)) * (radius_z / (numbercores * pertrusionpercore)) / lengthpercore;
  rog_xy[index] += (radius_xy / (numbercores * pertrusionpercore)) * (radius_xy / (numbercores * pertrusionpercore)) / lengthpercore;

  //        stdoutlog = fopen(stdoutlogname, "a");
  //        fprintf(stdoutlog, "rog[index]: %f\n",rog[index]/visits[index]);
  //        fclose(stdoutlog);

  double ssum[numbercores * pertrusionpercore] = {
    0
  };
  double tort = 0;

  int x1, y1, z1;
  int x2, y2, z2;

  for (int i = (0); i < numbercores * pertrusionpercore; i++) // performs the rotations
  {
    int s_vector[3 * lengthpercore] = {
      0
    };
    double s_avg[3] = {
      0
    };

    x1 = s_dx[i * lengthpercore + 0] - s_dx[i * lengthpercore + 0 + 1];
    x2 = s_dx[i * lengthpercore + 0] - s_dx[i * lengthpercore + 0 + 2];

    y1 = s_dy[i * lengthpercore + 0] - s_dy[i * lengthpercore + 0 + 1];
    y2 = s_dy[i * lengthpercore + 0] - s_dy[i * lengthpercore + 0 + 2];

    z1 = s_dz[i * lengthpercore + 0] - s_dz[i * lengthpercore + 0 + 1];
    z2 = s_dz[i * lengthpercore + 0] - s_dz[i * lengthpercore + 0 + 2];

    s_vector[3 * 0 + 0] = y1 * z2 - z1 * y2;
    s_vector[3 * 0 + 1] = z1 * x2 - x1 * z2;
    s_vector[3 * 0 + 2] = x1 * y2 - y1 * x2;

    for (int j = 1; j < lengthpercore - 2; j++) // performs the rotations
    {
      x1 = s_dx[i * lengthpercore + j] - s_dx[i * lengthpercore + j + 1];
      x2 = s_dx[i * lengthpercore + j] - s_dx[i * lengthpercore + j + 2];

      y1 = s_dy[i * lengthpercore + j] - s_dy[i * lengthpercore + j + 1];
      y2 = s_dy[i * lengthpercore + j] - s_dy[i * lengthpercore + j + 2];

      z1 = s_dz[i * lengthpercore + j] - s_dz[i * lengthpercore + j + 1];
      z2 = s_dz[i * lengthpercore + j] - s_dz[i * lengthpercore + j + 2];

      s_vector[3 * j + 0] = y1 * z2 - z1 * y2 + s_vector[3 * (j - 1) + 0];
      s_vector[3 * j + 1] = z1 * x2 - x1 * z2 + s_vector[3 * (j - 1) + 1];
      s_vector[3 * j + 2] = x1 * y2 - y1 * x2 + s_vector[3 * (j - 1) + 2];
    }
    for (int j = 0; j < lengthpercore - 2; j++) // performs the rotations
    {
      s_avg[0] += s_vector[3 * j + 0];
      s_avg[1] += s_vector[3 * j + 1];
      s_avg[2] += s_vector[3 * j + 2];
    }
    s_avg[0] /= ((double)((lengthpercore - 2)));
    s_avg[1] /= ((double)((lengthpercore - 2)));
    s_avg[2] /= ((double)((lengthpercore - 2)));

    for (int j = 0; j < lengthpercore - 2; j++) // performs the rotations
    {
      ssum[i] = (((double)(s_vector[3 * (j)])) - s_avg[0]) * (((double)(s_vector[3 * (j)])) - s_avg[0]) + (((double)(s_vector[3 * (j) + 1])) - s_avg[1]) * (((double)(s_vector[3 * (j) + 1])) - s_avg[1]) + (((double)(s_vector[3 * (j) + 2])) - s_avg[2]) * (((double)(s_vector[3 * (j) + 2])) - s_avg[2]);
    }

  }

  for (int i = 0; i < numbercores * pertrusionpercore; i++) {
    tort += sqrt(ssum[i] / (lengthpercore - 2));
  }

  tortuosity[index] += tort / (numbercores * pertrusionpercore);

  for (int i = (0); i < numbercores * pertrusionpercore * lengthpercore; i++) // performs the rotations
  {
    density_profile[index * L1dim_S_z + (poly_lattice_coordinates[3 * i + 2])]++;
  }

  //int fakeDPS[L1dim_S_z] = {0};
  int s_array[3] = {
    0
  };
  for (int i = (0); i < no_solvent_sites; i++) // performs the rotations
  {
    general_coord(s_array, solvent_loc[i]);
    density_profile_solvent[index * L1dim_S_z + (s_array[2])]++;
    //      fakeDPS[(s_array[2])]++;
  }

  // output_lattice(index,(radius/(numbercores*pertrusionpercore))*(radius/(numbercores*pertrusionpercore))/lengthpercore,tort/(numbercores*pertrusionpercore));

  //
  //    int summed_diff[3] = {0,0,0};
  //        int arm = 0;
  //        sprintf(filename, "Chain_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
  //        if ((file = fopen(filename, "w")) == NULL)
  //        {
  //            stdoutlog = fopen(stdoutlogname, "a");
  //            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
  //            fclose(stdoutlog);
  //            MPI_Abort(MPI_COMM_WORLD, 1);
  //        }
  //        else
  //        {
  //            for (int i = 0; i < numbercores*pertrusionpercore*lengthpercore; i++)
  //            {
  //                summed_diff[0]=poly_lattice_coordinates[3 * (arm*lengthpercore)];
  //                summed_diff[1]=poly_lattice_coordinates[3 * (arm*lengthpercore) + 1];
  //                summed_diff[2]=poly_lattice_coordinates[3 * (arm*lengthpercore) + 2];
  //
  //                summed_diff[0]+=dx[i];
  //                summed_diff[1]+=dy[i];
  //                summed_diff[2]+=dz[i];
  //
  //                fprintf(file, "%i\t%i\t%i\n",  summed_diff[0], summed_diff[1], summed_diff[2]);
  //
  //                if((i+1)%lengthpercore==0)
  //                {
  //                    arm++;
  //                }
  //            }
  //        }
  //        fclose(file);
  //
  //
  //            sprintf(filename, "Solvent_CCords_Proc%i_%i_%i.txt",myid,numbercores,lengthpercore);
  //            if ((file = fopen(filename, "w")) == NULL)
  //            {
  //                stdoutlog = fopen(stdoutlogname, "a");
  //                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
  //                fclose(stdoutlog);
  //                MPI_Abort(MPI_COMM_WORLD, 1);
  //            }
  //            else
  //            {
  //                for (int i = 0; i < no_solvent_sites; i++)
  //                {
  //                        general_coord(s_array,solvent_loc[i]);
  //
  //                    fprintf(file, "%i\t%i\t%i\n",  s_array[0], s_array[1], s_array[2]);
  //
  //                }
  //            }
  //            fclose(file);
  //
  //        if(myid==6)MPI_Abort(MPI_COMM_WORLD, 1);

}

void restart(int( & energie), int( & eold), int( & iteration), double( & lnf), int( & swap_every), double( & sweep)) {
  int stop = 0;
  sprintf(resetbuffer, "SysStatus.L%i.proc%04i", L1dim_S_xy, myid); // Preliminary output file that allows for process inspection
  FILE * rfptr;
  long long phldr;
  int Ewidth = Emaxindex - Eminindex;
  int fakeEmin1 = Eminindex;
  int fakeEmin2 = Eminindex;
  int fakeEmin3 = Eminindex;
  int fakeEmin4 = Eminindex;
  int fakeEmin5 = Eminindex;
  int fakeEmin6 = 0;
  int fakeEmin7 = ((Eminindex - 1) * L1dim_S_z) + L1dim_S_z;
  int fakeEmin8 = ((Eminindex - 1) * L1dim_S_z) + L1dim_S_z;
  int fakeEmin9 = ((Eminindex - 1) * numbercores) + numbercores;
  int fakeEmin10 = Eminindex;
  int fakeEmin11 = Eminindex;
  int fakeEmin12 = Eminindex;
  int chain_increment = 0;

  int counttt = 0;
  int otherrr = 0;
  rfptr = fopen(resetbuffer, "r");

  int resetiter = 1;

  while (fscanf(rfptr, "%lld\n", & phldr) > 0) {
    if (resetiter == 1) {
      // ratiovalue =0.9956;
    }
    if (resetiter == 2) {
      energie = phldr;
    }
    if (resetiter == 3) {
      eold = phldr;
    }
    if (resetiter == 4) {
      sweep = phldr;
    }
    if (resetiter == 5) {
      swap_every = phldr;
    }
    if (resetiter == 6) {
      lnf = ((double) phldr) / pow(10, 15);
    }
    if (resetiter == 7) {
      iteration = phldr;
    }
    if (resetiter >= (8) && resetiter < (8 + Ewidth + 1)) {
      HE[fakeEmin1] = ((int) phldr);
      fakeEmin1++;
    }
    if (resetiter >= (8 + Ewidth + 1) && resetiter < ((8 + Ewidth + 1) + Ewidth + 1)) {
      lngE[fakeEmin2] = ((double) phldr) / pow(10, 9);
      fakeEmin2++;
    }
    if (resetiter >= ((8 + Ewidth + 1) + Ewidth + 1) && resetiter < (((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1)) {
      rog[fakeEmin3] = ((double) phldr) / pow(10, 6);
      fakeEmin3++;
    }
    if (resetiter >= (((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) && resetiter < ((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1)) {
      tortuosity[fakeEmin4] = ((double) phldr) / pow(10, 6);
      fakeEmin4++;
    }
    if (resetiter >= ((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) && resetiter < (((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1)) {
      visits[fakeEmin5] = ((long long) phldr);
      fakeEmin5++;
    }
    if (resetiter >= (((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) && resetiter < ((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore))) {
      poly_lattice_indexes[chain_increment] = ((int) phldr);
      chain_increment++;
    }
    if (resetiter >= ((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) &&
      resetiter < ((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) {
      solvent_loc[fakeEmin6] = ((int) phldr);
      fakeEmin6++;
    }
    if (resetiter >= ((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites &&
      resetiter < (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z)) {
      density_profile[fakeEmin7] = ((double) phldr) / pow(10, 9);
      fakeEmin7++;
      counttt++;
    }
    if (resetiter >= (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) &&
      resetiter < (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z)) {
      density_profile_solvent[fakeEmin8] = ((double) phldr) / pow(10, 9);
      fakeEmin8++;
    }
    if (resetiter >= (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z) &&
      resetiter < (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (numbercores)) {
      indiv_rog[fakeEmin9] = ((double) phldr) / pow(10, 6);
      fakeEmin9++;
    }
    if (resetiter >= (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (numbercores) &&
      resetiter < (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (numbercores) + Ewidth + 1) {
      rog_z[fakeEmin10] = ((double) phldr) / pow(10, 6);
      fakeEmin10++;
    }
    if (resetiter >= (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (numbercores) + Ewidth + 1 &&
      resetiter < (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + (numbercores * lengthpercore * pertrusionpercore)) + no_solvent_sites) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (L1dim_S_z) + (Ewidth + 1) * (numbercores) + Ewidth + 1 + Ewidth + 1) {
      rog_xy[fakeEmin11] = ((double) phldr) / pow(10, 6);
      fakeEmin11++;
    }

    resetiter++;
  }
  fclose(rfptr);

  reset_indexes[0] = 0 - 1;
  reset_indexes[1] = no_solvent_sites;
  poly_latticepoint_clear(wl_solvent_loc, reset_indexes[0], reset_indexes[1]);
  poly_latticepoint_clear(ph_solvent_loc, reset_indexes[0], reset_indexes[1]);
  poly_latticepoint_clear(solvent_loc, reset_indexes[0], reset_indexes[1]);

  reset_indexes[0] = 0 - 1;
  reset_indexes[1] = numbercores * lengthpercore * pertrusionpercore;
  poly_latticepoint_clear(wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
  poly_latticepoint_clear(ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  poly_latticepoint_clear(poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);

  for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) {
    poly_coord_decompose(poly_lattice_indexes[i], i);
    //        stdoutlog = fopen(stdoutlogname, "a");
    //        fprintf(stdoutlog, "plc:(%i,%i,%i)\n",poly_lattice_coordinates[3*i],poly_lattice_coordinates[3*i+1],poly_lattice_coordinates[3*i+2]);
    //        fclose(stdoutlog);
  }

  for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) {
    wl_pseudo_chain[i] = poly_lattice_indexes[i];
    wl_pseudo_chain_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
    wl_pseudo_chain_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
    wl_pseudo_chain_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];

    ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
    ph_poly_lattice_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
    ph_poly_lattice_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
    ph_poly_lattice_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];
  }
  for (int i = 0; i < no_solvent_sites; i++) {
    wl_solvent_loc[i] = solvent_loc[i];
    ph_solvent_loc[i] = solvent_loc[i];
  }

  //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  //lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  //poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

  //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
  poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

  reset_indexes[0] = 0 - 1;
  reset_indexes[1] = no_solvent_sites;

  solvent_reset(ph_solvent_loc, solvent_loc);

  int eenergie = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc);
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

  if (eenergie != energie) {
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "2nd Estart out of range!!!  energie %i  eeenergie: %i\n", energie, eenergie);
    fclose(stdoutlog);

    stop = 1;
  } else {
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "2nd OK\n");
    fclose(stdoutlog);

  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (stop) MPI_Abort(MPI_COMM_WORLD, 1);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char * argv[]) {

  /*time(&timenow);
  std::cout << ctime(&timenow) << std::endl;*/

  clock_t tStart = clock();

  // set up MPI and MPI_COMM_WORLD
  MPI_Init( & argc, & argv);
  MPI_Comm_size(MPI_COMM_WORLD, & numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, & myid);

  // check command line arguments
  sprintf(stdoutlogname, "Error%04i.log", myid);
  if (argc < 5) {
    if (myid == 0) {
      stdoutlog = fopen(stdoutlogname, "w");
      fprintf(stdoutlog, "Unexpected number of command line arguments!\n");
      fprintf(stdoutlog, "Expect 4 arguments, %d were provided.\n", argc - 1);
      fprintf(stdoutlog, "Syntax: ./WLpotts_mpi [arg1] [arg2] [arg3] [arg4] \n\n");
      fprintf(stdoutlog, "Please provide the following command line arguments:\n");
      fprintf(stdoutlog, "1. Overlap between consecutive windows. [double, 0 <= overlap <= 1]\n");
      fprintf(stdoutlog, "2. Number of walkers per energy subwindow. [integer]\n");
      fprintf(stdoutlog, "3. Number of Monte Carlo steps between replica exchange. [integer]\n");
      fprintf(stdoutlog, "4. Random number seed. [integer]\n\n");
      fclose(stdoutlog);

      printf("ERROR: Unexpected number of command line arguments!\n");
      printf("       Expect 4 arguments, %d were provided.\n", argc - 1);
      printf("Syntax: ./WLpotts_mpi [arg1] [arg2] [arg3] [arg4] \n\n");
      printf("Please provide the following command line arguments:\n");
      printf("1. Overlap between consecutive windows. [double, 0 <= overlap <= 1]\n");
      printf("2. Number of walkers per energy subwindow. [integer]\n");
      printf("3. Number of Monte Carlo steps between replica exchange. [integer]\n");
      printf("4. Random number seed. [integer]\n\n");

    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // number of walkers per energy window
  multiple = atoi(argv[2]);

  // at the moment, the code works only for an _odd_ number of energy windows
  // (to make the RE in windows at the outside consistent)
  if ((numprocs / multiple) % 2 == 0) {
    if (myid == 0) {
      stdoutlog = fopen(stdoutlogname, "a");
      fprintf(stdoutlog, "ERROR: Even number of energy windows (%d) requested. Please request an odd number of energy windows.\n\n", numprocs / multiple);
      fclose(stdoutlog);

      printf("ERROR: Even number of energy windows (%d) requested. Please request an odd number of energy windows.\n\n", numprocs / multiple);
    }

    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // initiate random numbers
  rseed = atoi(argv[4]);
  rseed = time(NULL);
  //rseed = 1650860932;
  srand(rseed + myid);

  int restart_check = atoi(argv[5]);

  int swap_every = atoi(argv[3]); // after this number of sweeps try conformations swap
  int swap_every_init = swap_every;

  // init log file
  sprintf(stdoutlogname, "Proc%04i.sim%i.log", myid, rseed);

  // set local energy range
  ret_status = find_local_energy_range(Eglobalmin, Eglobalmax, atof(argv[1]), numprocs);

  MPI_Barrier(MPI_COMM_WORLD);

  if (ret_status != 0) {
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc. %3i: find_local_energy_range() returned %i\n", myid, ret_status);
    fclose(stdoutlog);

    MPI_Abort(MPI_COMM_WORLD, 1); // something went wrong in find_local_energy_range()
  }

  init_hists(); // moved above init_lattice() for calculation considerations
  init_neighbors();
  init_lattice(Emin, Emax); // 0 - random; 1 - all equal

  // calculate energy for the first time
  int eold, energie;
  energie = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc);
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Proc %3i: energy at start=%i\n", myid, energie);
  fclose(stdoutlog);
  energie -= Eglobalmin; // shift to positive values to use it as array index

  Estartindex = energie;

  int stop = 0;
  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Proc %3i: Parameter: Eminindex -- Emaxindex: %i -- %i; Estartindex=%i\t", myid, Eminindex, Emaxindex, Estartindex);
  fclose(stdoutlog);

  if ((Estartindex > Emaxindex) || (Estartindex < Eminindex)) {
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Estart out of range!!!\n");
    fclose(stdoutlog);

    stop = 1;
  } else {
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "OK\n");
    fclose(stdoutlog);

  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (stop) MPI_Abort(MPI_COMM_WORLD, 1);

  //Teststop
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   MPI_Abort(MPI_COMM_WORLD,1);

  // create new groups and communicators for each energy range window
  stdoutlog = fopen(stdoutlogname, "a");
  MPI_Comm_group(MPI_COMM_WORLD, & world); // get the group of processes in MPI_COMM_WORLD (i.e. all)
  int * ranks;
  ranks = (int * ) malloc(2 * multiple * sizeof(int));
  mpi_local_group = (MPI_Group * ) malloc(((numprocs / multiple) - 1) * sizeof(MPI_Group));
  mpi_local_comm = (MPI_Comm * ) malloc(((numprocs / multiple) - 1) * sizeof(MPI_Comm));

  for (int i = 0; i < ((numprocs / multiple) - 1); i++) // i is counter for the energy range windows
  {
    for (int j = 0; j < 2 * multiple; j++) {
      ranks[j] = i * multiple + j; // contains the ranks of processes in MPI_COMM_WORLD which should get into local group
      if (myid == 0) {
        fprintf(stdoutlog, "Proc %3i: %i will be part of communicator/group %i\n", myid, ranks[j], i);
      }
    }
    MPI_Group_incl(world, 2 * multiple, ranks, & mpi_local_group[i]); // create local group
    MPI_Comm_create(MPI_COMM_WORLD, mpi_local_group[i], & mpi_local_comm[i]); // create communicator for that group
  }
  fclose(stdoutlog);
  free(ranks);

  // get my local id (in my local communicators)
  if (myid < numprocs - multiple) {
    comm_id = 2 * (myid / (2 * multiple));
    MPI_Comm_rank(mpi_local_comm[comm_id], & mylocalid[0]);
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: I am part of communicator/group %i with local_id[0]=%i\n", myid, comm_id, mylocalid[0]);
    fclose(stdoutlog);

  } else {
    mylocalid[0] = INT_MAX; // just to give it a value
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: got local_id[0]=%i\n", myid, mylocalid[0]);
    fclose(stdoutlog);
  }

  if (myid >= multiple) {
    comm_id = 2 * ((myid - multiple) / (2 * multiple)) + 1;
    MPI_Comm_rank(mpi_local_comm[comm_id], & mylocalid[1]);
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: I am part of communicator/group %i with local_id[1]=%i\n", myid, comm_id, mylocalid[1]);
    fclose(stdoutlog);
  } else {
    mylocalid[1] = INT_MAX; // just to give it a value
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: got local_id[1]=%i\n", myid, mylocalid[1]);
    fclose(stdoutlog);
  }
  //fclose(stdoutlog);

  // log 'path'data (the same data is written, just differently sorted: with respect to the sample id and with respect to process id)
  sprintf(wanderlogname, "wanderlog_rseed%i_sample_%i.dat", rseed, latticepoint_S[numberspins_S + 2]);
  wanderlog = fopen(wanderlogname, "w");
  time( & timenow);
  fprintf(wanderlog, "#sweep\t\tproc-id\tconf-id\tenergy\t\ttime\n");
  fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", 0.0, myid, latticepoint_S[numberspins_S + 2], energie, ctime( & timenow));
  fclose(wanderlog);

  sprintf(wanderlogname, "wanderlog_rseed%i_walker_%i.dat", rseed, myid);
  wanderlog = fopen(wanderlogname, "w");
  time( & timenow);
  fprintf(wanderlog, "#sweep\t\tproc-id\tconf-id\tenergy\t\ttime\n");
  fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", 0.0, myid, latticepoint_S[numberspins_S + 2], energie, ctime( & timenow));
  fclose(wanderlog);

  //start simulation
  double wsk, dice; // wsk: probability
  int wiggle;
  int wiggletwo; // ID of spin to be updated (only single-spin update implemented)

  // the WL parameter should eventually go into the init file, too
  double countdown = 2;
  double lnf = 1.0; // my modification factor
  double lnf_slowest = lnf; // modification factor of slowest walker in my window
  long long recombine_counter = 0;
  long long recombine_counter_c = 0; // modification factor of slowest walker in my window
  //double lnfmin=log(1.000000001);
  double lnfmin = log(1.0000000001); // terminal modification factor
  double sweep = 0; // counter for MC sweeps
  int flat; // 0 - histogram not flat; 1 - flat
  int iteration = 1; // WL iteration counter
  int iter = 1;
  double check_flatness_every = 5000; // in number of sweeps
  int backup;
  int backuptwo;

  eold = energie;

  for (int i = 0; i <= hist_size; i++) HE[i] = 0; //init H(E)

  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Proc %3i: Start WL iteration78654\n", myid);
  fclose(stdoutlog);

  if (restart_check) {
    restart(energie, eold, iteration, lnf, swap_every, sweep);
  }

  int swtch;
  int found = 0;

  int prelim_counter = 0;

  stdoutlog = fopen(stdoutlogname, "a");
  fprintf(stdoutlog, "Proc %3i: Start WL iteration\n", myid);
  fclose(stdoutlog);

  for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) {
    wl_pseudo_chain[i] = poly_lattice_indexes[i];
    wl_pseudo_chain_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
    wl_pseudo_chain_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
    wl_pseudo_chain_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];
  }
  for (int i = 0; i < no_solvent_sites; i++) {
    wl_solvent_loc[i] = solvent_loc[i];
  }

  long long restartiter = 0;
  int eee;
  // Main Wang-Landau routine
  while (lnf_slowest > lnfmin) {

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Proc %3i: Start WL iteration\n", myid);
    fclose(stdoutlog);
    if (restartiter % 200 == 0) {
      //            for (int i = 0; i < numberspins_S; i++)
      //            {
      //                latticepoint_S[i] = 0;
      //            }
      //            for (int i = 0; i < numberspins_L; i++)
      //            {
      //                latticepoint_L[i] = 0;
      //            }
      //             int array[3] = { 0,0,0 };
      //            for (int i = 0; i < numberspins_S; i++)
      //            {
      //                general_coord(array, i);
      //
      //                //stdoutlog = fopen(stdoutlogname, "a");
      //                //fprintf(stdoutlog, "array: %i \t%i\t%i \n",array[0],array[1],array[2]);
      //                //fclose(stdoutlog);
      //
      //                if(array[2] == plane_z_max && latticepoint_S[i] == 0)
      //                {
      //                    latticepoint_S[i] = -1;
      //                }
      //            }
      //
      //            reset_indexes[0] = 0-1;
      //            reset_indexes[1]=numbercores * lengthpercore * pertrusionpercore;
      //            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
      //                poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
      //
      //
      //
      //                reset_indexes[0] = 0-1;
      //                reset_indexes[1]=no_solvent_sites;
      //
      //                solvent_reset(ph_solvent_loc,solvent_loc);
      //
      //
      eye();
      eee = total_energy(poly_lattice_coordinates, poly_lattice_indexes, solvent_loc);
      if (eee - Eglobalmin != energie) {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "energy different  eee:%i  energue:%i\n", eee, energie);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      sprintf(resetbuffer, "SysStatus.L%i.proc%04i", L1dim_S_xy, myid); // Preliminary output file that allows for process inspection

      FILE * filee;

      if ((filee = fopen(resetbuffer, "w")) == NULL) {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", resetbuffer);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
      } else {

        fprintf(filee, "%lld\n", (long long)(.99 * pow(10, 17)));
        fprintf(filee, "%i\n", energie);
        fprintf(filee, "%i\n", eold);
        fprintf(filee, "%lld\n", (long long)(sweep));
        fprintf(filee, "%i\n", swap_every);
        fprintf(filee, "%lld\n", (long long)(lnf * pow(10, 15)));
        fprintf(filee, "%i\n", iteration);

        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long) HE[i]);
        }

        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(lngE[i] * pow(10, 9)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(rog[i] * pow(10, 6)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(tortuosity[i] * pow(10, 6)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(visits[i]));
        }
        for (int i = 0; i < numbercores * pertrusionpercore * lengthpercore; i++) {
          fprintf(filee, "%lld\n", (long long) poly_lattice_indexes[i]);
        }
        for (int i = 0; i < no_solvent_sites; i++) {
          fprintf(filee, "%lld\n", (long long) solvent_loc[i]);
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          for (int j = 0; j < L1dim_S_z; j++) {
            fprintf(filee, "%lld\n", (long long)(density_profile[i * L1dim_S_z + j] * pow(10, 9)));
          }
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          for (int j = 0; j < L1dim_S_z; j++) {
            fprintf(filee, "%lld\n", (long long)(density_profile_solvent[i * L1dim_S_z + j] * pow(10, 9)));
          }
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          for (int i = 0; i < numbercores; i++) {
            fprintf(filee, "%lld\n", (long long)(indiv_rog[i * numbercores] * pow(10, 6)));
          }
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(rog_z[i] * pow(10, 6)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(rog_xy[i] * pow(10, 6)));
        }

      }
      MPI_Barrier(MPI_COMM_WORLD);
      sprintf(resetbuffer, "SysStatus.L%i.proc%04i", L1dim_S_xy, myid); // Preliminary output file that allows for process inspection

      if ((filee = fopen(resetbuffer, "w")) == NULL) {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Can not open file %s. Exit.\n", resetbuffer);
        fclose(stdoutlog);
        MPI_Abort(MPI_COMM_WORLD, 1);
      } else {

        fprintf(filee, "%lld\n", (long long)(.99 * pow(10, 17)));
        fprintf(filee, "%i\n", energie);
        fprintf(filee, "%i\n", eold);
        fprintf(filee, "%lld\n", (long long)(sweep));
        fprintf(filee, "%i\n", swap_every);
        fprintf(filee, "%lld\n", (long long)(lnf * pow(10, 15)));
        fprintf(filee, "%i\n", iteration);

        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long) HE[i]);
        }

        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(lngE[i] * pow(10, 9)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(rog[i] * pow(10, 6)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(tortuosity[i] * pow(10, 6)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(visits[i]));
        }
        for (int i = 0; i < numbercores * pertrusionpercore * lengthpercore; i++) {
          fprintf(filee, "%lld\n", (long long) poly_lattice_indexes[i]);
        }
        for (int i = 0; i < no_solvent_sites; i++) {
          fprintf(filee, "%lld\n", (long long) solvent_loc[i]);
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          for (int j = 0; j < L1dim_S_z; j++) {
            fprintf(filee, "%lld\n", (long long)(density_profile[i * L1dim_S_z + j] * pow(10, 9)));
          }
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          for (int j = 0; j < L1dim_S_z; j++) {
            fprintf(filee, "%lld\n", (long long)(density_profile_solvent[i * L1dim_S_z + j] * pow(10, 9)));
          }
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          for (int i = 0; i < numbercores; i++) {
            fprintf(filee, "%lld\n", (long long)(indiv_rog[i * numbercores] * pow(10, 6)));
          }
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(rog_z[i] * pow(10, 6)));
        }
        for (int i = Eminindex; i <= Emaxindex; i++) {
          fprintf(filee, "%lld\n", (long long)(rog_xy[i] * pow(10, 6)));
        }

      }
      restartiter = 0;
    }
    restartiter++;

    for (int k = 0; k < check_flatness_every; k++) {
      for (int i = 0; i < (numbercores * lengthpercore * pertrusionpercore) + no_solvent_sites; i++) // this does 1 MC sweep
      {

        //                int eeee = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
        //                stdoutlog = fopen(stdoutlogname, "a");
        //                fprintf(stdoutlog, "eeeb: %i\n\t",eeee);
        //                fclose(stdoutlog);

        energie = eold + propose_update(); // calculate energy of updated configuration

        if (MoveProposal == 0 || ((energie > Emaxindex) || (energie < Eminindex))) // boundary check
        {
          //                    stdoutlog = fopen(stdoutlogname, "a");
          //                    fprintf(stdoutlog, "10000");
          //                            fclose(stdoutlog);
          //
          //                    stdoutlog = fopen(stdoutlogname, "a");
          //                    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
          //                    {
          //                        fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
          //                    }
          //                    for (int i = 0; i < no_solvent_sites; i++)
          //                    {
          //                        fprintf(stdoutlog, "%4i\t\t%4i\t\t%i\n", solvent_loc[i],ph_solvent_loc[i],wl_solvent_loc[i]);
          //                    }

          //                    fclose(stdoutlog);
          if (MoveProposal == 1) {
            if (poly_solvent_mov == 0) {
              //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
              lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
              poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
              //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
              lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
              poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            }
            if (poly_solvent_mov == 1) {
              solvent_reset(solvent_loc, wl_solvent_loc);
              //                        stdoutlog = fopen(stdoutlogname, "a");
              //                        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
              //                        {
              //                            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
              //                        }
              //                        for (int i = 0; i < no_solvent_sites; i++)
              //                        {
              //                            fprintf(stdoutlog, "%4i\t\t%4i\t\t%i\n", solvent_loc[i],ph_solvent_loc[i],wl_solvent_loc[i]);
              //                        }
              //
              //                        fclose(stdoutlog);
              solvent_reset(ph_solvent_loc, solvent_loc);

            }
          }
          //                    stdoutlog = fopen(stdoutlogname, "a");
          //                    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
          //                    {
          //                        fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
          //                    }
          //                    for (int i = 0; i < no_solvent_sites; i++)
          //                    {
          //                        fprintf(stdoutlog, "%4i\t\t%4i\t\t%i\n", solvent_loc[i],ph_solvent_loc[i],wl_solvent_loc[i]);
          //                    }
          //
          //                    fclose(stdoutlog);
          energie = eold;
          //lngE[energie] = 1; // the real lngE of the system is set to a value of 1 if visitation occured at a specific energy (can be subtracted later) to force visitation in the real wl process
        } else // calculate acceptance propbability
        {

          // the regular wl process is kept using a place holder pseudolngE to force the system to not move to move to higher energies
          dice = (1.0 * rand() / (RAND_MAX + 1.0)); // roll the dice
          wsk = exp(lngE[eold - eold % 10] - lngE[energie - energie % 10]); // calculate acceptance probability
          if (dice > wsk) // reject
          {
            //                        stdoutlog = fopen(stdoutlogname, "a");
            //                        fprintf(stdoutlog, "20000");
            //                                fclose(stdoutlog);
            if (poly_solvent_mov == 0) {
              //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
              lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
              poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
              //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
              lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
              poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            }
            if (poly_solvent_mov == 1) {
              solvent_reset(solvent_loc, wl_solvent_loc);
              solvent_reset(ph_solvent_loc, solvent_loc);

            }

            energie = eold;

          } else {
            eold = energie; // accept

            //sysparam(energie-energie%10);
            //                        stdoutlog = fopen(stdoutlogname, "a");
            //                        fprintf(stdoutlog, "30000");
            //                                fclose(stdoutlog);
            if (poly_solvent_mov == 0) {
              //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
              lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
              poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            }
            if (poly_solvent_mov == 1) {
              solvent_reset(wl_solvent_loc, solvent_loc);

            }

            if (energie < minimum_en_config) {
              printed = 0;
              minimum_en_config = energie;
            }
            if (myid == 0 && energie <= minimum_en_config && printed == 0) {
              printed = 1;

              sprintf(filename, "Minimum_Config_Poly.txt");
              if ((file = fopen(filename, "w")) == NULL) {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD, 1);
              } else {
                for (int i = 0; i < numbercores * pertrusionpercore * lengthpercore; i++) {
                  fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
                }
              }
              fclose(file);

              sprintf(filename, "Minimum_Config_Solvent.txt");
              if ((file = fopen(filename, "w")) == NULL) {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD, 1);
              } else {
                for (int i = 0; i < no_solvent_sites; i++) {
                  fprintf(file, "%i\t%i\n", i, solvent_loc[i]);
                }
              }
              fclose(file);

              stdoutlog = fopen(stdoutlogname, "a");
              fprintf(stdoutlog, "Printed energie %i\n", energie);
              fclose(stdoutlog);

              eye();

            }
            if (energie > maximum_en_config) {
              max_printed = 0;
              maximum_en_config = energie;
            }
            if (myid == numprocs - 1 && energie >= maximum_en_config && max_printed == 0) {
              max_printed = 1;

              sprintf(filename, "Maximum_Config_Poly.txt");
              if ((file = fopen(filename, "w")) == NULL) {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD, 1);
              } else {
                for (int i = 0; i < numbercores * pertrusionpercore * lengthpercore; i++) {
                  fprintf(file, "%i\t%i\n", i, poly_lattice_indexes[i]);
                }
              }
              fclose(file);

              sprintf(filename, "Maximum_Config_Solvent.txt");
              if ((file = fopen(filename, "w")) == NULL) {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD, 1);
              } else {
                for (int i = 0; i < no_solvent_sites; i++) {
                  fprintf(file, "%i\t%i\n", i, solvent_loc[i]);
                }
              }
              fclose(file);

              stdoutlog = fopen(stdoutlogname, "a");
              fprintf(stdoutlog, "Printed energie %i\n", energie);
              fclose(stdoutlog);

              eye();
            }

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
        }
        // update histograms
        lngE[energie - energie % 10] += lnf;
        HE[energie - energie % 10]++;

        //                int eee = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
        //                stdoutlog = fopen(stdoutlogname, "a");
        //                fprintf(stdoutlog, "eee: %i   energie:%i  s_or_p:%i   r0:%i   r1:%i\n",eee,energie,poly_solvent_mov,reset_indexes[0],reset_indexes[1]);
        //                fclose(stdoutlog);
        //
        //                stdoutlog = fopen(stdoutlogname, "a");
        //                for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
        //                {
        //                    fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
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
      } // end 1 sweep

      sweep++; // sweep counter

      swap_every--; // RE counter (countdown)
      if (swap_every == 0) // ignition! (time to try to exchange configurations)
      {
        swap_every = swap_every_init; // reset countdown clock

        if (replica_exchange(sweep / swap_every, eold) == 1) // if configs were exchanged ('accept'-branch)
        {

          energie = switched;
          eold = energie;

          reset_indexes[0] = 0 - 1;
          reset_indexes[1] = no_solvent_sites;
          poly_latticepoint_clear(wl_solvent_loc, reset_indexes[0], reset_indexes[1]);

          reset_indexes[0] = 0 - 1;
          reset_indexes[1] = numbercores * lengthpercore * pertrusionpercore;
          poly_latticepoint_clear(wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);

          for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) {
            wl_pseudo_chain[i] = poly_lattice_indexes[i];
            wl_pseudo_chain_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
            wl_pseudo_chain_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
            wl_pseudo_chain_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];

            ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
            ph_poly_lattice_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
            ph_poly_lattice_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
            ph_poly_lattice_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];
          }
          for (int i = 0; i < no_solvent_sites; i++) {
            wl_solvent_loc[i] = solvent_loc[i];
            ph_solvent_loc[i] = solvent_loc[i];
          }

          //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
          //lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
          //poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

          //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
          lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
          poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);

          reset_indexes[0] = 0 - 1;
          reset_indexes[1] = no_solvent_sites;

          solvent_reset(ph_solvent_loc, solvent_loc);
          stdoutlog = fopen(stdoutlogname, "a");
          fprintf(stdoutlog, "Proc %3i: %i iteration, %e sweeps, Replica Exchange Success\n Energy Calculted: %i \t Recievced Index Energy: %i\n", myid, iteration, sweep, energie, switched);
          fclose(stdoutlog);

        }
        // update histograms (independently of whether RE happened or not)
        lngE[energie - energie % 10] += lnf;
        HE[energie - energie % 10]++;

      }
    } // end computation block

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Prelim\n");
    fclose(stdoutlog);

    sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim_S_xy, q, myid, iteration); // Preliminary output file that allows for process inspection
    if ((file = fopen(filename, "w")) == NULL) {
      stdoutlog = fopen(stdoutlogname, "a");
      fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
      fclose(stdoutlog);
      MPI_Abort(MPI_COMM_WORLD, 1);
    } else {
      for (int i = Eminindex; i <= Emaxindex; i++) {
        if (lngE[i] > 0.5) {
          fprintf(file, "%i\t%i\t%e\t%e\t%e\t%f\t%f\n", i, (i + (Eglobalmin)), 0.0, lngE[i], HE[i], rog[i] / visits[i], tortuosity[i] / visits[i]);
        }
      }
      fclose(file);
    }

    // check flatness
    //eye();
    flat = histflat(Eminindex, Emaxindex, 0.965);

    if (flat == 0) // histogram not flat, do nothing except maybe some logging
    {

      stdoutlog = fopen(stdoutlogname, "a"); // Prints out effictively flatness progress for inspection
      fprintf(stdoutlog, "Proc %3i: %i iteration, %e sweeps, Histogram noch nicht flach\nValue of Tried Minimum Histogram Value: %f    Value of Tried Ratio: %f\n", myid, iteration, sweep, flatmin, flatratio);

      fclose(stdoutlog);
      //phrepv = flatmin; // holds the minimum value of the histogram and checks change and prints out the current lattice if the value has reoccured frequently
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
    } else // histograms of all walkers are 'flat'
    {
      if (lnf > lnfmin) // as long as modification factor is still larger than terminal value
      {
        sprintf(filename, "L%iq%i.HE.proc%04i.iter%i", L1dim_S_xy, q, myid, iteration);
        if ((file = fopen(filename, "w")) == NULL) {
          stdoutlog = fopen(stdoutlogname, "a");
          fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
          fclose(stdoutlog);
          MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
          for (int i = Eminindex; i <= Emaxindex; i++) {
            if (HE[i] > 0.5) fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + (Eglobalmin)), 0.0, lngE[i], HE[i]);
          }
          fclose(file);
        }
      }

      // decrease modification factor
      // canonical method (reduce by sqrt(2)) implemented
      lnf /= 2.0;

      for (int i = 0; i < hist_size; i++) HE[i] = 0; //reset H(E)
      iteration++; // iteration counter

      if (merge_hists == 1) // merge g(E) estimators from multiple walkers in the same energy window
      {
        stdoutlog = fopen(stdoutlogname, "a");
        if (myid % multiple == 0) // 'root' in energy window, receive individual g(E) and send merged g(E)
        {
          for (int i = 1; i < multiple; i++) {
            MPI_Recv( & lngE_buf[0], hist_size, MPI_DOUBLE, myid + i, 77, MPI_COMM_WORLD, & status); // get other dens. of states
            fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, myid + i);
            for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j]; // sum up for average
          }
          for (int j = 0; j < hist_size; j++) lngE[j] /= (double) multiple; // normalize
          for (int i = 1; i < multiple; i++) {
            MPI_Send( & lngE[0], hist_size, MPI_DOUBLE, myid + i, 99, MPI_COMM_WORLD);
            fprintf(stdoutlog, "Proc %i: Sent merged lngE to Proc. %i\n", myid, myid + i);
          }
        } else // send individual g(E) and receive merged g(E)
        {
          MPI_Send( & lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 77, MPI_COMM_WORLD);
          fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, myid - (myid % multiple));
          MPI_Recv( & lngE_buf[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 99, MPI_COMM_WORLD, & status);
          fprintf(stdoutlog, "Proc %i: Received merged lngE from Proc. %i\n", myid, myid - (myid % multiple));
          for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j]; // replace individual lngE (could be done directly, yes)
        }
        fclose(stdoutlog);
      }
    }

    recombine_counter++;
    // check progress from all other windows
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce( & lnf, & lnf_slowest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce( & recombine_counter, & recombine_counter_c, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

    if (lnf_slowest <= (((double)(10000000000 / countdown)) * .0000000001) || recombine_counter_c > 20) {
      MPI_Barrier(MPI_COMM_WORLD);

      stdoutlog = fopen(stdoutlogname, "a");
      fprintf(stdoutlog, "access start\n");
      fclose(stdoutlog);
      accessiblelevels();
      partial_init_hists();
      stdoutlog = fopen(stdoutlogname, "a"); // Prints out effictively flatness progress for inspection
      fprintf(stdoutlog, "\tLeft Access\n");
      fclose(stdoutlog);

      stdoutlog = fopen(stdoutlogname, "a");
      fprintf(stdoutlog, "recombine\n");
      fclose(stdoutlog);
      recombine((((double)(10000000000 / countdown)) * .0000000001));
      stdoutlog = fopen(stdoutlogname, "a");
      fprintf(stdoutlog, "left recombine\n");
      fclose(stdoutlog);
      if (lnf_slowest <= (((double)(10000000000 / countdown)) * .0000000001)) {
        countdown *= 2;
      }
      recombine_counter = 0;
      partial_init_hists();
    }

    // just some logging
    if (flat == 1) {
      stdoutlog = fopen(stdoutlogname, "a");
      fprintf(stdoutlog, "Proc %3i: Start %i iteration, %e sweeps total so far, lnf now %e (lnfmin=%.2e, lnf_slowest=%e)\n", myid, iteration, sweep, lnf, lnfmin, lnf_slowest);
      fprintf(stdoutlog, "Proc %3i: tryleft: %i, exchangeleft %i (Akzeptanzleft:%.2lf) <--> tryright: %i, exchangeright %i (Akzeptanzright:%.2lf)\n", myid, tryleft, exchangeleft, (double) exchangeleft / (double) tryleft, tryright, exchangeright, (double) exchangeright / (double) tryright);
      fclose(stdoutlog);
    }

  } // end while(lnf_slowest>lnfmin) -> this terminates the simulation

  //   free(neighbor);
  //   free(latticepoint);

  // normalize results
  double norm = lngE[Eminindex] - log(q);
  for (int i = 0; i <= hist_size; i++) lngE[i] -= norm;

  // write final data to file
  sprintf(filename, "L%iq%i.proc%04i.lngE", L1dim_S_xy, q, myid);
  if ((file = fopen(filename, "w")) != NULL) {
    for (int i = Eminindex; i <= Emaxindex; i++)
      fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + Eglobalmin), 0.0, lngE[i], HE[i]);
    fclose(file);
  }

  //   free(HE);
  //   free(lngE);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  /*time(&timenow);
  std::cout << ctime(&timenow) << std::endl;*/

  printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

  return (0); // you did it!
}
