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
#include <cstdio>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include "Globals.h"
#include "Energy.h"
#include "Init.h"
#include "Mov.h"
#include "Mov_Pivot.h"
#include "Mov_Pull.h"
#include "Mov_Rebridge.h"
#include "Reorg.h"
#include "EXWL.h"

int switched = INT_MAX;

int iterr;	// Apart of an added debug output. Effictively just a counter for the number of repeated histogram values when flatness is checked
int repv = 0;	// Holds the value of the (potentially) repeated minvalue of the histogram
int phrepv;	// Holds the value of the previous minimum histogram value to compare with the prior variable

int localgroup;
int local_dup_group;
int num_dup_proc;
int dup_headproc;
// type of boundary condition: 0 - periodic; 1 - Braskamp Kunz

int *latticepoint;	// list containing values of all spins
int *neighbor;	// list containing indices of neighbors for all spins

//Value buffers for recombining function
double *HE;	// energy histogram
double *lngE;	// ln g(E) estimator
double *lngE_buf;	// ln g(E) estimator buffer for exchange
double *pseudolngE;
double *real_lngE;
double *real_lngE_buf;
double *microT;
double *microT_buf;

//Polymer chain constants
int *poly_chain;	// Polymer chain indexholder
int *pseudo_poly_chain;	// Pseudo polymer to hold previous indexes if change is rejected

int *doublepseudo_poly_chain;	// Pseudo polymer to hold previous indexes if change is rejected

int *pull_poly_chain;	// Pseudo polymer to hold previous indexes if change is rejected
int *wl_pseudo_chain;	// records an image of the polymer chain in case the movement is rejected
int *chain_sequence;	// stores the polymer sequence so it can be reset for all operations specifically rebridging
int *pivot_indexes;

int *digital_eye;	// holds the pivot indexes
int rebridge_index;

double *centeromass;
int *s_vector;
int *s_vector1;
int *s_vector2;
double *tortuosity;
double *rog;
long long * visits;
double *tortuosity_buf;
double *rog_buf;
long long * visits_buf;
double *s_avg;
double *endtoend;
double *endtoend_buf;

//polymer rebridging accomadations
int *latticepolymer;	// list containing indices of neighbors for all spins

int rseed;	// seed for random number generator
int energy, localenergy;

double Emin, Emax;	// doubles as calculation of boundaries in Energy is not in 'int'
int Eminindex, Emaxindex, Estartindex;	// local boundaries as index 

// MPI; set up with local communicators for replica exchange (RE)
// Needed so that two processes can talk w/o depending on / bothering the others
int numprocs, myid, multiple, comm_id;
// each process belongs to two local groups,
// one to communicate to left neighbor and one to communicate to right neighbor
// each process has different loca IDs in different communicatore, in general 
int mylocalid[2];	// id in local communicators 
MPI_Comm * mpi_local_comm;
MPI_Group * mpi_local_group;
MPI_Group world;
MPI_Status status;
int merge_hists = 1;	// flag whether or not to merge histograms between iterations for multiple walkers on the same energy range

// to keep track of exchange statistics
int tryleft, tryright, exchangeleft, exchangeright;

// File handlers for I/O
FILE * file;
FILE * stdoutlog;
FILE * wanderlog;
char filename[50];
char resetbuffer[50];
char stdoutlogname[128];
char wanderlogname[128];
int ret_status;

double flatratio;
double flatmin;

int happened;
int happened2;
int happened3 = -1;

int cf_value;

int MoveProposal = 0;

int elocvor;
int elocnach;
int r_picker;
int r_picke;
int compl_check;
//int buffer=0;
int hundred;

int axis;
int cos_value;
int sin_value;

int *pivot_indexes_r;
int xr;
int yr;
int zr;
int seq_plhdr;
int check_plhdr;
int pivotflag;

int xx;
int yy;
int zz;

int search_limit = 0;
int search_checks = 0;
int *configuration_search_r;
int *configuration_search_t;
// to track execution times
time_t timenow, timestart, timeend;

int pull_loc;
int pull_hort;

int pull_relax_loc;
int pull_relax_hort;

int wiggle_ht;

int pivot_loc;

int re_ht_si;
int re_ht_ht;

double *enindex;	// hold values for the radius of gyration at specfic energy levels
double *enindex2;	// hold values for the tortuosity at specfic energy levels
double *enindex3;	// hold values for the tortuosity at specfic energy levels
double *lowchecks;
double *highchecks;
int *checks;
double *checks_countdown;
long double unused1;
long double unused2;

int working = 0;

void eye();

void init_check()
{
	enindex = (double*) malloc(hist_size* sizeof(double));
	enindex2 = (double*) malloc(hist_size* sizeof(double));
	enindex3 = (double*) malloc(hist_size* sizeof(double));
	lowchecks = (double*) malloc(hist_size* sizeof(double));
	highchecks = (double*) malloc(hist_size* sizeof(double));
	checks_countdown = (double*) malloc(hist_size* sizeof(double));
	checks = (int*) malloc(hist_size* sizeof(int));

	for (int i = 0; i < hist_size; i++)
	{
		enindex[i] = 0;
		enindex2[i] = 0;
		enindex3[i] = 0;
		lowchecks[i] = 0;
		highchecks[i] = 2;
		checks_countdown[i] = 1;
		checks[i] = 0;
	}

	FILE * fptr;
	int value;
	double value2;
	double value3;
	double value4;
	int iii = 0;
	fptr = fopen("Combined_58mer_3_T_2.txt", "r");
	while (fscanf(fptr, "%i\t%lf\t%lf\t%lf\n", &value, &value2, &value3, &value4) > 0)
	{
		enindex[-value] = value2;
		enindex2[-value] = value3;
		enindex3[-value] = value4;
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", value, value2, value3, value4, enindex[-value], enindex2[-value], enindex3[-value], lngE[-value]);
		fclose(stdoutlog);
		iii++;
	}

	fclose(fptr);

}

/*
int poly_en()
{
	//energy computed for the entire chain length the head and tail are treated differently as the are disconnected (Need to matbe right solo function but no real need currently)
	int en = 0;
	int minus = chain_length-1;
	for (int i = 0; i < (chain_length); i++) {
		for (int j = 0; j < numberneighbors; j++) { 	switch (i) { 	case 0:
				if (latticepoint[neighbor[numberneighbors *poly_chain[i] + j]] == latticepoint[poly_chain[i]] && neighbor[numberneighbors *poly_chain[i] + j] != poly_chain[i+1] && latticepoint[poly_chain[i]]=1)
					en--;

				//stdoutlog = fopen(stdoutlogname, "a");
				//fprintf(stdoutlog, "\n");
				//fprintf(stdoutlog, "1  chain_length:%i  preceeding:%i  succeeding:%i    \n",latticepoint[neighbor[numberneighbors*poly_chain[i]+j]], neighbor[numberneighbors*poly_chain[i]+j],poly_chain[i+1]);

				//fprintf(stdoutlog, "\n");
				//fclose(stdoutlog);

				break;
			case minus:
				if (latticepoint[neighbor[numberneighbors *poly_chain[i] + j]] == latticepoint[poly_chain[i]] && neighbor[numberneighbors *poly_chain[i] + j] != poly_chain[i-1] && latticepoint[poly_chain[i]]==1)
					en--;
				break;
			default:
				if (latticepoint[neighbor[numberneighbors *poly_chain[i] + j]] == latticepoint[poly_chain[i]] && neighbor[numberneighbors *poly_chain[i] + j] != poly_chain[i+1] && neighbor[numberneighbors *poly_chain[i] + j] != poly_chain[i-1] && latticepoint[poly_chain[i]]==1)
					en--;
				break;
			}
		}
	}

	return en / 2;
}

*/

void sysparam(int index)
{
	//eye();
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".01 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	centeromass[0] = 0;
	centeromass[1] = 0;
	centeromass[2] = 0;
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".02 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	s_avg[0] = 0;
	s_avg[1] = 0;
	s_avg[2] = 0;
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".03 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	for (int i = 0; i < chain_length; i++)	// sets all points in lattice to default values
	{
		s_vector[3 *i] = 0;	//x-axis
		s_vector[3 *i + 1] = 0;	// y-axis
		s_vector[3 *i + 2] = 0;	//z-axis
		s_vector1[3 *i] = 0;	//x-axis
		s_vector1[3 *i + 1] = 0;	// y-axis
		s_vector1[3 *i + 2] = 0;	//z-axis
		//        stdoutlog = fopen(stdoutlogname, "a");
		//
		//        fprintf(stdoutlog, "%i  ",poly_chain[i]);
		//
		//        fclose(stdoutlog);
	}

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".04 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	pivot_extract(poly_chain[0], 0, 1);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".05 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".1 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;

	int pivot_increm = 0;
	for (int y = (0); y < chain_length; y++)	// performs the rotations
	{
		centeromass[0] += ((pivot_indexes[3 *pivot_increm]));	//- ((double) (pivot_indexes[3 *(pivot_increm-1)]));
		centeromass[1] += ((pivot_indexes[3 *pivot_increm + 1]));	//- ((double) (pivot_indexes[3 *(pivot_increm-1) + 1]));
		centeromass[2] += ((pivot_indexes[3 *pivot_increm + 2]));	//- ((double) (pivot_indexes[3 *(pivot_increm-1) + 2]));

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "1 x:%i  y:%i  z:%i    comx:%f   comy:%f   comz:%f",pivot_indexes[3 *pivot_increm],pivot_indexes[3 *pivot_increm+1],pivot_indexes[3 *pivot_increm+2],centeromass[0] ,centeromass[1] ,centeromass[2]);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);

		pivot_increm++;
	}

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".2 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	centeromass[0] /= ((double)(chain_length));
	centeromass[1] /= ((double)(chain_length));
	centeromass[2] /= ((double)(chain_length));

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;

	double rsum = 0.0;

	for (int y = (0 + 1); y < chain_length; y++)	// performs the rotations
	{
		rsum += (((double)(pivot_indexes[3 *y])) - centeromass[0]) *(((double)(pivot_indexes[3 *y])) - centeromass[0]) + (((double)(pivot_indexes[3 *y + 1])) - centeromass[1]) *(((double)(pivot_indexes[3 *y + 1])) - centeromass[1]) + (((double)(pivot_indexes[3 *y + 2])) - centeromass[2]) *(((double)(pivot_indexes[3 *y + 2])) - centeromass[2]);

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "2  rsum:%f x:%i  y:%i  z:%i    comx:%f   comy:%f   comz:%f",rsum,pivot_indexes[3 *y],pivot_indexes[3 *y+1],pivot_indexes[3 *y+2],centeromass[0] ,centeromass[1] ,centeromass[2]);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
	}

	rsum += centeromass[0] *centeromass[0] + centeromass[1] *centeromass[1] + centeromass[2] *centeromass[2];
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".3 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	rog[index] += sqrt(rsum / ((double) chain_length));

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;

	if (index + Eglobalmin == -400000)
	{
		int param_holder = (int)(sqrt(rsum / chain_length) *pow(10, 8));
		int similarity = 0;
		int i = 0;
		while (similarity == 0 && i <= search_limit)
		{
			if (configuration_search_r[i] == param_holder)
			{
				similarity = 1;
			}

			i++;
		}

		if (configuration_search_r[search_limit] == 0 && similarity == 0)
		{
			configuration_search_r[search_limit] = param_holder;
			search_limit++;
		}

		if (search_checks % 60 == 0)
		{
			sprintf(filename, "Configs_%i_%i.txt", myid, search_checks);	// Preliminary output file that allows for process inspection
			if ((file = fopen(filename, "w")) == NULL)
			{
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
				fclose(stdoutlog);
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			else
			{
			 	//fprintf(file, "%i\n", search_limit);
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(file, "%i\t%i\t%i\n", pivot_indexes[3 *i + 0], pivot_indexes[3 *i + 1], pivot_indexes[3 *i + 2]);
				}

				fclose(file);

				/*
				fprintf(file, "%i\n", search_limit);
				for (int i = 0; i <= search_limit; i++)
				{
					fprintf(file, "%i\n", configuration_search_r[i]);
				}

				fclose(file);
				 */
			}
		}

		search_checks++;
	}

	if (rog[index] < 0)
	{
		rsum = 0.0;
		for (int y = (0 + 1); y < chain_length; y++)	// performs the rotations
		{
			rsum += (((double)(pivot_indexes[3 *y])) - centeromass[0]) *(((double)(pivot_indexes[3 *y])) - centeromass[0]) + (((double)(pivot_indexes[3 *y + 1])) - centeromass[1]) *(((double)(pivot_indexes[3 *y + 1])) - centeromass[1]) + (((double)(pivot_indexes[3 *y + 2])) - centeromass[2]) *(((double)(pivot_indexes[3 *y + 2])) - centeromass[2]);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "2  rsum:%f x:%i  y:%i  z:%i    comx:%f   comy:%f   comz:%f", rsum, pivot_indexes[3 *y], pivot_indexes[3 *y + 1], pivot_indexes[3 *y + 2], centeromass[0], centeromass[1], centeromass[2]);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
		}

		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "3  rog:%f ",rog[index]);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".4 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;
	pivot_increm = 0;
	int iii;
	//int iii;
	int x1;
	int y1;
	int z1;

	int x2;
	int y2;
	int z2;
	for (int y = (0); y < chain_length - 2; y++)	// performs the rotations
	{
		iii = 0;
		s_vector[3 *y + 0] = 0;
		s_vector[3 *y + 1] = 0;
		s_vector[3 *y + 2] = 0;
		s_vector1[3 *y + 0] = 0;
		s_vector1[3 *y + 1] = 0;
		s_vector1[3 *y + 2] = 0;
		/*
                while(iii<=y)
                {
                    x1=pivot_indexes[3 *iii]-pivot_indexes[3 *(iii+1)];
                     y1=pivot_indexes[3 *iii + 1]-pivot_indexes[3 *(iii+1) + 1];
                     z1=pivot_indexes[3 *iii + 2] - pivot_indexes[3 *(iii+1) + 2];

                     x2=pivot_indexes[3 *iii] - pivot_indexes[3 *(iii+2)];
                     y2=pivot_indexes[3 *iii + 1] - pivot_indexes[3 *(iii+2) + 1];
                     z2=pivot_indexes[3 *iii + 2] - pivot_indexes[3 *(iii+2) + 2];

                    s_vector[3 *y + 0] += y1*z2-z1*y2;
                    s_vector[3 *y + 1] += z1*x2-x1*z2;
                    s_vector[3 *y + 2] += x1*y2-y1*x2;

       	//            stdoutlog = fopen(stdoutlogname, "a");
       	//            fprintf(stdoutlog, "\n");
       	//            fprintf(stdoutlog, "3.5 x1:%i  y1:%i  z1:%i   x2:%i  y2:%i  z2:%i",x1,y1,z1,x2,y2,z2);
       	//            fprintf(stdoutlog, "\n");
       	//            fclose(stdoutlog);

                    iii++;
                }

                */

		if (y == 0)
		{
			x1 = pivot_indexes[3 * 0] - pivot_indexes[3 *(0 + 1)];
			y1 = pivot_indexes[3 * 0 + 1] - pivot_indexes[3 *(0 + 1) + 1];
			z1 = pivot_indexes[3 * 0 + 2] - pivot_indexes[3 *(0 + 1) + 2];

			x2 = pivot_indexes[3 * 0] - pivot_indexes[3 *(0 + 2)];
			y2 = pivot_indexes[3 * 0 + 1] - pivot_indexes[3 *(0 + 2) + 1];
			z2 = pivot_indexes[3 * 0 + 2] - pivot_indexes[3 *(0 + 2) + 2];

			s_vector[3 * 0 + 0] = y1 *z2 - z1 * y2;
			s_vector[3 * 0 + 1] = z1 *x2 - x1 * z2;
			s_vector[3 * 0 + 2] = x1 *y2 - y1 * x2;

			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "3.5 x1:%i  y1:%i  z1:%i   x2:%i  y2:%i  z2:%i",x1,y1,z1,x2,y2,z2);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);

			//iii++;
		}
		else
		{
			x1 = pivot_indexes[3 *y] - pivot_indexes[3 *(y + 1)];
			y1 = pivot_indexes[3 *y + 1] - pivot_indexes[3 *(y + 1) + 1];
			z1 = pivot_indexes[3 *y + 2] - pivot_indexes[3 *(y + 1) + 2];

			x2 = pivot_indexes[3 *y] - pivot_indexes[3 *(y + 2)];
			y2 = pivot_indexes[3 *y + 1] - pivot_indexes[3 *(y + 2) + 1];
			z2 = pivot_indexes[3 *y + 2] - pivot_indexes[3 *(y + 2) + 2];

			s_vector[3 *y + 0] += y1 *z2 - z1 *y2 + s_vector[3 *(y - 1) + 0];
			s_vector[3 *y + 1] += z1 *x2 - x1 *z2 + s_vector[3 *(y - 1) + 1];
			s_vector[3 *y + 2] += x1 *y2 - y1 *x2 + s_vector[3 *(y - 1) + 2];

			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "3.5 x1:%i  y1:%i  z1:%i   x2:%i  y2:%i  z2:%i",x1,y1,z1,x2,y2,z2);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
		}

		pivot_increm++;
		/*
		if(s_vector[3*y]!=s_vector1[3*y] || s_vector[3*y+1]!=s_vector1[3*y+1]  || s_vector[3*y+2]!=s_vector1[3*y+2])
		{
		    stdoutlog = fopen(stdoutlogname, "a");
		    fprintf(stdoutlog, "real2nd\n");
		    fprintf(stdoutlog, "y=%i sx:%i  sy:%i  sz:%i\n",y,s_vector1[3 *y + 0],s_vector1[3 *y + 1],s_vector1[3 *y + 2]);
		    fprintf(stdoutlog, "y=%i sx:%i  sy:%i  sz:%i\n",y,y1*z2-z1*y2,z1*x2-x1*z2,x1*y2-y1*x2);
		    fprintf(stdoutlog, "\n");
		    fclose(stdoutlog);
		}

		 */

	}

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".5 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	pivot_increm = 0;
	s_avg[0] = 0.0;
	s_avg[1] = 0.0;
	s_avg[2] = 0.0;
	for (int y = (0); y < chain_length - 2; y++)	// performs the rotations
	{
		pivot_increm++;
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " pivot middle h_t:%i   angleif:%i  angleiff:%f angle:%f  axis:%i   cos:%f   sin:%f pivot_indexes[3 *pivot_increm]:%f pivot_indexes[3 *pivot_increm+1]:%f  pivot_indexes[3 *pivot_increm+2]: %f  give_up? : %f",h_t,angleif,angleiff,angle,axis ,cos_value,sin_value, pivot_indexes[3 *pivot_increm],pivot_indexes[3 *pivot_increm+1],pivot_indexes[3 *pivot_increm+2] ,pivot_indexes[3 *pivot_increm] *sin_value);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		s_avg[0] += ((s_vector[3 *pivot_increm + 0]));
		s_avg[1] += ((s_vector[3 *pivot_increm + 1]));
		s_avg[2] += ((s_vector[3 *pivot_increm + 2]));

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "3.5 sax:%f  sx:%i  say:%f   sy:%i  saz:%f  sz:%i",s_avg[0],s_vector[3*y+0],s_avg[1],s_vector[3*y+1],s_avg[2],s_vector[3*y+2]);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
	}

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".6 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "5 sax:%f  say:%f  saz:%f",s_avg[0],s_avg[1],s_avg[2]);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	s_avg[0] /= ((double)(chain_length - 2));
	s_avg[1] /= ((double)(chain_length - 2));
	s_avg[2] /= ((double)(chain_length - 2));
	double ssum = 0.0;
	for (int y = (0); y < chain_length - 2; y++)	// performs the rotations
	{
		ssum += (((double)(s_vector[3 *y])) - s_avg[0]) *(((double)(s_vector[3 *y])) - s_avg[0]) + (((double)(s_vector[3 *y + 1])) - s_avg[1]) *(((double)(s_vector[3 *y + 1])) - s_avg[1]) + (((double)(s_vector[3 *y + 2])) - s_avg[2]) *(((double)(s_vector[3 *y + 2])) - s_avg[2]);
	}

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;
	endtoend[index] += (sqrt((pivot_indexes[3 * 0] - pivot_indexes[3 *(chain_length - 1)]) *(pivot_indexes[3 * 0] - pivot_indexes[3 *(chain_length - 1)]) + (pivot_indexes[(3 *0) + 1] - pivot_indexes[3 *(chain_length - 1) + 1]) *(pivot_indexes[(3 *0) + 1] - pivot_indexes[3 *(chain_length - 1) + 1]) + (pivot_indexes[(3 *0) + 2] - pivot_indexes[3 *(chain_length - 1) + 2]) *(pivot_indexes[(3 *0) + 2] - pivot_indexes[3 *(chain_length - 1) + 2]))) / ((double) chain_length);

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;
	double ete = (sqrt((pivot_indexes[3 * 0] - pivot_indexes[3 *(chain_length - 1)]) *(pivot_indexes[3 * 0] - pivot_indexes[3 *(chain_length - 1)]) + (pivot_indexes[(3 *0) + 1] - pivot_indexes[3 *(chain_length - 1) + 1]) *(pivot_indexes[(3 *0) + 1] - pivot_indexes[3 *(chain_length - 1) + 1]) + (pivot_indexes[(3 *0) + 2] - pivot_indexes[3 *(chain_length - 1) + 2]) *(pivot_indexes[(3 *0) + 2] - pivot_indexes[3 *(chain_length - 1) + 2]))) / ((double) chain_length);

	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;
	/*
	if(index ==42)
	{
	    stdoutlog = fopen(stdoutlogname, "a");
	           fprintf(stdoutlog, "\n");
	           fprintf(stdoutlog, "\n");
	            fprintf(stdoutlog, "Index 42");
	           fprintf(stdoutlog, "\n");
	           fprintf(stdoutlog, "rog:%f  tort:%f  ete:%f",sqrt(rsum /((double) chain_length))/((double) chain_length),sqrt(ssum/(chain_length-2)),ete);
	           fprintf(stdoutlog, "\n");
	           fclose(stdoutlog);
	}

	*/

	if (0 != 0 && enindex[index] != 0)
	{
		if (checks[index] < 8 && sqrt(rsum / ((double) chain_length)) / ((double) chain_length) > (lowchecks[index] *((double) enindex[index])) &&

			sqrt(ssum / (chain_length - 2)) > ((lowchecks[index]) *((double) enindex2[index]))

			&&
			ete > ((lowchecks[index]) *((double) enindex3[index])))

		{

			if (sqrt(rsum / ((double) chain_length)) / ((double) chain_length) < ((highchecks[index]) *((double) enindex[index])) &&
				sqrt(ssum / (chain_length - 2)) < ((highchecks[index]) *((double) enindex2[index])) &&
				ete < ((highchecks[index]) *((double) enindex3[index])))
			{
				lowchecks[index] = ((lowchecks[index] + (.5 / checks_countdown[index])));
				highchecks[index] = ((highchecks[index] - (.5 / checks_countdown[index])));
				checks[index]++;

				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				   fprintf(stdoutlog, "index: %i \t lowchecks: %f \t highchecks: %f \t lc: %f \t hc: %f \t ccoutdown: %f ",index,lowchecks[index]+(.5/checks_countdown[index]),highchecks[index]-(.5/checks_countdown[index]),lowchecks[index],highchecks[index],checks_countdown[index]);
				    fprintf(stdoutlog, "\n");
				    fclose(stdoutlog);
				*/
				sprintf(filename, "En%i_%i_checks%i_CData_58mer.txt", -index, myid, checks[index]);	// Preliminary output file that allows for process inspection
				if ((file = fopen(filename, "w")) == NULL)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
					fclose(stdoutlog);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
				{
					for (int i = 0; i < chain_length; i++)
					{
						fprintf(file, "%i\n%i\n%i\n%i\n", pivot_indexes[3 *i + 0], pivot_indexes[3 *i + 1], pivot_indexes[3 *i + 2], chain_sequence[i]);
					}

					fclose(file);
				}

				sprintf(filename, "En%i_%i_CCord.txt", index, myid);	// Preliminary output file that allows for process inspection
				if ((file = fopen(filename, "w")) == NULL)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
					fclose(stdoutlog);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
				{
					fprintf(file, "checks=%i\n", checks[index]);
					for (int i = 0; i < chain_length; i++)
					{
						fprintf(file, "%i\t%i\t%i\n", pivot_indexes[3 *i + 0], pivot_indexes[3 *i + 1], pivot_indexes[3 *i + 2]);
					}

					fclose(file);
				}

				checks_countdown[index] *= 2;
				// sysparam(index);
			}
		}

		//.99609375
		//1.00390625
		//0.96875
		//1.03125
		if (checks[index] >= 8 && sqrt(rsum / ((double) chain_length)) / ((double) chain_length) > (lowchecks[index] *((double) enindex[index])) &&

			sqrt(ssum / (chain_length - 2)) > ((lowchecks[index]) *((double) enindex2[index]))

			&&
			ete > ((.99609375) *((double) enindex3[index])))

		{

			if (sqrt(rsum / ((double) chain_length)) / ((double) chain_length) < ((highchecks[index]) *((double) enindex[index])) &&
				sqrt(ssum / (chain_length - 2)) < ((highchecks[index]) *((double) enindex2[index])) &&
				ete < (1.00390625 *((double) enindex3[index])))
			{
				lowchecks[index] = ((lowchecks[index] + (.5 / checks_countdown[index])));
				highchecks[index] = ((highchecks[index] - (.5 / checks_countdown[index])));
				checks[index]++;

				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				   fprintf(stdoutlog, "index: %i \t lowchecks: %f \t highchecks: %f \t lc: %f \t hc: %f \t ccoutdown: %f ",index,lowchecks[index]+(.5/checks_countdown[index]),highchecks[index]-(.5/checks_countdown[index]),lowchecks[index],highchecks[index],checks_countdown[index]);
				    fprintf(stdoutlog, "\n");
				    fclose(stdoutlog);
				*/
				sprintf(filename, "En%i_%i_checks%i_CData_58mer.txt", -index, myid, checks[index]);	// Preliminary output file that allows for process inspection
				if ((file = fopen(filename, "w")) == NULL)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
					fclose(stdoutlog);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
				{
					for (int i = 0; i < chain_length; i++)
					{
						fprintf(file, "%i\n%i\n%i\n%i\n", pivot_indexes[3 *i + 0], pivot_indexes[3 *i + 1], pivot_indexes[3 *i + 2], chain_sequence[i]);
					}

					fclose(file);
				}

				sprintf(filename, "En%i_%i_CCord.txt", index, myid);	// Preliminary output file that allows for process inspection
				if ((file = fopen(filename, "w")) == NULL)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
					fclose(stdoutlog);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
				{
					fprintf(file, "checks=%i\n", checks[index]);
					for (int i = 0; i < chain_length; i++)
					{
						fprintf(file, "%i\t%i\t%i\n", pivot_indexes[3 *i + 0], pivot_indexes[3 *i + 1], pivot_indexes[3 *i + 2]);
					}

					fclose(file);
				}

				checks_countdown[index] *= 2;
				//  sysparam(index);
			}
		}
	}

	tortuosity[index] += sqrt(ssum / (chain_length - 2));
	visits[index]++;

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, ".8 poly_en():%i",poly_en());
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
}

void eye()
{
	for (int i = 0; i < chain_length; i++)
	{
		digital_eye[i] = -10000000;
	}

	int foundd = 0;
	for (int ii = 0; ii < chain_length; ii++)
	{
		for (int iii = 0; iii < numberneighbors; iii++)
		{
			if (latticepoint[neighbor[numberneighbors *poly_chain[ii] + iii]] > 0)
			{
				for (int i = 0; i < chain_length; i++)
				{
					if (neighbor[numberneighbors *poly_chain[ii] + iii] == poly_chain[i])
					{
						foundd = 1;
					}
				}

				if (foundd == 0)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "eye 1 Abort   happened:%i happened2:%i  offender:%i  check:%i \n\n", happened, happened2, neighbor[numberneighbors *poly_chain[ii] + iii], latticepoint[neighbor[numberneighbors *poly_chain[ii] + iii]]);
					for (int q = 0; q < chain_length; q++)
					{
						fprintf(stdoutlog, "%4i ", poly_chain[q]);
					}

					fclose(stdoutlog);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}
		}
	}

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, " 936 (%i %i %i %i %i %i)",neighbor[numberneighbors *936 + 0],neighbor[numberneighbors *936 + 1],neighbor[numberneighbors *936 + 2],neighbor[numberneighbors *936 + 3],neighbor[numberneighbors *936 + 4],neighbor[numberneighbors *936 + 5]);
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, " 869 (%i %i %i %i %i %i)",neighbor[numberneighbors *869 + 0],neighbor[numberneighbors *869 + 1],neighbor[numberneighbors *869 + 2],neighbor[numberneighbors *869 + 3],neighbor[numberneighbors *869 + 4],neighbor[numberneighbors *869 + 5]);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	for (int i = 0; i < chain_length; i++)
	{
		latticepoint[poly_chain[i]] = 0;
	}

	for (int i = 0; i < chain_length; i++)
	{
		if (latticepoint[poly_chain[i]] != 0)
		{
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "%i real eye 3 Abort", happened3);

			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}

			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		latticepoint[poly_chain[i]] = chain_sequence[i];
	}

	for (int i = 0; i < chain_length; i++)
	{
		digital_eye[i] = poly_chain[i];
		for (int ii = 0; ii < chain_length; ii++)
		{
			int neighbor_check = 2;
			if (ii > 0 && ii != chain_length - 1)
			{
				neighbor_check = 0;
				for (int iii = 0; iii < numberneighbors; iii++)
				{
					if (poly_chain[ii - 1] == neighbor[numberneighbors *poly_chain[ii] + iii])
					{
						neighbor_check++;
					}

					if (poly_chain[ii + 1] == neighbor[numberneighbors *poly_chain[ii] + iii])
					{
						neighbor_check++;
					}
				}

				if (neighbor_check != 2)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "%i real eye 2 Abort", happened3);

					fprintf(stdoutlog, "\n");
					for (int i = 0; i < chain_length; i++)
					{
						fprintf(stdoutlog, "%4i", poly_chain[i]);
					}

					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);

					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}

			if (digital_eye[i] != -10000000 && i != ii && neighbor_check != 2)
			{
				if (digital_eye[i] == digital_eye[ii] || digital_eye[i] < 0)
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "eye 2 Abort");
					fclose(stdoutlog);

					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}
		}
	}
}

void init_lattice(double emin, double emax)	// Changes made to bctype check and made latticepoint[numbersign equal to 0
{
	int e, r;

	latticepoint = (int*) malloc((numberspins + 2 + 1) *sizeof(int));
	// Note: we reserve space for 3 extra 'spins':
	// 2 extra to store fix values used for certain boundary conditions
	// 1 extra to carry an replica-id flag

	// find a fast way to create valid initial configurations
	// start at either Eglobalmin or Eglobalmax and change spins randomly
	// until energy is in the middle third of local energy range
	// global maximum of g(E) id at E/N=-0.2 

	//initialize lattice

	for (int i = 0; i < numberspins; i++)
	{
		latticepoint[i] = 0;
	}

	init_poly_chain();

	int populator = 1;
	for (int i = 0; i < chain_length; i++)
	{
		//latticepoint[poly_chain[i]] = populator;
		//populator++;
		latticepoint[poly_chain[i]] = chain_sequence[i];
	}

	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < numberspins; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", latticepoint[i]);
	//        if ((i + 1) % L1dim == 0)
	//            fprintf(stdoutlog, "\n");
	//
	//    }

	fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < numberspins; i++)
	//    {
	//        fprintf(stdoutlog, "%i (%i %i %i %i %i %i)",i,neighbor[numberneighbors *i + 0],neighbor[numberneighbors *i + 1],neighbor[numberneighbors *i + 2],neighbor[numberneighbors *i + 3],neighbor[numberneighbors *i + 4],neighbor[numberneighbors *i + 5]);
	//        if ((i + 1) % L1dim == 0)
	//            fprintf(stdoutlog, "\n");
	//
	//    }

	//fclose(stdoutlog);

	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);

	lattice_polymer();	//polymer function instancing moved before inside reference lattice
	//pull_1(); pull_relaxation();
	//slither_f();
	//slither_b();
	//cornerflip();
	//h_t_wiggle();
	//rebridging();
	//rebridging_h_t();

	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "cornerflip\n");
	//    for (int i = 0; i < numberspins; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", latticepoint[i]);
	//        if ((i + 1) % L1dim == 0)
	//            fprintf(stdoutlog, "\n");
	//
	//    }

	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "%i\n", poly_en());
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);

	/*for (int t; t < 10000; t++) {
		poly_mov();

		if(t%1000==0){
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		for (int i = 0; i < numberspins; i++)
		{
			fprintf(stdoutlog, "%4i", latticepoint[i]);
			if ((i + 1) % L1dim == 0)
				fprintf(stdoutlog, "\n");
		}

		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "%i\n",poly_en());
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		}
	}

	*/

	//    for (int i = 0; i < numberspins; i++)
	//    {
	//        latticepoint[i] = 0;
	//    }

	//populator = 1;
	for (int i = 0; i < chain_length; i++)
	{
		//latticepoint[poly_chain[i]] = populator;
		//populator++;
		latticepoint[poly_chain[i]] = chain_sequence[i];
	}

	e = poly_en();

	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "Proc. %i: Initialized lattice with energy e=%i, create setup with %lf<e<%lf\n", myid, e, (emin + (emax - emin) / 3), (emin + 2 *(emax - emin) / 3));
	fclose(stdoutlog);

	//for (int i = 0; i < numberspins; i++)
	//{
	//	latticepoint[i] = (i + count) % (q - 1);
	//	if ((i + 1) % L1dim == 0)
	//	{
	//		int rhld = rand() % (q - 1);
	//		while (rhld == count)
	//		{
	//			rhld = rand() % (q - 1);
	//			count = rhld;
	//		}

	//	}

	//}

	/*
	int count = 0;

	for (int i = 0; i < numberspins; i++)
	{
		latticepoint[i] = (i + count) % (q - 1);
		if ((i + 1) % L1dim == 0)
			count++;
	}

	*/
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < numberspins; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", latticepoint[i]);
	//        if ((i + 1) % L1dim == 0)
	//            fprintf(stdoutlog, "\n");
	//
	//    }

	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "%i", poly_en());
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);

	// run to a valid energy for the energy window

	// as long as energy is outside the middle third of local energy range
	long increm = 0;

	while (((e < (emin + (emax - emin) / 7)) || (e > (emin + 6 *(emax - emin) / 7))))
	{
		poly_mov();
		e = poly_en();
		//eye();
		increm = increm + 1;
		if (increm % 70000 == 0)
		{
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "increment: %li", increm);

			//stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < numberspins; i++)
			//            {            	//                fprintf(stdoutlog, "%4i", latticepoint[i]);
			//                if ((i + 1) % L1dim == 0)
			//                    fprintf(stdoutlog, "\n");
			//
			//            }

			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "%i\n", poly_en());
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i\t", poly_chain[i]);
			}
		}
	}

	if (bctype == 1)
	{
		latticepoint[numberspins] = 1;	// Braskamp Kunz has fixed spin values at boundaries
		latticepoint[numberspins + 1] = -1;
	}

	// use last spin to store replica-id; this is just a marker to
	// allows us to keep track of replicas as they are exchanged
	latticepoint[numberspins + 2] = myid;
	// Pseudo WL Process is started for all processes to explore the energy levels
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\nPseudo WL Proces has started for process: %i\n", myid);
	fclose(stdoutlog);
	//if (myid == 0 || myid == 1 || myid == 2) lngE[4] = 1;
	//init_check();
	for (int i = Eminindex; i <= Emaxindex; i++)
	{
		visits[i] = 0;
		visits_buf[i] = 0;
		rog[i] = 0.0;
		rog_buf[i] = 0.0;
		tortuosity[i] = 0.0;
		tortuosity_buf[i] = 0.0;
		endtoend[i] = 0.0;
		endtoend_buf[i] = 0.0;
	}

	pseudowl();
	//lngE[4]

	stdoutlog = fopen(stdoutlogname, "a");
	if (local_dup_group == localgroup && merge_hists == 1)	// merge g(E) estimators from multiple walkers in the same energy window
	{
		if (localgroup == 6) lngE[44] = 1;
		fprintf(stdoutlog, "ugig");
		if (myid == localgroup)	// 'root' in energy window, receive individual g(E) and send merged g(E)
		{
			fprintf(stdoutlog, "upig");
			for (int i = (dup_headproc + 1); i < numprocs; i++)
			{
				MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, i, 76, MPI_COMM_WORLD, &status);	// get other dens. of states
				fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, i);
				for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j];	// sum up for average
			}

			for (int i = (dup_headproc + 1); i < numprocs; i++)
			{
				MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, i, 98, MPI_COMM_WORLD);
				fprintf(stdoutlog, "Proc %i: Sent combined lngE to Proc. %i\n", myid, i);
			}
		}
		else	// send individual g(E) and receive merged g(E)
		{
			fprintf(stdoutlog, "ufig");
			MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, dup_headproc, 76, MPI_COMM_WORLD);
			fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, dup_headproc);
			MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, dup_headproc, 98, MPI_COMM_WORLD, &status);
			fprintf(stdoutlog, "Proc %i: Received combined lngE from Proc. %i\n", myid, dup_headproc);
			for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];	// replace individual lngE (could be done directly, yes)
		}

		fclose(stdoutlog);
	}

	if (multiple > 1 && merge_hists == 1)	// merge g(E) estimators from multiple walkers in the same energy window
	{
		stdoutlog = fopen(stdoutlogname, "a");
		if (myid % multiple == 0)	// 'root' in energy window, receive individual g(E) and send merged g(E)
		{
			for (int i = 1; i < multiple; i++)
			{
				MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid + i, 76, MPI_COMM_WORLD, &status);	// get other dens. of states
				fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, myid + i);
				for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j];	// sum up for average
			}

			for (int i = 1; i < multiple; i++)
			{
				MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid + i, 98, MPI_COMM_WORLD);
				fprintf(stdoutlog, "Proc %i: Sent combined lngE to Proc. %i\n", myid, myid + i);
			}
		}
		else	// send individual g(E) and receive merged g(E)
		{
			MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 76, MPI_COMM_WORLD);
			fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, myid - (myid % multiple));
			MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 98, MPI_COMM_WORLD, &status);
			fprintf(stdoutlog, "Proc %i: Received combined lngE from Proc. %i\n", myid, myid - (myid % multiple));
			for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];	// replace individual lngE (could be done directly, yes)
		}

		fclose(stdoutlog);
	}

	//Print lattice to logfile
	stdoutlog = fopen(stdoutlogname, "a");
	for (int i = 0; i < numberspins; i++)
	{
		fprintf(stdoutlog, "%4i", latticepoint[i]);
		if ((i + 1) % L1dim == 0)
			fprintf(stdoutlog, "\n");
	}

	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "exit %i", poly_en());
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);
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

	if (fopen("Ewindows.dat", "r") == NULL)
	{
		fprintf(stdoutlog, "\nProc %i: Can't find file Ewindows.dat. Will calculate equal-size windows with overlap %lf\n", myid, overlap);
		double Ewidth = (Eglobmax - Eglobmin) / (1.0 + ((double)(N / multiple) - 1.0) *(1.0 - overlap));

		Emin = Eglobmin + (double)(myid / multiple) *(1.0 - overlap) *Ewidth;
		Emax = Emin + Ewidth;

		Eminindex = floor(Emin + -Eglobalmin);
		Emaxindex = ceil(Emax + -Eglobalmin);

		time(&timenow);
		fprintf(stdoutlog, "Proc %3i: Parameter: Eglobmin -- Eglobmax: %lf -- %lf; overlap=%i percent, Ewindowwidth=%lf, Emin -- Emax: %d -- %d;  %s", myid, Eglobmin, Eglobmax, (int)(overlap *100.0), Ewidth, Eminindex, Emaxindex, ctime(&timenow));
		//fprintf(stdoutlog, "Proc %3i: Parameter: Eglobmin -- Eglobmax: %lf -- %lf; %s", myid, Emin, Emax, ctime(&timenow));
	}
	else	// rewrite of prior file reading system to be more portable
	{
		FILE * fptr;
		int value;
		int valuee;
		int valueee;
		int valueeee;

		fptr = fopen("Ewindows_R.dat", "r");

		while (fscanf(fptr, "%i,%i,%i,%i\n", &value, &valuee, &valueee, &valueeee) > 0)
		{
			if (value == myid / multiple)
			{
				Eminindex = valuee;
				Emaxindex = valueee;
				localgroup = valueeee;
				fprintf(stdoutlog, "%i %i %i %i\n", value, Eminindex, Emaxindex, localgroup);
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
int replica_exchange(int swap_direction, int index_akt)
{
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < numberspins; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", latticepoint[i]);
	//        if ((i + 1) % L1dim == 0)
	//            fprintf(stdoutlog, "\n");
	//
	//    }

	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "%i\n", poly_en());
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);

	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "1\n");
	fclose(stdoutlog);

	int i_new;	// histogram index of my configuration
	int Ecur;	// current energy

	// frac refers to local exchange probability
	// wk is combined exchange probability
	double myfrac, otherfrac, randx, wk;

	int change = 0;	// boolean: 0 for not exchanging, 1 for exchanging
	int swap_partner = -1;	// id of swap partner (receive_buffer)

	swap_direction = swap_direction % 2;	// comes actually as number of swap attempt

	// everyone has to find its swap-partner

	int *pairs;	// array containing the partners of each process (send_buffer)
	pairs = (int*) malloc(2 *multiple* sizeof(int));

	if (mylocalid[swap_direction] == 0)	// 'head-node' in the energy window determines pairs of flippartners
	{
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "2\n");
		fclose(stdoutlog);

		int choose_from = multiple;	// number of free partners in higher window of communicator
		int select;	// storage for random number

		int *libre;	// list of free partners from higher window in communicator
		libre = (int*) malloc(multiple* sizeof(int));

		for (int i = 0; i < multiple; i++) libre[i] = multiple + i;	// initialise

		// idea: processes from the lower window choose someone from the higher window at random
		// of course, the chosen walker can't have an exchange partner yet
		for (int i = 0; i < multiple; i++)	// loop over processes in the lower window
		{
			select = rand() % choose_from;
			pairs[i] = libre[select];
			pairs[libre[select]] = i;	// the 'vice-versa pair'
			// update list
			choose_from--;
			for (int j = select; j < choose_from; j++)
				libre[j] = libre[j + 1];
		}

		//       stdoutlog=fopen(stdoutlogname,"a");
		//       fprintf(stdoutlog,"Proc %3i: Drew the following swap partners:\n",myid);
		//       for (int i=0;i<2*multiple;i++)
		// 	fprintf(stdoutlog,"Proc %3i: %i -- %i (local ids in communicator)\n",myid,i,pairs[i]);
		//       fclose(stdoutlog);

		free(libre);
	}

	// at this point, every walker has a swap partner assigned, now they must be communicated
	if ((swap_direction == 0) && (myid < (numprocs - multiple)))	// the walkers from the last node should not swap
	{
		comm_id = 2 *(myid / (2 *multiple));	// ! all integer, the '/' is a (div) ! Not the same as myid/multiple !
		MPI_Scatter(pairs, 1, MPI_INT, &swap_partner, 1, MPI_INT, 0, mpi_local_comm[comm_id]);
	}

	if ((swap_direction == 1) && (myid >= multiple))	// the walkers from the zero-node should not swap
	{
		comm_id = ((myid - multiple) / (2 *multiple)) *2 + 1;	// ! all integer, the '/' is a (div) ! See above
		MPI_Scatter(pairs, 1, MPI_INT, &swap_partner, 1, MPI_INT, 0, mpi_local_comm[comm_id]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	free(pairs);

	if (swap_partner != -1)	// i.e. if there is a swap-partner for me (if I am at a boundary, I might not have a swap partner this time)
	{
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "3\n");
		fclose(stdoutlog);

		// statistics
		if (swap_partner > mylocalid[swap_direction]) tryright++;
		else tryleft++;

		// safety cross check
		Ecur = poly_en();
		//eye();
		if (Ecur + (-Eglobalmin) != index_akt)
		{
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "Proc %3i, replica_exchange(): Something is wrong here! Received index=%i, calculated index=%i totalenergy=%i. Abort.\n", myid, index_akt, Ecur + (-Eglobalmin), poly_en());
			fclose(stdoutlog);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		i_new = index_akt;

		// get histogram index from my swap partner
		MPI_Sendrecv_replace(&i_new, 1, MPI_INT, swap_partner, 1, swap_partner, 1, mpi_local_comm[comm_id], &status);

		if ((i_new > Emaxindex) || (i_new < Eminindex))	// energyranges must overlap!
		{
			myfrac = -1.0;
		}
		else
		{
			// calculate my part of the exchange probability
			myfrac = exp(lngE[index_akt] - lngE[i_new]);	// g(myE)/g(otherE)
		}

		if (mylocalid[swap_direction] < multiple)	// I am receiver and calculator
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
		else	// I just send my part of exchange probability and await decision
		{
			MPI_Send(&myfrac, 1, MPI_DOUBLE, swap_partner, 2, mpi_local_comm[comm_id]);
			MPI_Recv(&change, 1, MPI_INT, swap_partner, 3, mpi_local_comm[comm_id], &status);
		}

		// if decision was made to exchange configurations
		if (change == 1)
		{
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "4\n");
			fclose(stdoutlog);

			// exchange spin conformations (incl. the 3 'special' lattice sites)
			//MPI_Sendrecv_replace(&latticepoint[0], numberspins + 2 + 1, MPI_INT, swap_partner, 1, swap_partner, 1, mpi_local_comm[comm_id], &status);

			MPI_Sendrecv_replace(&poly_chain[0], chain_length, MPI_INT, swap_partner, 1, swap_partner, 1, mpi_local_comm[comm_id], &status);

			switched = i_new;

			// statistics
			if (swap_partner > mylocalid[swap_direction]) exchangeright++;
			else exchangeleft++;
		}
	}

	return (change);
	// returns whether or not configs were actually exchanged
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

	int restart = atoi(argv[7]);
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
	//rseed = 1619117003;
	srand(rseed + myid);

	int swap_every = atoi(argv[3]);	// after this number of sweeps try conformations swap
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

		MPI_Abort(MPI_COMM_WORLD, 1);	// something went wrong in find_local_energy_range()
	}

	init_hists();	// moved above init_lattice() for calculation considerations
	init_neighbors();

	MPI_Barrier(MPI_COMM_WORLD);
	//if (stop) MPI_Abort(MPI_COMM_WORLD, 1);

	//Teststop
	//   MPI_Barrier(MPI_COMM_WORLD);
	//   MPI_Abort(MPI_COMM_WORLD,1);

	// create new groups and communicators for each energy range window
	stdoutlog = fopen(stdoutlogname, "a");
	MPI_Comm_group(MPI_COMM_WORLD, &world);	// get the group of processes in MPI_COMM_WORLD (i.e. all)
	int *ranks;
	ranks = (int*) malloc(2 *multiple* sizeof(int));
	mpi_local_group = (MPI_Group*) malloc(((numprocs / multiple) - 1) *sizeof(MPI_Group));
	mpi_local_comm = (MPI_Comm*) malloc(((numprocs / multiple) - 1) *sizeof(MPI_Comm));
	fprintf(stdoutlog, "ready");
	fclose(stdoutlog);
	for (int i = 0; i < ((numprocs / multiple) - 1); i++)	// i is counter for the energy range windows
	{
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "set");
		fclose(stdoutlog);
		for (int j = 0; j < 2 * multiple; j++)
		{
			ranks[j] = i *multiple + j;	// contains the ranks of processes in MPI_COMM_WORLD which should get into local group
			if (myid == 0)
			{
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "Proc %3i: %i will be part of communicator/group %i\n", myid, ranks[j], i);
				fclose(stdoutlog);
			}
		}

		MPI_Group_incl(world, 2 *multiple, ranks, &mpi_local_group[i]);	// create local group
		MPI_Comm_create(MPI_COMM_WORLD, mpi_local_group[i], &mpi_local_comm[i]);	// create communicator for that group
	}

	free(ranks);

	stdoutlog = fopen(stdoutlogname, "a");
	// get my local id (in my local communicators)
	if (myid < numprocs - multiple)
	{
		comm_id = 2 *(myid / (2 *multiple));
		MPI_Comm_rank(mpi_local_comm[comm_id], &mylocalid[0]);
		fprintf(stdoutlog, "Proc %3i: I am part of communicator/group %i with local_id[0]=%i\n", myid, comm_id, mylocalid[0]);
	}
	else
	{
		mylocalid[0] = INT_MAX;	// just to give it a value
		fprintf(stdoutlog, "Proc %3i: got local_id[0]=%i\n", myid, mylocalid[0]);
	}

	fclose(stdoutlog);

	stdoutlog = fopen(stdoutlogname, "a");
	if (myid >= multiple)
	{
		comm_id = 2 *((myid - multiple) / (2 *multiple)) + 1;
		MPI_Comm_rank(mpi_local_comm[comm_id], &mylocalid[1]);
		fprintf(stdoutlog, "Proc %3i: I am part of communicator/group %i with local_id[1]=%i\n", myid, comm_id, mylocalid[1]);
	}
	else
	{
		mylocalid[1] = INT_MAX;	// just to give it a value
		fprintf(stdoutlog, "Proc %3i: got local_id[1]=%i\n", myid, mylocalid[1]);
	}

	fclose(stdoutlog);

	init_lattice(Emin, Emax);	// 0 - random; 1 - all equal

	// calculate energy for the first time
	int eold, energie;
	energie = poly_en();
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "Proc %3i: energy at start=%i\n", myid, energie);
	fclose(stdoutlog);
	energie = poly_en() - Eglobalmin;	// shift to positive values to use it as array index

	Estartindex = energie;

	int stop = 0;
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "Proc %3i: Parameter: Eminindex -- Emaxindex: %i -- %i; Estartindex=%i\t", myid, Eminindex, Emaxindex, Estartindex);
	if ((Estartindex > Emaxindex) || (Estartindex < Eminindex))
	{
		fprintf(stdoutlog, "Estart out of range!!!\n");
		stop = 0;	// not needed
	}
	else
	{
		fprintf(stdoutlog, "OK\n");
	}

	fclose(stdoutlog);

	MPI_Barrier(MPI_COMM_WORLD);
	if (stop) MPI_Abort(MPI_COMM_WORLD, 1);

	// log 'path'data (the same data is written, just differently sorted: with respect to the sample id and with respect to process id)
	sprintf(wanderlogname, "wanderlog_rseed%i_sample_%i.dat", rseed, latticepoint[numberspins + 2]);
	wanderlog = fopen(wanderlogname, "w");
	time(&timenow);
	fprintf(wanderlog, "#sweep\t\tproc-id\tconf-id\tenergy\t\ttime\n");
	fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", 0.0, myid, latticepoint[numberspins + 2], energie, ctime(&timenow));
	fclose(wanderlog);

	sprintf(wanderlogname, "wanderlog_rseed%i_walker_%i.dat", rseed, myid);
	wanderlog = fopen(wanderlogname, "w");
	time(&timenow);
	fprintf(wanderlog, "#sweep\t\tproc-id\tconf-id\tenergy\t\ttime\n");
	fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", 0.0, myid, latticepoint[numberspins + 2], energie, ctime(&timenow));
	fclose(wanderlog);

	//start simulation
	double wsk, dice;	// wsk: probability  
	int wiggle;
	int wiggletwo;	// ID of spin to be updated (only single-spin update implemented)

	// the WL parameter should eventually go into the init file, too
	double countdown = 2;
	double lnf = 1.0;	// my modification factor
	double lnf_slowest = lnf;	// modification factor of slowest walker in my window
	//double lnfmin=log(1.000000001);
	double lnfmin = log(1.0000000001);	// terminal modification factor
	double sweep = 0;	// counter for MC sweeps
	int flat;	// 0 - histogram not flat; 1 - flat
	int iteration = 1;	// WL iteration counter
	int iter = 1;
	double check_flatness_every = 10000;	// in number of sweeps
	int backup;
	int backuptwo;

	eold = energie;

	if (restart == 1)
	{
		int Ewidth = Emaxindex - Eminindex;

		for (int u = 0; u < chain_length; u++)
		{
			wl_pseudo_chain[u] = poly_chain[u];
		}

		sprintf(resetbuffer, "SysStatus.L%i.proc%04i", L1dim, myid);	// Preliminary output file that allows for process inspection
		FILE * rfptr;
		long long phldr;

		int fakeEmin1 = Eminindex;
		int fakeEmin2 = Eminindex;
		int fakeEmin3 = Eminindex;
		int fakeEmin4 = Eminindex;
		int fakeEmin5 = Eminindex;
		int fakeEmin6 = Eminindex;
		int fakeEmin7 = Eminindex;
		int chain_increment = 0;
		rfptr = fopen(resetbuffer, "r");

		int resetiter = 1;

		while (fscanf(rfptr, "%lld\n", &phldr) > 0)
		{
			if (resetiter == 1)
			{
			 	//ratiovalue =0.9956;
			}

			if (resetiter == 2)
			{
				energie = phldr;
			}

			if (resetiter == 3)
			{
				eold = phldr;
			}

			if (resetiter == 4)
			{
				sweep = phldr;
			}

			if (resetiter == 5)
			{
				swap_every = phldr;
			}

			if (resetiter == 6)
			{
			 	//lnf = ((double)phldr) / pow(10, 12);
			}

			if (resetiter == 7)
			{
				iteration = phldr;
			}

			if (resetiter >= (8) && resetiter < (8 + Ewidth + 1))
			{
				HE[fakeEmin1] = ((int) phldr);

				fakeEmin1++;
			}

			if (resetiter >= (8 + Ewidth + 1) && resetiter < ((8 + Ewidth + 1) + Ewidth + 1))
			{
				lngE[fakeEmin2] = ((double) phldr) / pow(10, 9);
				fakeEmin2++;
			}

			if (resetiter >= ((8 + Ewidth + 1) + Ewidth + 1) && resetiter < (((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1))
			{
				rog[fakeEmin2] = ((double) phldr) / pow(10, 6);
				fakeEmin3++;
			}

			if (resetiter >= (((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) && resetiter < ((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1))
			{
				tortuosity[fakeEmin2] = ((double) phldr) / pow(10, 6);
				fakeEmin4++;
			}

			if (resetiter >= ((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) && resetiter < (((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1))
			{
				visits[fakeEmin2] = ((long long) phldr);
				fakeEmin5++;
			}

			if (resetiter >= (((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) && resetiter < ((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + chain_length))
			{
				poly_chain[chain_increment] = ((int) phldr);
				chain_increment++;
			}

			if (resetiter >= ((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + chain_length) && resetiter < (((((((8 + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + Ewidth + 1) + chain_length) + Ewidth + 1))
			{
				endtoend[chain_increment] = ((double) phldr) / pow(10, 6);
				fakeEmin6++;
			}

			resetiter++;

			fclose(rfptr);
		}

		double counter_iteration = iteration;
		while (counter_iteration != 1)
		{
			lnf /= 2.;
			counter_iteration -= 1;
		}

		lp_reorg_non_rebrid(wl_pseudo_chain, poly_chain, 0, chain_length);	// rejection needs to change the indexs back to the former storage value
		int pop = 1;
		lattice_reorg_non_rebrid(wl_pseudo_chain, poly_chain, pop, 0, chain_length);
		// makes a copy of the lattice if the change is rejected
		for (int u = 0; u < chain_length; u++)
		{
			wl_pseudo_chain[u] = poly_chain[u];
			latticepoint[poly_chain[u]] = chain_sequence[u];
		}
	}

	int swtch;
	int found = 0;

	int prelim_counter = 0;
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "Proc %3i: Start WL iteration\n", myid);
	fclose(stdoutlog);

	long long restartiter = 0;

	for (int i = 0; i <= -Eglobalmin; i++) HE[i] = 0;	//init H(E)

	// Main Wang-Landau routine
	while (lnf_slowest > lnfmin)
	{
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "Proc %3i: Start WL iteration\n", myid);
		fclose(stdoutlog);
		if (restartiter % 500 == 0)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			sprintf(resetbuffer, "SysStatus.L%i.proc%04i", L1dim, myid);	// Preliminary output file that allows for process inspection

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
				fprintf(filee, "%lld\n", (long long)(.99* pow(10, 17)));
				fprintf(filee, "%i\n", energie);
				fprintf(filee, "%i\n", eold);
				fprintf(filee, "%lld\n", (long long)(sweep));
				fprintf(filee, "%i\n", swap_every);
				fprintf(filee, "%lld\n", (long long)(lnf* pow(10, 12)));
				fprintf(filee, "%i\n", iteration);

				for (int i = Eminindex; i <= Emaxindex; i++)
				{
					fprintf(filee, "%lld\n", (long long) HE[i]);
				}

				for (int i = Eminindex; i <= Emaxindex; i++)
				{
					fprintf(filee, "%lld\n", (long long)(lngE[i] *pow(10, 9)));
				}

				for (int i = Eminindex; i <= Emaxindex; i++)
				{
					fprintf(filee, "%lld\n", (long long)(rog[i] *pow(10, 6)));
				}

				for (int i = Eminindex; i <= Emaxindex; i++)
				{
					fprintf(filee, "%lld\n", (long long)(tortuosity[i] *pow(10, 6)));
				}

				for (int i = Eminindex; i <= Emaxindex; i++)
				{
					fprintf(filee, "%lld\n", (long long)(visits[i]));
				}

				for (int i = 0; i < (chain_length); i++)
				{
					fprintf(filee, "%lld\n", (long long) poly_chain[i]);
				}

				for (int i = Eminindex; i <= Emaxindex; i++)
				{
					fprintf(filee, "%lld\n", (long long)(endtoend[i] *pow(10, 6)));
				}

				fclose(filee);
			}

			restartiter = 0;
		}

		restartiter++;

		for (int k = 0; k < check_flatness_every; k++)
		{
			for (int i = 0; i < L1dim; i++)	// this does 1 MC sweep
			{
				energie = eold + propose_update(eold);	// calculate energy of updated configuration
				//eye();
				//if(MoveProposal == 1)
				//{ 	/*stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "energie: %i  actual: %i",energie,poly_en()-Eglobalmin);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				/*
				if(energie!=poly_en()-Eglobalmin)
				{
					swap_every = swap_every_init;	// reset countdown clock
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					for (int i = 0; i < numberspins; i++)
					{
						fprintf(stdoutlog, "%4i", latticepoint[i]);
						if ((i + 1) % L1dim == 0)
							fprintf(stdoutlog, "\n");
					}

					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "en: %i, poly_en()-Eglobmin: %i , energie:%i  ,eold : %i, happened: %i  happened_2: %i ",poly_en(),poly_en()-Eglobalmin,energie,eold,happened,happened2);
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					for (int i = 0; i < chain_length; i++)
					{
						fprintf(stdoutlog, "%4i", poly_chain[i]);
					}

					fclose(stdoutlog);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}

				 */

				//eye();
				//poly_mov();
				//energie = poly_en()-Eglobalmin;

				// reject because you got out of bounds of energy window
				if ((energie > Emaxindex) || (energie < Eminindex))
				{
					lp_reorg_non_rebrid(poly_chain, wl_pseudo_chain, 0, chain_length);	// rejection needs to change the indexs back to the former storage value
					int pop = 1;
					lattice_reorg_non_rebrid(poly_chain, wl_pseudo_chain, pop, 0, chain_length);
					// makes a copy of the lattice if the change is rejected
					for (int u = 0; u < chain_length; u++)
					{
						poly_chain[u] = wl_pseudo_chain[u];
						latticepoint[poly_chain[u]] = chain_sequence[u];
					}

					//                    for (int i = 0; i < numberspins; i++)
					//                    {                     		//                        latticepoint[i] = 0;
					//                    }

					pop = 1;
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "reject wl");
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					energie = eold;
				}
				else	// calculate acceptance propbability
				{
					dice = (1.0* rand() / (RAND_MAX + 1.0));	// roll the dice
					wsk = exp(lngE[eold] - lngE[energie]);	// calculate acceptance probability
					if (dice > wsk || MoveProposal == 0)	// reject
					{
						if (MoveProposal > 0)
						{
							lp_reorg_non_rebrid(poly_chain, wl_pseudo_chain, 0, chain_length);	// rejection needs to change the indexs back to the former storage value
							int pop = 1;
							lattice_reorg_non_rebrid(poly_chain, wl_pseudo_chain, pop, 0, chain_length);

							for (int u = 0; u < chain_length; u++)
							{
								poly_chain[u] = wl_pseudo_chain[u];
								latticepoint[poly_chain[u]] = chain_sequence[u];
							}

							//                        for (int i = 0; i < numberspins; i++)
							//                        {                         				//                            latticepoint[i] = 0;
							//                        }

							//sysparam(energie);
						}

						energie = eold;
					}
					else
					{
						sysparam(energie);
						eold = energie;	// accept
						for (int u = 0; u < chain_length; u++)
						{
							wl_pseudo_chain[u] = poly_chain[u];
						}
					}
				}

				// update histograms
				//}

				lngE[energie] += lnf;
				HE[energie]++;
			}	// end 1 sweep

			sweep++;	// sweep counter

			swap_every--;	// RE counter (countdown)
			if (swap_every == 0)	// ignition! (time to try to exchange configurations)
			{
				swap_every = swap_every_init;	// reset countdown clock
				//stdoutlog = fopen(stdoutlogname, "a");
				//fprintf(stdoutlog, "\n");
				//                for (int i = 0; i < numberspins; i++)
				//                {                 	//                    fprintf(stdoutlog, "%4i", latticepoint[i]);
				//                    if ((i + 1) % L1dim == 0)
				//                        fprintf(stdoutlog, "\n");
				//
				//                }

				//fprintf(stdoutlog, "\n");
				//fprintf(stdoutlog, "\n");
				//fprintf(stdoutlog, "en: %i, poly_en()-Eglobmin: %i , energie:%i  ,eold : %i  ",poly_en(),poly_en()-Eglobalmin,energie,eold);
				//fprintf(stdoutlog, "\n");
				//fprintf(stdoutlog, "\n");
				//fclose(stdoutlog);
				if (replica_exchange(sweep / swap_every, eold) == 1)	// if configs were exchanged ('accept'-branch)
				{
				 		//stdoutlog = fopen(stdoutlogname, "a");
					//fprintf(stdoutlog, "\n");
					//for (int i = 0; i < chain_length; i++)
					//{ 		//    fprintf(stdoutlog, "%4i", wl_pseudo_chain[i]);
					//}

					//fprintf(stdoutlog, "\n");
					//fprintf(stdoutlog, "\n");
					//for (int i = 0; i < chain_length; i++)
					//{ 		//    fprintf(stdoutlog, "%4i", poly_chain[i]);
					//}

					//fclose(stdoutlog);

					lp_reorg_non_rebrid(wl_pseudo_chain, poly_chain, 0, chain_length);
					int pop = 1;
					lattice_reorg_non_rebrid(wl_pseudo_chain, poly_chain, pop, 0, chain_length);
					lp_reorg_non_rebrid(wl_pseudo_chain, poly_chain, 0, chain_length);
					//int pop = 1;
					lattice_reorg_non_rebrid(wl_pseudo_chain, poly_chain, pop, 0, chain_length);

					//                    for (int i = 0; i < numberspins; i++)
					//                    {                     		//                        latticepoint[i] = 0;
					//                    }

					pop = 1;
					for (int i = 0; i < chain_length; i++)
					{
						latticepoint[poly_chain[i]] = pop;
						pop++;
						latticepoint[poly_chain[i]] = chain_sequence[i];
					}

					for (int u = 0; u < chain_length; u++)
					{
						wl_pseudo_chain[u] = poly_chain[u];
					}

					//stdoutlog = fopen(stdoutlogname, "a");
					//fprintf(stdoutlog, "\n");
					//for (int i = 0; i < numberspins; i++)
					//{ 		//    fprintf(stdoutlog, "%4i", latticepoint[i]);
					//    if ((i + 1) % L1dim == 0)
					//        fprintf(stdoutlog, "\n");

					//}

					//fprintf(stdoutlog, "\n");
					//fprintf(stdoutlog, "%i\n", poly_en());
					//fprintf(stdoutlog, "\n");
					//fclose(stdoutlog);
					eye();

					if ((poly_en() - Eglobalmin) != switched) MPI_Abort(MPI_COMM_WORLD, 1);

					energie = switched;
					eold = energie;
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "Proc %3i: %i iteration, %e sweeps, Replica Exchange Success\n Energy Calculted: %i \t Recievced Index Energy: %i\n", myid, iteration, sweep, energie, switched);
					fclose(stdoutlog);
				}

				// update histograms (independently of whether RE happened or not)
				lngE[energie] += lnf;
				HE[energie]++;

				// record position (process) of the actual sample
				sprintf(wanderlogname, "wanderlog_rseed%i_sample_%i.dat", rseed, latticepoint[numberspins + 2]);
				wanderlog = fopen(wanderlogname, "a");
				time(&timenow);
				fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", sweep, myid, latticepoint[numberspins + 2], eold, ctime(&timenow));
				fclose(wanderlog);

				// record actual sample id for each walker
				sprintf(wanderlogname, "wanderlog_rseed%i_walker_%i.dat", rseed, myid);
				wanderlog = fopen(wanderlogname, "a");
				time(&timenow);
				fprintf(wanderlog, "%e\t%i\t%i\t%i\t%s", sweep, myid, latticepoint[numberspins + 2], eold, ctime(&timenow));
				fclose(wanderlog);
			}
		}	// end computation block

		prelim_counter++;

		if (prelim_counter % 10 == 0)
		{
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim, q, myid, iteration);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim, q, myid, iteration);	// Preliminary output file that allows for process inspection
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
					fprintf(file, "%i\t%i\t%e\t%e\t%e\t%d\t%e\t%e\t%e\n", i, (i + (Eglobalmin)), 0.0, lngE[i], HE[i], 2, rog[i] / chain_length / visits[i], tortuosity[i] / visits[i], endtoend[i] / visits[i]);
				}

				fclose(file);
			}

			//fclose(file);
		}

		// check flatness
		flat = histflat(Eminindex, Emaxindex, 0.992);
		if (flat == 0)	// histogram not flat, do nothing except maybe some logging
		{
			stdoutlog = fopen(stdoutlogname, "a");	// Prints out effictively flatness progress for inspection
			fprintf(stdoutlog, "Proc %3i: %i iteration, %e sweeps, Histogram noch nicht flach\nValue of Tried Minimum Histogram Value: %f    Value of Tried Ratio: %f\n", myid, iteration, sweep, flatmin, flatratio);

			phrepv = flatmin;	// holds the minimum value of the histogram and checks change and prints out the current lattice if the value has reoccured frequently
			if (repv == phrepv)
			{
				iterr++;
				repv = phrepv;
			}
			else
			{
				repv = phrepv;
			}

			//            if (iterr == 20)
			//            {            	//                iterr = 0;
			//
			//                fprintf(stdoutlog, "\n");
			//
			//                for (int i = 0; i < numberspins; i++)
			//                {                	//                    fprintf(stdoutlog, "%4i", latticepoint[i]);
			//                    if ((i + 1) % L1dim == 0)
			//                        fprintf(stdoutlog, "\n");
			//
			//                }

			//
			//                fprintf(stdoutlog, "\n Sytem Energy Index: %i \n", energie);
			//
			//                fprintf(stdoutlog, "\n");
			//            }

			fclose(stdoutlog);
		}
		else	// histograms of all walkers are 'flat'
		{
			if (lnf > lnfmin)	// as long as modification factor is still larger than terminal value
			{
				sprintf(filename, "L%iq%i.HE.proc%04i.iter%i", L1dim, q, myid, iteration);
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

				//fclose(file);
			}

			// decrease modification factor
			// canonical method (reduce by sqrt(2)) implemented
			lnf /= 2.0;

			for (int i = 0; i <= -Eglobalmin; i++) HE[i] = 0;	//reset H(E)
			iteration++;	// iteration counter

			if (merge_hists == 1)	// merge g(E) estimators from multiple walkers in the same energy window
			{
				stdoutlog = fopen(stdoutlogname, "a");
				if (myid % multiple == 0)	// 'root' in energy window, receive individual g(E) and send merged g(E)
				{
					for (int i = 1; i < multiple; i++)
					{
						MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid + i, 77, MPI_COMM_WORLD, &status);	// get other dens. of states
						fprintf(stdoutlog, "Proc %i: Received lngE from Proc. %i\n", myid, myid + i);
						for (int j = 0; j < hist_size; j++) lngE[j] += lngE_buf[j];	// sum up for average
					}

					for (int j = 0; j < hist_size; j++) lngE[j] /= (double) multiple;	// normalize
					for (int i = 1; i < multiple; i++)
					{
						MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid + i, 99, MPI_COMM_WORLD);
						fprintf(stdoutlog, "Proc %i: Sent merged lngE to Proc. %i\n", myid, myid + i);
					}
				}
				else	// send individual g(E) and receive merged g(E)
				{
					MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 77, MPI_COMM_WORLD);
					fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, myid - (myid % multiple));
					MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, myid - (myid % multiple), 99, MPI_COMM_WORLD, &status);
					fprintf(stdoutlog, "Proc %i: Received merged lngE from Proc. %i\n", myid, myid - (myid % multiple));
					for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j];	// replace individual lngE (could be done directly, yes)
				}

				fclose(stdoutlog);
			}
		}

		// check progress from all other windows
		MPI_Allreduce(&lnf, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		if (lnf_slowest <= (((double)(10000000000 / countdown)) *.0000000001))
		{
			recombine((((double)(10000000000 / countdown)) *.0000000001));
			countdown *= 2;
			partial_init_hists();
		}

		// just some logging
		if (flat == 1)
		{
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "Proc %3i: Start %i iteration, %e sweeps total so far, lnf now %e (lnfmin=%.2e, lnf_slowest=%e)\n", myid, iteration, sweep, lnf, lnfmin, lnf_slowest);
			fprintf(stdoutlog, "Proc %3i: tryleft: %i, exchangeleft %i (Akzeptanzleft:%.2lf) < --> tryright: %i, exchangeright %i (Akzeptanzright:%.2lf)\n", myid, tryleft, exchangeleft, (double) exchangeleft / (double) tryleft, tryright, exchangeright, (double) exchangeright / (double) tryright);
			fclose(stdoutlog);
		}
	}	// end while(lnf_slowest>lnfmin) -> this terminates the simulation

	//   free(neighbor);
	//   free(latticepoint);

	// normalize results
	double norm = lngE[Eminindex] - log(q);
	for (int i = 0; i <= 2 *(-Eglobalmin); i++) lngE[i] -= norm;

	// write final data to file
	sprintf(filename, "L%iq%i.proc%04i.lngE", L1dim, q, myid);
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

	return (0);	// you did it!
}
