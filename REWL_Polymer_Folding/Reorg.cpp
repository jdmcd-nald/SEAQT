#include "Globals.h"

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


#include "Energy.h"
#include "Mov.h"
#include "Mov_Pivot.h"
#include "Mov_Pull.h"
#include "Mov_Rebridge.h"
#include "Reorg.h"
#include "EXWL.h"

int lp_reorg_non_rebrid(int* p_chain, int* c_chain, int start, int end) // reorganizes the lattice_polymer for non rebridging moves takes the input of the former chain and the updated chain
{
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "q\n");
	//fclose(stdoutlog);
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "\nhuh   poly_chain[1]: %i  pchain:%i\n",poly_chain[1],p_chain[1]);
	//fclose(stdoutlog);
	for (int i = start; i < end; i++)
	{
		latticepolymer[5 * p_chain[i]] = p_chain[i];     // lattice index
		latticepolymer[5 * p_chain[i] + 1] = -1; // next lattice point on chain (if not applicable -1 assigned)
		latticepolymer[5 * p_chain[i] + 2] = -1;
		latticepolymer[5 * p_chain[i] + 3] = -1;
		latticepolymer[5 * p_chain[i] + 4] = -1;
		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "qq:  %i     p_chain[i]:%i   p_chain[1]: %i\n",i,p_chain[i],p_chain[1]);
		//fclose(stdoutlog);
	}
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "qq\n");
	//fclose(stdoutlog);
	for (int i = start; i < end; i++)
	{
		latticepolymer[5 * c_chain[i]] = c_chain[i];     // lattice index
		latticepolymer[5 * c_chain[i] + 4] = i;
		if (i == 0)
		{
			latticepolymer[5 * c_chain[i] + 1] = c_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * c_chain[i] + 2] = -100; //previous lattice point on chain (if applicable -1 assigned)
		}
		if (i > 0 && i < (chain_length - 1))
		{
			latticepolymer[5 * c_chain[i] + 1] = c_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * c_chain[i] + 2] = c_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
		}
		if (i == (chain_length - 1))
		{
			latticepolymer[5 * c_chain[i] + 1] = -200; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * c_chain[i] + 2] = c_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
		}
	}
	return 0;
}

int lattice_reorg_non_rebrid(int* p_chain, int* c_chain, int fill, int start, int end) // reorganizes the latticepoints more efficient than previous scrubbing over entire lattice
{
	for (int i = start; i < end; i++)
	{
		latticepoint[p_chain[i]] = 0;     // lattice index
	}
	for (int i = start; i < end; i++)
	{
		//latticepoint[c_chain[i]] = fill;     // lattice index // need to replace with chain sequence value to model hydrophobic properties

		/*
		if(latticepoint[c_chain[i]]>0)
		{
			for (int ii = 0; ii < numberspins; ii++)
			{
				latticepoint[i] = 0;     // lattice index
			}
			for (int ii = 0; ii < chain_length; ii++)
			{
				latticepoint[p_chain[ii]] = chain_sequence[ii];     // lattice index
			}
			for (int y = 0; y < chain_length; y++)
			{
				poly_chain[y]=p_chain[y];
			}
			for (int ii = 0; ii < numberspins; ii++)
			{
				latticepolymer[5 * ii] = ii;     // lattice index
				latticepolymer[5 * ii+ 1] = -1; // next lattice point on chain (if not applicable -1 assigned)
				latticepolymer[5 * ii + 2] = -1;
				latticepolymer[5 * ii + 3] = -1;
			}

			for (int ii = 0; ii < chain_length; ii++)
			{
				latticepolymer[5 * p_chain[ii]] = p_chain[ii];     // lattice index

				if (ii == 0)
				{
					latticepolymer[5 * c_chain[ii] + 1] = p_chain[ii + 1]; // next lattice point on chain (if not applicable -1 assigned)
					latticepolymer[5 * c_chain[ii] + 2] = -100; //previous lattice point on chain (if applicable -1 assigned)
				}
				if (ii != 0 && i < (chain_length - 1))
				{
					latticepolymer[5 * c_chain[ii] + 1] = p_chain[ii + 1]; // next lattice point on chain (if not applicable -1 assigned)
					latticepolymer[5 * c_chain[ii] + 2] = p_chain[ii - 1]; //previous lattice point on chain (if applicable -1 assigned)
				}
				if (ii == (chain_length - 1))
				{
					latticepolymer[5 * c_chain[ii] + 1] = -200; // next lattice point on chain (if not applicable -1 assigned)
					latticepolymer[5 * c_chain[ii] + 2] = p_chain[ii - 1]; //previous lattice point on chain (if applicable -1 assigned)
				}
			}
			return 0;
		}
		*/

		//        if(latticepoint[c_chain[i]]>0)
		//        {
		//            for (int i = 0; i < chain_length; i++)
		//            {
		//                latticepoint[c_chain[i]] = 0;     // lattice index
		//            }
		//            for (int i = 0; i < chain_length; i++)
		//            {
		//                latticepoint[p_chain[i]] = chain_sequence[i];     // lattice index
		//            }
		//            for (int y = 0; y < chain_length; y++)
		//            {
		//                poly_chain[y]=p_chain[y];
		//            }
		//            return 0;
		//        }

		latticepoint[c_chain[i]] = chain_sequence[i];


		poly_chain[i] = c_chain[i];
	}
	return 0;
}

int reorg_wl(int start, int end) // reorganizes the lattice_polymer for non rebridging moves takes the input of the former chain and the updated chain
{
	for (int i = start; i < end; i++)
	{
		latticepolymer[5 * poly_chain[i]] = poly_chain[i];     // lattice index
		latticepolymer[5 * poly_chain[i] + 1] = -1; // next lattice point on chain (if not applicable -1 assigned)
		latticepolymer[5 * poly_chain[i] + 2] = -1;
		latticepolymer[5 * poly_chain[i] + 3] = -1;
		latticepolymer[5 * poly_chain[i] + 4] = -1;
		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "qq:  %i     p_chain[i]:%i   p_chain[1]: %i\n",i,p_chain[i],p_chain[1]);
		//fclose(stdoutlog);
		latticepoint[poly_chain[i]] = 0;     // lattice index
	}
	for (int i = start; i < end; i++)
	{
		latticepoint[wl_pseudo_chain[i]] = chain_sequence[i];
		poly_chain[i] = wl_pseudo_chain[i];

		latticepolymer[5 * wl_pseudo_chain[i]] = wl_pseudo_chain[i];     // lattice index
		latticepolymer[5 * wl_pseudo_chain[i] + 4] = i;
		if (i > 0 && i < (chain_length - 1))
		{
			latticepolymer[5 * wl_pseudo_chain[i] + 1] = wl_pseudo_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * wl_pseudo_chain[i] + 2] = wl_pseudo_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
		}
		else if (i == 0)
		{
			latticepolymer[5 * wl_pseudo_chain[i] + 1] = wl_pseudo_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * wl_pseudo_chain[i] + 2] = -100; //previous lattice point on chain (if applicable -1 assigned)
		}
		else if (i == (chain_length - 1))
		{
			latticepolymer[5 * wl_pseudo_chain[i] + 1] = -200; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * wl_pseudo_chain[i] + 2] = wl_pseudo_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
		}
	}
	return 0;
}

int lp_reorg_rebrid(int index, int h_t) //The reorganization program for the rebridgining operations would like for the tqo peograms to be joined but the rebridging is too different to be generalized
{
	if (h_t == 1) // assuming the argument index is the head
	{

		int starting_index = index; //starting index to check along the length of the chain
		int current_moving_index = index; // updating indexing for the moving check
		int prev_moving_index = -2000; // previous index to ensure the system doesnt count previous point
		int moved_check = 0;

		for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector to be updated
		{
			pseudo_poly_chain[y] = poly_chain[y];
		}

		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, " lp_reorg_rebrid");
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);

		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < chain_length; i++)
		//{
		//    fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < chain_length; i++)
		//{
		//    fprintf(stdoutlog, "%4i", poly_chain[i]);
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);

		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		for (int i = 0; i < chain_length; i++)
		{
			fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
		}
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/


		if (latticepolymer[5 * index + 1] >= 0)  // offsets the initial lattice count by one to start the sequence
		{
			prev_moving_index = current_moving_index;
			current_moving_index = latticepolymer[5 * index + 1];

			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " start 1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			latticepolymer[5 * prev_moving_index + 1] = current_moving_index;
			latticepolymer[5 * prev_moving_index + 2] = -100;
			latticepolymer[5 * prev_moving_index + 3] = -1;
			poly_chain[0] = prev_moving_index;

		}
		else
		{
			if (latticepolymer[5 * index + 2] >= 0)  // offsets the initial lattice count by one to start the sequence
			{
				prev_moving_index = current_moving_index;
				current_moving_index = latticepolymer[5 * index + 2];

				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " start 2 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/

				latticepolymer[5 * prev_moving_index + 1] = current_moving_index;
				latticepolymer[5 * prev_moving_index + 2] = -100;
				latticepolymer[5 * prev_moving_index + 3] = -1;
				poly_chain[0] = index;
			}
		}

		int i = 1;
		moved_check = 1;


		while (moved_check != 0) // increments through the chain flipping the lattice polymer points to be inline with expectations to reset the lattice
		{
			moved_check = 0;
			if (latticepolymer[5 * current_moving_index + 1] >= 0 && latticepolymer[5 * current_moving_index + 1] != prev_moving_index && moved_check == 0 && i != chain_length - 1) // checks if movement is possible
			{
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " 2lpre+1 i:%i  pmindex: %4i  cmindex: %i  lp0:%i   lp1:%i    lp2:%i",i, prev_moving_index,current_moving_index,latticepolymer[5 * current_moving_index + 0],latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2]);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/

				latticepolymer[5 * current_moving_index + 2] = prev_moving_index; // stated before previous index is updated to the curent index
				latticepolymer[5 * current_moving_index + 3] = -1;
				poly_chain[i] = current_moving_index;

				prev_moving_index = current_moving_index;
				current_moving_index = latticepolymer[5 * current_moving_index + 1];

				latticepolymer[5 * prev_moving_index + 1] = current_moving_index;
				moved_check = 1;

				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " 3lpre+1 i:%i  pmindex: %4i  cmindex: %i  lp0:%i   lp1:%i    lp2:%i", i, prev_moving_index,current_moving_index,latticepolymer[5 * current_moving_index + 0],latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2]);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
			}
			if (latticepolymer[5 * current_moving_index + 2] >= 0 && latticepolymer[5 * current_moving_index + 2] != prev_moving_index && moved_check == 0 && i != chain_length - 1)// checks if movement is possible
			{
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " 2lpre+2 i:%i  pmindex: %4i  cmindex: %i  lp0:%i   lp1:%i    lp2:%i",i, prev_moving_index,current_moving_index,latticepolymer[5 * current_moving_index + 0],latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2]);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/

				int storage;
				storage = latticepolymer[5 * current_moving_index + 2];
				latticepolymer[5 * current_moving_index + 2] = prev_moving_index;
				latticepolymer[5 * current_moving_index + 3] = -1;
				poly_chain[i] = current_moving_index;

				prev_moving_index = current_moving_index;
				current_moving_index = storage;

				latticepolymer[5 * prev_moving_index + 1] = current_moving_index;
				moved_check = 1;

				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " 3lpre+2  i:%i    pmindex: %4i  cmindex: %i  lp0:%i   lp1:%i    lp2:%i",i, prev_moving_index,current_moving_index,latticepolymer[5 * current_moving_index + 0],latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2]);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " 3.5lpre+2  lp0:%i   lp1:%i    lp2:%i",latticepolymer[5 * prev_moving_index + 0],latticepolymer[5 * prev_moving_index + 1],latticepolymer[5 * prev_moving_index + 2]);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
			}
			i++;
		}
		// sets the values for the tail as the system ends before it gets to the tail as it doesnt have a next lattice point

		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " 4lpre  i:%i    pmindex: %4i  cmindex: %i  lp0:%i   lp1:%i    lp2:%i",i, prev_moving_index,current_moving_index,latticepolymer[5 * current_moving_index + 0],latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2]);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/

		latticepolymer[5 * current_moving_index + 2] = prev_moving_index;
		latticepolymer[5 * current_moving_index + 3] = -1;
		latticepolymer[5 * current_moving_index + 1] = -200;
		poly_chain[chain_length - 1] = current_moving_index;

		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		for (int i = 0; i < chain_length; i++)
		{
			fprintf(stdoutlog, "%4i", poly_chain[i]);
		}
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/


	}
	return 0;
}
