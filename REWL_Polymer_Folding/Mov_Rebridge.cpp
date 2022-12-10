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
#include "Reorg.h"
#include "Mov_Rebridge.h"

int circle(int index, int next_index)
{
	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, " circle");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/
	int starting_index = index; //starting index to check the creation of the circle
	int current_moving_index = index; // updating indexing for the moving check
	latticepolymer[5 * current_moving_index + 3] = 1; //The kinks in the polymer could cause the created circle from rebridging to attempt a rebridge with itself this tags all seen points in a potential cirle to prevent this
	int prev_moving_index = -2000; // previous index to ensure the system doesnt count previous point in loop intial value is random, loop or not any formation from rebridging must gave 2 points
	int lp_check = -3000; //Initial lattice polymer directionm check may not need
	int lp_moving_dir = -4000; // Initial lattice polymer directionm check may not need
	int moved_check = 0; // checks if a move around the loop could be performed if not the segment is assumed to not be a loop and the next created segment is considered
	int circle_check = 0; // checks if a loop is formed may be redundant
	int tag = 1;
	// checks to see if movement can be perforemed in either available directions to set a value for prev_moving_index
	/*
	if (latticepolymer[5 * current_moving_index + 1] >= 0)
	{
		prev_moving_index = current_moving_index;
		current_moving_index = latticepolymer[5 * prev_moving_index + 1];
		latticepolymer[5 * current_moving_index + 3] = 1; // tags the indexes of the new loop
		lp_check = 1;
	}
	if (latticepolymer[5 * current_moving_index + 2] >= 0)
	{
		prev_moving_index = current_moving_index;
		current_moving_index = latticepolymer[5 * prev_moving_index + 2];
		latticepolymer[5 * current_moving_index + 3] = 1;
		lp_check = 2;
	}
	 */
	if (latticepolymer[5 * index + 1] >= 0)  // offsets the initial lattice count by one to start the sequence
	{
		prev_moving_index = current_moving_index;
		current_moving_index = latticepolymer[5 * index + 1];

		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " start 1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/

		latticepolymer[5 * prev_moving_index + 3] = tag;
		latticepolymer[5 * current_moving_index + 3] = tag;
		lp_check = 1;
	}
	else
	{
		if (latticepolymer[5 * index + 2] >= 0)  // offsets the initial lattice count by one to start the sequence
		{
			prev_moving_index = current_moving_index;
			current_moving_index = latticepolymer[5 * index + 2];

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " start 2 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			latticepolymer[5 * prev_moving_index + 3] = tag;
			latticepolymer[5 * current_moving_index + 3] = tag;

			lp_check = 1;
		}
	}

	int i = 1;
	moved_check = 1;
	//moves around loop and checks if the system ends up where it started
	while (current_moving_index != starting_index)
	{
		if (latticepolymer[5 * current_moving_index + 1] >= 0 && latticepolymer[5 * current_moving_index + 1] != prev_moving_index) // checks if movement is possible
		{

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 2+1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			prev_moving_index = current_moving_index;
			current_moving_index = latticepolymer[5 * current_moving_index + 1];
			latticepolymer[5 * current_moving_index + 3] = tag;
			moved_check = 1;

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 3+1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			circle_check = 1;
		}
		if (latticepolymer[5 * current_moving_index + 2] >= 0 && latticepolymer[5 * current_moving_index + 2] != prev_moving_index)// checks if movement is possible
		{
			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 3+1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			prev_moving_index = current_moving_index;
			current_moving_index = latticepolymer[5 * current_moving_index + 2];
			latticepolymer[5 * current_moving_index + 3] = tag;
			moved_check = 1;

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 3+2 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			circle_check = 1;
		}

		if (latticepolymer[5 * current_moving_index + 1] == -200 || latticepolymer[5 * current_moving_index + 2] == -100)// checks if movement is possible
		{
			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 4 next checknext: %4i  prevcheck: %i", latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2] );
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			current_moving_index = starting_index;
			circle_check = 0;

		}
	}
	int final_index;

	if (circle_check == 0)
	{
		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " circle222222222");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/
		starting_index = next_index; //starting index to check the creation of the circle
		current_moving_index = next_index; // updating indexing for the moving check
		latticepolymer[5 * current_moving_index + 3] = 1; //The kinks in the polymer could cause the created circle from rebridging to attempt a rebridge with itself this tags all seen points in a potential cirle to prevent this
		prev_moving_index = -2000; // previous index to ensure the system doesnt count previous point in loop intial value is random, loop or not any formation from rebridging must gave 2 points
		lp_check = -3000; //Initial lattice polymer directionm check may not need
		lp_moving_dir = -4000; // Initial lattice polymer directionm check may not need
		moved_check = 0; // checks if a move around the loop could be performed if not the segment is assumed to not be a loop and the next created segment is considered
		 // checks if a loop is formed may be redundant
		tag = 2;

		if (latticepolymer[5 * next_index + 1] >= 0)  // offsets the initial lattice count by one to start the sequence
		{
			prev_moving_index = current_moving_index;
			current_moving_index = latticepolymer[5 * next_index + 1];

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " start next 1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			latticepolymer[5 * prev_moving_index + 3] = tag;
			latticepolymer[5 * current_moving_index + 3] = tag;
			lp_check = 1;
		}
		else
		{
			if (latticepolymer[5 * next_index + 2] >= 0)  // offsets the initial lattice count by one to start the sequence
			{
				prev_moving_index = current_moving_index;
				current_moving_index = latticepolymer[5 * next_index + 2];

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " start next 2 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				latticepolymer[5 * prev_moving_index + 3] = tag;
				latticepolymer[5 * current_moving_index + 3] = tag;

				lp_check = 1;
			}
		}

		int i = 1;
		moved_check = 1;
		//moves around loop and checks if the system ends up where it started
		while (current_moving_index != starting_index)
		{
			if (latticepolymer[5 * current_moving_index + 1] >= 0 && latticepolymer[5 * current_moving_index + 1] != prev_moving_index) // checks if movement is possible
			{

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "next 2+1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				prev_moving_index = current_moving_index;
				current_moving_index = latticepolymer[5 * current_moving_index + 1];
				latticepolymer[5 * current_moving_index + 3] = tag;
				moved_check = 1;

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "next 3+1 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				circle_check = 1;

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "circle_check %i ",circle_check);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
			}
			if (latticepolymer[5 * current_moving_index + 2] >= 0 && latticepolymer[5 * current_moving_index + 2] != prev_moving_index)// checks if movement is possible
			{
				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "next 3+2 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				prev_moving_index = current_moving_index;
				current_moving_index = latticepolymer[5 * current_moving_index + 2];
				latticepolymer[5 * current_moving_index + 3] = tag;
				moved_check = 1;

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "next 3+2 pmindex: %4i  cmindex: %i", prev_moving_index,current_moving_index);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				circle_check = 1;

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "circle_check %i ",circle_check);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
			}

			if (latticepolymer[5 * current_moving_index + 1] == -200 || latticepolymer[5 * current_moving_index + 2] == -100)// checks if movement is possible
			{
				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " 4 next checknext: %4i  prevcheck: %i", latticepolymer[5 * current_moving_index + 1],latticepolymer[5 * current_moving_index + 2] );
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
				current_moving_index = starting_index;
				circle_check = 0;

			}
		}

		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "end circle_check %i ",circle_check);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/
	}
	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "5circle_check %i ",circle_check);
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/
	if (circle_check == 1)
	{


		final_index = current_moving_index;
		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "final_index 1: %i   2: %i",latticepolymer[5 * final_index + 1],latticepolymer[5 * final_index + 2]);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/
	}
	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	for (int i = 0; i < numberspins; i++)
	{
		fprintf(stdoutlog, "%4i", latticepolymer[5*i+3]);
		if ((i + 1) % L1dim == 0)
			fprintf(stdoutlog, "\n");
	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n\n");
	fclose(stdoutlog);*/
	rebridge_index = final_index;
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "happended: %i  tag: %i",9,tag);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", poly_chain[i]);
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+3]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "rebridge_index 1: %i   2: %i",latticepolymer[5 * rebridge_index + 1],latticepolymer[5 * rebridge_index + 2]);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/
	return tag;
}

int h_t_wiggle() // this funtion is for rebridging of bonds around the head or tail (Flipping of the bonds)
{
	int complete = 0; // complete check
	int i = 1;
	int random_nieghbour = rand() % numberneighbors;

	for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}

	int h_t = (rand() % (2)) + 1; // starts the indexing at a random location either the head or the tail

	while (i < 2) // checks
	{
		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "0  h&t:%i",h_t);
		fprintf(stdoutlog, "\n\n");
		fclose(stdoutlog);*/

		if (h_t == 1) // Starts looking from the head for any neighbours
		{
			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			int main_index = poly_chain[0];
			//for (int ii = 0; ii < numberneighbors; ii++)
			//{
			random_nieghbour++;
			random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value
			// checks to see if the index before the head has any neighbouring lattice

			while ((neighbor[numberneighbors * poly_chain[1] + random_nieghbour]) == poly_chain[0])
			{
				random_nieghbour = rand() % numberneighbors;
			}

			if ((latticepoint[neighbor[numberneighbors * poly_chain[1] + random_nieghbour]]) == 0)
			{

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "2h  random neighbor:%i",random_nieghbour);
				fprintf(stdoutlog, "\n\n");
				fclose(stdoutlog);*/

				//latticepolymer[5 * poly_chain[0] + 1] = -1; // must be stated before the poly_chain[0] is redefined
				//latticepolymer[5 * poly_chain[0] + 2] = -1;
				//latticepoint[poly_chain[0]] = 0; //need to reset the lattice point
				poly_chain[0] = neighbor[numberneighbors * poly_chain[1] + random_nieghbour];
				//latticepoint[poly_chain[0]] = 1;
				//latticepoint[poly_chain[0]] = chain_sequence[0];

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				//latticepolymer[5 * poly_chain[0] + 1] = poly_chain[1];
				//latticepolymer[5 * poly_chain[0] + 2] = -100;

				//latticepolymer[5 * poly_chain[1] + 2] = poly_chain[0]; //updates the record of the attached bonds
				complete = 1;

			}
			if (complete == 1)
			{
				int pop = 1; // Sets a population value
				//reorganize lattice polymer and lattice points
				wiggle_ht = 1;
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, 3);
				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, 0, 3);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    fprintf(stdoutlog, "\n");
				//                    fprintf(stdoutlog, "1 happended: %i",6);
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                    }
				//                    fclose(stdoutlog);
				return complete;
			}

			//random_nieghbour++;
		//}
			if (complete == 0)
			{
				return complete;
				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				random_nieghbour = rand() % numberneighbors;
				int real = chain_length;
				for (int ii = 0; ii < numberneighbors; ii++) // does the fail function which is the same as the head function and checks for neighbours and flips the bonds accordingly
				{
					random_nieghbour++;
					random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value

					if ((latticepoint[neighbor[numberneighbors * poly_chain[chain_length - 2] + random_nieghbour]]) == 0)
					{

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3h  random neighbor:%i",random_nieghbour);
						fprintf(stdoutlog, "\n\n");
						fclose(stdoutlog);*/

						//latticepolymer[5 * poly_chain[chain_length - 1] + 1] = -1;
						//latticepolymer[5 * poly_chain[chain_length - 1] + 2] = -1;
						//latticepoint[poly_chain[chain_length - 1]] = 0; //need to reset the lattice point
						poly_chain[chain_length - 1] = neighbor[numberneighbors * poly_chain[chain_length - 2] + random_nieghbour];
						//latticepoint[poly_chain[chain_length - 1]] = 1;
						//latticepoint[poly_chain[chain_length-1]] = chain_sequence[chain_length-1];

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "%4i", poly_chain[i]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						//latticepolymer[5 * poly_chain[chain_length - 1] + 1] = -200;
						//latticepolymer[5 * poly_chain[chain_length - 1] + 2] = poly_chain[chain_length - 2];

						//latticepolymer[5 * poly_chain[chain_length-2] + 1] = poly_chain[chain_length-1]; //updates the record of the attached bonds
						complete = 1;

					}
					if (complete == 1)
					{
						int pop = 1; // Sets a population value
						//reorganize lattice polymer and lattice points
						lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);
						lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, 0, chain_length);

						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        fprintf(stdoutlog, "\n");
						//                        fprintf(stdoutlog, "12 happended: %i",6);
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                        }
						//                        fclose(stdoutlog);
						return complete;
					}

					//random_nieghbour++;
				}
				if (complete == 0)
				{
					int pop = 1; // Sets a population value
					//reorganize lattice polymer and lattice points
					lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
					lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, pop, 0, chain_length);

					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "12f happended: %i",6);
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
					//                    }
					//                    fclose(stdoutlog);
					return complete;
				}
			}
			i = 2;
		}
		if (h_t == 2) // same function but beigns at the tail
		{
			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/
			int main_index = poly_chain[0];
			//for (int ii = 0; ii < numberneighbors; ii++) //
			//{
			random_nieghbour++;
			random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value
			int stored_index; // stores the index of the new head as it gets quickly overwritten

			while ((neighbor[numberneighbors * poly_chain[chain_length - 2] + random_nieghbour]) == poly_chain[chain_length - 1])
			{
				random_nieghbour = rand() % numberneighbors;
			}

			if ((latticepoint[neighbor[numberneighbors * poly_chain[chain_length - 2] + random_nieghbour]]) == 0)
			{

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "2t  random neighbor:%i",random_nieghbour);
				fprintf(stdoutlog, "\n\n");
				fclose(stdoutlog);*/

				//latticepolymer[5 * poly_chain[chain_length - 1] + 1] = -1;
				//latticepolymer[5 * poly_chain[chain_length - 1] + 2] = -1;
				//latticepoint[poly_chain[chain_length - 1]] = 0; //need to reset the lattice point
				poly_chain[chain_length - 1] = neighbor[numberneighbors * poly_chain[chain_length - 2] + random_nieghbour];
				//latticepoint[poly_chain[chain_length - 1]] = 1;
				//latticepoint[poly_chain[chain_length-1]] = chain_sequence[chain_length-1];

				//latticepolymer[5 * poly_chain[chain_length - 1] + 1] = -200;
				//latticepolymer[5 * poly_chain[chain_length - 1] + 2] = poly_chain[chain_length - 2];

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				//latticepolymer[5 * poly_chain[chain_length-2] + 1] = poly_chain[chain_length-1]; //updates the record of the attached bonds

				complete = 1;

			}
			if (complete == 1)
			{
				int pop = 1; // Sets a population value
				//reorganize lattice polymer and lattice points
				wiggle_ht = 2;
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, chain_length - 3, chain_length);
				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, chain_length - 3, chain_length);
				//lp_reorg_rebrid(poly_chain[0], 1);

//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "\n");
//                    fprintf(stdoutlog, "2 happended: %i",6);
//                    fprintf(stdoutlog, "\n");
//                    fclose(stdoutlog);
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    for (int i = 0; i < chain_length; i++)
//                    {
//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
//                    }
//                    fclose(stdoutlog);
				return complete;
			}
			//random_nieghbour++;
		//}
			if (complete == 0)
			{
				return complete;
				random_nieghbour = rand() % numberneighbors;
				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "%4i", poly_chain[i]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
				int real = chain_length;
				for (int ii = 0; ii < numberneighbors; ii++) //
				{
					random_nieghbour++;
					random_nieghbour = (random_nieghbour) % numberneighbors;// reset the direction value for overflow so I can start with any value
					//changes position of the head if possible very minimal change
					if ((latticepoint[neighbor[numberneighbors * poly_chain[1] + random_nieghbour]]) == 0)
					{
						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3t  random neighbor:%i",random_nieghbour);
						fprintf(stdoutlog, "\n\n");
						fclose(stdoutlog);*/

						//latticepolymer[5 * poly_chain[0] + 1] = -1;
						//latticepolymer[5 * poly_chain[0] + 2] = -1;
						//latticepoint[poly_chain[0]] = 0; //need to reset the lattice point
						poly_chain[0] = neighbor[numberneighbors * poly_chain[1] + random_nieghbour];
						//latticepoint[poly_chain[0]] = 1;
						//latticepoint[poly_chain[0]] = chain_sequence[0];

						//latticepolymer[5 * poly_chain[0] + 1] = poly_chain[1];
						//latticepolymer[5 * poly_chain[0] + 2] = -100;

						complete = 1;

						//latticepolymer[5 * poly_chain[1] + 2] = poly_chain[0]; //updates the record of the attached bonds
						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "%4i", poly_chain[i]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/
					}
					if (complete == 1)
					{
						int pop = 1; // Sets a population value
						//reorganize lattice polymer and lattice points
						lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);
						lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, 0, chain_length);

						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        fprintf(stdoutlog, "\n");
						//                        fprintf(stdoutlog, "21 happended: %i",6);
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                        }
						//                        fclose(stdoutlog);
						return complete;
					}

					//random_nieghbour++;
				}
				if (complete == 0)
				{

					int pop = 1; // Sets a population value
					//reorganize lattice polymer and lattice points
					lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
					lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, pop, 0, chain_length);

					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "21f happended: %i",6);
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
					//                    }
					//                    fclose(stdoutlog);
					return complete;
				}
			}
			i = 2;
		}
		i = 2;
	}
	int pop = 1; // Sets a population value
	//reorganize lattice polymer and lattice points
	lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
	lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, pop, 0, chain_length);
	return complete;
	return complete;
}

int rebridging_h_t() // this funtion is for rebridging of bonds around the head or tail (Flipping of the bonds)
{
	int complete = 0; // complete check
	int i = 1;
	int random_nieghbour = (rand()) % numberneighbors;
	int h_t = (rand() % (2)) + 1; // starts the indexing at a random location either the head or the tail
	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "rebridging_h_t start  h_t:%i",h_t);
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/

	for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "starting happended: %i",8);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", poly_chain[i]);
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
	//
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	while (i < 2) // checks
	{
		if (h_t == 1) // Starts looking from the head for any neighbours
		{
			int main_index = poly_chain[0];
			//for (int ii = 0; ii < numberneighbors; ii++)
			//{
			random_nieghbour++;
			random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value
			// checks to see if the head has any neighbouring lattice and ensures that the interaction is not in sequence
			int stored_index; // stores the index of the new head as it gets quickly overwritten

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "1st ht1  polymer_chain[random neigh]:%i  latticepoint:%i", neighbor[numberneighbors * poly_chain[0] + random_nieghbour],latticepoint[neighbor[numberneighbors * poly_chain[0] + random_nieghbour]]);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			while ((neighbor[numberneighbors * poly_chain[0] + random_nieghbour]) == poly_chain[1]) {
				random_nieghbour = (rand()) % numberneighbors;
			}

			if ((latticepoint[neighbor[numberneighbors * poly_chain[0] + random_nieghbour]]) > 0 && (neighbor[numberneighbors * poly_chain[0] + random_nieghbour]) != poly_chain[1])
			{
				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "ht1  nieghboring latticepoint:%i   previous neighbor index/stored index:%i   polymer_chain[0]:%i",neighbor[numberneighbors * poly_chain[0] + random_nieghbour],latticepolymer[5 * neighbor[numberneighbors * poly_chain[0] + random_nieghbour] + 2], poly_chain[0]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				stored_index = latticepolymer[5 * neighbor[numberneighbors * poly_chain[0] + random_nieghbour] + 2];
				latticepolymer[5 * neighbor[numberneighbors * poly_chain[0] + random_nieghbour] + 2] = poly_chain[0];
				latticepolymer[5 * poly_chain[0] + 2] = neighbor[numberneighbors * poly_chain[0] + random_nieghbour]; // updates the previous interaction

				latticepolymer[5 * stored_index + 1] = -100;
				complete = 1;

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
			}
			if (complete == 1)
			{
				re_ht_si= latticepolymer[5 * stored_index + 4];
				re_ht_ht = h_t;
				lp_reorg_rebrid(latticepolymer[5 * stored_index], 1);

				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);

				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    fprintf(stdoutlog, "\n");
				//                    fprintf(stdoutlog, "1 happended: %i  %i",8,stored_index);
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
				//
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
				//
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
				//
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				return complete;
			}

			//random_nieghbour++;
		//}
			if (complete == 0)
			{
				return complete;

				int real = chain_length;
				for (int ii = 0; ii < numberneighbors; ii++) // does the fail function which is the same as the head function and checks for neighbours and flips the bonds accordingly
				{
					random_nieghbour++;
					random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value
					int stored_index; // stores the index of the new head as it gets quickly overwritten

					/* stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "1st ht12  polymer_chain[random neigh]:%i", neighbor[numberneighbors * poly_chain[chain_length-1] + random_nieghbour]);
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);*/


					if ((latticepoint[neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour]]) > 0 && (neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour]) != poly_chain[chain_length - 2])
					{

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "ht12  nieghboring latticepoint:%i   previous neighbor index/stored index:%i   polymer_chain[0]:%i",neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour],latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour] + 2], poly_chain[chain_length - 1]);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						stored_index = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour] + 1];
						latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour] + 1] = poly_chain[chain_length - 1];
						latticepolymer[5 * poly_chain[chain_length - 1] + 1] = neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour];

						latticepolymer[5 * stored_index + 2] = -200; // the previous index on the new tail is changed because its new
						complete = 1;

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/
					}
					if (complete == 1)
					{
						lp_reorg_rebrid(poly_chain[0], 1);

						lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
						lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);

						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        fprintf(stdoutlog, "\n");
						//                        fprintf(stdoutlog, "12 happended: %i    %i",8,stored_index);
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                        }
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
						//
						//                        }
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
						//
						//                        }
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						return complete;
					}

					//random_nieghbour++;
				}
				if (complete == 0)
				{
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "impp1 happended: %i",8);
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
					//                    }
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
					//
					//                    }
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
					//
					//                    }
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					return complete;
				}
			}
			i = 2;
		}
		if (h_t == 2) // same function but beigns at the tail
		{
			int main_index = poly_chain[0];
			//for (int ii = 0; ii < numberneighbors; ii++) //
			//{
			random_nieghbour++;
			random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value
			int stored_index; // stores the index of the new head as it gets quickly overwritten

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "1st ht2  polymer_chain[random neigh]:%i", neighbor[numberneighbors * poly_chain[chain_length-1] + random_nieghbour]);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/


			while ((neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour]) == poly_chain[chain_length - 2]) {
				random_nieghbour = (rand()) % numberneighbors;
			}

			if ((latticepoint[neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour]]) > 0 && (neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour]) != poly_chain[chain_length - 2])
			{

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "ht2  nieghboring latticepoint:%i   previous neighbor index/stored index:%i   polymer_chain[0]:%i",neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour],latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour] + 2], poly_chain[chain_length - 1]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				stored_index = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour] + 1];
				latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour] + 1] = poly_chain[chain_length - 1];
				latticepolymer[5 * poly_chain[chain_length - 1] + 1] = neighbor[numberneighbors * poly_chain[chain_length - 1] + random_nieghbour];

				latticepolymer[5 * stored_index + 2] = -200; // the previous index on the new tail is changed because its new
				complete = 1;

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				for (int i = 0; i < chain_length; i++)
				{
					fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
				}
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
			}
			if (complete == 1)
			{
				re_ht_si = latticepolymer[5 * stored_index + 4];
				re_ht_ht = h_t;
				lp_reorg_rebrid(poly_chain[0], 1);

				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);

				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    fprintf(stdoutlog, "\n");
				//                    fprintf(stdoutlog, "2 happended: %i    %i",8,stored_index);
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
				//
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    for (int i = 0; i < chain_length; i++)
				//                    {
				//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
				//
				//                    }
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				return complete;
			}
			// random_nieghbour++;
		//}
			if (complete == 0)
			{
				return complete;

				int real = chain_length;
				for (int ii = 0; ii < numberneighbors; ii++) //
				{
					random_nieghbour++;
					random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value
					int stored_index; // stores the index of the new head as it gets quickly overwritten

					/* stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "1st ht21  polymer_chain[random neigh]:%i", neighbor[numberneighbors * poly_chain[0] + random_nieghbour]);
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);*/


					if ((latticepoint[neighbor[numberneighbors * poly_chain[0] + random_nieghbour]]) > 0 && (neighbor[numberneighbors * poly_chain[0] + random_nieghbour]) != poly_chain[1])
					{
						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "ht21  nieghboring latticepoint:%i   previous neighbor index/stored index:%i   polymer_chain[0]:%i",neighbor[numberneighbors * poly_chain[0] + random_nieghbour],latticepolymer[5 * neighbor[numberneighbors * poly_chain[0] + random_nieghbour] + 2], poly_chain[0]);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						stored_index = latticepolymer[5 * neighbor[numberneighbors * poly_chain[0] + random_nieghbour] + 2];
						latticepolymer[5 * neighbor[numberneighbors * poly_chain[0] + random_nieghbour] + 2] = poly_chain[0];
						latticepolymer[5 * poly_chain[0] + 2] = neighbor[numberneighbors * poly_chain[0] + random_nieghbour]; // updates the previous interaction

						latticepolymer[5 * stored_index + 1] = -100;
						complete = 1;

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/
					}
					if (complete == 1)
					{
						lp_reorg_rebrid(latticepolymer[5 * stored_index], 1);

						lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
						lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);

						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        fprintf(stdoutlog, "\n");
						//                        fprintf(stdoutlog, "21 happended: %i    %i",8,stored_index);
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                        }
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
						//
						//                        }
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
						//
						//                        }
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						return complete;
					}

					//random_nieghbour++;
				}
				if (complete == 0)
				{
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "impp2 happended: %i",8);
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
					//                    }
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
					//
					//                    }
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    for (int i = 0; i < chain_length; i++)
					//                    {
					//                        fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
					//
					//                    }
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
					return complete;
				}
			}
			i = 2;
		}
	}
	return complete;
}






int rebridging_2(int index, int tagged) // Reconnects the segment which has formed a loop
{
	// The two cases that need to be accounted for and I cant
	//*******************
	//********|------|h**
	//********|******||**
	//********|--x---||**
	//***t_______x____|**
	//*******************

	//******************
	//***|------------|*
	//***|**|------|h*|*
	//***|**|******||*|*
	//***|**|--x---||*|*
	//***|_____x_t**|_|*
	//******************

	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, " rebridge_22222");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/

	int current_index = index; // function utilizes previous function for indexing from circle() so current and previous indexes are stated
	int previous_index = -1908;
	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, " rebridge_22222");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/
	// previous indexes are set these two will actually loop back on each other but its fine as I have no way of knowing a priori
	/*if (latticepolymer[5 * current_index + 1] >= 0)
	{
		previous_index = current_index;
		current_index = latticepolymer[5 * current_index + 1];
	}
	if (latticepolymer[5 * current_index + 2] >= 0)
	{
		previous_index = current_index;
		current_index = latticepolymer[5 * current_index + 2];
	}*/

	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "index 1: %i   2: %i",latticepolymer[5 * index + 1],latticepolymer[5 * index + 2]);
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/

	if (latticepolymer[5 * index + 1] >= 0)  // offsets the initial lattice count by one to start the sequence
	{
		previous_index = current_index;
		current_index = latticepolymer[5 * index + 1];

		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " start 1 pmindex: %4i  cmindex: %i", previous_index,current_index);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/
	}
	else
	{
		if (latticepolymer[5 * index + 2] >= 0)  // offsets the initial lattice count by one to start the sequence
		{
			previous_index = current_index;
			current_index = latticepolymer[5 * index + 2];

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " start 2 pmindex: %4i  cmindex: %i", previous_index,current_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/
		}
	}


	int random_offset = rand() % (chain_length - 1); // sets the starting marker for the rejoining at a random site
	int i = 0;
	while (i < random_offset) // simply to accomodate for offset effectively so the rebridging does not occur directly next to the previous broken bond site
	{
		if (latticepolymer[5 * current_index + 1] >= 0 && latticepolymer[5 * current_index + 1] != previous_index) // checks if movement is possible
		{
			previous_index = current_index;
			current_index = latticepolymer[5 * current_index + 1];

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " random_offset 1 i:%i pmindex: %4i  cmindex: %i",i, previous_index,current_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/
		}
		if (latticepolymer[5 * current_index + 2] >= 0 && latticepolymer[5 * current_index + 2] != previous_index)// checks if movement is possible
		{
			previous_index = current_index;
			current_index = latticepolymer[5 * current_index + 2];

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " random_offset 2 i:%i pmindex: %4i  cmindex: %i",i, previous_index,current_index);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/
		}
		i++;
	}

	int reattached = 0; // reattached check

	while (reattached == 0) // checks
	{
		//checks to ensure indexes are nnot the tail and head

		int sorround_check_prev = 100;
		int sorround_check = 100; // unoccupiable values that store the neighbor index assignment

		int random_nieghbour = (rand() % numberneighbors) + 1;  // Chooses a random starting direction

		int lp_index = 100;
		int lp_index_2 = 100;



		//Turn into While loop
		//for (int ii = 0; ii < numberneighbors; ii++) {
			//random_nieghbour++;
		random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value

		if (neighbor[numberneighbors * previous_index + random_nieghbour] != current_index) //checks the neighbour for the first point to ensure its not the other index
		{

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "2  previous_index: %i   current_index: %i  neigh: %i", previous_index, current_index, neighbor[numberneighbors * previous_index + random_nieghbour]);
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/


			if (latticepoint[neighbor[numberneighbors * previous_index + random_nieghbour]] > 0 && latticepoint[neighbor[numberneighbors * current_index + random_nieghbour]] > 0) //ensures both indexes sorround adjacent populated sites
			{

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "3  neighprev: %4i  neighcurr: %i", neighbor[numberneighbors * previous_index + random_nieghbour], neighbor[numberneighbors * current_index + random_nieghbour]);
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				//checks for if difference in the tag between the current points of intrest and the neighbouring adjacent points for attempted rebridging to ensure the are different so that the loop does not rebridge with itself
				if (latticepolymer[5 * neighbor[numberneighbors * previous_index + random_nieghbour] + 3] != tagged && latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + 3] != tagged)
				{

					/* stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "4  neighprev: %4i  neighcurr: %i", neighbor[numberneighbors * previous_index + random_nieghbour], neighbor[numberneighbors * current_index + random_nieghbour]);
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);*/

					if (latticepolymer[5 * neighbor[numberneighbors * previous_index + random_nieghbour] + 1] == neighbor[numberneighbors * current_index + random_nieghbour]) // performs a check to see if the point is the connected to the other stores the lattice polymer index if so and updates the other lattice polymer site
					{
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "happended: %i  rb2 1",9);
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+3]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
													/* stdoutlog = fopen(stdoutlogname, "a");
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "5a");
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "\n");
													fclose(stdoutlog);*/
													/* stdoutlog = fopen(stdoutlogname, "a");
													fprintf(stdoutlog, "\n");
													for (int i = 0; i < chain_length; i++)
													{
														fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
													}
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "\n");
													fclose(stdoutlog);*/

						lp_index = 1; // stores the index of the lattice polymer indexes
						lp_index_2 = 2;

						latticepolymer[5 * (neighbor[numberneighbors * previous_index + random_nieghbour]) + lp_index] = previous_index; // changes the lattice polymer index for all of the indexes involved
						latticepolymer[5 * previous_index + 1] = latticepolymer[5 * neighbor[numberneighbors * previous_index + random_nieghbour] + 0];

						latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + lp_index_2] = current_index;
						latticepolymer[5 * current_index + 2] = latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + 0];

						reattached = 1;

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "happended: %i  rb2 12",9);
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+3]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
					}
					if (latticepolymer[5 * neighbor[numberneighbors * latticepolymer[5 * previous_index] + random_nieghbour] + 2] == neighbor[numberneighbors * latticepolymer[5 * current_index] + random_nieghbour]) // performs a check to see if the point is the connected to the other stores the lattice polymer index if so and updates the other lattice polymer site
					{

						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "happended: %i  rb2 2",9);
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+3]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);

													/* stdoutlog = fopen(stdoutlogname, "a");
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "5b");
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "\n");
													fclose(stdoutlog);*/
													/* stdoutlog = fopen(stdoutlogname, "a");
													fprintf(stdoutlog, "\n");
													for (int i = 0; i < chain_length; i++)
													{
														fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
													}
													fprintf(stdoutlog, "\n");
													fprintf(stdoutlog, "\n");
													fclose(stdoutlog);*/

													//stores the lattice polymer indexes
						lp_index = 2;
						lp_index_2 = 1;

						//update the lattice polymer connections
						latticepolymer[5 * neighbor[numberneighbors * previous_index + random_nieghbour] + lp_index] = latticepolymer[5 * previous_index];
						latticepolymer[5 * previous_index + 1] = latticepolymer[5 * neighbor[numberneighbors * previous_index + random_nieghbour] + 0];

						latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + lp_index_2] = current_index;
						latticepolymer[5 * current_index + 2] = latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + 0];

						reattached = 1;

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < chain_length; i++)
						{
							fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "happended: %i  rb2 21",9);
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            fprintf(stdoutlog, "\n");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
						//                            stdoutlog = fopen(stdoutlogname, "a");
						//                            for (int i = 0; i < chain_length; i++)
						//                            {
						//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+3]);
						//
						//                            }
						//                            fprintf(stdoutlog, "\n");
						//                            fprintf(stdoutlog, "\n");
						//                            fclose(stdoutlog);
					}
					if (reattached == 1)
					{
						return 7;
					}
				}
			}
		}
		// itterates the index if none is chosen
   //}
   // moves the index as it cannot move along a set line
		int moved_check = 0;
		if (latticepolymer[5 * current_index + 1] >= 0 && latticepolymer[5 * current_index + 1] != previous_index && moved_check == 0) // checks if movement is possible
		{
			previous_index = current_index;
			current_index = latticepolymer[5 * current_index + 1];
			moved_check = 1;
		}
		if (latticepolymer[5 * current_index + 2] >= 0 && latticepolymer[5 * current_index + 2] != previous_index && moved_check == 0)// checks if movement is possible
		{
			previous_index = current_index;
			current_index = latticepolymer[5 * current_index + 2];
			moved_check = 1;
		}

		if (reattached == 0)
		{
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, 0, chain_length);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
			return 0;
		}


	}
	return 0;
}

int rebridging()
{
	int complete = 0; // complete check
	int plane_check; // checks to see if adjacent variables are on the same plane
	int current_index;
	int holder = -1;
	int i = 1;

	int random_brid = rand() % (chain_length); // starts the indexing at a random location along the chain
	//must be two indexes away from the end as the end is coiled to make the chain pull
	int prev_index = poly_chain[random_brid]; // stores the value one of random pull
	int prev_brid = random_brid;
	random_brid = (random_brid + 1) % chain_length; // stores the value of 1 + random pull


	/* stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, " rebridge main");
	fprintf(stdoutlog, "\n");
	for (int i = 0; i < chain_length; i++)
	{
		fprintf(stdoutlog, "%4i", poly_chain[i]);
	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);
	 */

	for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}

	while (i < chain_length) // checks
	{
		random_brid = (random_brid) % chain_length;
		prev_index = poly_chain[random_brid];//resets value of
		prev_brid = random_brid;
		random_brid++;
		random_brid = (random_brid) % chain_length;
		current_index = poly_chain[random_brid];


		while (random_brid == 0)
		{
			random_brid = rand() % (chain_length);
			random_brid = (random_brid) % chain_length;
			prev_index = poly_chain[random_brid];//resets value of
			prev_brid = random_brid;
			random_brid++;
			random_brid = (random_brid) % chain_length;
			current_index = poly_chain[random_brid];
		}

		// stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "0  pindex:%i   cindex:%i    random_pull:%i",prev_index,current_index,random_pull);
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);*/

		//checks to ensure indexes are not the tail and head
		if (random_brid != 0)
		{
			///* stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\n");
			//fprintf(stdoutlog, "1  pindex:%i   cindex:%i",prev_index,current_index);
			//fprintf(stdoutlog, "\n");
			//fclose(stdoutlog);*/*

			int sorround_check_prev = 100;
			int sorround_check = 100; // unoccupiable values that store the neighbor index assignment

			int random_nieghbour = (rand() % numberneighbors);  // Chooses a random starting direction

			// set values for the lattice polymer index that need to be changed for the the found rebridging segment
			int lp_index = 100;
			int lp_index_2 = 100;

			// Given the new knowledge of formable structures I need a check so only one configuration requires the tagging from circle check
			int circle_check = -1;

			//Turn into While loop
			for (int ii = 0; ii < numberneighbors; ii++)
			{
				random_nieghbour++;
				random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value

				while (neighbor[numberneighbors * prev_index + random_nieghbour] == current_index)
				{
					random_nieghbour = (rand()) % numberneighbors;
				}

				if (neighbor[numberneighbors * prev_index + random_nieghbour] != current_index) //checks the neighbour for the first point to ensure its not the other index
				{
					/* stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "2  prev_brid: %4i  random_brid: %i   previous_index: %i   current_index: %i  neigh: %i", prev_brid , random_brid, prev_index, current_index, neighbor[numberneighbors * prev_index + random_nieghbour]);
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);*/

					if (latticepoint[neighbor[numberneighbors * prev_index + random_nieghbour]] > 0 && latticepoint[neighbor[numberneighbors * current_index + random_nieghbour]] > 0) //ensures both indexes sorround adjacent populated sites
					{
						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3  neighprev: %4i  neighcurr: %i", neighbor[numberneighbors * prev_index + random_nieghbour], neighbor[numberneighbors * current_index + random_nieghbour]);
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						if (latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour] + 1] == neighbor[numberneighbors * current_index + random_nieghbour] &&
							(latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour]] != prev_index) &&
							(latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour]] != current_index) &&
							poly_chain[random_brid + 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[random_brid + 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[random_brid - 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[random_brid - 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[random_brid] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[random_brid] != neighbor[numberneighbors * prev_index + random_nieghbour] &&

							poly_chain[prev_brid + 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[prev_brid + 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[prev_brid - 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[prev_brid - 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[prev_brid] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[prev_brid] != neighbor[numberneighbors * prev_index + random_nieghbour]) // performs a check to see if the point is the connected to the other stores the lattice polymer index if so and updates the other lattice polymer site
						{

							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "1 happended: %i",9);
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);

														/* stdoutlog = fopen(stdoutlogname, "a");
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "4a");
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "\n");
														fclose(stdoutlog);*/
														/* stdoutlog = fopen(stdoutlogname, "a");
														fprintf(stdoutlog, "\n");
														for (int i = 0; i < chain_length; i++)
														{
															fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
														}
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "\n");
														fclose(stdoutlog);*/

							lp_index = 1; // stores the index of the lattice polymer indexes
							lp_index_2 = 2;

							// changes the lattice polymer index for all of the indexes involved
							// the lp values are in place to ensure that the correct lattice polymer value for the connecting index
							// In this case the index which matches with pindex is the fursthest monomer in the chain so lpin
							//  **|-----|**
							//  **|-x-**|**
							//  *---x---|**
							//Change from before prev_brid will always references the same location
							//This configuration will never form a circle
							latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour] + lp_index] = prev_index;
							latticepolymer[5 * prev_index + 1] = latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour]];

							latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + lp_index_2] = current_index;
							latticepolymer[5 * current_index + 2] = latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour]];

							circle_check = 0;

							/* stdoutlog = fopen(stdoutlogname, "a");
							fprintf(stdoutlog, "\n");
							for (int i = 0; i < chain_length; i++)
							{
								fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
							}
							fprintf(stdoutlog, "\n");
							fprintf(stdoutlog, "\n");
							fclose(stdoutlog);*/

							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "12 happended: %i",9);
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);

						}
						if (latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour] + 2] == neighbor[numberneighbors * current_index + random_nieghbour] &&
							(latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour]] != prev_index) &&
							(latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour]] != current_index) &&
							poly_chain[random_brid + 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[random_brid + 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[random_brid - 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[random_brid - 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[random_brid] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[random_brid] != neighbor[numberneighbors * prev_index + random_nieghbour] &&

							poly_chain[prev_brid + 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[prev_brid + 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[prev_brid - 1] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[prev_brid - 1] != neighbor[numberneighbors * prev_index + random_nieghbour] &&
							poly_chain[prev_brid] != neighbor[numberneighbors * current_index + random_nieghbour] &&
							poly_chain[prev_brid] != neighbor[numberneighbors * prev_index + random_nieghbour])


							// performs a check to see if the point is the connected to the other stores the lattice polymer index if so and updates the other lattice polymer site
						{

							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "2 happended: %i",9);
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//
														/* stdoutlog = fopen(stdoutlogname, "a");
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "4b");
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "\n");
														fclose(stdoutlog);*/
														/* stdoutlog = fopen(stdoutlogname, "a");
														fprintf(stdoutlog, "\n");
														for (int i = 0; i < chain_length; i++)
														{
															fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
														}
														fprintf(stdoutlog, "\n");
														fprintf(stdoutlog, "\n");
														fclose(stdoutlog);*/

														//stores the lattice polymer indexes
							lp_index = 2;
							lp_index_2 = 1;

							//update the lattice polymer connections
							//  **|--|******
							//  **|**|-x-***
							//  **|----x--**
							// This configuration will always form a circle
							// Additionally now the loop coordinates are exact and dont need to changed
							latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour] + lp_index] = prev_index;
							latticepolymer[5 * prev_index + 1] = latticepolymer[5 * neighbor[numberneighbors * prev_index + random_nieghbour]];

							latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + lp_index_2] = current_index;
							latticepolymer[5 * current_index + 2] = latticepolymer[5 * neighbor[numberneighbors * current_index + random_nieghbour] + 0];

							circle_check = 1;

							/* stdoutlog = fopen(stdoutlogname, "a");
							fprintf(stdoutlog, "\n");
							for (int i = 0; i < chain_length; i++)
							{
								fprintf(stdoutlog, "lp0:%4i   lp1:%i     lp2:%i\n", latticepolymer[5*poly_chain[i]],latticepolymer[5*poly_chain[i]+1],latticepolymer[5*poly_chain[i]+2]);
							}
							fprintf(stdoutlog, "\n");
							fprintf(stdoutlog, "\n");
							fclose(stdoutlog);*/

							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "22 happended: %i",9);
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            fprintf(stdoutlog, "\n");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i", poly_chain[i]);
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);
							//
							//                            stdoutlog = fopen(stdoutlogname, "a");
							//                            for (int i = 0; i < chain_length; i++)
							//                            {
							//                                fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
							//
							//                            }
							//                            fprintf(stdoutlog, "\n");
							//                            fprintf(stdoutlog, "\n");
							//                            fclose(stdoutlog);

						}
						if (circle_check >= 0)
						{
							if (circle_check == 0) // If the final matrix is still continuous
							{
								//polymer chain reorganizing function goes here

								lp_reorg_rebrid(poly_chain[0], 1);
								lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length); // resets the poly chain lattice points
								lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length); // reset the plymer chain
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "\n");
//                                fprintf(stdoutlog, "happended: %i  reroute",9);
//                                fprintf(stdoutlog, "\n");
//                                fclose(stdoutlog);
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "\n");
//                                for (int i = 0; i < chain_length; i++)
//                                {
//                                    fprintf(stdoutlog, "%4i", poly_chain[i]);
//                                }
//                                fprintf(stdoutlog, "\n");
//                                fclose(stdoutlog);
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                for (int i = 0; i < chain_length; i++)
//                                {
//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
//
//                                }
//                                fprintf(stdoutlog, "\n");
//                                fclose(stdoutlog);
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                for (int i = 0; i < chain_length; i++)
//                                {
//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
//
//                                }
//                                fprintf(stdoutlog, "\n");
//                                fclose(stdoutlog);
//
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                for (int i = 0; i < chain_length; i++)
//                                {
//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
//
//                                }
//                                fprintf(stdoutlog, "\n");
//                                fprintf(stdoutlog, "\n");
//                                fclose(stdoutlog);

								complete = 1;
								return complete;
							}
							if (circle_check == 1) // if the segment which contains p_index is the circle
							{
								int completion_check;

								int tag;
								tag = circle_check = circle(current_index, prev_index); // tags the new loop
								//rebridging_2(prev_index, circle_check); // attaches the loop
								completion_check = rebridging_2(rebridge_index, tag);
								if (completion_check == 0)
								{
									lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
									lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, 0, chain_length);
									return 0;
								}
								//rebridging_2(prev_index);
								lp_reorg_rebrid(poly_chain[0], 1);
								lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length); // resets the poly chain lattice points
								lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length); // reset the plymer chain
								complete = 1;
								//                                stdoutlog = fopen(stdoutlogname, "a");
								//                                fprintf(stdoutlog, "\n");
								//                                fprintf(stdoutlog, "happended: %i  rb2 f",9);
								//                                fprintf(stdoutlog, "\n");
								//                                fclose(stdoutlog);
								//                                stdoutlog = fopen(stdoutlogname, "a");
								//                                fprintf(stdoutlog, "\n");
								//                                for (int i = 0; i < chain_length; i++)
								//                                {
								//                                    fprintf(stdoutlog, "%4i", poly_chain[i]);
								//                                }
								//                                fprintf(stdoutlog, "\n");
								//                                fclose(stdoutlog);
								//                                stdoutlog = fopen(stdoutlogname, "a");
								//                                for (int i = 0; i < chain_length; i++)
								//                                {
								//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]]);
								//
								//                                }
								//                                fprintf(stdoutlog, "\n");
								//                                fclose(stdoutlog);
								//                                stdoutlog = fopen(stdoutlogname, "a");
								//                                for (int i = 0; i < chain_length; i++)
								//                                {
								//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+1]);
								//
								//                                }
								//                                fprintf(stdoutlog, "\n");
								//                                fclose(stdoutlog);
								//                                stdoutlog = fopen(stdoutlogname, "a");
								//                                for (int i = 0; i < chain_length; i++)
								//                                {
								//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+2]);
								//
								//                                }
								//                                fprintf(stdoutlog, "\n");
								//                                fprintf(stdoutlog, "\n");
								//                                fclose(stdoutlog);
								//                                stdoutlog = fopen(stdoutlogname, "a");
								//                                for (int i = 0; i < chain_length; i++)
								//                                {
								//                                    fprintf(stdoutlog, "%4i",latticepolymer[5*poly_chain[i]+3]);
								//
								//                                }
								//                                fprintf(stdoutlog, "\n");
								//                                fprintf(stdoutlog, "\n");
								//                                fclose(stdoutlog);
								return complete;
							}
						}
						sorround_check = random_nieghbour;
					}
				}
				//random_nieghbour++; // itterates the index if none is chosen
				return complete;
			}

			// stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\n");
			//fprintf(stdoutlog, "3  sorround_check:%i  check1:%i    check2:%i",sorround_check, neighbor[4*p_index+sorround_check], neighbor[4*c_index+sorround_check]);
			//fprintf(stdoutlog, "\n");
			//fclose(stdoutlog);
		}
		return complete;
		i++;
	}
	return complete;
}
