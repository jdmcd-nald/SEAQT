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
#include "Mov_Other.h"

int cornerflip(int f) //First of the movement options for the polymer when atleast one change in direction occurs in a kink the corner is allowed to move to a connected adjecent empty lattice point
{
	int o_iii = 0;
	int flip_index = -1;
	int empty_index = -1;
	int chain_location = -1;
	int random_cf = rand() % (chain_length - 2);//Picks a random index to start the search for pull
	random_cf = (random_cf) % (chain_length - 2);

	for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}

	for (int ii = 0; ii < chain_length - 2; ii++)
	{

		random_cf = (random_cf) % (chain_length - 2);
		cf_value = random_cf;
		if (f >= 0)
		{
			random_cf = f;
		}

		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "corner flip 1  chain_length:%i  preceeding:%i  succeeding:%i    \n",random_cf, poly_chain[random_cf],poly_chain[random_cf+2]);
		fprintf(stdoutlog, "1 left if 1:%i  right if 1:%i left if 2:%i   right if two:%i ",poly_chain[random_cf]%L1dim, poly_chain[random_cf+2]%L1dim, ((poly_chain[random_cf]-(poly_chain[random_cf]%L1dim))),((poly_chain[random_cf+2]-(poly_chain[random_cf+2]%L1dim))))  ;
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/

		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		for (int i = 0; i < numberspins; i++)
		{
			fprintf(stdoutlog, "%4i", latticepoint[i]);
			if ((i + 1) % L1dim == 0)
				fprintf(stdoutlog, "\n");
		}
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/

		/*
		int cr_index=random_cf+1;
		int p_diff = poly_chain[random_cf] - poly_chain[cr_index];
		int c_diff = poly_chain[cr_index] - poly_chain[random_cf];

		if(abs(p_diff) == 900||abs(c_diff) == 900)
		{
			if(abs(p_diff) == 900) p_diff=100;
			if(abs(c_diff) == 900) c_diff=100;
		}
		if(abs(p_diff) == 90||abs(c_diff) == 90)
		{
			if(abs(p_diff) == 90) p_diff=10;
			if(abs(c_diff) == 90) c_diff=10;
		}
		if(abs(p_diff) == 9||abs(c_diff) == 9)
		{
			if(abs(p_diff) == 9) p_diff=1;
			if(abs(c_diff) == 9) c_diff=1;
		}
		*/



		int equal_check = 0;

		for (int i = 0; i < numberneighbors; i++)
		{
			for (int iii = 0; iii < numberneighbors; iii++)
			{
				if (neighbor[numberneighbors * poly_chain[random_cf] + iii] == neighbor[numberneighbors * poly_chain[random_cf + 2] + i])
				{
					equal_check++;
				}
			}
		}

		o_iii = 0;
		// If statemet simplified from before now only checks if there is a three chain link in up/down dir in 3d ( also checks top and bottom connections on the planeand checks if there is are coinciding rows
		if (equal_check == 2)//(poly_chain[random_cf] % L1dim != poly_chain[random_cf + 2] % L1dim) && (poly_chain[random_cf] - poly_chain[random_cf] % L1dim) != (poly_chain[random_cf + 2] - poly_chain[random_cf + 2] % L1dim))
		{
			if (neighbor[(numberneighbors * poly_chain[random_cf]) + rand() % numberneighbors] != poly_chain[random_cf + 1])
			{
				return 0;
			}


			//cycles through neighbors list for sorrounding sites
			for (int iii = 0; iii < numberneighbors * numberneighbors; iii++)
			{
				//neighbor list for polymer adjacent to the corner
				if (iii != 0 && iii % numberneighbors == 0)
				{
					o_iii++;
				}

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "1  random_cf: %i  random_cf+1: %i  randsom_cf+2: %i  random_cf index: %i\n  leftif: %i  rightif: %i\n",poly_chain[random_cf],poly_chain[random_cf+1],poly_chain[random_cf+2],random_cf,neighbor[numberneighbors*poly_chain[random_cf]+(iii%numberneighbors)],neighbor[numberneighbors*poly_chain[random_cf+2]+o_iii]);

				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				//checks for coinciding latticepoint between the points on each side of the flip
				if (neighbor[numberneighbors * poly_chain[random_cf] + (iii % numberneighbors)] == neighbor[numberneighbors * poly_chain[random_cf + 2] + o_iii])
				{

					/* stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "2  iii:%i  ipercent:%i  o_iii:%i    \n",iii, iii%numberneighbors,o_iii);
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);*/

					//stores the index and the chain monomener *(last one unneeded)
					if (poly_chain[random_cf + 1] == neighbor[numberneighbors * poly_chain[random_cf + 2] + o_iii])
					{
						flip_index = poly_chain[random_cf + 1];
						chain_location = (random_cf + 1);

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3flip  flipindex:%i  chainlocation:%i   \n",flip_index, chain_location);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

					}

					if (latticepoint[neighbor[numberneighbors * poly_chain[random_cf] + (iii % numberneighbors)]] == 0)
					{
						empty_index = neighbor[numberneighbors * poly_chain[random_cf] + (iii % numberneighbors)];

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3empty  emptyindex:%i  check:%i  check2: %i   check3:%i\n",empty_index,neighbor[numberneighbors*poly_chain[random_cf]+(iii%numberneighbors)], neighbor[numberneighbors*poly_chain[random_cf+2]+o_iii],latticepoint[neighbor[numberneighbors*poly_chain[random_cf]+(iii%numberneighbors)]]);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/
					}

					if (empty_index >= 0 && flip_index >= 0 & chain_location >= 0) // if corner flipped performed
					{


						//corner flip latticepoints updated
						//latticepoint[empty_index] = latticepoint[flip_index];
						//latticepoint[flip_index] = 0;

						//corner flip latticepoints updated
						//latticepolymer[5 * empty_index + 1] = latticepolymer[5 * flip_index + 1];
						//latticepolymer[5 * empty_index + 2] = latticepolymer[5 * flip_index + 2];
						//latticepolymer[5 * flip_index + 1] = -1;
						//latticepolymer[5 * flip_index + 2] = -1;

						poly_chain[chain_location] = empty_index;

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "4finished empty_index: %i emptyindexlatticepoint: %i  latticepolymer empty: %i  latticepolymerflip: %i\n",empty_index,latticepoint[empty_index] ,latticepolymer[5 * empty_index + 1] , latticepolymer[5 * flip_index + 1]);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < numberspins; i++)
						{
							fprintf(stdoutlog, "%4i", latticepoint[i]);
							if ((i + 1) % L1dim == 0)
								fprintf(stdoutlog, "\n");
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n\n");
						fclose(stdoutlog);*/

						int pop = 1; // Sets a population value
						//reorganize lattice polymer and lattice points
						lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, random_cf, random_cf + 3);
						lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, random_cf, random_cf + 3);

						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        fprintf(stdoutlog, "\n");
						//                        fprintf(stdoutlog, "happended: %i",4);
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                        }
						//                        fclose(stdoutlog);


						return 1;
					}
				}
			}
		}

		if (empty_index >= 0 && flip_index >= 0 & chain_location >= 0) // if corner flipped performed
		{


			//corner flip latticepoints updated
			//latticepoint[empty_index] = latticepoint[flip_index];
			//latticepoint[flip_index] = 0;

			//corner flip latticepoints updated
			//latticepolymer[5 * empty_index + 1] = latticepolymer[5 * flip_index + 1];
			//latticepolymer[5 * empty_index + 2] = latticepolymer[5 * flip_index + 2];
			//latticepolymer[5 * flip_index + 1] = -1;
			//latticepolymer[5 * flip_index + 2] = -1;

			//updating adjacent monomers
			//latticepolymer[5 * poly_chain[chain_location+1] + 2] = empty_index;
			//latticepolymer[5 * poly_chain[chain_location-1] + 1] = empty_index;

			poly_chain[chain_location] = empty_index;

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "4finished empty_index: %i emptyindexlatticepoint: %i  latticepolymer empty: %i  latticepolymerflip: %i\n",empty_index,latticepoint[empty_index] ,latticepolymer[5 * empty_index + 1] , latticepolymer[5 * flip_index + 1]);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/

			lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, random_cf, random_cf + 3);
			lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, random_cf, random_cf + 3);

			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "happended: %ib",4);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", poly_chain[i]);
			//            }
			//            fclose(stdoutlog);


			return 1;
		}
		return 0;
		if (f >= 0)
		{
			return 0;
		}
		flip_index = -1;
		chain_location = -1;
		empty_index = -1;
		random_cf++;
		random_cf = (random_cf) % (chain_length - 2);
	}

	int pop = 1; // Sets a population value
	//reorganize lattice polymer and lattice points
	lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
	lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, pop, 0, chain_length);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "happended: %i  f",4);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", poly_chain[i]);
	//    }
	//    fclose(stdoutlog);


	return 0;
}

int doublecornerflip(int f) //First of the movement options for the polymer when atleast one change in direction occurs in a kink the corner is allowed to move to a connected adjecent empty lattice point
{
	int o_iii = 0;
	int flip_index = -1;
	int empty_index = -1;
	int chain_location = -1;
	int random_cf = rand() % (chain_length - 2);//Picks a random index to start the search for pull
	random_cf = (random_cf) % (chain_length - 2);

	for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
	{
		doublepseudo_poly_chain[y] = pseudo_poly_chain[y];
	}

	for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}

	int lorr = rand() % 2 + 1;

	int check = 0;

	for (int ii = 0; ii < chain_length - 2; ii++)
	{
		if (cf_value == 0)
		{
			random_cf = cf_value + 1;
			if (check == 1)
			{
				return 0;
			}
		}

		if (cf_value == chain_length - 3)
		{
			random_cf = cf_value - 1;
			if (check == 1)
			{
				return 0;
			}
		}

		if (cf_value > 0 && cf_value < chain_length - 3)
		{
			if (lorr == 1 && check == 0)  random_cf = cf_value - 1;

			if (lorr == 2 && check == 1)  random_cf = cf_value - 1;

			if (lorr == 2 && check == 0)  random_cf = cf_value + 1;

			if (lorr == 1 && check == 1)  random_cf = cf_value + 1;

			if (check == 2)
			{
				return 0;
			}
		}


		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "corner flip 1  chain_length:%i  preceeding:%i  succeeding:%i    \n",random_cf, poly_chain[random_cf],poly_chain[random_cf+2]);
		fprintf(stdoutlog, "1 left if 1:%i  right if 1:%i left if 2:%i   right if two:%i ",poly_chain[random_cf]%L1dim, poly_chain[random_cf+2]%L1dim, ((poly_chain[random_cf]-(poly_chain[random_cf]%L1dim))),((poly_chain[random_cf+2]-(poly_chain[random_cf+2]%L1dim))))  ;
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/

		/* stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		for (int i = 0; i < numberspins; i++)
		{
			fprintf(stdoutlog, "%4i", latticepoint[i]);
			if ((i + 1) % L1dim == 0)
				fprintf(stdoutlog, "\n");
		}
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/

		/*
		int cr_index=random_cf+1;
		int p_diff = poly_chain[random_cf] - poly_chain[cr_index];
		int c_diff = poly_chain[cr_index] - poly_chain[random_cf];

		if(abs(p_diff) == 900||abs(c_diff) == 900)
		{
			if(abs(p_diff) == 900) p_diff=100;
			if(abs(c_diff) == 900) c_diff=100;
		}
		if(abs(p_diff) == 90||abs(c_diff) == 90)
		{
			if(abs(p_diff) == 90) p_diff=10;
			if(abs(c_diff) == 90) c_diff=10;
		}
		if(abs(p_diff) == 9||abs(c_diff) == 9)
		{
			if(abs(p_diff) == 9) p_diff=1;
			if(abs(c_diff) == 9) c_diff=1;
		}
		*/
		int equal_check = 0;

		for (int i = 0; i < numberneighbors; i++)
		{
			for (int iii = 0; iii < numberneighbors; iii++)
			{
				if (neighbor[numberneighbors * poly_chain[random_cf] + iii] == neighbor[numberneighbors * poly_chain[random_cf + 2] + i])
				{
					equal_check++;
				}
			}
		}

		o_iii = 0;
		// If statemet simplified from before now only checks if there is a three chain link in up/down dir in 3d ( also checks top and bottom connections on the planeand checks if there is are coinciding rows
		if (equal_check == 2)//(poly_chain[random_cf] % L1dim != poly_chain[random_cf + 2] % L1dim) && (poly_chain[random_cf] - poly_chain[random_cf] % L1dim) != (poly_chain[random_cf + 2] - poly_chain[random_cf + 2] % L1dim))
		{
			//cycles through neighbors list for sorrounding sites
			for (int iii = 0; iii < numberneighbors * numberneighbors; iii++)
			{
				//neighbor list for polymer adjacent to the corner
				if (iii != 0 && iii % numberneighbors == 0)
				{
					o_iii++;
				}

				/* stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "1  random_cf: %i  random_cf+1: %i  randsom_cf+2: %i  random_cf index: %i\n  leftif: %i  rightif: %i\n",poly_chain[random_cf],poly_chain[random_cf+1],poly_chain[random_cf+2],random_cf,neighbor[numberneighbors*poly_chain[random_cf]+(iii%numberneighbors)],neighbor[numberneighbors*poly_chain[random_cf+2]+o_iii]);

				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				//checks for coinciding latticepoint between the points on each side of the flip
				if (neighbor[numberneighbors * poly_chain[random_cf] + (iii % numberneighbors)] == neighbor[numberneighbors * poly_chain[random_cf + 2] + o_iii])
				{

					/* stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "\n");
					fprintf(stdoutlog, "2  iii:%i  ipercent:%i  o_iii:%i    \n",iii, iii%numberneighbors,o_iii);
					fprintf(stdoutlog, "\n");
					fclose(stdoutlog);*/

					//stores the index and the chain monomener *(last one unneeded)
					if (poly_chain[random_cf + 1] == neighbor[numberneighbors * poly_chain[random_cf + 2] + o_iii])
					{
						flip_index = poly_chain[random_cf + 1];
						chain_location = (random_cf + 1);

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3flip  flipindex:%i  chainlocation:%i   \n",flip_index, chain_location);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

					}

					if (latticepoint[neighbor[numberneighbors * poly_chain[random_cf] + (iii % numberneighbors)]] == 0)
					{
						empty_index = neighbor[numberneighbors * poly_chain[random_cf] + (iii % numberneighbors)];

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "3empty  emptyindex:%i  check:%i  check2: %i   check3:%i\n",empty_index,neighbor[numberneighbors*poly_chain[random_cf]+(iii%numberneighbors)], neighbor[numberneighbors*poly_chain[random_cf+2]+o_iii],latticepoint[neighbor[numberneighbors*poly_chain[random_cf]+(iii%numberneighbors)]]);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/
					}

					if (empty_index >= 0 && flip_index >= 0 & chain_location >= 0) // if corner flipped performed
					{


						//corner flip latticepoints updated
						//latticepoint[empty_index] = latticepoint[flip_index];
						//latticepoint[flip_index] = 0;

						//corner flip latticepoints updated
						//latticepolymer[5 * empty_index + 1] = latticepolymer[5 * flip_index + 1];
						//latticepolymer[5 * empty_index + 2] = latticepolymer[5 * flip_index + 2];
						//latticepolymer[5 * flip_index + 1] = -1;
						//latticepolymer[5 * flip_index + 2] = -1;

						poly_chain[chain_location] = empty_index;

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "4finished empty_index: %i emptyindexlatticepoint: %i  latticepolymer empty: %i  latticepolymerflip: %i\n",empty_index,latticepoint[empty_index] ,latticepolymer[5 * empty_index + 1] , latticepolymer[5 * flip_index + 1]);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

						/* stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						for (int i = 0; i < numberspins; i++)
						{
							fprintf(stdoutlog, "%4i", latticepoint[i]);
							if ((i + 1) % L1dim == 0)
								fprintf(stdoutlog, "\n");
						}
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "\n\n");
						fclose(stdoutlog);*/

						int pop = 1; // Sets a population value
						//reorganize lattice polymer and lattice points
						cf_value = random_cf;
						lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, random_cf, random_cf + 3);
						lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, random_cf, random_cf + 3);

						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        fprintf(stdoutlog, "\n");
						//                        fprintf(stdoutlog, "happended: %i",4);
						//                        fprintf(stdoutlog, "\n");
						//                        fclose(stdoutlog);
						//                        stdoutlog = fopen(stdoutlogname, "a");
						//                        for (int i = 0; i < chain_length; i++)
						//                        {
						//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
						//                        }
						//                        fclose(stdoutlog);


						return 1;
					}
				}
			}
		}

		if (empty_index >= 0 && flip_index >= 0 & chain_location >= 0) // if corner flipped performed
		{


			//corner flip latticepoints updated
			//latticepoint[empty_index] = latticepoint[flip_index];
			//latticepoint[flip_index] = 0;

			//corner flip latticepoints updated
			//latticepolymer[5 * empty_index + 1] = latticepolymer[5 * flip_index + 1];
			//latticepolymer[5 * empty_index + 2] = latticepolymer[5 * flip_index + 2];
			//latticepolymer[5 * flip_index + 1] = -1;
			//latticepolymer[5 * flip_index + 2] = -1;

			//updating adjacent monomers
			//latticepolymer[5 * poly_chain[chain_location+1] + 2] = empty_index;
			//latticepolymer[5 * poly_chain[chain_location-1] + 1] = empty_index;

			poly_chain[chain_location] = empty_index;

			/* stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "4finished empty_index: %i emptyindexlatticepoint: %i  latticepolymer empty: %i  latticepolymerflip: %i\n",empty_index,latticepoint[empty_index] ,latticepolymer[5 * empty_index + 1] , latticepolymer[5 * flip_index + 1]);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/
			cf_value = random_cf;
			lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, random_cf, random_cf + 3);
			lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, random_cf, random_cf + 3);

			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "happended: %ib",4);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", poly_chain[i]);
			//            }
			//            fclose(stdoutlog);


			return 1;
		}

		check++;
		flip_index = -1;
		chain_location = -1;
		empty_index = -1;
		random_cf++;
		random_cf = (random_cf) % (chain_length - 2);
	}

	int pop = 1; // Sets a population value
	//reorganize lattice polymer and lattice points
	lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);
	lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, pop, 0, chain_length);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "happended: %i  f",4);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", poly_chain[i]);
	//    }
	//    fclose(stdoutlog);


	return 0;
}

int slither_f() //ovement of the chain from the head
{
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "\n");
	//fprintf(stdoutlog, "2");
	//fprintf(stdoutlog, "\n");
	//fclose(stdoutlog);

	int complete = 0; // check for if a movement action has been completed

	int sorround_check_prev = 100;
	int sorround_check = 100; // stores a random impossible value for sorround check if no values are available
	int random_slither = rand() % numberneighbors; // Picks a random offset direction to check movement direction options
	//for (int ii = 0; ii < numberneighbors; ii++)
	//{

	random_slither = (random_slither) % numberneighbors;//Changes the value of the randome direction with each movement
	//stores the value for random_slither if the checked index is empty

	while (neighbor[numberneighbors * poly_chain[0] + random_slither] == poly_chain[1])
	{
		random_slither = rand() % numberneighbors;
	}

	if (latticepoint[neighbor[numberneighbors * poly_chain[0] + random_slither]] == 0)
	{
		sorround_check = random_slither;
	}
	//random_slither++;
//}

//stdoutlog = fopen(stdoutlogname, "a");
//fprintf(stdoutlog, "\n");
//fprintf(stdoutlog, "3  sorround_check:%i  check1:%i   ",sorround_check, neighbor[4*poly_chain[0]+sorround_check]);
//fprintf(stdoutlog, "\n");
//fclose(stdoutlog);

	if (sorround_check != 100) // If the head of the chain is able to move
	{
		for (int y = 0; y < chain_length; y++) // stores the index values of the chain to a fake storage vector
		{
			pseudo_poly_chain[y] = poly_chain[y];
		}

		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < chain_length; i++)
		//{
		//    fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "same check");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < chain_length; i++)
		//{
		//    fprintf(stdoutlog, "%4i", poly_chain[i]);
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);

		poly_chain[0] = neighbor[numberneighbors * poly_chain[0] + sorround_check]; // updates the head of the chain to the new storage


		for (int iii = 1; iii < (chain_length); iii++)  // Updates the indexes of the remainder of the chain using an offset of the stored pseudo_poly_chain
		{
			poly_chain[iii] = pseudo_poly_chain[iii - 1];
		}

		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < chain_length; i++)
		//{
		//    fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "same check");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < chain_length; i++)
		//{
		//    fprintf(stdoutlog, "%4i", poly_chain[i]);
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);


		int pop = 1; // Sets a population value
		//reorganize lattice polymer and lattice points
		lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, 0, chain_length);
		lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);
		//for (int i = 0; i < numberspins; i++) //Clears all the indexes of the system
		//{
		//    latticepoint[i] = 0;
		//
		//}

		//for (int i = 0; i < chain_length; i++)
		//{
		   //latticepoint[poly_chain[i]] = pop;
		   //pop++;
		//   latticepoint[poly_chain[i]] = chain_sequence[i];
		//}

		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < numberspins; i++)
		//{
		//    fprintf(stdoutlog, "%4i", latticepoint[i]);
		//    if ((i + 1) % L1dim == 0)
		//        fprintf(stdoutlog, "\n");
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);

		complete = 1;


		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i",2);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);


		return 1;
	}
	else
	{

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i  f",2);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);


		return complete;
	}
}

int slither_b() //Same as Slither_f but works backwards
{
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "\n");
	//fprintf(stdoutlog, "2");
	//fprintf(stdoutlog, "\n");
	//fclose(stdoutlog);

	int complete = 0;

	int sorround_check_prev = 100;
	int sorround_check = 100;
	int random_slither = rand() % numberneighbors;
	//for (int ii = 0; ii < numberneighbors; ii++) {

	random_slither = (random_slither) % numberneighbors;

	while (neighbor[numberneighbors * poly_chain[chain_length - 1] + random_slither] == poly_chain[chain_length - 2])
	{
		random_slither = rand() % numberneighbors;
	}

	if (latticepoint[neighbor[numberneighbors * poly_chain[chain_length - 1] + random_slither]] == 0)
	{
		sorround_check = random_slither;
	}
	random_slither++;
	//}

	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "\n");
	//fprintf(stdoutlog, "3  sorround_check:%i  check1:%i   ",sorround_check, neighbor[4*poly_chain[0]+sorround_check]);
	//fprintf(stdoutlog, "\n");
	//fclose(stdoutlog);

	if (sorround_check != 100)
	{
		for (int y = 0; y < chain_length; y++)
		{
			pseudo_poly_chain[y] = poly_chain[y];
		}

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
		//fclose(stdoutlog);

		poly_chain[chain_length - 1] = neighbor[numberneighbors * poly_chain[chain_length - 1] + sorround_check];


		for (int iii = 1; iii < (chain_length); iii++)
		{
			poly_chain[chain_length - 1 - iii] = pseudo_poly_chain[chain_length - 1 - iii + 1];
		}

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
		//fclose(stdoutlog);
		int pop = 1;
		lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop, 0, chain_length);
		lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);

		//for (int i = 0; i < numberspins; i++)
		//{
		//    latticepoint[i] = 0;
		//}
		//pop = 1;
		//for (int i = 0; i < chain_length; i++)
		//{
			//latticepoint[poly_chain[i]] = pop;
			//pop++;
		//    latticepoint[poly_chain[i]] = chain_sequence[i];
		//}
		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//for (int i = 0; i < numberspins; i++)
		//{
		//    fprintf(stdoutlog, "%4i", latticepoint[i]);
		//    if ((i + 1) % L1dim == 0)
		//        fprintf(stdoutlog, "\n");
		//}
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);

		complete = 1;


		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i",3);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);


		return complete;

	}
	else
	{

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i  f",3);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);


		return complete;
	}
}

