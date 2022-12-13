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
#include "Mov_Pull.h"

int pull_2(int p_index, int c_index, int q, int c_chainlocation, int ht)
{
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "2");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);

	int sorround_check_prev = 100;
	int sorround_check = 100; // unoccupiable values that store the neighbor index assignment

	int random_nieghbour = (rand() % numberneighbors) + 1;  // Chooses a random starting direction

	//Turn into While loop
	//for (int ii = 0; ii < numberneighbors; ii++)
	//{
	random_nieghbour = (random_nieghbour) % numberneighbors; // reset the direction value for overflow so I can start with any value

	while (neighbor[numberneighbors * p_index + random_nieghbour] == c_index)
	{
		random_nieghbour = (rand() % numberneighbors);
	}

	if (latticepoint[neighbor[numberneighbors * p_index + random_nieghbour]] != latticepoint[c_index]) //checks the neighbour for the first point to ensure its not the other index
	{
		if (latticepoint[neighbor[numberneighbors * p_index + random_nieghbour]] == 0 && latticepoint[neighbor[numberneighbors * c_index + random_nieghbour]] == 0) //ensures both indexes sorround emptyspace
		{
			sorround_check = random_nieghbour;
		}
	}
	random_nieghbour++; // itterates the index if none is chosen
//}

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "\n");
//    fprintf(stdoutlog, "3  sorround_check:%i  check1:%i    check2:%i  p_index: %i  c_index: %i",sorround_check, neighbor[4*p_index+sorround_check], neighbor[4*c_index+sorround_check],p_index,c_index);
//    fprintf(stdoutlog, "\n");
//    fclose(stdoutlog);

	if (sorround_check != 100)
	{
		// stores a copy of the current chain
		if (q == 0)
		{
			for (int y = 0; y < chain_length; y++)
			{
				pseudo_poly_chain[y] = poly_chain[y];
			}
		}
		for (int y = 0; y < chain_length; y++)
		{
			pull_poly_chain[y] = poly_chain[y];
		}
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i    %i   %i    %i   %i   %i    %i    %i   rwhejt",1,q,p_index,c_index,latticepoint[neighbor[numberneighbors * p_index + sorround_check]],latticepoint[neighbor[numberneighbors * c_index + sorround_check]],neighbor[numberneighbors * p_index + sorround_check],neighbor[numberneighbors * c_index + sorround_check]);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
		//        }
		//        fprintf(stdoutlog, "\n");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", pull_poly_chain[i]);
		//        }
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
				/*stdoutlog = fopen(stdoutlogname, "a");
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
				fclose(stdoutlog);*/
		int hort;
		if (q == 0)
		{
			hort = (rand() % 3) + 1;
		}
		if (q == 1)
		{
			hort = (rand() % 2) + 1;
		}
		if (ht == 1)
		{
			hort = 2;
		}
		if (ht == 2)
		{
			hort = 1;
		}
		if (ht == 0)
		{
			if (c_chainlocation > chain_length - 3)
			{
				hort = 2;
			}
			if (c_chainlocation < 3)
			{
				hort = 1;
			}
		}

		if (q == 0)
		{
			if (hort == 2)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				int stall = 0; // UNUSED
				for (int iii = 0; iii < (chain_length - 2); iii++)
				{
					stall = iii;
					poly_chain[iii + ad] = pseudo_poly_chain[iii]; //Potential offset for poly chain
					if (p_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
						poly_chain[iii + 1] = neighbor[p_index * numberneighbors + sorround_check];
						poly_chain[iii + 2] = neighbor[c_index * numberneighbors + sorround_check];
					}
				}
				pull_loc = c_chainlocation;
				pull_hort = hort;
				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, c_chainlocation - 1, chain_length);
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, c_chainlocation - 1, chain_length);
			}
			if (hort == 1)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				int stall = 0; // UNUSED
				for (int iii = chain_length - 1; iii > 1; iii--)
				{
					stall = iii;
					poly_chain[iii - ad] = pseudo_poly_chain[iii]; //Potential offset for poly chain
					if (c_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
						poly_chain[iii - 2] = neighbor[p_index * numberneighbors + sorround_check];
						poly_chain[iii - 1] = neighbor[c_index * numberneighbors + sorround_check];
					}
				}
				pull_loc = c_chainlocation;
				pull_hort = hort;
				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, c_chainlocation + 3);
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, c_chainlocation + 3);
			}

			if (hort == 3)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				int stall = 0; // UNUSED
				for (int iii = 0; iii < (chain_length - 2); iii++)
				{
					stall = iii;
					poly_chain[iii + ad] = pseudo_poly_chain[iii + 1]; //Potential offset for poly chain
					if (p_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
						poly_chain[iii + 1] = neighbor[p_index * numberneighbors + sorround_check];
						poly_chain[iii + 2] = neighbor[c_index * numberneighbors + sorround_check];
					}
				}
				pull_loc = c_chainlocation;
				pull_hort = hort;
				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
				lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);
			}
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//             fprintf(stdoutlog, "regular %i %i %i  %i\n",hort,happened3,c_index,p_index);
			//             fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < chain_length; i++)
			//                    {
			//                        fprintf(stdoutlog, "%4i", pull_poly_chain[i]);
			//                   }
			//             fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < chain_length; i++)
			//                    {
			//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
			//                   }
			//             fprintf(stdoutlog, "\n");
			//             fclose(stdoutlog);
		}

		if (q == 1)
		{
			if (hort == 2)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				int stall = 0; // UNUSED
				for (int iii = 0; iii < (chain_length - 2); iii++)
				{
					stall = iii;
					poly_chain[iii + ad] = pull_poly_chain[iii]; //Potential offset for poly chain

//                  stdoutlog = fopen(stdoutlogname, "a");
//                  fprintf(stdoutlog, "%4i : %i", poly_chain[iii+ad],iii+ad);
//                  fclose(stdoutlog);

					if (p_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
						poly_chain[iii + 1] = neighbor[p_index * numberneighbors + sorround_check];
						poly_chain[iii + 2] = neighbor[c_index * numberneighbors + sorround_check];
					}
				}
				//lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1,0,chain_length);
				//lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain,0,chain_length);
			}
			if (hort == 1)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				int stall = 0; // UNUSED
				for (int iii = chain_length - 1; iii > 1; iii--)
				{
					stall = iii;
					poly_chain[iii - ad] = pull_poly_chain[iii]; //Potential offset for poly chain

//                  stdoutlog = fopen(stdoutlogname, "a");
//                  fprintf(stdoutlog, "%4i : %i", poly_chain[iii+ad],iii+ad);
//                  fclose(stdoutlog);

					if (c_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
						poly_chain[iii - 2] = neighbor[p_index * numberneighbors + sorround_check];
						poly_chain[iii - 1] = neighbor[c_index * numberneighbors + sorround_check];
					}
				}
				//lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1,0,chain_length);
				//lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain,0,chain_length);

			}
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//             fprintf(stdoutlog, "relax %i  %i  %i   %i\n",hort,happened3,c_index,p_index);
			//             fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < chain_length; i++)
			//                    {
			//                        fprintf(stdoutlog, "%4i", pull_poly_chain[i]);
			//                   }
			//             fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < chain_length; i++)
			//                    {
			//                        fprintf(stdoutlog, "%4i", poly_chain[i]);
			//                   }
			//             fprintf(stdoutlog, "\n");
			//             fclose(stdoutlog);
		}

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i immedialey after",1);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);

				/*stdoutlog = fopen(stdoutlogname, "a");
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
				fclose(stdoutlog);*/

				//resets the system lattice and repolpoulate the points of the chain
		int pop = 1;

		//for (int i = 0; i < chain_length; i++)
		//{
		//    latticepoint[pseudo_poly_chain[i]] = 0;
		//}
		//pop = 1;
		//for (int i = 0; i < chain_length; i++)
		//{
		//    latticepoint[poly_chain[i]] = pop;
		//  pop++;
		//    latticepoint[poly_chain[i]] = chain_sequence[i];
		//}
		/*stdoutlog = fopen(stdoutlogname, "a");
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


		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i",1);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);


		return 5;
	}
	else
	{

		//        stdoutlog = fopen(stdoutlogname, "a");
		//        fprintf(stdoutlog, "\n");
		//        fprintf(stdoutlog, "happended: %i f",1);
		//        fprintf(stdoutlog, "\n");
		//        fclose(stdoutlog);
		//        stdoutlog = fopen(stdoutlogname, "a");
		//        for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//        }
		//        fclose(stdoutlog);


		return 0;
	}
}

int pull_1(int q, int ht) //the chain forms a kink along an empty plane along random points crucial to avoid bottlenecking
{

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "pull_1(%i)",q);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    
	int complete = 0; // complete check
	int plane_check; // checks to see if adjacent variables are on the same plane
	int current_index = -1;
	int holder = -1; //stores the return value for the pull2() function
	int i = 1;

	int random_pull = rand() % (chain_length); // starts the indexing at a random location along the chain
	//must be two indexes away from the end as the end is coiled to make the chain pull
	int prev_index = poly_chain[random_pull]; // stores the value one of random pull
	int prev_pull = random_pull;
	random_pull = (random_pull + 1) % (chain_length); // stores the value of 1 + random pull +1 so it starts away from prev_pull

	int relax_check = 0;

	while (i < (chain_length)) // cheks
	{
		relax_check = 0;

		if (ht == 0)
		{
			prev_pull = (random_pull) % (chain_length);
			prev_index = poly_chain[prev_pull];//resets value of
			// prev_pull = random_pull;
			random_pull = random_pull + 1;
			random_pull = (random_pull) % (chain_length);
			current_index = poly_chain[random_pull];
		}
		if (ht == 1)
		{
			prev_pull = (random_pull) % (chain_length - 2);
			prev_index = poly_chain[prev_pull];//resets value of
			// prev_pull = random_pull;
			random_pull = random_pull + 1;
			random_pull = (random_pull) % (chain_length - 2);
			current_index = poly_chain[random_pull];
		}
		if (ht == 2)
		{
			prev_pull = (random_pull) % (chain_length);
			prev_index = poly_chain[prev_pull];//resets value of
			// prev_pull = random_pull;
			random_pull = random_pull + 1;
			random_pull = (random_pull) % (chain_length);
			current_index = poly_chain[random_pull];
			if (random_pull < 3)
			{
				relax_check = 2;
			}
		}

		while ((random_pull == 0 && (prev_pull + 1) != random_pull && relax_check != 0))
		{
			random_pull = rand() % (chain_length);

			relax_check = 0;

			if (ht == 0)
			{
				prev_pull = (random_pull) % (chain_length);
				prev_index = poly_chain[prev_pull];//resets value of
				// prev_pull = random_pull;
				random_pull = random_pull + 1;
				random_pull = (random_pull) % (chain_length);
				current_index = poly_chain[random_pull];
			}
			if (ht == 1)
			{
				prev_pull = (random_pull) % (chain_length - 2);
				prev_index = poly_chain[prev_pull];//resets value of
				// prev_pull = random_pull;
				random_pull = random_pull + 1;
				random_pull = (random_pull) % (chain_length - 2);
				current_index = poly_chain[random_pull];
			}
			if (ht == 2)
			{
				prev_pull = (random_pull) % (chain_length);
				prev_index = poly_chain[prev_pull];//resets value of
				// prev_pull = random_pull;
				random_pull = random_pull + 1;
				random_pull = (random_pull) % (chain_length);
				current_index = poly_chain[random_pull];
				if (random_pull < 3)
				{
					relax_check = 2;
				}
			}
		}

		/*stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "0  pindex:%i   cindex:%i    random_pull:%i",prev_index,current_index,random_pull);
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/

		//checks to ensure indexes are nnot the tail and head
		if (random_pull != 0 && (prev_pull + 1) == random_pull && relax_check == 0)
		{
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "1  pindex:%i   cindex:%i",prev_index,current_index);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);

			holder = pull_2(prev_index, current_index, q, random_pull, ht);//completes the actual pull movement
			if (holder == 5)
			{
				complete = 7;
				return complete;
			}
			//return complete;
		}
		i++;
	}
	return complete;
}

int pull_relaxation(int k) //this will accomate both relaxation and repulling of the polymer chain
{

	/*stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "cornerflip\n");
	for (int i = 0; i < numberspins; i++)
	{
		fprintf(stdoutlog, "%4i", latticepoint[i]);
		if ((i + 1) % L1dim == 0)
			fprintf(stdoutlog, "\n");

	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "%i\n", poly_en());
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);*/

	int complete = 0; // completion check
	int current_index;
	int i = 1;
	int random_pull = rand() % (chain_length - 3); // starts the indexing at a random location along the chain
	//must be two indexes away from the end as the end is coiled to make the chain pull
	int prev_index = poly_chain[random_pull + 3]; // stores the value one of random pull
	int prev_pull = random_pull + 3;
	random_pull = (random_pull + 1) % (chain_length - 3); // stores the value of 1 + random pull

	int h_t = (rand() % (2)) + 1;

	int hort = (rand() % (2)) + 1;

	h_t = k;

	for (int y = 0; y < chain_length; y++)
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}

	while (i < (chain_length - 3)) // cheks
	{
		random_pull++;
		random_pull = (random_pull) % (chain_length - 3); // gives chain poistions
		current_index = poly_chain[random_pull];
		prev_pull = (random_pull + 3) % (chain_length);
		prev_index = poly_chain[((random_pull + 3) % (chain_length))];//resets value of

		while ((prev_pull - 3) != random_pull)
		{
			random_pull = rand() % (chain_length - 3); // starts the indexing at a random location along the chain
			random_pull = (random_pull) % (chain_length - 3); // gives chain poistions
			current_index = poly_chain[random_pull];
			prev_pull = (random_pull + 3) % (chain_length);
			prev_index = poly_chain[((random_pull + 3) % (chain_length))];//resets value of
		}

		/*stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "-2  pindex:%i   cindex:%i    random_pull:%i",prev_index,current_index,random_pull);
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);*/
		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "cin: %i  prev_in:%i",current_index,prev_index);
		//fprintf(stdoutlog, "\n");
		//fclose(stdoutlog);
		//checks to ensure indexes are always
		int neighborcheck = 0;
		for (int i = 0; i < numberneighbors; i++)
		{
			//stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\n");
			//fprintf(stdoutlog, "i: %i   n:%i",i,neighbor[numberneighbors*poly_chain[current_index] + i]);
			//fprintf(stdoutlog, "\n");
			//fclose(stdoutlog);
			if (neighbor[numberneighbors * poly_chain[random_pull] + i] == prev_index && poly_chain[random_pull + 3] == prev_index)
			{
				neighborcheck = 1;
			}
		}

		if ((prev_pull > (random_pull)) && ((prev_pull - 3) == random_pull) && neighborcheck == 1)//((current_index - prev_index) == ((L1dim) * L1dim * L1dim)) || ((current_index - prev_index) == (-(L1dim) * L1dim))))
		{
			if (random_pull == 0)
			{
				hort = 1;
			}
			if (prev_pull == chain_length - 1)
			{
				hort = 2;
			}

			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "happended:%i   %i   %i",5,current_index,prev_index);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", poly_chain[i]);
			//            }
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);


						/*stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "-1  cindex - pindex: %i ",current_index-prev_index);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

			for (int y = 0; y < chain_length; y++)
			{
				pseudo_poly_chain[y] = poly_chain[y];
			}

			/*stdoutlog = fopen(stdoutlogname, "a");
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
			fclose(stdoutlog);*/

			if (hort == 1)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				for (int iii = 0; iii < (chain_length - 2); iii++)
				{
					poly_chain[iii] = pseudo_poly_chain[iii + ad]; //Potential offset for poly chain
					if (current_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
					}

					//                stdoutlog = fopen(stdoutlogname, "a");
					//                fprintf(stdoutlog, "%4i", poly_chain[iii]);
					//                fclose(stdoutlog);

				}
				poly_chain[chain_length - 2] = current_index;
				poly_chain[chain_length - 1] = prev_index;

				//                stdoutlog = fopen(stdoutlogname, "a");
				//                fprintf(stdoutlog, "\n");
				//                fprintf(stdoutlog, "\n");
				//                 fprintf(stdoutlog, "really relax INIT %i  %i  %i   %i\n",hort,happened3,current_index,prev_index);
				//                 fprintf(stdoutlog, "\n");
				//                for (int i = 0; i < chain_length; i++)
				//                        {
				//                            fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
				//                       }
				//                 fprintf(stdoutlog, "\n");
				//                for (int i = 0; i < chain_length; i++)
				//                        {
				//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                       }
				//                 fprintf(stdoutlog, "\n");
				//                 fprintf(stdoutlog, "\n");
				//                 fclose(stdoutlog);


				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
			}

			if (hort == 2)
			{
				int ad = 0; // adjustment value for pull move to offset the poly_chain assignment and truncate the final two indexes when performed
				for (int iii = chain_length - 1; iii > (1); iii--)
				{
					poly_chain[iii] = pseudo_poly_chain[iii - ad]; //Potential offset for poly chain
					if (prev_index == poly_chain[iii] && ad == 0) // if statement the adds the two new indexs and offsets the poly chain indexes
					{
						ad = 2;
					}

					//                  stdoutlog = fopen(stdoutlogname, "a");
					//                  fprintf(stdoutlog, "%4i", poly_chain[iii]);
					//                  fclose(stdoutlog);

				}
				poly_chain[1] = current_index;
				poly_chain[0] = prev_index;

				//                stdoutlog = fopen(stdoutlogname, "a");
				//                fprintf(stdoutlog, "\n");
				//                fprintf(stdoutlog, "\n");
				//                 fprintf(stdoutlog, "relax really INIT %i  %i  %i   %i\n",hort,happened3,current_index,prev_index);
				//                 fprintf(stdoutlog, "\n");
				//                for (int i = 0; i < chain_length; i++)
				//                        {
				//                            fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
				//                       }
				//                 fprintf(stdoutlog, "\n");
				//                for (int i = 0; i < chain_length; i++)
				//                        {
				//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                       }
				//                 fprintf(stdoutlog, "\n");
				//                 fprintf(stdoutlog, "\n");
				//                 fclose(stdoutlog);

				lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
			}
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "relax INIT %i  %i  %i   %i\n",hort,happened3,current_index,prev_index);
			//            fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			//            }
			//            fprintf(stdoutlog, "\n");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", poly_chain[i]);
			//            }
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "happended: %i ihhi",5);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", poly_chain[i]);
			//            }
			//            fclose(stdoutlog);

						/*stdoutlog = fopen(stdoutlogname, "a");
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
						fclose(stdoutlog);*/

						/*stdoutlog = fopen(stdoutlogname, "a");
						fprintf(stdoutlog, "\n");
						fprintf(stdoutlog, "h_t: %i",h_t);
						fprintf(stdoutlog, "\n");
						fclose(stdoutlog);*/

			if (h_t == 1)
			{


				/*stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "pull commence");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/


				//                stdoutlog = fopen(stdoutlogname, "a");
				//                fprintf(stdoutlog, "\n");
				//                fprintf(stdoutlog, "happended: %i  ajooj",5);
				//                fprintf(stdoutlog, "\n");
				//                fclose(stdoutlog);
				//                stdoutlog = fopen(stdoutlogname, "a");
				//                for (int i = 0; i < chain_length; i++)
				//                {
				//                    fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                }
				//                fclose(stdoutlog);


				happened2 = 0;
				complete = pull_1(1, hort);

				if (complete == 0) {
					if (hort == 2)
					{
						lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, 0, prev_pull + 1);
						//lp_reorg_non_rebrid(poly_chain,pseudo_poly_chain,0,prev_pull+1);
					}
					if (hort == 1)
					{
						lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, random_pull, chain_length);
						//lp_reorg_non_rebrid(poly_chain,pseudo_poly_chain,random_pull,chain_length);
					}

					happened2 = 1;
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "happended: %i  afailinm0i",5);
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
				}
				else
				{
					happened2 = 3;
					lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, chain_length);
					lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, chain_length);
					//                    stdoutlog = fopen(stdoutlogname, "a");
					//                    fprintf(stdoutlog, "\n");
					//                    fprintf(stdoutlog, "happended: %i  1success0j0ni",5);
					//                    fprintf(stdoutlog, "\n");
					//                    fclose(stdoutlog);
				}

				/*stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "pull ends");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/
				return complete;
			}

			if (h_t == 2)
			{
				/*stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "relaxation commence");
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);*/

				int o_iii = rand() % numberneighbors;
				int anchor = poly_chain[chain_length - 3];
				//for (int iii = 0; iii < numberneighbors * numberneighbors; iii++)
				//{
					//neighbor list for polymer adjacent to the corner
					//if (iii != 0 && iii % numberneighbors == 0)
					//{
					//    o_iii++;
					//}

				int act = rand() % numberneighbors;

				while (hort == 2 && neighbor[numberneighbors * neighbor[numberneighbors * poly_chain[2] + o_iii] + act] == poly_chain[2])
				{
					o_iii = rand() % numberneighbors;
					act = rand() % numberneighbors;
				}
				while (hort == 1 && neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act] == anchor)
				{
					o_iii = rand() % numberneighbors;
					act = rand() % numberneighbors;
				}
				if (hort == 2 && latticepoint[neighbor[numberneighbors * poly_chain[2] + o_iii]] == 0 && latticepoint[neighbor[numberneighbors * neighbor[numberneighbors * poly_chain[2] + o_iii] + act]] == 0)
				{
					//                        stdoutlog = fopen(stdoutlogname, "a");
					//                        fprintf(stdoutlog, "\n");
					//                        fprintf(stdoutlog, "happended: %i  b089jh9    %i    %i    %i     %i %i   %i    %i",5,neighbor[numberneighbors * anchor + o_iii],neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act],latticepoint[neighbor[numberneighbors * anchor + o_iii]],latticepoint[neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + iii]], o_iii,act,anchor);
					//                        fprintf(stdoutlog, "\n");
					//                        fclose(stdoutlog);
					//                        stdoutlog = fopen(stdoutlogname, "a");
					//                        fprintf(stdoutlog, "\n");
					//                        fprintf(stdoutlog, "First index: %i last index: %i",neighbor[numberneighbors * poly_chain[2] + o_iii],neighbor[numberneighbors * neighbor[numberneighbors * poly_chain[2] + o_iii] + iii]);
					//                        fprintf(stdoutlog, "\n");
					//                        fclose(stdoutlog);
					happened2 = 4;
					poly_chain[1] = neighbor[numberneighbors * poly_chain[2] + o_iii];
					poly_chain[0] = neighbor[numberneighbors * neighbor[numberneighbors * poly_chain[2] + o_iii] + act];
					complete = 1;
					pull_relax_loc = prev_pull;
					pull_relax_hort = 2;
					lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, 0, prev_pull + 1);
					lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 0, prev_pull + 1);


					//                        stdoutlog = fopen(stdoutlogname, "a");
					//                        fprintf(stdoutlog, "\n");
					//                        fprintf(stdoutlog, "happended: %i  b089jh9    %i    %i    %i     %i %i   %i    %i",5,neighbor[numberneighbors * anchor + o_iii],neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + iii],latticepoint[neighbor[numberneighbors * anchor + o_iii]],latticepoint[neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act]], o_iii,act,anchor);
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

				if (hort == 1 && latticepoint[neighbor[numberneighbors * anchor + o_iii]] == 0 && latticepoint[neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act]] == 0)
				{
					//                        stdoutlog = fopen(stdoutlogname, "a");
					//                        fprintf(stdoutlog, "\n");
					//                        fprintf(stdoutlog, "happended: %i  b089jh9    %i    %i    %i     %i %i   %i    %i",5,neighbor[numberneighbors * anchor + o_iii],neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act],latticepoint[neighbor[numberneighbors * anchor + o_iii]],latticepoint[neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + iii]], o_iii,act,anchor);
					//                        fprintf(stdoutlog, "\n");
					//                        fclose(stdoutlog);
											/*stdoutlog = fopen(stdoutlogname, "a");
											fprintf(stdoutlog, "\n");
											fprintf(stdoutlog, "First index: %i last index: %i",neighbor[numberneighbors * anchor + o_iii],neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + iii]);
											fprintf(stdoutlog, "\n");
											fclose(stdoutlog);*/
					happened2 = 4;
					poly_chain[chain_length - 2] = neighbor[numberneighbors * anchor + o_iii];
					poly_chain[chain_length - 1] = neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act];
					complete = 1;
					pull_relax_loc = random_pull;
					pull_relax_hort = 1;
					lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, random_pull, chain_length);
					lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, random_pull, chain_length);


					//                        stdoutlog = fopen(stdoutlogname, "a");
					//                        fprintf(stdoutlog, "\n");
					//                        fprintf(stdoutlog, "happended: %i  b089jh9    %i    %i    %i     %i %i   %i    %i",5,neighbor[numberneighbors * anchor + o_iii],neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + iii],latticepoint[neighbor[numberneighbors * anchor + o_iii]],latticepoint[neighbor[numberneighbors * neighbor[numberneighbors * anchor + o_iii] + act]], o_iii,act,anchor);
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

				//}
				if (hort == 2)
				{
					lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, 0, prev_pull + 1);
					//lp_reorg_non_rebrid(poly_chain,pseudo_poly_chain,0,prev_pull+1);
				}
				if (hort == 1)
				{
					lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, random_pull, chain_length);
					//lp_reorg_non_rebrid(poly_chain,pseudo_poly_chain,random_pull,chain_length);
				}
				happened2 = 5;



				//                stdoutlog = fopen(stdoutlogname, "a");
				//                fprintf(stdoutlog, "\n");
				//                fprintf(stdoutlog, "happended: %i  fb08ino",5);
				//                fprintf(stdoutlog, "\n");
				//                fclose(stdoutlog);
				//                stdoutlog = fopen(stdoutlogname, "a");
				//                for (int i = 0; i < chain_length; i++)
				//                {
				//                    fprintf(stdoutlog, "%4i", poly_chain[i]);
				//                }
				//                fclose(stdoutlog);


				return complete;
			}
			/*stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "1  pindex:%i   cindex:%i",prev_index,current_index);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);*/


			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "happended: %i unknown-9jion",5);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);

		}
		i++;
		return complete;
	}

	happened2 = 6;
	lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, 0, chain_length);
	lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 0, chain_length);

	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "happended: %i  fj0ionjih",5);
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%4i", poly_chain[i]);
	//    }
	//    fclose(stdoutlog);


	return complete;
}

