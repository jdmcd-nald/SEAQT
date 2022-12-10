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
#include "Reorg.h"
#include "EXWL.h"


void pseudowl() // A Fake WL Function used to explore the energy landscape before the main WL operation to better enesure all lnge values of the system have been explored
{
	int maxe = Eglobalmin; //holds the index of the maximum energy found by the system used to reset the for loop counter if a higher energy is found

	int qhold = q; // Holds the value of q from before (could be better coded but as far as I know would require significant retooling of some functions)
	
	double lnf = 1.0;

	double wsk, dice; // wsk: probability  
	int wiggle;
	int wiggletwo;

	// terminal modification factor             
    double check_flatness_every = 50;   // in number of sweeps changed to a large number to encourage visitation of all levels (reduces time in the long run)
	int backup;
	int backuptwo;

	int eold, energie;
	energie = poly_en();
	energie -= Eglobalmin;
	eold = energie;

	int swtch;
	int found = 0;

	long long iqq = 0;

	for (int i = 0; i <= -Eglobalmin; i++) HE[i] = 0; //init H(E)

	for (int u = 0; u < chain_length; u++) {
		wl_pseudo_chain[u] = poly_chain[u];
	}
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
	for (int k = 0; k < check_flatness_every; k++)
	{
		for (int i = 0; i < L1dim; i++) // this does 1 MC sweep
		{


			energie = eold + propose_update(eold); // calculate energy of updated configuration
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "\n");
//            fprintf(stdoutlog, "energie: %i  &  poly_en():%i",energie,poly_en());
//            fprintf(stdoutlog, "\n");
//            fclose(stdoutlog);


			//eye();
			/*stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "energie: %i  actual: %i",energie,poly_en()-Eglobalmin);
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			//          stdoutlog = fopen(stdoutlogname, "a");
			//          fprintf(stdoutlog, "\n");
			//          fprintf(stdoutlog, "energie: %i  &  poly_en():%i",energie,poly_en());
			//          fprintf(stdoutlog, "\n");
			//          fclose(stdoutlog);
						/*
						if(energie!=poly_en()-Eglobalmin)
						{
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
							fprintf(stdoutlog, "en: %i, poly_en()-Eglobmin: %i  , energie:%i  ,eold : %i, happened: %i  happened_2: %i ",poly_en(),poly_en()-Eglobalmin,energie,eold,happened,happened2);
							fprintf(stdoutlog, "\n");
							fprintf(stdoutlog, "\n");
							fprintf(stdoutlog, "\n");
							fprintf(stdoutlog, "\n");
							for (int i = 0; i < chain_length; i++)
							{
								fprintf(stdoutlog, "%4i\t", wl_pseudo_chain[i]);
							}
							for (int i = 0; i < chain_length; i++)
							{
								fprintf(stdoutlog, "%4i\t", poly_chain[i]);
							}
							fclose(stdoutlog);
							MPI_Abort(MPI_COMM_WORLD, 1);
						}

						eye();
						 */
						 //poly_mov();
						 //energie = poly_en()-Eglobalmin;
						 // reject because you got out of bounds of energy window

			if ((energie > Emaxindex) || (energie < Eminindex)) // boundary check
			{
				lp_reorg_non_rebrid(poly_chain, wl_pseudo_chain, 0, chain_length); // rejection needs to change the indexs back to the former storage value
                int pop = 1;
                lattice_reorg_non_rebrid(poly_chain, wl_pseudo_chain, pop, 0, chain_length);
                // makes a copy of the lattice if the change is rejected
                for (int u = 0; u < chain_length; u++) {
                    poly_chain[u] = wl_pseudo_chain[u];
                    latticepoint[poly_chain[u]] = chain_sequence[u];
                }
				//                    stdoutlog = fopen(stdoutlogname, "a");
				//                    fprintf(stdoutlog, "\n");
				//                    fprintf(stdoutlog, "reject wl");
				//                    fprintf(stdoutlog, "\n");
				//                    fclose(stdoutlog);
				energie = eold;
				lngE[eold] = 1; // the real lngE of the system is set to a value of 1 if visitation occured at a specific energy (can be subtracted later) to force visitation in the real wl process
			}
			else // calculate acceptance propbability
			{
				lngE[energie] = 1;
				if (energie > maxe) // checks to see if new energy is higher than max e
				{
					maxe = energie;
					k = 0; // reset for loop iterative
				}

				// the regular wl process is kept using a place holder pseudolngE to force the system to not move to move to higher energies
				dice = (1.0 * rand() / (RAND_MAX + 1.0)); // roll the dice
				wsk = exp(pseudolngE[eold] - pseudolngE[energie]); // calculate acceptance probability
				if (dice > wsk || MoveProposal == 0)  // reject
				{
					if (MoveProposal > 0)
					{
						lp_reorg_non_rebrid(poly_chain, wl_pseudo_chain, 0, chain_length); // rejection needs to change the indexs back to the former storage value
                        int pop = 1;
                        lattice_reorg_non_rebrid(poly_chain, wl_pseudo_chain, pop, 0, chain_length);
                        // makes a copy of the lattice if the change is rejected
                        for (int u = 0; u < chain_length; u++) {
                            poly_chain[u] = wl_pseudo_chain[u];
                            latticepoint[poly_chain[u]] = chain_sequence[u];
                        }
						//sysparam(energie);
					}

					energie = eold;
				}
				else {
					//eye();
					eold = energie; // accept
					//sysparam(energie);
					for (int u = 0; u < chain_length; u++) {
						wl_pseudo_chain[u] = poly_chain[u];
					}
				}

				//if ((energie > eold))
				//{
				//	latticepoint[wiggle] = backup;
				//	latticepoint[wiggletwo] = backuptwo;
				//	energie = eold;
				//}
				//if ((energie < eold))
				//{
				//	k = 0;
				//	eold = energie;
				//}
			}

			pseudolngE[energie] += lnf;
			lngE[energie] = 1;


			if (iqq % (4 * 800000) == 0)
			{
				//if (myid == 0 || myid == 1 || myid == 2) { lngE[4] = 1; }
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, "iqq: %lli", iqq);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim, q, myid, -1); // Preliminary output file that allows for process inspection
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
						fprintf(file, "%i\t%i\t%lld\t%e\t%e\t%i\t%f\t%f\n", i, (i + (Eglobalmin)), iqq, pseudolngE[i], lngE[i], k, rog[i] / visits[i], tortuosity[i] / visits[i]);
					}
					fclose(file);
				}
			}
			iqq++;
		}

	}


	// The occucpied placeholder lnge and the values of the real lnge are out putted into a prelim file
	sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim, q, myid, 0);
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
			fprintf(file, "%i\t%i\t%e\t%e\t%e\n", i, (i + (Eglobalmin)), 0.0, lngE[i], pseudolngE[i]);
		}
		fclose(file);
	}
}

int histflat(int imin, int imax, double ratio)
{
	// checks flatness of histograms for _all_ walkers in energy window
	int myflat, otherflat;
	int flatproc, ioffset;
	int flatprocr, ioffsetr;
	int flatness_crit = 1;

    int ddflat;
    
	int multimyflat = 0; // bool value of wheter in of the concurrent walkers in a window are flat

	int merge_crit = 1; // merge_hist criteria

	// check own flatness first
	if (flatness_crit == 0)       // Zhou criterion
	{
		// NOT AVAILABLE YET!
		// TODO: CALCULATE histmin FOR ZHOU CRITERION
		//    myflat=1;
		//    for (int x=imin;x<=imax;x++)
		//    if (HE[x]<histmin) myflat=0;
	}
	else if (flatness_crit == 1)  // "Original" percentage criterion
	{
		myflat = 1;
		double minval;
		minval = HE[imin];        // take GS hits as reference


		/*if (myid == 0)
		{
			minval = HE[10];
		}*/

		minval = DBL_MAX; // Ridiculously large arbitrary value
		for (int x = imin; x <= imax; x++)
		{
			if ((lngE[x] > 0.5 || HE[x] > 0.5) && HE[x] <= minval)
			{
				minval = HE[x];
			}
		}

		double average = 0.0;
		double count = 0.0;

		for (int x = imin; x <= imax; x++)
		{
			/*if (((x > 5) || (x == 4)) && (HE[x] < minval))
				minval = HE[x];*/
				//(I am not sure right now why I included the first condition ...)

			if (lngE[x] > 0.5 || HE[x] > 0.5)
			{
				average += HE[x];
				count++;
			}
		}
		average /= count;

		flatratio = (ratio * average);
		flatmin = (minval);

		if (minval < (ratio * average))
			myflat = 0;
	}

	// now talk to all the other walkers in the energy window
	// (! this whole thing can be reduced to an MPI_Allreduce once there 
	// are separate communicators for energy windows !)
    if (local_dup_group == localgroup)
    {
        if (merge_crit == 1 && merge_hists == 1)                   // check flatness of other walkers in window
        {
            if (myid == dup_headproc)             // 'root' in energy window, receive individual flatnesses
            {
                if (myflat == 1) //Main node multimyflat (checks for the flat process in an energy window
                {
                    flatproc = myid; // id of flat process
                    ioffset = dup_headproc; // offset from the main energy window node
                    multimyflat = 1;
                }

                for (int i = (dup_headproc + 1); i < numprocs; i++)
                {
                    MPI_Recv(&otherflat, 1, MPI_INT, i, 66, MPI_COMM_WORLD, &status);

                    if (otherflat == 1) //sets the value of the two variable based on information recieved from the process each one is communicating with
                    {
                        flatproc = i;
                        ioffset = myid - dup_headproc;
                        multimyflat = 1;
                    }

                    myflat *= otherflat;        // all have to be '1' to give '1' in the end (manual '&&')
                }
                for (int i = (dup_headproc + 1); i < numprocs; i++) // and let everybody know
                {
                    MPI_Send(&myflat, 1, MPI_INT, i, 88, MPI_COMM_WORLD);
                    MPI_Send(&multimyflat, 1, MPI_INT, i, 86, MPI_COMM_WORLD);
                    if (multimyflat == 1) // if multimyflat is found to be equal to one flatproc and ioffset
                    {
                        MPI_Send(&flatproc, 1, MPI_INT, i, 90, MPI_COMM_WORLD);
                        MPI_Send(&ioffset, 1, MPI_INT, i, 92, MPI_COMM_WORLD);
                    }
                }
            }
            else                                // send individual flatness and receive merged status variables
            {
                MPI_Send(&myflat, 1, MPI_INT, dup_headproc, 66, MPI_COMM_WORLD);
                MPI_Recv(&otherflat, 1, MPI_INT, dup_headproc, 88, MPI_COMM_WORLD, &status);
                MPI_Recv(&multimyflat, 1, MPI_INT, dup_headproc, 86, MPI_COMM_WORLD, &status);
                if (multimyflat == 1)
                {
                    MPI_Recv(&flatprocr, 1, MPI_INT, dup_headproc, 90, MPI_COMM_WORLD, &status);
                    flatproc = flatprocr;
                    MPI_Recv(&ioffsetr, 1, MPI_INT, dup_headproc, 92, MPI_COMM_WORLD, &status);
                    ioffset = ioffsetr;
                }
                myflat = otherflat;             // replace individual flatness by merged
            }


            if (multimyflat == 1)
            {
                if (myid != flatproc) // non flat process recieve merged flat process density of states and the myflat status flagged is called to perform the lnge merging procedures in the wl routine
                {
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "Proc %3i: recievied flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\t flatproc: %3i\n", myid, flatproc, ioffset, multimyflat, flatproc);
                    fclose(stdoutlog);
                    MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, flatproc, 17, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j]; // overrides density of states of non flat processes
                    myflat = 1;
                }
                else // flat process sends out its density of states values
                {
                    for (int i = (dup_headproc); i < numprocs; i++)
                    {
                        if (myid == flatproc && i != flatproc)
                        {
                            stdoutlog = fopen(stdoutlogname, "a");
                            fprintf(stdoutlog, "Proc %3i: sent flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\t flatproc: %3i\t %i\n", myid, myid - (myid % multiple) + i, ioffset, multimyflat, flatproc,i);
                            fclose(stdoutlog);
                            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, i, 17, MPI_COMM_WORLD);
                            myflat = 1;
                        }
                    }
                }
            }
        }
    }
    
	if (multiple > 1)
	{
		if (merge_crit == 0 && merge_hists == 1)                   // check flatness of other walkers in window
		{
			if (myid % multiple == 0)             // 'root' in energy window, receive individual flatnesses
			{
				for (int i = 1; i < multiple; i++)
				{
					MPI_Recv(&otherflat, 1, MPI_INT, myid + i, 66, MPI_COMM_WORLD, &status);
					myflat *= otherflat;        // all have to be '1' to give '1' in the end (manual '&&')
				}
				for (int i = 1; i < multiple; i++)  // and let everybody know
				{
					MPI_Send(&myflat, 1, MPI_INT, myid + i, 88, MPI_COMM_WORLD);
				}
			}
			else                                // send individual flatness and receive 'merged' flatness
			{
				MPI_Send(&myflat, 1, MPI_INT, myid - (myid % multiple), 66, MPI_COMM_WORLD);
				MPI_Recv(&otherflat, 1, MPI_INT, myid - (myid % multiple), 88, MPI_COMM_WORLD, &status);
				myflat = otherflat;             // replace individual flatness by merged
			}
		}

		if (merge_crit == 1 && merge_hists == 1)                   // check flatness of other walkers in window
		{
			if (myid % multiple == 0)             // 'root' in energy window, receive individual flatnesses
			{
				if (myflat == 1) //Main node multimyflat (checks for the flat process in an energy window
				{
					flatproc = myid; // id of flat process
					ioffset = (myid - (myid % multiple)) - myid; // offset from the main energy window node
					multimyflat = 1;
				}

				for (int i = 1; i < multiple; i++)
				{
					MPI_Recv(&otherflat, 1, MPI_INT, myid + i, 66, MPI_COMM_WORLD, &status);

					if (otherflat == 1) //sets the value of the two variable based on information recieved from the process each one is communicating with
					{
						flatproc = (myid + i);
						ioffset = (myid - (myid % multiple)) - myid;
						multimyflat = 1;
					}

					myflat *= otherflat;        // all have to be '1' to give '1' in the end (manual '&&')
				}
				for (int i = 1; i < multiple; i++)  // and let everybody know
				{
					MPI_Send(&myflat, 1, MPI_INT, myid + i, 88, MPI_COMM_WORLD);
					MPI_Send(&multimyflat, 1, MPI_INT, myid + i, 86, MPI_COMM_WORLD);
					if (multimyflat == 1) // if multimyflat is found to be equal to one flatproc and ioffset 
					{
						MPI_Send(&flatproc, 1, MPI_INT, myid + i, 90, MPI_COMM_WORLD);
						MPI_Send(&ioffset, 1, MPI_INT, myid + i, 92, MPI_COMM_WORLD);
					}
				}
			}
			else                                // send individual flatness and receive merged status variables
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
				myflat = otherflat;             // replace individual flatness by merged
			}


			if (multimyflat == 1)
			{
				if (myid != flatproc) // non flat process recieve merged flat process density of states and the myflat status flagged is called to perform the lnge merging procedures in the wl routine
				{
					stdoutlog = fopen(stdoutlogname, "a");
					fprintf(stdoutlog, "Proc %3i: recievied flat DOS from Proc: %3i \t offset: %3i \t multimyflat: %3i\n", myid, flatproc, ioffset, multimyflat);
					fclose(stdoutlog);
					MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, flatproc, 77, MPI_COMM_WORLD, &status);
					for (int j = 0; j < hist_size; j++) lngE[j] = lngE_buf[j]; // overrides density of states of non flat processes
					myflat = 1;
				}
				else // flat process sends out its density of states values
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


void recombine(double countd)
{

	stdoutlog = fopen(stdoutlogname, "a");

	double init_dos = 382.86414;

	int init_check = 0; int init_check_2 = 0;
	double rec_lnge; double rec_real_lnge;
	int rec_en; int rec_en_other;
	int rec_en_2; int rec_en_other_2;
	double microT_compare;

	if (myid == 0) // 'root' in energy window, receive individual g(E) and send merged g(E)
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
		sprintf(filename, "TestL%iq%i.HE.proc%04i.iter0", L1dim, q, myid);
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
			MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, ii, 76, MPI_COMM_WORLD, &status); // get other dens. of states
			fprintf(stdoutlog, "\nProc %i: Received lngE from Proc. %i\n", myid, ii);

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
					fprintf(stdoutlog, "\n%i\n", i);
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

		sprintf(filename, "Test2L%iq%i.HE.proc%04i.iter0", L1dim, q, myid);
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
				if (real_lngE[i] > 0.5) fprintf(file, "%i\t%f\n", (i + (Eglobalmin)) + (L1dim * L1dim) * 2, real_lngE[i]);
			}
			fclose(file);
		}

		for (int ii = 1; ii < numprocs; ii++)
        {
            MPI_Recv(&tortuosity_buf[0], hist_size, MPI_DOUBLE, ii, 83, MPI_COMM_WORLD, &status);
            MPI_Recv(&rog_buf[0], hist_size, MPI_DOUBLE, ii, 84, MPI_COMM_WORLD, &status);
            MPI_Recv(&visits_buf[0], hist_size, MPI_LONG_LONG, ii, 85, MPI_COMM_WORLD, &status);
            MPI_Recv(&endtoend_buf[0], hist_size, MPI_DOUBLE, ii, 87, MPI_COMM_WORLD, &status);

            for (int i = Eminindex; i <= Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.5)
                {
                    tortuosity[i] += tortuosity_buf[i];
                    rog[i] += rog_buf[i];
                    visits[i] += visits_buf[i];
                    endtoend[i] += endtoend_buf[i];
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
                if (real_lngE[i] > 0.5) fprintf(file, "%i\t%f\t%f\t%f\t%lli\n", i, (float)rog[i] / ((double)chain_length) / ((double)visits[i]), tortuosity[i] / ((double)visits[i]), endtoend[i] / ((double)visits[i]),visits[i]);
            }
            fclose(file);
        }

	}
	else // send individual g(E) and receive merged g(E)
	{
		fprintf(stdoutlog, "ufig");
		MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, 0, 76, MPI_COMM_WORLD);
		fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, 0);

		MPI_Send(&tortuosity[0], hist_size, MPI_DOUBLE, 0, 83, MPI_COMM_WORLD);
		MPI_Send(&rog[0], hist_size, MPI_DOUBLE, 0, 84, MPI_COMM_WORLD);
		MPI_Send(&visits[0], hist_size, MPI_LONG_LONG, 0, 85, MPI_COMM_WORLD);
        MPI_Send(&endtoend[0], hist_size, MPI_LONG_LONG, 0, 87, MPI_COMM_WORLD);
	}
	fclose(stdoutlog);
}

