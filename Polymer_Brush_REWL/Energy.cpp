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
#include "Mov.h"
#include "Energy.h"

int poly_en()
{
	//energy computed for the entire chain length the head and tail are treated differently as the are disconnected (Need to matbe right solo function but no real need currently)
	int en = 0;
	for (int i = 0; i < (chain_length); i++) {
		for (int j = 0; j < numberneighbors; j++) {
			switch (i) {
                    default:
                    if (latticepoint[poly_chain[i]] == 1 && latticepoint[neighbor[numberneighbors * poly_chain[i] + j]] == 1 && (neighbor[numberneighbors * poly_chain[i] + j] != poly_chain[i + 1] && neighbor[numberneighbors * poly_chain[i] + j] != poly_chain[i - 1]))
                        en++;
                    break;
			case 0:
				if (latticepoint[poly_chain[i]] == 1 && latticepoint[neighbor[numberneighbors * poly_chain[i] + j]] == 1 && (neighbor[numberneighbors * poly_chain[i] + j] != poly_chain[i + 1]))
					en++;

				//stdoutlog = fopen(stdoutlogname, "a");
				//fprintf(stdoutlog, "\n");
				//fprintf(stdoutlog, "1  chain_length:%i  preceeding:%i  succeeding:%i    \n",latticepoint[neighbor[numberneighbors*poly_chain[i]+j]], neighbor[numberneighbors*poly_chain[i]+j],poly_chain[i+1]);

				//fprintf(stdoutlog, "\n");
				//fclose(stdoutlog);

				break;
			case 57:
				if (latticepoint[poly_chain[i]] == 1 && latticepoint[neighbor[numberneighbors * poly_chain[i] + j]] == 1 && (neighbor[numberneighbors * poly_chain[i] + j] != poly_chain[i - 1]))
					en++;
				break;
			
			}
		}
	}
	/*
	if(en/2==-38)
	{
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		for (int i = 0; i < chain_length; i++)
		{
			fprintf(stdoutlog, "%i ", poly_chain[i]);
		}
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
	}*/
	return ((en / 2)+Eglobalmin);
}

int pseudo_poly_en()
{
	//energy computed for the entire chain length the head and tail are treated differently as the are disconnected (Need to matbe right solo function but no real need currently)
	int en = 0;
	for (int i = 0; i < (chain_length); i++) {
		for (int j = 0; j < numberneighbors; j++) {
			switch (i) {
			case 0:
				if (latticepoint[wl_pseudo_chain[i]] == 1 && latticepoint[neighbor[numberneighbors * wl_pseudo_chain[i] + j]] == 1 && (neighbor[numberneighbors * wl_pseudo_chain[i] + j] != wl_pseudo_chain[i + 1]))
					en--;

				//stdoutlog = fopen(stdoutlogname, "a");
				//fprintf(stdoutlog, "\n");
				//fprintf(stdoutlog, "1  chain_length:%i  preceeding:%i  succeeding:%i    \n",latticepoint[neighbor[numberneighbors*poly_chain[i]+j]], neighbor[numberneighbors*poly_chain[i]+j],poly_chain[i+1]);

				//fprintf(stdoutlog, "\n");
				//fclose(stdoutlog);

				break;
			case 47:
				if (latticepoint[wl_pseudo_chain[i]] == 1 && latticepoint[neighbor[numberneighbors * wl_pseudo_chain[i] + j]] == 1 && (neighbor[numberneighbors * wl_pseudo_chain[i] + j] != wl_pseudo_chain[i - 1]))
					en--;
				break;
			default:
				if (latticepoint[wl_pseudo_chain[i]] == 1 && latticepoint[neighbor[numberneighbors * wl_pseudo_chain[i] + j]] == 1 && (neighbor[numberneighbors * wl_pseudo_chain[i] + j] != wl_pseudo_chain[i + 1] && neighbor[numberneighbors * wl_pseudo_chain[i] + j] != wl_pseudo_chain[i - 1]))
					en--;
				break;
			}
		}
	}

	//if(en/2==-25)
	//{
	//    stdoutlog = fopen(stdoutlogname, "a");
	 //   fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    for (int i = 0; i < chain_length; i++)
	//    {
	//        fprintf(stdoutlog, "%i ", poly_chain[i]);
	//    }
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "\n");
	//    fclose(stdoutlog);
	//}
	return en / 2;
}



int propose_update(int yyy)    // returns energy change if a spin would be updated
{
	// ! this could be coded better: this function also changes the actual spin !
	// ! if update is not accepted, this has to be undone explicitly !

	poly_mov();   // spin _is_ changed here
	if (MoveProposal > 0)
	{
		elocvor = yyy + Eglobalmin;//pseudo_poly_en();
		elocnach = poly_en();
	}
	else
	{
		return 0;
	}
	return(elocnach - elocvor);
}

