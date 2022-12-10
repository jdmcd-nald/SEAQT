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
#include "Mov_Pivot.h"

int init_pivot() // instances the vector that will hold the pivot indexes
{


	for (int i = 0; i < chain_length; i++)   // sets all points in lattice to default values
	{
		pivot_indexes[3 * i] = -1000;     //x-axis
		pivot_indexes[3 * i + 1] = -1000;         // y-axis
		pivot_indexes[3 * i + 2] = -1000; //z-axis
	}
	return 0;
}

int pivot_rotations(int place)
{
	if (axis == 0) // abot the x-axis
	{
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " double check pivot middle h_t:%i   angleif:%i  angleiff:%f angle:%f  axis:%i   cos:%f   sin:%f pivot_indexes[3 * pivot_increm]:%f pivot_indexes[3 * pivot_increm+1]:%f  pivot_indexes[3 * pivot_increm+2]: %f  give_up? : %f",h_t,angleif,angleiff,angle,axis ,cos_value,sin_value, pivot_indexes[3 * pivot_increm],pivot_indexes[3 * pivot_increm+1],pivot_indexes[3 * pivot_increm+2] ,pivot_indexes[3 * pivot_increm] * sin_value);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		xx = pivot_indexes[3 * place];
		yy = pivot_indexes[3 * place + 1];
		zz = pivot_indexes[3 * place + 2];

		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " axis x    place: %i",place);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " %i  =   %i    +0     +0",xx * 1 + yy * 0 + zz * 0,pivot_indexes[3 * place] * 1);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " %i  =   %i    +%i    +%i",xx * 0 + yy * cos_value + zz * sin_value,pivot_indexes[3 * place] * 0,pivot_indexes[3 * place + 1] * cos_value,pivot_indexes[3 * place + 2] * sin_value);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " %i  =   %i    +%i    +%i",xx * 0 + yy * -sin_value + zz * cos_value,pivot_indexes[3 * place] * 0,pivot_indexes[3 * place + 1] * -sin_value,pivot_indexes[3 * place + 2] * cos_value);
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		pivot_indexes_r[3 * place] = xx;
		pivot_indexes_r[3 * place + 1] = yy * cos_value + zz * sin_value;
		pivot_indexes_r[3 * place + 2] = yy * -sin_value + zz * cos_value;



	}
	if (axis == 1) // about the y-axis
	{
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " double check pivot middle h_t:%i   angleif:%i  angleiff:%f angle:%i  axis:%i   cos:%i   sin:%i pivot_indexes[3 * place]:%i pivot_indexes[3 * place+1]:%i  pivot_indexes[3 * place+2]: %i  give_up? : %i",h_t,angleif,angleiff,angle,axis ,cos_value,sin_value, pivot_indexes[3 * place],pivot_indexes[3 * place+1],pivot_indexes[3 * place+2] ,pivot_indexes[3 * place] * sin_value);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		xx = pivot_indexes[3 * place];
		yy = pivot_indexes[3 * place + 1];
		zz = pivot_indexes[3 * place + 2];
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " axis y    place: %i",place);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "  %i  =   %i    +%i    +%i",xx * cos_value + yy * 0 + zz * -sin_value,pivot_indexes[3 * place] * cos_value,pivot_indexes[3 * place + 1] * 0 ,pivot_indexes[3 * place + 2] * -sin_value);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "  %i  =   %i    +%i    +%i",xx * 0 + yy * 1 + zz * 0,pivot_indexes[3 * place] * 0,pivot_indexes[3 * place + 1] * 1 ,pivot_indexes[3 * place + 2] * 0);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "  %i  =   %i    +%i    +%i",xx * sin_value + (yy * 0) + (zz * cos_value),pivot_indexes[3 * place] * sin_value , pivot_indexes[3 * place + 1] * 0 , pivot_indexes[3 * place + 2] * cos_value);
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		pivot_indexes_r[3 * place] = xx * cos_value + zz * -sin_value;
		pivot_indexes_r[3 * place + 1] = yy;
		pivot_indexes_r[3 * place + 2] = xx * sin_value + (zz * cos_value);




	}
	if (axis == 2) // abot the z-axis
	{
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " double check pivot middle h_t:%i   angleif:%i  angleiff:%i angle:%i  axis:%i   cos:%i   sin:%i pivot_indexes[3 * place]:%i pivot_indexes[3 * place+1]:%i  pivot_indexes[3 * place+2]: %i  give_up? : %i",h_t,angleif,angleiff,angle,axis ,cos_value,sin_value, pivot_indexes[3 * place],pivot_indexes[3 * place+1],pivot_indexes[3 * place+2] ,pivot_indexes[3 * place] * sin_value);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/

		xx = pivot_indexes[3 * place];
		yy = pivot_indexes[3 * place + 1];
		zz = pivot_indexes[3 * place + 2];

		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " axis z    place: %i",place);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog,"  %i  =   %i    +%i    +%i",xx * cos_value + yy * sin_value + zz * 0, pivot_indexes[3 * place] * cos_value,pivot_indexes[3 * place + 1] * sin_value,pivot_indexes[3 * place + 2] * 0);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "  %i  =   %i    +%i    +%i",xx * -sin_value + yy * cos_value + zz * 0,pivot_indexes[3 * place] * -sin_value , pivot_indexes[3 * place + 1] * cos_value,pivot_indexes[3 * place + 2] * 0);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog,"  %i  =   %i    +%i    +%i",xx * 0 + yy * 0 + zz * 1,pivot_indexes[3 * place] * 0,pivot_indexes[3 * place + 1] * 0, pivot_indexes[3 * place + 2] * 1);
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		pivot_indexes_r[3 * place] = xx * cos_value + yy * sin_value;
		pivot_indexes_r[3 * place + 1] = xx * -sin_value + yy * cos_value;
		pivot_indexes_r[3 * place + 2] = zz;




	}
	return 0;
}

int pivot_assignment(int place, int chain_y, int chain_pivot)
{
	int complete = 0;
	if (pivot_indexes_r[3 * place] - pivot_indexes_r[3 * (place - 1)] == 1) // abot the x-axis
	{

		if (latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 1] + 4] > pivotflag)
		{
			pivotflag = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 1] + 4];
		}

		if (latticepoint[neighbor[numberneighbors * poly_chain[chain_y - 1] + 1]] > 0 && latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 1] + 4] <= chain_pivot)
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " fail 5");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			 */
			 /*
			 stdoutlog = fopen(stdoutlogname, "a");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "same check");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "\n");
			 fclose(stdoutlog);
			 */
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, chain_pivot - 1, pivotflag + 1);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, chain_pivot - 1, pivotflag + 1);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			complete = 0;
			for (int ii = 0; ii < chain_length; ii++)
			{
				poly_chain[ii] = pseudo_poly_chain[ii];
			}
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "pivot indeces\n");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 1]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 2]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			return complete;

		}
		poly_chain[chain_y] = neighbor[numberneighbors * poly_chain[chain_y - 1] + 1];
		//seq_plhdr = latticepolymer[5*poly_chain[chain_y]+4];
		//latticepolymer[5*poly_chain[chain_y]+5]=1;
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, " 1 %i: (%i %i %i %i %i %i) %i ",poly_chain[chain_y - 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 0],neighbor[numberneighbors * poly_chain[chain_y - 1] + 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 2],neighbor[numberneighbors * poly_chain[chain_y - 1] + 3],neighbor[numberneighbors * poly_chain[chain_y - 1] + 4],neighbor[numberneighbors * poly_chain[chain_y - 1] + 5],poly_chain[chain_y]);
		fclose(stdoutlog);
		*/
		latticepoint[poly_chain[chain_y]] = chain_sequence[chain_y];
		latticepolymer[5 * poly_chain[chain_y] + 4] = chain_y;
	}
	if (pivot_indexes_r[3 * place] - pivot_indexes_r[3 * (place - 1)] == -1) // abot the x-axis
	{

		if (latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 4] + 4] > pivotflag)
		{
			pivotflag = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 4] + 4];
		}

		if (latticepoint[neighbor[numberneighbors * poly_chain[chain_y - 1] + 4]] > 0 && latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 4] + 4] <= chain_pivot)
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " fail 5");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			 */
			 /*
			 stdoutlog = fopen(stdoutlogname, "a");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "same check");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "\n");
			 fclose(stdoutlog);
			 */
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, chain_pivot - 1, pivotflag + 1);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, chain_pivot - 1, pivotflag + 1);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			complete = 0;
			for (int ii = 0; ii < chain_length; ii++)
			{
				poly_chain[ii] = pseudo_poly_chain[ii];
			}
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "pivot indeces\n");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 1]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 2]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			return complete;
		}
		poly_chain[chain_y] = neighbor[numberneighbors * poly_chain[chain_y - 1] + 4];
		//seq_plhdr = latticepolymer[5*poly_chain[chain_y]+4];
		//latticepolymer[5*poly_chain[chain_y]+5]=1;
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, " 4 %i: (%i %i %i %i %i %i) %i ",poly_chain[chain_y - 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 0],neighbor[numberneighbors * poly_chain[chain_y - 1] + 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 2],neighbor[numberneighbors * poly_chain[chain_y - 1] + 3],neighbor[numberneighbors * poly_chain[chain_y - 1] + 4],neighbor[numberneighbors * poly_chain[chain_y - 1] + 5],poly_chain[chain_y]);
		fclose(stdoutlog);
		*/
		latticepoint[poly_chain[chain_y]] = chain_sequence[chain_y];
		latticepolymer[5 * poly_chain[chain_y] + 4] = chain_y;
	}

	if (pivot_indexes_r[3 * place + 1] - pivot_indexes_r[3 * (place - 1) + 1] == 1) // abot the y-axis
	{

		if (latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 0] + 4] > pivotflag)
		{
			pivotflag = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 0] + 4];
		}

		if (latticepoint[neighbor[numberneighbors * poly_chain[chain_y - 1] + 0]] > 0 && latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 0] + 4] <= chain_pivot)
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " fail 5");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			 */
			 /*
			 stdoutlog = fopen(stdoutlogname, "a");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "same check");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "\n");
			 fclose(stdoutlog);
			 */
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, chain_pivot - 1, pivotflag + 1);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, chain_pivot - 1, pivotflag + 1);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			complete = 0;
			for (int ii = 0; ii < chain_length; ii++)
			{
				poly_chain[ii] = pseudo_poly_chain[ii];
			}
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "pivot indeces\n");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 1]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 2]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			return complete;
		}
		poly_chain[chain_y] = neighbor[numberneighbors * poly_chain[chain_y - 1] + 0];
		//seq_plhdr = latticepolymer[5*poly_chain[chain_y]+4];
		//latticepolymer[5*poly_chain[chain_y]+5]=1;
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, " 0 %i: (%i %i %i %i %i %i) %i ",poly_chain[chain_y - 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 0],neighbor[numberneighbors * poly_chain[chain_y - 1] + 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 2],neighbor[numberneighbors * poly_chain[chain_y - 1] + 3],neighbor[numberneighbors * poly_chain[chain_y - 1] + 4],neighbor[numberneighbors * poly_chain[chain_y - 1] + 5],poly_chain[chain_y]);
		fclose(stdoutlog);
		*/
		latticepoint[poly_chain[chain_y]] = chain_sequence[chain_y];
		latticepolymer[5 * poly_chain[chain_y] + 4] = chain_y;
	}
	if (pivot_indexes_r[3 * place + 1] - pivot_indexes_r[3 * (place - 1) + 1] == -1) // abot the y-axis
	{

		if (latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 3] + 4] > pivotflag)
		{
			pivotflag = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 3] + 4];
		}

		if (latticepoint[neighbor[numberneighbors * poly_chain[chain_y - 1] + 3]] > 0 && latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 3] + 4] <= chain_pivot)
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " fail 5");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			 */
			 /*
			 stdoutlog = fopen(stdoutlogname, "a");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "same check");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "\n");
			 fclose(stdoutlog);
			 */
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, chain_pivot - 1, pivotflag + 1);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, chain_pivot - 1, pivotflag + 1);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			complete = 0;
			for (int ii = 0; ii < chain_length; ii++)
			{
				poly_chain[ii] = pseudo_poly_chain[ii];
			}
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "pivot indeces\n");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 1]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 2]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/


			return complete;
		}

		poly_chain[chain_y] = neighbor[numberneighbors * poly_chain[chain_y - 1] + 3];
		//seq_plhdr = latticepolymer[5*poly_chain[chain_y]+4];
		/*
		latticepolymer[5*poly_chain[chain_y]+5]=1;
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, " 3 %i: (%i %i %i %i %i %i) %i ",poly_chain[chain_y - 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 0],neighbor[numberneighbors * poly_chain[chain_y - 1] + 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 2],neighbor[numberneighbors * poly_chain[chain_y - 1] + 3],neighbor[numberneighbors * poly_chain[chain_y - 1] + 4],neighbor[numberneighbors * poly_chain[chain_y - 1] + 5],poly_chain[chain_y]);
		fclose(stdoutlog);
		*/
		latticepoint[poly_chain[chain_y]] = chain_sequence[chain_y];
		latticepolymer[5 * poly_chain[chain_y] + 4] = chain_y;
	}

	if (pivot_indexes_r[3 * place + 2] - pivot_indexes_r[3 * (place - 1) + 2] == 1) // abot the z-axis
	{

		if (latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 2] + 4] > pivotflag)
		{
			pivotflag = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 2] + 4];
		}

		if (latticepoint[neighbor[numberneighbors * poly_chain[chain_y - 1] + 2]] > 0 && latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 2] + 4] <= chain_pivot)
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " fail 5");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			 */
			 /*
			 stdoutlog = fopen(stdoutlogname, "a");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "same check");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "\n");
			 fclose(stdoutlog);
			 */
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, chain_pivot - 1, pivotflag + 1);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, chain_pivot - 1, pivotflag + 1);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			complete = 0;
			for (int ii = 0; ii < chain_length; ii++)
			{
				poly_chain[ii] = pseudo_poly_chain[ii];
			}
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "pivot indeces\n");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 1]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 2]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			return complete;
		}
		poly_chain[chain_y] = neighbor[numberneighbors * poly_chain[chain_y - 1] + 2];
		//seq_plhdr = latticepolymer[5*poly_chain[chain_y]+4];
		//latticepolymer[5*poly_chain[chain_y]+5]=1;
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, " 2 %i: (%i %i %i %i %i %i) %i ",poly_chain[chain_y - 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 0],neighbor[numberneighbors * poly_chain[chain_y - 1] + 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 2],neighbor[numberneighbors * poly_chain[chain_y - 1] + 3],neighbor[numberneighbors * poly_chain[chain_y - 1] + 4],neighbor[numberneighbors * poly_chain[chain_y - 1] + 5],poly_chain[chain_y]);
		fclose(stdoutlog);
		*/
		latticepoint[poly_chain[chain_y]] = chain_sequence[chain_y];
		latticepolymer[5 * poly_chain[chain_y] + 4] = chain_y;
	}
	if (pivot_indexes_r[3 * place + 2] - pivot_indexes_r[3 * (place - 1) + 2] == -1) // abot the z-axis
	{

		if (latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 5] + 4] > pivotflag)
		{
			pivotflag = latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 5] + 4];
		}

		if (latticepoint[neighbor[numberneighbors * poly_chain[chain_y - 1] + 5]] > 0 && latticepolymer[5 * neighbor[numberneighbors * poly_chain[chain_y - 1] + 5] + 4] <= chain_pivot)
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " fail 5");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			 */
			 /*
			 stdoutlog = fopen(stdoutlogname, "a");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "same check");
			 fprintf(stdoutlog, "\n");
			 for (int i = 0; i < chain_length; i++)
			 {
				 fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			 }
			 fprintf(stdoutlog, "\n");
			 fprintf(stdoutlog, "\n");
			 fclose(stdoutlog);
			 */
			lattice_reorg_non_rebrid(poly_chain, pseudo_poly_chain, 1, chain_pivot - 1, pivotflag + 1);
			lp_reorg_non_rebrid(poly_chain, pseudo_poly_chain, chain_pivot - 1, pivotflag + 1);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);

			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "same check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			complete = 0;
			for (int ii = 0; ii < chain_length; ii++)
			{
				poly_chain[ii] = pseudo_poly_chain[ii];
			}
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "pivot indeces\n");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 1]);
			}
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			for (int j = 0; j < chain_length; j++)
			{
				fprintf(stdoutlog, "%4i",pivot_indexes_r[3 * j + 2]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			return 0;
		}
		poly_chain[chain_y] = neighbor[numberneighbors * poly_chain[chain_y - 1] + 5];
		//seq_plhdr = latticepolymer[5*poly_chain[chain_y]+4];
		//latticepolymer[5*poly_chain[chain_y]+5]=1;
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, " 5 %i: (%i %i %i %i %i %i) %i ",poly_chain[chain_y - 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 0],neighbor[numberneighbors * poly_chain[chain_y - 1] + 1],neighbor[numberneighbors * poly_chain[chain_y - 1] + 2],neighbor[numberneighbors * poly_chain[chain_y - 1] + 3],neighbor[numberneighbors * poly_chain[chain_y - 1] + 4],neighbor[numberneighbors * poly_chain[chain_y - 1] + 5],poly_chain[chain_y]);
		fclose(stdoutlog);
		*/
		latticepoint[poly_chain[chain_y]] = chain_sequence[chain_y];
		latticepolymer[5 * poly_chain[chain_y] + 4] = chain_y;
	}
	return 1;
}

int pivot_extract(int reference_index, int chain_reference, int type) // will extract the relative indexes for the pivot operation reference_index is the relative pivot index and chain reference is the poly chain referecnce
{
	for (int y = 0; y < chain_length; y++)
	{
		pseudo_poly_chain[y] = poly_chain[y];
	}
	/*
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	for (int i = 0; i < chain_length; i++)
	{
		fprintf(stdoutlog, "%4i\t", poly_chain[i]);
	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "pivot extract");
	fprintf(stdoutlog, "\n");
	for (int i = 0; i < chain_length; i++)
	{
		fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);
	*/
	//if(type==0)
	//{
	//    for (int yy = (chain_reference + 1); yy< chain_length; yy++)
	//    {
	//        latticepoint[pseudo_poly_chain[yy]] = 0;
	//    }
	//}

	//int difference;
	int pivot_increm = 0;
	pivot_indexes[3 * 0] = 0;
	pivot_indexes[3 * 0 + 1] = 0;
	pivot_indexes[3 * 0 + 2] = 0;

	pivot_indexes_r[3 * pivot_increm] = 0;
	pivot_indexes_r[3 * pivot_increm + 1] = 0;
	pivot_indexes_r[3 * pivot_increm + 2] = 0;
	int holder;
	int neighbor_holder; // hold the neighbor index for the point from the last
	int iii = 0;
	pivotflag = 1;
	for (int yy = (chain_reference + 1); yy < chain_length; yy++)
	{
		pivot_increm++;
		//difference = (poly_chain[yy]-pseudo_poly_chain[yy - 1]);
		//check_plhdr=latticepoint[pseudo_poly_chain[yy]];
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " pseudo_poly_chain[yy - 1]:%i   chain_reference:%i   difference:%i   pivot_increm:%i  yy:%i",pseudo_poly_chain[yy - 1],poly_chain[yy],difference,pivot_increm,yy);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		neighbor_holder = -1;
		iii = 0;
		while (neighbor_holder == -1)
		{
			if (neighbor[numberneighbors * pseudo_poly_chain[yy - 1] + iii] == poly_chain[yy])
			{
				neighbor_holder = iii;
			}
			iii++;
		}


        
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " neighborholder: %i",neighbor_holder);
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/
		if (yy > pivotflag)
		{
			pivotflag = yy;
		}

        pivot_indexes[3 * 0] = 0;
        pivot_indexes[3 * 0 + 1] = 0;
        pivot_indexes[3 * 0 + 2] = 0;
        pivot_indexes[3 * pivot_increm] = 0;
        pivot_indexes[3 * pivot_increm + 1] = 0;
        pivot_indexes[3 * pivot_increm + 2] = 0;
        
		if (neighbor_holder == 0) //forward neighbor check
		{
			pivot_indexes[3 * pivot_increm] = (pivot_indexes[3 * (pivot_increm - 1)]);
			pivot_indexes[3 * pivot_increm + 1] = (pivot_indexes[3 * (pivot_increm - 1) + 1]) + 1;
			pivot_indexes[3 * pivot_increm + 2] = (pivot_indexes[3 * (pivot_increm - 1) + 2]);


			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 0");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm],(pivot_indexes[3 * (pivot_increm - 1)]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 1",pivot_indexes[3 * pivot_increm+ 1],(pivot_indexes[3 * (pivot_increm - 1)+ 1]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 2],(pivot_indexes[3 * (pivot_increm - 1)+ 2]));
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			if (type == 0)
			{
				pivot_rotations(pivot_increm);
				holder = pivot_assignment(pivot_increm, yy, chain_reference);
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " latticepolymer[5*pseudo_poly_chain[yy]+4]:%i   check_plhdr:%i   um:%i    uh:%i    latticepolymer[5*pseudo_poly_chain[yy]+5]: %i  poly yy: %i",latticepolymer[5*pseudo_poly_chain[yy]+4],latticepoint[pseudo_poly_chain[yy]],latticepolymer[5*poly_chain[yy]+4],seq_plhdr,latticepolymer[5*pseudo_poly_chain[yy]+5],latticepoint[poly_chain[yy]]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				if (holder == 0) {
					return 0;
				}


			}
		}
		if (neighbor_holder == 1) // right neighbor check
		{
			pivot_indexes[3 * pivot_increm] = (pivot_indexes[3 * (pivot_increm - 1)]) + 1;
			pivot_indexes[3 * pivot_increm + 1] = (pivot_indexes[3 * (pivot_increm - 1) + 1]);
			pivot_indexes[3 * pivot_increm + 2] = (pivot_indexes[3 * (pivot_increm - 1) + 2]);

			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 1");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 1",pivot_indexes[3 * pivot_increm],(pivot_indexes[3 * (pivot_increm - 1)]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 1],(pivot_indexes[3 * (pivot_increm - 1)+ 1]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 2],(pivot_indexes[3 * (pivot_increm - 1)+ 2]));
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			if (type == 0)
			{

				//seq_plhdr = latticepolymer[5*pseudo_poly_chain[yy]+4];
				pivot_rotations(pivot_increm);
				holder = pivot_assignment(pivot_increm, yy, chain_reference);
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " latticepolymer[5*pseudo_poly_chain[yy]+4]:%i   check_plhdr:%i   um:%i    uh:%i    latticepolymer[5*pseudo_poly_chain[yy]+5]: %i  poly yy: %i",latticepolymer[5*pseudo_poly_chain[yy]+4],latticepoint[pseudo_poly_chain[yy]],latticepolymer[5*poly_chain[yy]+4],seq_plhdr,latticepolymer[5*pseudo_poly_chain[yy]+5],latticepoint[poly_chain[yy]]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				if (holder == 0) {
					return 0;
				}



			}
		}
		if (neighbor_holder == 2) // above neighbor check
		{
			pivot_indexes[3 * pivot_increm] = (pivot_indexes[3 * (pivot_increm - 1)]);
			pivot_indexes[3 * pivot_increm + 1] = (pivot_indexes[3 * (pivot_increm - 1) + 1]);
			pivot_indexes[3 * pivot_increm + 2] = (pivot_indexes[3 * (pivot_increm - 1) + 2]) + 1;

			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 2");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm],(pivot_indexes[3 * (pivot_increm - 1)]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 1],(pivot_indexes[3 * (pivot_increm - 1)+ 1]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 1",pivot_indexes[3 * pivot_increm+ 2],(pivot_indexes[3 * (pivot_increm - 1)+ 2]));
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			if (type == 0)
			{

				//seq_plhdr = latticepolymer[5*pseudo_poly_chain[yy]+4];
				pivot_rotations(pivot_increm);
				holder = pivot_assignment(pivot_increm, yy, chain_reference);
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " latticepolymer[5*pseudo_poly_chain[yy]+4]:%i   check_plhdr:%i   um:%i    uh:%i    latticepolymer[5*pseudo_poly_chain[yy]+5]: %i  poly yy: %i",latticepolymer[5*pseudo_poly_chain[yy]+4],latticepoint[pseudo_poly_chain[yy]],latticepolymer[5*poly_chain[yy]+4],seq_plhdr,latticepolymer[5*pseudo_poly_chain[yy]+5],latticepoint[poly_chain[yy]]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				if (holder == 0) {
					return 0;
				}


			}
		}
		if (neighbor_holder == 3) //behind neighbor check
		{
			pivot_indexes[3 * pivot_increm] = (pivot_indexes[3 * (pivot_increm - 1)]);
			pivot_indexes[3 * pivot_increm + 1] = (pivot_indexes[3 * (pivot_increm - 1) + 1]) - 1;
			pivot_indexes[3 * pivot_increm + 2] = (pivot_indexes[3 * (pivot_increm - 1) + 2]);

			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 3");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm],(pivot_indexes[3 * (pivot_increm - 1)]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   - 1",pivot_indexes[3 * pivot_increm+ 1],(pivot_indexes[3 * (pivot_increm - 1)+ 1]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 2],(pivot_indexes[3 * (pivot_increm - 1)+ 2]));
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			if (type == 0)
			{

				//seq_plhdr = latticepolymer[5*pseudo_poly_chain[yy]+4];
				pivot_rotations(pivot_increm);
				holder = pivot_assignment(pivot_increm, yy, chain_reference);
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " latticepolymer[5*pseudo_poly_chain[yy]+4]:%i   check_plhdr:%i   um:%i    uh:%i    latticepolymer[5*pseudo_poly_chain[yy]+5]: %i  poly yy: %i",latticepolymer[5*pseudo_poly_chain[yy]+4],latticepoint[pseudo_poly_chain[yy]],latticepolymer[5*poly_chain[yy]+4],seq_plhdr,latticepolymer[5*pseudo_poly_chain[yy]+5],latticepoint[poly_chain[yy]]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				if (holder == 0) {
					return 0;
				}


			}
		}
		if (neighbor_holder == 4) // left neighbor check
		{
			pivot_indexes[3 * pivot_increm] = (pivot_indexes[3 * (pivot_increm - 1)]) - 1;
			pivot_indexes[3 * pivot_increm + 1] = (pivot_indexes[3 * (pivot_increm - 1) + 1]);
			pivot_indexes[3 * pivot_increm + 2] = (pivot_indexes[3 * (pivot_increm - 1) + 2]);

			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 4");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   - 1",pivot_indexes[3 * pivot_increm],(pivot_indexes[3 * (pivot_increm - 1)]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 1],(pivot_indexes[3 * (pivot_increm - 1)+ 1]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 2],(pivot_indexes[3 * (pivot_increm - 1)+ 2]));
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			if (type == 0)
			{

				//seq_plhdr = latticepolymer[5*pseudo_poly_chain[yy]+4];
				pivot_rotations(pivot_increm);

				holder = pivot_assignment(pivot_increm, yy, chain_reference);
				/*stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " latticepolymer[5*pseudo_poly_chain[yy]+4]:%i   check_plhdr:%i   um:%i    uh:%i    latticepolymer[5*pseudo_poly_chain[yy]+5]: %i  poly yy: %i",latticepolymer[5*pseudo_poly_chain[yy]+4],latticepoint[pseudo_poly_chain[yy]],latticepolymer[5*poly_chain[yy]+4],seq_plhdr,latticepolymer[5*pseudo_poly_chain[yy]+5],latticepoint[poly_chain[yy]]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				if (holder == 0) {
					return 0;
				}


			}
		}
		if (neighbor_holder == 5) // below neighbor check
		{
			pivot_indexes[3 * pivot_increm] = (pivot_indexes[3 * (pivot_increm - 1)]);
			pivot_indexes[3 * pivot_increm + 1] = (pivot_indexes[3 * (pivot_increm - 1) + 1]);
			pivot_indexes[3 * pivot_increm + 2] = (pivot_indexes[3 * (pivot_increm - 1) + 2]) - 1;

			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " 5");
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm],(pivot_indexes[3 * (pivot_increm - 1)]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   + 0",pivot_indexes[3 * pivot_increm+ 1],(pivot_indexes[3 * (pivot_increm - 1)+ 1]));
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, " %i  =  %i   - 1",pivot_indexes[3 * pivot_increm+ 2],(pivot_indexes[3 * (pivot_increm - 1)+ 2]));
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			if (type == 0)
			{

				//seq_plhdr = latticepolymer[5*pseudo_poly_chain[yy]+4];
				pivot_rotations(pivot_increm);
				holder = pivot_assignment(pivot_increm, yy, chain_reference);
				/*
				stdoutlog = fopen(stdoutlogname, "a");
				fprintf(stdoutlog, "\n");
				fprintf(stdoutlog, " latticepolymer[5*pseudo_poly_chain[yy]+4]:%i   check_plhdr:%i   um:%i    uh:%i    latticepolymer[5*pseudo_poly_chain[yy]+5]: %i  poly yy: %i",latticepolymer[5*pseudo_poly_chain[yy]+4],latticepoint[pseudo_poly_chain[yy]],latticepolymer[5*poly_chain[yy]+4],seq_plhdr,latticepolymer[5*pseudo_poly_chain[yy]+5],latticepoint[poly_chain[yy]]);
				fprintf(stdoutlog, "\n");
				fclose(stdoutlog);
				*/
				if (holder == 0) {
					return 0;
				}


			}
		}
	}
	return 1;
}


int pivot_movement() // will perform the main pivot operation
{
	int complete = 0; // completion check
	int current_index;
	int i = 1;
	int random_pull = rand() % (chain_length - 1); // starts the indexing at a random location along the chain
	//must be two indexes away from the end as the end is coiled to make the chain pull
	int prev_index = poly_chain[random_pull + 1]; // stores the value one of random pull
	int prev_pull = random_pull + 1;
	random_pull = (random_pull + 1) % (chain_length - 1); // stores the value of 1 + random pull

	int h_t = (rand() % (2)) + 1;
	int angleif = (rand() % (3)) + 1;


	if (angleif == 1) // 90
	{
		cos_value = 0;
		sin_value = 1;
	}
	if (angleif == 2)//180
	{
		cos_value = -1;
		sin_value = 0;
	}
	if (angleif == 3) // 270
	{
		cos_value = 0;
		sin_value = -1;
	}

	axis = (rand() % (3));
	/*
	stdoutlog = fopen(stdoutlogname, "a");
	fprintf(stdoutlog, "\n");
	for (int i = 0; i < chain_length; i++)
	{
		fprintf(stdoutlog, "%4i\t", poly_chain[i]);
	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "pivot start");
	fprintf(stdoutlog, "\n");
	for (int i = 0; i < chain_length; i++)
	{
		fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
	}
	fprintf(stdoutlog, "\n");
	fprintf(stdoutlog, "\n");
	fclose(stdoutlog);
	*/
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "\n");

	//for (int j = 0; j < chain_length; j++)
	//{
	//    fprintf(stdoutlog, "%i     \t",pseudo_poly_chain[j]);

	//}
	//fprintf(stdoutlog, "\n");
	//fclose(stdoutlog);
	//stdoutlog = fopen(stdoutlogname, "a");
	//fprintf(stdoutlog, "\n");
	//fprintf(stdoutlog, " pivot start h_t:%i   angleif:%i  angleiff:%f angle:%f  axis:%i",h_t,angleif,angleiff,angle,axis);
	//fprintf(stdoutlog, "\n");
	//fprintf(stdoutlog, "\n");
	//fclose(stdoutlog);

	while (i < (chain_length - 2)) // cheks
	{
		while (random_pull >= (chain_length - 3) || random_pull <= (2))
		{
			random_pull = rand() % (chain_length);
		}
		//random_pull++;
		random_pull = (random_pull) % (chain_length - 1); // gives chain poistions
		current_index = poly_chain[random_pull];
		//prev_index = poly_chain[(random_pull+2) % (chain_length - 2)]; // stores the value one of random pull
		//prev_pull = (random_pull+2) % (chain_length - 2);

		//if(random_pull<=chain_length-2)
		//{
		/*
		stdoutlog = fopen(stdoutlogname, "a");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, " call pivot extract");
		fprintf(stdoutlog, "\n");
		fprintf(stdoutlog, "\n");
		fclose(stdoutlog);
		*/

		complete = pivot_extract(current_index, random_pull, 0);

		if (complete == 1)
		{
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            fprintf(stdoutlog, "\n");
			//            fprintf(stdoutlog, "1f happened:%i",7);
			//            fprintf(stdoutlog, "\n");
			//            fclose(stdoutlog);
			//            stdoutlog = fopen(stdoutlogname, "a");
			//            for (int i = 0; i < chain_length; i++)
			//            {
			//                fprintf(stdoutlog, "%4i", poly_chain[i]);
			//            }
			//            fclose(stdoutlog);
			//stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\n");

			//for (int j = 0; j < chain_length; j++)
			//{
			//    fprintf(stdoutlog, "%i    \t",pseudo_poly_chain[j]);

			//}
			//fprintf(stdoutlog, "\n");
			//fclose(stdoutlog);
			//stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\ncc   poly_chain[1]: %i  random_pull:%i\n",poly_chain[1],random_pull);
			//fclose(stdoutlog);
			int pop = 1;
			//lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, pop,0,chain_length);
			//stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\ncc   pchain[1]: %i  random_pull:%i\n",poly_chain[1],random_pull);
			//fclose(stdoutlog);

			//stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "\n");

			//for (int j = 0; j < chain_length; j++)
			//{
			//    fprintf(stdoutlog, "%i    \t",pseudo_poly_chain[j]);

			//}
			//fprintf(stdoutlog, "\n");
			//fclose(stdoutlog);

			/*stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "final  b check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "final  b check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			//return complete;
			pivot_loc = random_pull - 1;
			lattice_reorg_non_rebrid(pseudo_poly_chain, poly_chain, 1, random_pull - 1, chain_length);
			lp_reorg_non_rebrid(pseudo_poly_chain, poly_chain, random_pull - 1, chain_length);
			//stdoutlog = fopen(stdoutlogname, "a");
			//fprintf(stdoutlog, "ccc");
			//fclose(stdoutlog);

		//            for (int y = 0; y < chain_length; y++)
		//            {
		//                pseudo_poly_chain[y]=poly_chain[y] ;
		//            }
		//            stdoutlog = fopen(stdoutlogname, "a");
		//            fprintf(stdoutlog, "\n");
		//            fprintf(stdoutlog, "happended: %i  , axis:%i  , random_pull:  %i,  angle:%f",7,axis,random_pull,angle);
		//            fprintf(stdoutlog, "\n");
		//            fclose(stdoutlog);
		//            stdoutlog = fopen(stdoutlogname, "a");
		//            for (int i = 0; i < chain_length; i++)
		//            {
		//                fprintf(stdoutlog, "%4i", poly_chain[i]);
		//            }
		//            fclose(stdoutlog);
		//            stdoutlog = fopen(stdoutlogname, "a");
		//            fprintf(stdoutlog, "\n");
		//            fprintf(stdoutlog, "pivot indeces\n");
		//            for (int j = 0; j < chain_length; j++)
		//            {
		//                fprintf(stdoutlog, "%4i",pivot_indexes[3 * j]);
		//
		//            }
		//            fprintf(stdoutlog, "\n");
		//            fclose(stdoutlog);
		//            stdoutlog = fopen(stdoutlogname, "a");
		//            for (int j = 0; j < chain_length; j++)
		//            {
		//                fprintf(stdoutlog, "%4i",pivot_indexes[3 * j + 1]);
		//
		//            }
		//            fprintf(stdoutlog, "\n");
		//            fclose(stdoutlog);
		//            stdoutlog = fopen(stdoutlogname, "a");
		//            for (int j = 0; j < chain_length; j++)
		//            {
		//                fprintf(stdoutlog, "%4i",pivot_indexes[3 * j + 2]);
		//
		//            }
		//            fprintf(stdoutlog, "\n");
		//            fprintf(stdoutlog, "\n");
		//            fclose(stdoutlog);

					//stdoutlog = fopen(stdoutlogname, "a");
					//fprintf(stdoutlog, "c");
					//fclose(stdoutlog);
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "final check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "final check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/
			return complete;
		}
		else
		{
			/*
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "final check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", pseudo_poly_chain[i]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			stdoutlog = fopen(stdoutlogname, "a");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "final check");
			fprintf(stdoutlog, "\n");
			for (int i = 0; i < chain_length; i++)
			{
				fprintf(stdoutlog, "%4i", latticepoint[pseudo_poly_chain[i]]);
			}
			fprintf(stdoutlog, "\n");
			fprintf(stdoutlog, "\n");
			fclose(stdoutlog);
			*/

			return 0;
		}
		//}
		for (int y = 0; y < chain_length; y++)
		{
			latticepoint[poly_chain[y] = pseudo_poly_chain[y]] = 0;
			latticepoint[poly_chain[y] = pseudo_poly_chain[y]] = chain_sequence[y];
		}
		return 0;
	}
	for (int y = 0; y < chain_length; y++)
	{
		poly_chain[y] = pseudo_poly_chain[y];
	}
	//    stdoutlog = fopen(stdoutlogname, "a");
	//    fprintf(stdoutlog, "\n");
	//    fprintf(stdoutlog, "unknown happended: %i",7);
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

