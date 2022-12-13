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
#include "Mov_Pull.h"
#include "Mov_Rebridge.h"
#include "Mov_Other.h"
#include "Mov.h"

void poly_mov()
{
	//r_picker;
	compl_check = 0;
	//int buffer=0;
	//hundred;
	//while(compl_check==0)
	//{
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "\n");
//        fprintf(stdoutlog, "eye check   %i\n",happened3);
//        for (int i = 0; i < chain_length; i++)
//                {
//                    fprintf(stdoutlog, "%4i", poly_chain[i]);
//               }
//         fprintf(stdoutlog, "\n");
//        fprintf(stdoutlog, "\n");
//         fclose(stdoutlog);
		//eye();
	hundred = rand() % 100;
	if (hundred < 75)
	{
		r_picker = (rand() % 6) + 2;
		r_picke = (rand() % chain_length);
		if (r_picke == 0 || r_picke == chain_length - 1)
		{
			r_picker = (rand() % 3) + 6;
		}
		else {
			r_picker = (rand() % 5) + 1;
		}
	}
	if (hundred >= 75 && hundred < 98)
	{
		r_picker = (rand() % 2) + 10;
		r_picke = (rand() % chain_length);
		if (r_picke == 0 || r_picke == chain_length - 1)
		{
			r_picker = 10;
		}
		else {
			r_picker = 11;
		}
	}
	if (hundred >= 98 && hundred < 100)
	{
		r_picker = 9;
	}
	happened = r_picker;
	int i;
	switch (r_picker) {
	case 1:

		compl_check = pull_relaxation(1);


		MoveProposal = compl_check;
		happened3 = 1;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 2:
		compl_check = pull_1(0, 0);

		MoveProposal = compl_check;
		happened3 = 2;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 3:
		compl_check = cornerflip(-1);

		MoveProposal = compl_check;
		happened3 = 3;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 4:
		i = 0;
		compl_check = cornerflip(-1);
		if (compl_check != 0)
		{
			compl_check = doublecornerflip(-1);
			if (compl_check == 0)
			{
				lp_reorg_non_rebrid(poly_chain, doublepseudo_poly_chain, cf_value, cf_value + 3); // rejection needs to change the indexs back to the former storage value
				int pop = 1;
				lattice_reorg_non_rebrid(poly_chain, doublepseudo_poly_chain, 1, cf_value, cf_value + 3);  // makes a copy of the lattice if the change is rejected
			}
		}

		MoveProposal = compl_check;

		happened3 = 4;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 5:

		compl_check = pull_relaxation(2);
		MoveProposal = compl_check;
		happened3 = 5;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 6:
		compl_check = h_t_wiggle();
		MoveProposal = compl_check;
		happened3 = 6;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 7:
		compl_check = slither_f();
		MoveProposal = compl_check;
		happened3 = 7;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 8:
		compl_check = slither_b();
		MoveProposal = compl_check;
		happened3 = 8;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 9:
		compl_check = pivot_movement();
		MoveProposal = compl_check;
		happened3 = 9;

		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//       }
		// fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		// fclose(stdoutlog);

		break;
	case 10:
		compl_check = rebridging_h_t();
		MoveProposal = compl_check;
		happened3 = 10;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 11:

		compl_check = rebridging();
		MoveProposal = compl_check;
		happened3 = 11;

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	default:
		break;
	}
}

void poly_mov_reset()
{
	switch (happened3) {
	case 1:
		if (pull_hort == 2)
		{
			//reorg_wl(pull_loc-1, chain_length);
            reorg_wl(0, chain_length);
		}
		if (pull_hort == 1)
		{
			//reorg_wl (0, pull_loc + 3);
            reorg_wl(0, chain_length);
		}
		if (pull_hort == 3)
		{
			reorg_wl(0, chain_length);
		}
		break;
	case 2:
		reorg_wl(0, chain_length);
		break;
	case 3:
		reorg_wl(cf_value, cf_value+4);
        //    reorg_wl(0, chain_length);
		break;
	case 4:
		//reorg_wl(cf_value, cf_value + 3);
        reorg_wl(0, chain_length);
		break;
	case 5:
		if (pull_relax_hort == 2)
		{
			reorg_wl( 0, pull_relax_loc + 1);
            //reorg_wl(0, chain_length);
		}

		if (pull_relax_hort == 1)
		{
			reorg_wl( pull_relax_loc, chain_length);
            //reorg_wl(0, chain_length);
		}
		break;
	case 6:
		if (wiggle_ht == 1)
		{
			reorg_wl(0, 3);
            //reorg_wl(0, chain_length);
		}

		if (wiggle_ht == 2)
		{
			reorg_wl(chain_length - 3, chain_length);
            //reorg_wl(0, chain_length);
		}
		break;
	case 7:
		reorg_wl(0, chain_length);
        //    reorg_wl(0, chain_length);

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 8:
		reorg_wl(0, chain_length);
        //    reorg_wl(0, chain_length);

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 9:
		reorg_wl(pivot_loc - 2, chain_length);
        //    reorg_wl(0, chain_length);
		//stdoutlog = fopen(stdoutlogname, "a");
		//fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//for (int i = 0; i < chain_length; i++)
		//        {
		//            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//       }
		// fprintf(stdoutlog, "\n");
		//fprintf(stdoutlog, "\n");
		// fclose(stdoutlog);

		break;
	case 10:
		if (re_ht_ht == 1)
		{
			reorg_wl(0, re_ht_si+2);
            //reorg_wl(0, chain_length);
		}

		if (re_ht_ht == 2)
		{
			//reorg_wl(re_ht_si-1, chain_length);
            reorg_wl(0, chain_length);
		}

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	case 11:

		reorg_wl(0, chain_length);

		//                stdoutlog = fopen(stdoutlogname, "a");
		//                fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "eye check  2 %i\n",happened3);
		//                for (int i = 0; i < chain_length; i++)
		//                        {
		//                            fprintf(stdoutlog, "%4i", poly_chain[i]);
		//                       }
		//                 fprintf(stdoutlog, "\n");
		//                fprintf(stdoutlog, "\n");
		//                 fclose(stdoutlog);
						//eye();
		break;
	default:
		break;
	}
}
