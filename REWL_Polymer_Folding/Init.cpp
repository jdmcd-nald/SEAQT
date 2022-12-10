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
#include "Init.h"

void init_poly_chain() // initialize polymer chain indexes and initial values
{
	poly_chain = (int*)malloc(chain_length * sizeof(int)); 
	pseudo_poly_chain = (int*)malloc(chain_length * sizeof(int));
	doublepseudo_poly_chain = (int*)malloc(chain_length * sizeof(int));
	pull_poly_chain = (int*)malloc(chain_length * sizeof(int));
	wl_pseudo_chain = (int*)malloc(chain_length * sizeof(int));
	chain_sequence = (int*)malloc(chain_length * sizeof(int));
	digital_eye = (int*)malloc(chain_length * sizeof(int));
	pivot_indexes = (int*)malloc(((chain_length * 3) + 1) * sizeof(int));

	pivot_indexes_r = (int*)malloc(((chain_length * 3) + 1) * sizeof(int));

	centeromass = (double*)malloc(3 + 1 * sizeof(double));
	tortuosity = (double*)malloc(hist_size * sizeof(double));
	tortuosity_buf = (double*)malloc(hist_size * sizeof(double));
	rog = (double*)malloc(hist_size * sizeof(double));
	rog_buf = (double*)malloc(hist_size * sizeof(double));
	visits = (long long*)malloc(hist_size * sizeof(long long));
	visits_buf = (long long*)malloc(hist_size * sizeof(long long));
	s_avg = (double*)malloc(3 + 1 * sizeof(double));
	s_vector = (int*)malloc(((chain_length * 3) + 1) * sizeof(int));

    endtoend = (double*)malloc(hist_size * sizeof(double));
    endtoend_buf = (double*)malloc(hist_size * sizeof(double));
    
s_vector1 = (int*)malloc(((chain_length * 3) + 1) * sizeof(int));
s_vector2 = (int*)malloc(((chain_length * 3) + 1) * sizeof(int));

	configuration_search_r = (int*)malloc(20000000 * sizeof(int));
	configuration_search_t = (int*)malloc(20000000 * sizeof(int));

	for (int y = 0; y < 20000000; y++)
	{
		configuration_search_r[y] = 0;
	}

	for (int y = 0; y < hist_size; y++)
	{
		visits[y] = 0;
		visits_buf[y] = 0;
		rog[y] = 0.0;
		rog_buf[y] = 0.0;
		tortuosity[y] = 0.0;
		tortuosity_buf[y] = 0.0;
        endtoend[y] = 0;
        endtoend_buf[y] = 0;
	}

	//Listing of polymer chain MUST BE CONNECTED IN 3D SPACE indexes are basded on vector neighbor listing
	// From Top Down
	//Ones place controls left and right movement 0 farthest left L1dim-1 farthest right
	//Tens place controls vertical movement 0 upmost and L1dim-1 farthest down
	//Hundreds place controls depth layer 0 the highest layer and L1dim the lowest layer

	//14 mer .1
	/*
	poly_chain[0] = 1;
	poly_chain[1] = 2;
	poly_chain[2] = 3;
	poly_chain[3] = 4;
	poly_chain[4] = 5;
	poly_chain[5] = 6;
	poly_chain[6] = 7;
	poly_chain[7] = 8;
	poly_chain[8] = 9;
	poly_chain[9] = 10;
	poly_chain[10] = 11;
	poly_chain[11] = 12;
	poly_chain[12] = 13;
	poly_chain[13] = 14;
	chain_sequence[0] = 1;
	chain_sequence[1] = 2;
	chain_sequence[2] = 1;
	chain_sequence[3] = 2;
	chain_sequence[4] = 1;
	chain_sequence[5] = 1;
	chain_sequence[6] = 2;
	chain_sequence[7] = 1;
	chain_sequence[8] = 2;
	chain_sequence[9] = 1;
	chain_sequence[10] =1;
	chain_sequence[11] =2;
	chain_sequence[12] =2;
	chain_sequence[13] = 1;
     */


	//67 mer
		/*
		poly_chain[0]=112553;
		poly_chain[1]=112622;
		poly_chain[2]=112691;
		poly_chain[3]=112692;
		poly_chain[4]=112693;
		poly_chain[5]=117454;
		poly_chain[6]=117523;
		poly_chain[7]=112762;
		poly_chain[8]=112761;
		poly_chain[9]=108000;
		poly_chain[10]=103239;
		poly_chain[11]=103240;
		poly_chain[12]=108001;
		poly_chain[13]=107932;
		poly_chain[14]=107931;
		poly_chain[15]=103170;
		poly_chain[16]=103101;
		poly_chain[17]=103102;
		poly_chain[18]=107863;
		poly_chain[19]=107864;
		poly_chain[20]=107933;
		poly_chain[21]=112694;
		poly_chain[22]=112625;
		poly_chain[23]=112624;
		poly_chain[24]=117385;
		poly_chain[25]=117386;
		poly_chain[26]=117455;
		poly_chain[27]=117456;
		poly_chain[28]=117525;
		poly_chain[29]=112764;
		poly_chain[30]=112695;
		poly_chain[31]=107934;
		poly_chain[32]=107935;
		poly_chain[33]=108004;
		poly_chain[34]=108003;
		poly_chain[35]=108002;
		poly_chain[36]=108071;
		poly_chain[37]=108070;
		poly_chain[38]=112831;
		poly_chain[39]=112832;
		poly_chain[40]=112763;
		poly_chain[41]=117524;
		poly_chain[42]=117593;
		poly_chain[43]=117592;
		poly_chain[44]=117661;
		poly_chain[45]=117660;
		poly_chain[46]=117591;
		poly_chain[47]=122352;
		poly_chain[48]=122353;
		poly_chain[49]=127114;
		poly_chain[50]=127115;
		poly_chain[51]=127046;
		poly_chain[52]=122285;
		poly_chain[53]=122286;
		poly_chain[54]=122217;
		poly_chain[55]=122216;
		poly_chain[56]=122215;
		poly_chain[57]=122284;
		poly_chain[58]=122283;
		poly_chain[59]=117522;
		poly_chain[60]=117453;
		poly_chain[61]=117452;
		poly_chain[62]=117383;
		poly_chain[63]=117384;
		poly_chain[64]=112623;
		poly_chain[65]=107862;
		poly_chain[66]=107861;
		*/

		//64 mer
		/*
		poly_chain[0] = 38714;
		poly_chain[1] =  33953;
		poly_chain[2] =  34022;
		poly_chain[3] = 34023;
		poly_chain[4] = 33954;
		poly_chain[5] =  29193;
		poly_chain[6] =  29262;
		poly_chain[7] = 29261;
		poly_chain[8] = 29192;
		poly_chain[9] = 29123;
		poly_chain[10] = 29054;
		poly_chain[11] =  24293;
		poly_chain[12] =  24362;
		poly_chain[13] =  24361;
		poly_chain[14] =  24430;
		poly_chain[15] = 29191;
		poly_chain[16] =  29260;
		poly_chain[17] =  24499;
		poly_chain[18] =  24498;
		poly_chain[19] = 24567;
		poly_chain[20] =  29328;
		poly_chain[21] = 29329;
		poly_chain[22] = 24568;
		poly_chain[23] =  24569;
		poly_chain[24] =  29330;
		poly_chain[25] = 29331;
		poly_chain[26] =  24570;
		poly_chain[27] =  24571;
		poly_chain[28] = 19810;
		poly_chain[29] = 19741;
		poly_chain[30] =  19740;
		poly_chain[31] =  14979;
		poly_chain[32] = 14910;
		poly_chain[33] = 19671;
		poly_chain[34] = 19672;
		poly_chain[35] = 24433;
		poly_chain[36] = 24502;
		poly_chain[37] = 24501;
		poly_chain[38] = 24432;
		poly_chain[39] = 24431;
		poly_chain[40] = 24500;
		poly_chain[41] = 19739;
		poly_chain[42] = 14978;
		poly_chain[43] = 15047;
		poly_chain[44] = 19808;
		poly_chain[45] = 19877;
		poly_chain[46] = 19876;
		poly_chain[47] = 24637;
		poly_chain[48] = 24636;
		poly_chain[49] = 19875;
		poly_chain[50] = 19806;
		poly_chain[51] = 19807;
		poly_chain[52] = 15046;
		poly_chain[53] = 14977;
		poly_chain[54] = 19738;
		poly_chain[55] = 19669;
		poly_chain[56] = 19670;
		poly_chain[57] = 19601;
		poly_chain[58] = 19602;
		poly_chain[59] = 19603;
		poly_chain[60] = 24364;
		poly_chain[61] = 24363;
		poly_chain[62] = 29124;
		poly_chain[63] = 29125;
		*/
    
	poly_chain[0] = 181794;
	poly_chain[1] = 181795;
	poly_chain[2] = 178195;
	poly_chain[3] = 178255;
	poly_chain[4] = 181855;
	poly_chain[5] = 181856;
	poly_chain[6] = 181857;
	poly_chain[7] = 181797;
	poly_chain[8] = 181796;
	poly_chain[9] = 185396;
	poly_chain[10] = 185397;
	poly_chain[11] = 188997;
	poly_chain[12] = 188996;
	poly_chain[13] = 192596;
	poly_chain[14] = 192656;
	poly_chain[15] = 189056;
	poly_chain[16] = 189116;
	poly_chain[17] = 189115;
	poly_chain[18] = 185515;
	poly_chain[19] = 185516;
	poly_chain[20] = 185456;
	poly_chain[21] = 185455;
	poly_chain[22] = 189055;
	poly_chain[23] = 192655;
	poly_chain[24] = 192595;
	poly_chain[25] = 192594;
	poly_chain[26] = 188994;
	poly_chain[27] = 185394;
	poly_chain[28] = 185395;
	poly_chain[29] = 188995;
	poly_chain[30] = 188935;
	poly_chain[31] = 188936;
	poly_chain[32] = 185336;
	poly_chain[33] = 185335;
	poly_chain[34] = 181735;
	poly_chain[35] = 181675;
	poly_chain[36] = 181676;
	poly_chain[37] = 181736;
	poly_chain[38] = 178136;
	poly_chain[39] = 178196;
	poly_chain[40] = 174596;
	poly_chain[41] = 170996;
	poly_chain[42] = 171056;
	poly_chain[43] = 174656;
	poly_chain[44] = 178256;
	poly_chain[45] = 178257;
	poly_chain[46] = 174657;
	poly_chain[47] = 174717;
	poly_chain[48] = 174716;
	poly_chain[49] = 174715;
	poly_chain[50] = 178315;
	poly_chain[51] = 178316;
	poly_chain[52] = 178317;
	poly_chain[53] = 181917;
	poly_chain[54] = 181916;
	poly_chain[55] = 181976;
	poly_chain[56] = 181975;
	poly_chain[57] = 181915;





	//   poly_chain[48] = 49;
	//   poly_chain[49] = 50;
	//   poly_chain[50] = 51;
	//   poly_chain[51] = 52;
	//   poly_chain[52] = 53;
	//   poly_chain[53] = 54;
	//   poly_chain[54] = 55;
	//   poly_chain[55] = 56;
	//   poly_chain[56] = 57;
	//   poly_chain[57] = 58;
	//   poly_chain[58] = 59;
	//   poly_chain[59] = 60;
	//   poly_chain[60] = 61;
	//   poly_chain[61] = 62;
	//   poly_chain[62] = 63;
	//   poly_chain[63] = 64;
	//   poly_chain[64] = 65;
	//   poly_chain[65] = 66;
	//   poly_chain[66] = 67;


		/*
		chain_sequence[0] =  2;
		chain_sequence[1] =  1;
		chain_sequence[2] =  2;
		chain_sequence[3] =  1;
		chain_sequence[4] =  1;
		chain_sequence[5] =  2;
		chain_sequence[6] =  1;
		chain_sequence[7] =  1;
		chain_sequence[8] =  2;
		chain_sequence[9] =  1;
		chain_sequence[10] = 2;
		chain_sequence[11] = 2;
		chain_sequence[12] = 1;
		chain_sequence[13] = 1;
		chain_sequence[14] = 1;
		chain_sequence[15] = 2;
		chain_sequence[16] = 2;
		chain_sequence[17] = 2;
		chain_sequence[18] = 1;
		chain_sequence[19] = 2;
		chain_sequence[20] = 1;
		chain_sequence[21] = 1;
		chain_sequence[22] = 2;
		chain_sequence[23] = 1;
		chain_sequence[24] = 1;
		chain_sequence[25] = 2;
		chain_sequence[26] = 1;
		chain_sequence[27] = 2;
		chain_sequence[28] = 2;
		chain_sequence[29] = 1;
		chain_sequence[30] = 1;
		chain_sequence[31] = 1;
		chain_sequence[32] = 2;
		chain_sequence[33] = 2;
		chain_sequence[34] = 2;
		chain_sequence[35] = 1;
		chain_sequence[36] = 2;
		chain_sequence[37] = 1;
		chain_sequence[38] = 1;
		chain_sequence[39] = 2;
		chain_sequence[40] = 1;
		chain_sequence[41] = 1;
		chain_sequence[42] = 2;
		chain_sequence[43] = 1;
		chain_sequence[44] = 2;
		chain_sequence[45] = 2;
		chain_sequence[46] = 1;
		chain_sequence[47] = 1;
		chain_sequence[48] = 1;
		chain_sequence[49] = 2;
		chain_sequence[50] = 2;
		chain_sequence[51] = 2;
		chain_sequence[52] = 1;
		chain_sequence[53] = 2;
		chain_sequence[54] = 1;
		chain_sequence[55] = 1;
		chain_sequence[56] = 2;
		chain_sequence[57] = 1;
		chain_sequence[58] = 1;
		chain_sequence[59] = 2;
		chain_sequence[60] = 1;
		chain_sequence[61] = 2;
		chain_sequence[62] = 2;
		chain_sequence[63] = 1;
		chain_sequence[64] = 1;
		chain_sequence[65] = 1;
		chain_sequence[66] = 2;
		*/


	chain_sequence[0] = 1;
	chain_sequence[1] = 2;
	chain_sequence[2] = 1;
	chain_sequence[3] = 1;
	chain_sequence[4] = 2;
	chain_sequence[5] = 2;
	chain_sequence[6] = 1;
	chain_sequence[7] = 1;
	chain_sequence[8] = 1;
	chain_sequence[9] = 1;
	chain_sequence[10] = 2;
	chain_sequence[11] = 1;
	chain_sequence[12] = 1;
	chain_sequence[13] = 1;
	chain_sequence[14] = 2;
	chain_sequence[15] = 2;
	chain_sequence[16] = 1;
	chain_sequence[17] = 1;
	chain_sequence[18] = 2;
	chain_sequence[19] = 2;
	chain_sequence[20] = 1;
	chain_sequence[21] = 2;
	chain_sequence[22] = 1;
	chain_sequence[23] = 1;
	chain_sequence[24] = 1;
	chain_sequence[25] = 2;
	chain_sequence[26] = 1;
	chain_sequence[27] = 2;
	chain_sequence[28] = 1;
	chain_sequence[29] = 1;
	chain_sequence[30] = 2;
	chain_sequence[31] = 2;
	chain_sequence[32] = 1;
	chain_sequence[33] = 1;
	chain_sequence[34] = 2;
	chain_sequence[35] = 2;
	chain_sequence[36] = 2;
	chain_sequence[37] = 1;
	chain_sequence[38] = 2;
	chain_sequence[39] = 2;
	chain_sequence[40] = 2;
	chain_sequence[41] = 2;
	chain_sequence[42] = 2;
	chain_sequence[43] = 2;
	chain_sequence[44] = 2;
	chain_sequence[45] = 2;
	chain_sequence[46] = 1;
	chain_sequence[47] = 1;
	chain_sequence[48] = 2;
	chain_sequence[49] = 2;
	chain_sequence[50] = 2;
	chain_sequence[51] = 2;
	chain_sequence[52] = 2;
	chain_sequence[53] = 2;
	chain_sequence[54] = 2;
	chain_sequence[55] = 2;
	chain_sequence[56] = 1;
	chain_sequence[57] = 1;



	chain_sequence[0] = 2;
	chain_sequence[1] = 1;
	chain_sequence[2] = 2;
	chain_sequence[3] = 1;
	chain_sequence[4] = 1;
	chain_sequence[5] = 1;
	chain_sequence[6] = 2;
	chain_sequence[7] = 1;
	chain_sequence[8] = 1;
	chain_sequence[9] = 1;
	chain_sequence[10] = 2;
	chain_sequence[11] = 2;
	chain_sequence[12] = 1;
	chain_sequence[13] = 1;
	chain_sequence[14] = 2;
	chain_sequence[15] = 1;
	chain_sequence[16] = 2;
	chain_sequence[17] = 1;
	chain_sequence[18] = 1;
	chain_sequence[19] = 2;
	chain_sequence[20] = 1;
	chain_sequence[21] = 1;
	chain_sequence[22] = 1;
	chain_sequence[23] = 2;
	chain_sequence[24] = 1;
	chain_sequence[25] = 2;
	chain_sequence[26] = 1;
	chain_sequence[27] = 2;
	chain_sequence[28] = 1;
	chain_sequence[29] = 1;
	chain_sequence[30] = 2;
	chain_sequence[31] = 2;
	chain_sequence[32] = 1;
	chain_sequence[33] = 1;
	chain_sequence[34] = 1;
	chain_sequence[35] = 2;
	chain_sequence[36] = 2;
	chain_sequence[37] = 1;
	chain_sequence[38] = 2;
	chain_sequence[39] = 1;
	chain_sequence[40] = 2;
	chain_sequence[41] = 2;
	chain_sequence[42] = 2;
	chain_sequence[43] = 2;
	chain_sequence[44] = 1;
	chain_sequence[45] = 2;
	chain_sequence[46] = 2;
	chain_sequence[47] = 1;
	chain_sequence[48] = 2;
	chain_sequence[49] = 2;
	chain_sequence[50] = 1;
	chain_sequence[51] = 1;
	chain_sequence[52] = 2;
	chain_sequence[53] = 2;
	chain_sequence[54] = 1;
	chain_sequence[55] = 2;
	chain_sequence[56] = 2;
	chain_sequence[57] = 1;



	chain_sequence[0] = 2;
	chain_sequence[1] = 1;
	chain_sequence[2] = 2;
	chain_sequence[3] = 1;
	chain_sequence[4] = 1;
	chain_sequence[5] = 1;
	chain_sequence[6] = 2;
	chain_sequence[7] = 1;
	chain_sequence[8] = 1;
	chain_sequence[9] = 1;
	chain_sequence[10] = 2;
	chain_sequence[11] = 2;
	chain_sequence[12] = 1;
	chain_sequence[13] = 1;
	chain_sequence[14] = 2;
	chain_sequence[15] = 1;
	chain_sequence[16] = 2;
	chain_sequence[17] = 1;
	chain_sequence[18] = 1;
	chain_sequence[19] = 2;
	chain_sequence[20] = 1;
	chain_sequence[21] = 1;
	chain_sequence[22] = 1;
	chain_sequence[23] = 2;
	chain_sequence[24] = 1;
	chain_sequence[25] = 2;
	chain_sequence[26] = 1;
	chain_sequence[27] = 2;
	chain_sequence[28] = 1;
	chain_sequence[29] = 1;
	chain_sequence[30] = 2;
	chain_sequence[31] = 2;
	chain_sequence[32] = 1;
	chain_sequence[33] = 1;
	chain_sequence[34] = 1;
	chain_sequence[35] = 2;
	chain_sequence[36] = 2;
	chain_sequence[37] = 1;
	chain_sequence[38] = 2;
	chain_sequence[39] = 1;
	chain_sequence[40] = 2;
	chain_sequence[41] = 2;
	chain_sequence[42] = 2;
	chain_sequence[43] = 2;
	chain_sequence[44] = 1;
	chain_sequence[45] = 2;
	chain_sequence[46] = 2;
	chain_sequence[47] = 1;
	chain_sequence[48] = 2;
	chain_sequence[49] = 2;
	chain_sequence[50] = 1;
	chain_sequence[51] = 1;
	chain_sequence[52] = 2;
	chain_sequence[53] = 2;
	chain_sequence[54] = 1;
	chain_sequence[55] = 2;
	chain_sequence[56] = 2;
	chain_sequence[57] = 1;

	//    chain_sequence[48] = 2;
	//    chain_sequence[49] = 1;
	//    chain_sequence[50] = 2;
	//    chain_sequence[51] = 1;
	//    chain_sequence[52] = 2;
	//    chain_sequence[53] = 2;
	//    chain_sequence[54] = 1;
	//    chain_sequence[55] = 1;
	//    chain_sequence[56] = 1;
	//    chain_sequence[57] = 2;
	//    chain_sequence[58] = 1;
	//    chain_sequence[59] = 1;
	//    chain_sequence[60] = 2;
	//    chain_sequence[61] = 1;
	//    chain_sequence[62] = 1;
	//    chain_sequence[63] = 2;


	//    chain_sequence[14] = 2;
	//    chain_sequence[15] = 2;
	//    chain_sequence[16] = 2;
	//    chain_sequence[17] = 2;
	//    chain_sequence[18] = 1;
	//    chain_sequence[19] = 1;


		/*poly_chain[0]=91;
		poly_chain[1]=92;
		poly_chain[2]=93;
		poly_chain[3]=94;
		poly_chain[4]=95;
		poly_chain[5]=96;
		poly_chain[6]=97;
		poly_chain[7]=98;
		poly_chain[8]=99;
		poly_chain[9]=9;
		poly_chain[10]=109;
		poly_chain[11]=100;
		poly_chain[12]=200;
		poly_chain[13]=210;
		poly_chain[14]=219;
		poly_chain[15]=218;
		poly_chain[16]=217;
		poly_chain[17]=216;
		poly_chain[18]=215;
		poly_chain[19]=214;*/

}



void init_neighbors() // create neighbor list  // for bc type on made
{
	// neighbor contains the index of the neighboring spin
	// for each spin there are four neighbors in this order: above, right, below, left
	neighbor = (int*)malloc(((numberspins * numberneighbors) + 1) * sizeof(int));

	for (int i = 0; i < numberspins; i++)    // in general
	{
		neighbor[6 * i] = i - L1dim;     // up
		neighbor[6 * i + 1] = i + 1;         // right
		neighbor[6 * i + 2] = i - L1dim * L1dim; //above
		neighbor[6 * i + 3] = i + L1dim;     // down
		neighbor[6 * i + 4] = i - 1;         // left
		neighbor[6 * i + 5] = i + L1dim * L1dim; //below
	};

	if (bctype == 0)   // periodic BC
		for (int i = 0; i < numberspins; i++)                  // now treat boundaries separately sets neighbors for all boundary points in system
		{
			if (i % (L1dim * L1dim) < L1dim)                                 // top row
				neighbor[6 * i] = i + (L1dim * (L1dim - 1));

			if ((i + 1) % L1dim == 0)                          // rightmost column
				neighbor[6 * i + 1] = i - (L1dim - 1);

			if (i < L1dim * L1dim)                                 // highest plane
				neighbor[6 * i + 2] = i + (L1dim * L1dim * (L1dim - 1));

			if ((i % (L1dim * L1dim)) > (L1dim * (L1dim - 1) - 1))            // bottom row
				neighbor[6 * i + 3] = i - (L1dim * (L1dim - 1));

			if (i % L1dim == 0)                              // leftmost column
				neighbor[6 * i + 4] = i + (L1dim - 1);

			if (i >= L1dim * L1dim * (L1dim - 1))                    // deepest row
				neighbor[6 * i + 5] = i - (L1dim * L1dim * (L1dim - 1));

		};

	if (bctype == 1) // Braskamp Kunz BC UNUSED
		for (int i = 0; i < numberspins; i++)                  // now treat boundaries separately
		{
			if ((i + 1) % L1dim == 0)                          // rightmost column
				neighbor[4 * i + 1] = i + 1 - L1dim;
			if (i % L1dim == 0)                              // leftmost column
				neighbor[4 * i + 3] = i - 1 + L1dim;
			if (i < L1dim)                                 // top row
				neighbor[4 * i] = numberspins;
			if (i > (numberspins - L1dim - 1))            // bottom row
				neighbor[4 * i + 2] = numberspins + (i % 2);
		};

	//print neighbor list (for manual debugging)
	//  for (int i=0; i<numberspins*numberneighbors; i++)
	//    {
	//      printf("%4d", neighbor[i]);
	//      if ((i+1)%numberneighbors == 0) printf("\n");
	//    };
}




int lattice_polymer()
{
	latticepolymer = (int*)malloc(((numberspins * 5) + 1) * sizeof(int));

	int poly_chain_site = 0; //tracks movement along polymer chain;

	for (int i = 0; i < numberspins; i++)   // sets all points in lattice to default values
	{
		latticepolymer[5 * i] = i;     // lattice index
		latticepolymer[5 * i + 1] = -1;         // next lattice point on chain (if not applicable -1 assigned)
		latticepolymer[5 * i + 2] = -1; //previous lattice point on chain (if applicable -1 assigned)
		latticepolymer[5 * i + 3] = -1; //stores a tag almost exclusivly utilized for rebridging to set tags on the created circle so it does not rebridge with itself
		latticepolymer[5 * i + 4] = -1; //stores a tag almost exclusivly utilized for rebridging to set tags on the created circle so it does not rebridge with itself

	}

	for (int i = 0; i < chain_length; i++)
	{
		latticepolymer[5 * poly_chain[i]] = poly_chain[i];
		latticepolymer[5 * poly_chain[i] + 4] = i;


		if (i == 0)
		{
			latticepolymer[5 * poly_chain[i] + 1] = poly_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * poly_chain[i] + 2] = -100; //previous lattice point on chain (if applicable -1 assigned)
		}
		if (i > 0 && i < (chain_length - 1))
		{
			latticepolymer[5 * poly_chain[i] + 1] = poly_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * poly_chain[i] + 2] = poly_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
		}
		if (i == (chain_length - 1))
		{
			latticepolymer[5 * poly_chain[i] + 1] = -200; // next lattice point on chain (if not applicable -1 assigned)
			latticepolymer[5 * poly_chain[i] + 2] = poly_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
		}
	}
	return 0;
}

void init_hists() // initialize histograms
{
	lngE = (double*)malloc(hist_size * sizeof(double));
	lngE_buf = (double*)malloc(hist_size * sizeof(double));
	real_lngE = (double*)malloc(hist_size * sizeof(double));
	real_lngE_buf = (double*)malloc(hist_size * sizeof(double));
	microT = (double*)malloc(hist_size * sizeof(double));
	microT_buf = (double*)malloc(hist_size * sizeof(double));
	pseudolngE = (double*)malloc(hist_size * sizeof(double)); // added for pseudo wl process
	HE = (double*)malloc(hist_size * sizeof(double));

	for (int i = 0; i < hist_size; i++)
	{
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
	for (int i = 0; i < hist_size; i++)
	{
		lngE_buf[i] = 0.0;
		real_lngE[i] = 0.0;
		real_lngE_buf[i] = 0.0;
		microT[i] = 0.0;
		microT_buf[i] = 0.0;
	}
}


