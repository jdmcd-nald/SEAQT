#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <random>

#include "Constants.h"
#include "General.h"
#include "Energy.h"

int total_energy()  // returns total energy of system
{
  int e = 0;
  int array[3];

  //    stdoutlog = fopen(stdoutlogname, "a");
  //    fprintf(stdoutlog, "Made it to totalenergy\n");
  //    fclose(stdoutlog);

  for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) {
    if (i % lengthpercore != 0)  // excludes the monomer in the wall
    {
      if (poly_lattice_coordinates[3 * i + 2] ==
          L1dim_S_z - 2)  // checks for distance from wall
      {
        e += 2 * en_array_wall[0][0][2];
      }
      if (poly_lattice_coordinates[3 * i + 2] ==
          L1dim_S_z - 3)  // checks for distance from wall
      {
        e += 2 * en_array_wall[0][0][3];
      }

      if (poly_point_distance(poly_lattice_coordinates, i - 1, i, array) <
          max_dist)  // checks the distance between the current monomer and the
                     // one below for bond energy alculation
      {
        e += 2 * en_array_bond[array[0]][array[1]][array[2]];  // abs are taken
                                                               // because
                                                               // array[x] may
                                                               // not be
                                                               // positive
      }

      // if((i+1)%lengthpercore !=0) // was intended to calculate bond energy in
      // the reverse but does not seem necessary as the energy in the previous
      // calculation has been doubled
      //{
      //    if(poly_point_distance(poly_lattice_coordinates,i,i+1,array)<=max_dist)
      //    {
      //    e+=en_array_bond[abs(array[0])][abs(array[1])][abs(array[2])];
      //    }
      //}

      for (int j = 0; j < numbercores * lengthpercore * pertrusionpercore;
           j++) {
        if (j % lengthpercore != 0) {
          if (poly_point_distance(poly_lattice_coordinates, j, i, array) <=
              max_en_dist)  // checks the distance between the current monomer
                            // and the one below for bond energy alculation
          {
            e += en_array[array[0]][array[1]][array[2]];  // abs are taken
                                                          // because array[x]
                                                          // may not be positive
          }
        }
      }

      for (int j = 0; j < no_solvent_sites;
           j++)  // calculatw=es energy for solvent interactions
      {
        point_distance(poly_lattice_indexes[i], solvent_loc[j], array);

        if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist) {
          e += 2 * en_array[array[0]][array[1]][array[2]];
        }
      }
      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "energy = %i\n",e);
      //            fclose(stdoutlog);
    }
  }
  for (int j = 0; j < no_solvent_sites;
       j++)  // calculatw=es energy for solvent interactions
  {
    for (int k = 0; k < no_solvent_sites;
         k++)  // calculatw=es energy for solvent interactions
    {
      if (j != k) {
        point_distance(solvent_loc[k], solvent_loc[j], array);
        //                stdoutlog = fopen(stdoutlogname, "a");
        //                fprintf(stdoutlog, "\t\t\tarray =
        //                {%i,%i,%i}\n",array[0],array[1],array[2]);
        //                fclose(stdoutlog);
        if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist) {
          e += en_array[array[0]][array[1]][array[2]];

          //                    stdoutlog = fopen(stdoutlogname, "a");
          //                    fprintf(stdoutlog, "\t i=%i jj=%i solventenergy
          //                    = %i   i = %i  p_l_i[i]: %i  s_l[j]:
          //                    %i\n",i,j,e,i,p_l_i[i],s_l[j]);
          //                    fclose(stdoutlog);
          //                    stdoutlog = fopen(stdoutlogname, "a");
          //                                   fprintf(stdoutlog, "\t\tarray =
          //                                   {%i,%i,%i}\n",array[0],array[1],array[2]);
          //                                   fclose(stdoutlog);
        }
      }
    }
  }

  return (e / 2);
}

int total_energy(int* p_l_c, int* p_l_i,
                 int* s_l)  // returns total energy of system
{
  int e = 0;
  int array[3];

  //    stdoutlog = fopen(stdoutlogname, "a");
  //    fprintf(stdoutlog, "Made it to totalenergy\n");
  //    fclose(stdoutlog);
  //
  for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) {
    if (i % lengthpercore != 0)  // excludes the monomer in the wall
    {
      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\t%i  %i\n\n",i,p_l_i[i]);
      //            fclose(stdoutlog);

      int reset_ind = i;

      if (p_l_c[3 * i + 2] == L1dim_S_z - 2)  // checks for distance from wall
      {
        e += 2 * en_array_wall[0][0][2];
      }

      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\tprelimenergy = %i\n",e);
      //            fclose(stdoutlog);

      if (p_l_c[3 * i + 2] == L1dim_S_z - 3)  // checks for distance from wall
      {
        e += 2 * en_array_wall[0][0][3];
      }

      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\tprelimenergy = %i\n",e);
      //            fclose(stdoutlog);

      if (poly_point_distance(p_l_c, i - 1, i, array) <
          max_dist)  // checks the distance between the current monomer and the
                     // one below for bond energy alculation
      {
        // stdoutlog = fopen(stdoutlogname, "a");
        //                fprintf(stdoutlog, "\t\tarray =
        //                {%i,%i,%i}\n",array[0],array[1],array[2]);
        //                fclose(stdoutlog);

        e += 2 * en_array_bond[array[0]][array[1]][array[2]];  // abs are taken
                                                               // because
                                                               // array[x] may
                                                               // not be
                                                               // positive
      }

      //            if((i+1)%lengthpercore !=0) // was intended to calculate
      //            bond energy in the reverse but does not seem necessary as
      //            the energy in the previous calculation has been doubled
      //            {
      //                if(poly_point_distance(poly_lattice_coordinates,i,i+1,array)<=max_dist)
      //                {
      //                e+=en_array_bond[abs(array[0])][abs(array[1])][abs(array[2])];
      //                }
      //            }

      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\tsprelimenergy = %i\n",e);
      //            fclose(stdoutlog);

      for (int j = 0; j < numbercores * lengthpercore * pertrusionpercore;
           j++) {
        if (j % lengthpercore != 0 && i != j) {
          if (poly_point_distance(p_l_c, j, i, array) <=
              max_en_dist)  // checks the distance between the current monomer
                            // and the one below for bond energy alculation
          {
            e += en_array_pp[array[0]][array[1]][array[2]];  // abs are taken
                                                             // because array[x]
                                                             // may not be
                                                             // positive

            //
            //                                        stdoutlog =
            //                                        fopen(stdoutlogname, "a");
            //                                        fprintf(stdoutlog, "\tpp
            //                                        i=%i jj=%i  i = %i
            //                                        p_l_i[i]: %i  s_l[j]:
            //                                        %i\n",i,j,i,p_l_i[i],p_l_i[j]);
            //                                        fclose(stdoutlog);
            //                                        stdoutlog =
            //                                        fopen(stdoutlogname, "a");
            //                                                       fprintf(stdoutlog,
            //                                                       "\t\tarray
            //                                                       =
            //                                                       {%i,%i,%i}\n",array[0],array[1],array[2]);
            //                                                       fclose(stdoutlog);
          }
        }
      }

      for (int j = 0; j < no_solvent_sites;
           j++)  // calculatw=es energy for solvent interactions
      {
        point_distance(p_l_i[i], s_l[j], array);

        //                stdoutlog = fopen(stdoutlogname, "a");
        //                fprintf(stdoutlog, "\t\t\tarray =
        //                {%i,%i,%i}\n",array[0],array[1],array[2]);
        //                fclose(stdoutlog);

        if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist) {
          e += 2 * en_array_ps[array[0]][array[1]][array[2]];

          // if(j==4)
          //{
          ////                    stdoutlog = fopen(stdoutlogname, "a");
          ////                    fprintf(stdoutlog, "\t i=%i jj=%i
          ///solventenergy = %i   i = %i  p_l_i[i]: %i  s_l[j]: %i
          ///%i\n",i,j,e,i,p_l_i[i],s_l[j],2*en_array_ps[array[0]][array[1]][array[2]]);
          ////                    fclose(stdoutlog);
          ////                    stdoutlog = fopen(stdoutlogname, "a");
          ////                                   fprintf(stdoutlog, "\t\tarray =
          ///{%i,%i,%i}\n",array[0],array[1],array[2]);
          ////                                   fclose(stdoutlog);
          // }
        }
      }
      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "energy = %i   i = %i  p_l_i[i]:
      //            %i\n",e,i,p_l_i[i]);
      //            fclose(stdoutlog);
    }
  }

  for (int j = 0; j < no_solvent_sites;
       j++)  // calculatw=es energy for solvent interactions
  {
    if (s_l[j] >=
        ((L1dim_S_z * L1dim_S_xy * L1dim_S_xy) -
         3 * (L1dim_S_xy * L1dim_S_xy)))  // checks for distance from wall
    {
      e += 2 * en_array_wall[0][0][2];
    } else {
      if (s_l[j] >=
          ((L1dim_S_z * L1dim_S_xy * L1dim_S_xy) -
           2 * (L1dim_S_xy * L1dim_S_xy)))  // checks for distance from wall
      {
        e += 2 * en_array_wall[0][0][3];
      }
    }
    //            stdoutlog = fopen(stdoutlogname, "a");
    //            fprintf(stdoutlog, "\tprelimenergy = %i\n",e);
    //            fclose(stdoutlog);

    for (int k = 0; k < no_solvent_sites;
         k++)  // calculatw=es energy for solvent interactions
    {
      point_distance(s_l[k], s_l[j], array);
      //                stdoutlog = fopen(stdoutlogname, "a");
      //                fprintf(stdoutlog, "\t\t\tarray =
      //                {%i,%i,%i}\n",array[0],array[1],array[2]);
      //                fclose(stdoutlog);
      if (dist_array[array[0]][array[1]][array[2]] <= max_en_dist) {
        e += en_array_ss[array[0]][array[1]][array[2]];
        // if(j==4 || k==4)
        //{
        //                        stdoutlog = fopen(stdoutlogname, "a");
        //                        fprintf(stdoutlog, "\t i=%i jj=%i
        //                        solventenergy = %i   i = %i  s_l[i]: %i
        //                        s_l[j]: %i
        //                        %i\n",j,k,e,j,s_l[j],s_l[k],en_array_ss[array[0]][array[1]][array[2]]);
        //                        fclose(stdoutlog);
        //                        stdoutlog = fopen(stdoutlogname, "a");
        //                                       fprintf(stdoutlog, "\t\tarray =
        //                                       {%i,%i,%i}\n",array[0],array[1],array[2]);
        //                                       fclose(stdoutlog);
        //                    }
      }
    }
  }

  return (e / 2);
}

int local_energy(int* p_l_c, int* p_l_i,
                 int* s_l)  // returns energy of a single spin
{
  int eloc = 0;
  int array[3];

  //    stdoutlog = fopen(stdoutlogname, "a");
  //    fprintf(stdoutlog, "\nMade it to local_energy\n");
  //    fprintf(stdoutlog, "chosen_mov:%i \t
  //    poly_solv:%i\n",chosen_mov,poly_solvent_mov);
  //    fprintf(stdoutlog, "reset_indexes[0]+1:%i \t
  //    reset_indexes[1]:%i\n",reset_indexes[0]+1,reset_indexes[1]);
  //    fclose(stdoutlog);

  if (poly_solvent_mov ==
      0)  // for polymwe movements only increments between moved monomer
  {
    for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++) {
      // if(i%lengthpercore !=0 )
      //{
      if (p_l_c[3 * i + 2] == L1dim_S_z - 2) {
        eloc += en_array_wall[0][0][2] + en_array_wall[0][0][2];
      }
      if (p_l_c[3 * i + 2] == L1dim_S_z - 3) {
        eloc += en_array_wall[0][0][3] + en_array_wall[0][0][3];
      }
      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\tprelimenergy = %i\n",eloc);
      //            fclose(stdoutlog);
      //                int x =
      //                reduced_lattice_optimization(p_l_i[i],p_l_i[i-1]);
      //
      //
      //                if(x<int_max_dist) // checks backwards
      //                {// energy calculation is doubled because in
      //                total_energy it is calculated twice for each bond
      //                    if(reset_indexes[0]+1 == i) // check put in place
      //                    {
      //
      //                        stdoutlog = fopen(stdoutlogname, "a");
      //                                    fprintf(stdoutlog, "\t bondenergy =
      //                                    %i x=%i p_l_i[i]: %i  s_l[j]:
      //                                    %i\n",int_en_array_bond[x]+int_en_array_bond[x],x,p_l_i[i],p_l_i[i-1]);
      //                                   fclose(stdoutlog);
      //                    }
      //                    else{
      //                        //eloc+=int_en_array_bond[x];
      //                    }
      //                }
      //                else
      //                {
      //                    stdoutlog = fopen(stdoutlogname, "a");
      //                                           fprintf(stdoutlog,
      //                                           "\tbchecking eye = %i\n x:%i
      //                                           i:%i  j:%i   plc:%i
      //                                           plc:%i\nplc {%i,%i,%i}
      //                                           plc{%i,%i,%i}",eloc,x,i,i-1,p_l_i[i],p_l_i[i-1],p_l_c[3*i],p_l_c[3*i+1],p_l_c[3*i+2],p_l_c[3*(i-1)],p_l_c[3*(i-1)+1],p_l_c[3*(i-1)+2]);
      //                                           fclose(stdoutlog);
      //
      //
      //                                           MPI_Abort(MPI_COMM_WORLD,1);
      //                }
      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\tprelimenergy = %i\n",eloc);
      //            fclose(stdoutlog);
      //                if((i+1)%lengthpercore !=0)
      //                {
      //                    x =
      //                    reduced_lattice_optimization(p_l_i[i],p_l_i[i+1]);
      //
      //                    if(x<int_max_dist)
      //                    {
      //                        if(reset_indexes[1]-1 == i) // check put in
      //                        place
      //                        {
      //                            eloc+=int_en_array_bond[x]+int_en_array_bond[x];
      //                        }
      //                        else{
      //                            eloc+=int_en_array_bond[x];
      //                        }
      //                    }
      //                    else
      //                    {
      //                        stdoutlog = fopen(stdoutlogname, "a");
      //                                               fprintf(stdoutlog,
      //                                               "\tfchecking eye = %i\n
      //                                               x:%i  i:%i  j:%i   plc:%i
      //                                               plc:%i\nplc {%i,%i,%i}
      //                                               plc{%i,%i,%i}",eloc,x,i,(i-1),p_l_i[i],p_l_i[(i-1)],p_l_c[3*i],p_l_c[3*i+1],p_l_c[3*i+2],p_l_c[3*(i-1)],p_l_c[3*(i-1)+1],p_l_c[3*(i-1)+2]);
      //                                               fclose(stdoutlog);
      //
      //
      //                                               MPI_Abort(MPI_COMM_WORLD,1);
      //                    }
      //                }
      //            stdoutlog = fopen(stdoutlogname, "a");
      //            fprintf(stdoutlog, "\tprelimenergy = %i\n",eloc);
      //            fclose(stdoutlog);

      eloc +=
          local_energy_op_2(p_l_i, s_l, i, reset_indexes[0], reset_indexes[1]);
      //                stdoutlog = fopen(stdoutlogname, "a");
      //                fprintf(stdoutlog, "%i 0 localenergy = %i
      //                %i\n",i,eloc,p_l_i[i]);
      //                fclose(stdoutlog);
    }
    //}
  }
  if (poly_solvent_mov == 1) {
    if (s_l[reset_indexes[0] + 1] >=
        ((L1dim_S_z * L1dim_S_xy * L1dim_S_xy) -
         3 * (L1dim_S_xy * L1dim_S_xy)))  // checks for distance from wall
    {
      eloc += 2 * en_array_wall[0][0][2];
    } else {
      if (s_l[reset_indexes[0] + 1] >=
          ((L1dim_S_z * L1dim_S_xy * L1dim_S_xy) -
           2 * (L1dim_S_xy * L1dim_S_xy)))  // checks for distance from wall
      {
        eloc += 2 * en_array_wall[0][0][3];
      }
    }

    eloc += local_energy_op_2(p_l_i, s_l, reset_indexes[0] + 1,
                              reset_indexes[0], reset_indexes[1]);

    //        for (int j = reset_indexes[0]+1; j < reset_indexes[1]; j++) //
    //        calculatw=es energy for solvent interactions
    //                    {
    //
    //                        for (int k = 0; k < no_solvent_sites; k++) //
    //                        calculatw=es energy for solvent interactions
    //                                    {
    //
    //
    //                        point_distance(s_l[k],s_l[j], array);
    ////                        stdoutlog = fopen(stdoutlogname, "a");
    ////                        fprintf(stdoutlog, "\t\t\tarray =
    ///{%i,%i,%i}\n",array[0],array[1],array[2]);
    ////                        fclose(stdoutlog);
    //
    //                        if
    //                        (dist_array[array[0]][array[1]][array[2]]<=max_en_dist)
    //                        {
    //                            //eloc+=2*;
    //
    //                            stdoutlog = fopen(stdoutlogname, "a");
    //                            fprintf(stdoutlog, "\t j=%i k=%i solventenergy
    //                            = %i   s_l[j]: %i  s_l[k]: %i    en:
    //                            %i\n",j,k,en_array[array[0]][array[1]][array[2]],s_l[j],s_l[k],en_array[array[0]][array[1]][array[2]]);
    //                            fclose(stdoutlog);
    //                            stdoutlog = fopen(stdoutlogname, "a");
    //                                           fprintf(stdoutlog, "\t\tarray =
    //                                           {%i,%i,%i}\n",array[0],array[1],array[2]);
    //                                           fclose(stdoutlog);
    //                        }
    //
    //                                        }
    //                        }
  }

  //    stdoutlog = fopen(stdoutlogname, "a");
  //    fprintf(stdoutlog, "\n");
  //    fclose(stdoutlog);

  return (eloc / 2);
}

int local_energy_op(int* p_l_i, int* s_l_i) {
  int eloc = 0;

  // int considered_spin = poly_solvent_mov;
  int i = 2;
  int j = 0;
  int k = 0;
  int c_x;
  int c_y;
  int c_z;

  int phh;
  int phhh;
  int phhhh;
  int phhhhh;

  if (poly_solvent_mov == 0) {
    phh = p_l_i[reset_indexes[0] + 1];
    phhh = p_l_i[(reset_indexes[0] + 1) - 1];
    phhhh = p_l_i[(reset_indexes[0] + 1) + 1];
    phhhhh = (reset_indexes[0] + 1 + 1) % lengthpercore;
    c_x = index_x[p_l_i[reset_indexes[0] + 1]];
    c_y = index_y[p_l_i[reset_indexes[0] + 1]];
    c_z = index_z[p_l_i[reset_indexes[0] + 1]];
  }
  if (poly_solvent_mov == 1) {
    c_x = index_x[s_l_i[reset_indexes[0] + 1]];
    c_y = index_y[s_l_i[reset_indexes[0] + 1]];
    c_z = index_z[s_l_i[reset_indexes[0] + 1]];
  }

  for (int a = 0; a < 2; a++) {
    i = 2;  //+(-4*(a==1))
    if (a == 1) {
      i = -2;
    }

    for (int b = 0; b < 5; b++) {
      j = b;  //+(-4*(b==3))+(-6*(b==4));
      if (b >= 3) {
        j = b - 5;
      }

      for (int c = 0; c < 3; c++) {
        k = c;  //+(-3*(c==2));
        if (c == 2) {
          k = c - 3;
        }

        int x = 0;
        int y = 1;
        int z = 2;

        for (int q = 0; q < 3; q++) {
          x = i * (q == 0) + j * (q == 1) + k * (q == 2);
          y = j * (q == 0) + k * (q == 1) + i * (q == 2);
          z = k * (q == 0) + i * (q == 1) + j * (q == 2);

          int checked_index =
              xyz_index[(z + c_z + 5) * (L1dim_S_xy + 10) * (L1dim_S_xy + 10) +
                        (y + c_y + 5) * (L1dim_S_xy + 10) + (x + c_x + 5)];

          //                    stdoutlog = fopen(stdoutlogname, "a");
          //                                            fprintf(stdoutlog,
          //                                            "\t\tbc %i    %i   %i
          //                                            \n",checked_index,poly_lattice_indexes[(reset_indexes[0]+1)-1],latticepoint_S[checked_index]);
          //                                            fclose(stdoutlog);

          if (poly_solvent_mov == 0) {
            int holdover = latticepoint_S[checked_index];
            if (holdover == 1 && phh != checked_index) {
              //                        int kkp[3]={c_x,c_y,c_z};
              //                        stdoutlog = fopen(stdoutlogname, "a");
              //                        fprintf(stdoutlog, "\t\tpc_index %i %i
              //                        checked index %i latticepoint_S %i
              //                        op_array = {%i,%i,%i}   c={%i,%i,%i}
              //                        xyz={%i,%i,%i} eloc %i    \t
              //                        %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
              //                        ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
              //                        fclose(stdoutlog);
              if (checked_index < ((L1dim_S_z - 1) * L1dim_S_xy * L1dim_S_xy)) {
                eloc += 2 * int_en_array[i * i + j * j + k * k];
              }

              if (checked_index == phhh) {
                eloc += 2 * int_en_array_bond[i * i + j * j + k * k];
                //                        int kkp[3]={c_x,c_y,c_z};
                //                                                stdoutlog =
                //                                                fopen(stdoutlogname,
                //                                                "a");
                //                                                fprintf(stdoutlog,
                //                                                "\t\tbc_index
                //                                                %i %i  checked
                //                                                index %i
                //                                                latticepoint_S
                //                                                %i op_array =
                //                                                {%i,%i,%i}
                //                                                c={%i,%i,%i}
                //                                                xyz={%i,%i,%i}
                //                                                eloc %i    \t
                //                                                %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
                //                                                ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array_bond[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
                //                                                fclose(stdoutlog);
              }
              if (phhhhh != 0 && checked_index == phhhh) {
                eloc += 2 * int_en_array_bond[i * i + j * j + k * k];
              }
            }

            // int kkp[3]={c_x,c_y,c_z};
            if (latticepoint_S[checked_index] == 2) {
              //                        stdoutlog = fopen(stdoutlogname, "a");
              //                        fprintf(stdoutlog, "\t\tc_index %i %i
              //                        checked index %i latticepoint_S %i
              //                        op_array = {%i,%i,%i}   c={%i,%i,%i}
              //                        xyz={%i,%i,%i} eloc %i\t
              //                        %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
              //                        ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
              //                        fclose(stdoutlog);

              eloc += 2 * int_en_array[i * i + j * j + k * k];
            }
          }
          if (poly_solvent_mov == 1) {
            if (latticepoint_S[checked_index] == 1) {
              //                        stdoutlog = fopen(stdoutlogname, "a");
              //                        fprintf(stdoutlog, "\t\tc_index %i %i
              //                        checked index %i latticepoint_S %i
              //                        op_array = {%i,%i,%i}   c={%i,%i,%i}
              //                        xyz={%i,%i,%i} eloc %i\t
              //                        %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
              //                        ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
              //                        fclose(stdoutlog);

              if (checked_index < ((L1dim_S_z - 1) * L1dim_S_xy * L1dim_S_xy)) {
                eloc += 2 * int_en_array[i * i + j * j + k * k];
              }
            }

            // int kkp[3]={c_x,c_y,c_z};
            if (latticepoint_S[checked_index] == 2 &&
                checked_index != ph_solvent_loc[reset_indexes[0] + 1]) {
              //                                                stdoutlog =
              //                                                fopen(stdoutlogname,
              //                                                "a");
              //                                                fprintf(stdoutlog,
              //                                                "\t\tc_index %i
              //                                                %i  checked
              //                                                index %i
              //                                                latticepoint_S
              //                                                %i op_array =
              //                                                {%i,%i,%i}
              //                                                c={%i,%i,%i}
              //                                                xyz={%i,%i,%i}
              //                                                eloc
              //                                                %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
              //                                                ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k]);
              //                                                fclose(stdoutlog);

              eloc += 2 * int_en_array[i * i + j * j + k * k];
            }
          }
        }
      }
    }
  }

  i = 3;
  j = 0;
  k = 0;
  int x = 0;
  int y = 1;
  int z = 2;
  for (int a = 0; a < 2; a++) {
    i = 3;  //+(-4*(a==1))
    if (a == 1) {
      i = -3;
    }
    for (int q = 0; q < 3; q++) {
      x = i * (q == 0) + j * (q == 1) + k * (q == 2);
      y = j * (q == 0) + k * (q == 1) + i * (q == 2);
      z = k * (q == 0) + i * (q == 1) + j * (q == 2);

      int checked_index =
          xyz_index[(z + c_z + 5) * (L1dim_S_xy + 10) * (L1dim_S_xy + 10) +
                    (y + c_y + 5) * (L1dim_S_xy + 10) + (x + c_x + 5)];
      if (poly_solvent_mov == 0) {
        int holdover = latticepoint_S[checked_index];
        if (holdover == 1 && phh != checked_index) {
          //                            int kkp[3]={c_x,c_y,c_z};
          //                            stdoutlog = fopen(stdoutlogname, "a");
          //                            fprintf(stdoutlog, "\t\t3pc_index %i %i
          //                            checked index %i latticepoint_S %i
          //                            op_array = {%i,%i,%i}   c={%i,%i,%i}
          //                            xyz={%i,%i,%i} eloc %i    \t
          //                            %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
          //                            ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
          //                            fclose(stdoutlog);
          if (checked_index < ((L1dim_S_z - 1) * L1dim_S_xy * L1dim_S_xy)) {
            eloc += 2 * int_en_array[i * i + j * j + k * k];
          }

          if (checked_index == phhh) {
            eloc += 2 * int_en_array_bond[i * i + j * j + k * k];
          }
          if ((phhhhh) != 0 && checked_index == phhhh) {
            eloc += 2 * int_en_array_bond[i * i + j * j + k * k];
          }
        }

        if (holdover == 2) {
          //              int kkp[3]={c_x,c_y,c_z};
          //              stdoutlog = fopen(stdoutlogname, "a");
          //              fprintf(stdoutlog, "\t\t3c_index %i %i  checked index
          //              %i latticepoint_S %i op_array = {%i,%i,%i}
          //              c={%i,%i,%i}  xyz={%i,%i,%i} eloc %i    \t
          //              %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
          //              ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
          //              fclose(stdoutlog);

          eloc += 2 * int_en_array[i * i + j * j + k * k];
        }
      }
      if (poly_solvent_mov == 1) {
        if (latticepoint_S[checked_index] == 1) {
          //                        stdoutlog = fopen(stdoutlogname, "a");
          //                        fprintf(stdoutlog, "\t\tc_index %i %i
          //                        checked index %i latticepoint_S %i op_array
          //                        = {%i,%i,%i}   c={%i,%i,%i}  xyz={%i,%i,%i}
          //                        eloc %i\t
          //                        %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
          //                        ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
          //                        fclose(stdoutlog);
          if (checked_index < ((L1dim_S_z - 1) * L1dim_S_xy * L1dim_S_xy)) {
            eloc += 2 * int_en_array[i * i + j * j + k * k];
          }
        }

        // int kkp[3]={c_x,c_y,c_z};
        if (latticepoint_S[checked_index] == 2 &&
            checked_index != ph_solvent_loc[reset_indexes[0] + 1]) {
          //                       stdoutlog = fopen(stdoutlogname, "a");
          //                       fprintf(stdoutlog, "\t\tc_index %i %i
          //                       checked index %i latticepoint_S %i op_array =
          //                       {%i,%i,%i}   c={%i,%i,%i}  xyz={%i,%i,%i}
          //                       eloc %i\t en: %i
          //                       \n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
          //                       ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],2*int_en_array[i*i+j*j+k*k]);
          //                       fclose(stdoutlog);
          //
          eloc += 2 * int_en_array[i * i + j * j + k * k];
        }
      }
    }
  }
  return eloc;
}

void exclusion(int(&arr)[120], int index) {
  int ph = index * numberneighbors_L;
  for (int i = 1; i < numberneighbors_L; i++) {
    arr[what_to_exclude[ph + i]] = 1;
  }
}

int local_energy_op_2(int* p_l_i, int* s_l_i, int reset_ind, int s, int e) {
  int eloc = 0;

  int phh;
  int phhh;
  int phhhh;
  int phhhhh;

  int index_of_concern;

  if (poly_solvent_mov == 0) {
    phh = p_l_i[reset_ind];
    phhh = p_l_i[(reset_ind)-1];
    phhhh = p_l_i[(reset_ind) + 1];
    phhhhh = (reset_ind + 1) % lengthpercore;
    index_of_concern = p_l_i[reset_ind];
  }
  if (poly_solvent_mov == 1) {
    index_of_concern = s_l_i[reset_indexes[0] + 1];
  }
  if (poly_solvent_mov == 1) {
    for (int i = 0; i < 120; i++) {
      //                stdoutlog = fopen(stdoutlogname, "a");
      //        fprintf(stdoutlog, "\t\tijk = {%i,%i,%i}  %i =
      //        local_en_indexes_op[index_of_concern*96+i] index_of_concern
      //        %i\n",x,y,z,checked_index,index_of_concern);
      //        fclose(stdoutlog);

      int checked_index = local_en_indexes_op[index_of_concern * 120 + i];

      int holder = latticepoint_S[checked_index];
      if (holder == 1) {
        //                stdoutlog = fopen(stdoutlogname, "a");
        //            fprintf(stdoutlog, "\t\tijk = {%i,%i,%i}  %i =
        //            local_en_indexes_op[index_of_concern*96+i]
        //            index_of_concern
        //            %i\n",x,y,z,checked_index,index_of_concern);
        //            fclose(stdoutlog);
        int sum = local_en_dist[i];
        //                        stdoutlog = fopen(stdoutlogname, "a");
        //                        fprintf(stdoutlog, "\t\tc_index %i %i  checked
        //                        index %i latticepoint_S %i op_array =
        //                        {%i,%i,%i}   c={%i,%i,%i}  xyz={%i,%i,%i} eloc
        //                        %i\t
        //                        %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
        //                        ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
        //                        fclose(stdoutlog);

        if (checked_index < ((L1dim_S_z - 1) * L1dim_S_xy * L1dim_S_xy)) {
          eloc += 2 * int_en_array_ps[sum];

          //                                                                stdoutlog
          //                                                                =
          //                                                                fopen(stdoutlogname,
          //                                                                "a");
          //                        fprintf(stdoutlog, "p_s\t\t\t
          //                        elocc:%i\n\n",2*int_en_array_ps[sum]);
          //                            fclose(stdoutlog);
        }
        // exclusion(exclude_indexes,i);
      }

      // int kkp[3]={c_x,c_y,c_z};
      if (holder == 2) {
        //                                                          stdoutlog =
        //                                                          fopen(stdoutlogname,
        //                                                          "a");
        //            fprintf(stdoutlog, "\t\tijk = {%i,%i,%i}  %i =
        //            local_en_indexes_op[index_of_concern*96+i]
        //            index_of_concern
        //            %i\n",x,y,z,checked_index,index_of_concern);
        //            fclose(stdoutlog);

        int sum = local_en_dist[i];
        //                                                stdoutlog =
        //                                                fopen(stdoutlogname,
        //                                                "a");
        //                                                fprintf(stdoutlog,
        //                                                "\t\tc_index %i %i
        //                                                checked index %i
        //                                                latticepoint_S %i
        //                                                op_array = {%i,%i,%i}
        //                                                c={%i,%i,%i}
        //                                                xyz={%i,%i,%i} eloc
        //                                                %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
        //                                                ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k]);
        //                                                fclose(stdoutlog);
        //                                                             stdoutlog
        //                                                             =
        //                                                             fopen(stdoutlogname,
        //                                                             "a");
        //                        fprintf(stdoutlog, "s_s\t\t\t
        //                        elocc:%i\n\n",2*int_en_array_ss[sum]);
        //                            fclose(stdoutlog);
        eloc += 2 * int_en_array_ss[sum];
        // exclusion(exclude_indexes,i);
      }
    }
  }

  if (poly_solvent_mov == 0) {
    for (int i = 0; i < 120; i++) {
      //        stdoutlog = fopen(stdoutlogname, "a");
      //        fprintf(stdoutlog, "\t\tijk = {%i,%i,%i}  %i =
      //        local_en_indexes_op[index_of_concern*96+i] index_of_concern
      //        %i\n",x,y,z,checked_index,index_of_concern);
      //        fclose(stdoutlog);

      int checked_index = local_en_indexes_op[index_of_concern * 120 + i];

      int holdover = latticepoint_S[checked_index];

      //                    stdoutlog = fopen(stdoutlogname, "a");
      //                    fprintf(stdoutlog, "\t\treset_ind: %i ijk =
      //                    {%i,%i,%i}  %i =
      //                    local_en_indexes_op[index_of_concern*96+i] %i
      //                    index_of_concern %i\n\t\t\t%i %i   %i
      //                    %i\n",reset_ind,x,y,z,holdover,checked_index,index_of_concern,checked_index<((L1dim_S_z-1)*L1dim_S_xy*L1dim_S_xy),phhh,phhhhh,phhhh);
      //                    fclose(stdoutlog);

      if (holdover == 1 && phh != checked_index) {
        int sum = local_en_dist[i];
        if (checked_index < ((L1dim_S_z - 1) * L1dim_S_xy * L1dim_S_xy)) {
          for (int i = reset_indexes[0] + 1; i < reset_indexes[1]; i++) {
            if (checked_index == p_l_i[i]) {
              eloc -= int_en_array_pp[sum];
            }
          }
          eloc += 2 * int_en_array_pp[sum];

          //                        stdoutlog = fopen(stdoutlogname, "a");
          //                        fprintf(stdoutlog, "p_p\t\t\t
          //                        elocc:%i\n\n",eloc);
          //                            fclose(stdoutlog);
        }

        if (checked_index == phhh) {
          eloc += int_en_array_bond[sum];
          if (phhh == p_l_i[s]) {
            eloc += int_en_array_bond[sum];
          }

          //                        stdoutlog = fopen(stdoutlogname, "a");
          //                        fprintf(stdoutlog, "bpb\t\t\t
          //                        elocc:%i\n\n",eloc);
          //                        fclose(stdoutlog);

          //                        int kkp[3]={c_x,c_y,c_z};
          //                                                stdoutlog =
          //                                                fopen(stdoutlogname,
          //                                                "a");
          //                                                fprintf(stdoutlog,
          //                                                "\t\tbc_index %i %i
          //                                                checked index %i
          //                                                latticepoint_S %i
          //                                                op_array =
          //                                                {%i,%i,%i}
          //                                                c={%i,%i,%i}
          //                                                xyz={%i,%i,%i} eloc
          //                                                %i    \t
          //                                                %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
          //                                                ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array_bond[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
          //                                                fclose(stdoutlog);
        }
        if (phhhhh != 0 && checked_index == phhhh) {
          eloc += int_en_array_bond[sum];
          if (phhhh == p_l_i[e]) {
            eloc += int_en_array_bond[sum];
          }

          //                        stdoutlog = fopen(stdoutlogname, "a");
          //                                                   fprintf(stdoutlog,
          //                                                   "fpb\t\t\t
          //                                                   elocc:%i\n\n",eloc);
          //                                                   fclose(stdoutlog);
          //
        }
        // exclusion(exclude_indexes,i);
      }

      // int kkp[3]={c_x,c_y,c_z};
      if (holdover == 2) {
        int sum = local_en_dist[i];
        //                        stdoutlog = fopen(stdoutlogname, "a");
        //                        fprintf(stdoutlog, "\t\tc_index %i %i  checked
        //                        index %i latticepoint_S %i op_array =
        //                        {%i,%i,%i}   c={%i,%i,%i}  xyz={%i,%i,%i} eloc
        //                        %i\t
        //                        %i\n",(c_z)*(L1dim_S_xy)*(L1dim_S_xy)+(c_y)*(L1dim_S_xy)+(c_x),general_index(kkp)
        //                        ,checked_index,latticepoint_S[checked_index],i,j,k,c_x,c_y,c_z,x,y,z,int_en_array[i*i+j*j+k*k],lattice_dist_optomize[polymer_optomize[general_index(kkp)]+checked_index]);
        //                        fclose(stdoutlog);

        eloc += 2 * int_en_array_ps[sum];
        // exclusion(exclude_indexes,i);

        //                    stdoutlog = fopen(stdoutlogname, "a");
        //                                               fprintf(stdoutlog,
        //                                               "sp\t\t\t
        //                                               elocc:%i\n\n",eloc);
        //                                               fclose(stdoutlog);
      }
    }
  }
  return eloc;
}
