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

#include "General.h"
#include "Energy.h"
#include "Mov.h"


int poly_mov()
{
    
    int mov_choice;
    //std::random_device rd;  // Will be used to obtain a seed for the random number engine
    
    mov_choice = (rand() % (100))+1;
    //mov_choice = 1;
    //int mov_completion = 0;
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to poly move() Move choice:%i\n",mov_choice);
//    fclose(stdoutlog);
    
    if (mov_choice < 30)
    {
     
        poly_solvent_mov = 0; //monomer mov
        chosen_mov =0; //translation
        MoveProposal = mov_translation();
        
        return MoveProposal;
    }
    if (mov_choice >= 30  && mov_choice <99)
    {
     
       chosen_mov =4; //solvent
        poly_solvent_mov = 1; //solvent mov
        MoveProposal = solvent_translation();
        return MoveProposal;
        
        return MoveProposal;
    }
    else
    {
        int sec_choice = (rand() % (4)) +  1;
        //sec_choice = 1;
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Had to make a second choice poly move() sec choice:%i\n",sec_choice);
//        fclose(stdoutlog);
        
        if (sec_choice == 1)
        {
            chosen_mov =1; //crankshaft
            poly_solvent_mov = 0; //monomer mov
            MoveProposal = mov_pivot();
            
            return MoveProposal;
        }
        if (sec_choice == 2)
        {
            chosen_mov =2; //pivot
            poly_solvent_mov = 0; //monomer mov
            MoveProposal = mov_crankshaft();
            
            return MoveProposal;
        }
        if (sec_choice == 3)
        {
            chosen_mov =3; //reflection
            poly_solvent_mov = 0; //monomer mov
            MoveProposal = mov_planar_reflection();
            
            return MoveProposal;
        }
        if (sec_choice == 4)
        {
            chosen_mov =4; //solvent
            poly_solvent_mov = 1; //solvent mov
            MoveProposal = solvent_translation();
            return MoveProposal;
        }

    }


    return 0;
}


int propose_update()    // returns energy change if a spin would be updated
{
    loc_en_before = 0;
    //        stdoutlog = fopen(stdoutlogname, "a");
    //        if(increm == 685)
    //        {
    //            fprintf(stdoutlog, "686 1 latticepoint 3233 %i\n",latticepoint_L[3233]);
    //        }
    //        fclose(stdoutlog);
            
            //int ep = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc); // compares the energy before and after a movement has been preformed
            
            poly_mov();
            
    //        stdoutlog = fopen(stdoutlogname, "a");
    //        fprintf(stdoutlog, "Left polymov()\n");
    //        fclose(stdoutlog);
            
            //int e = total_energy(ph_poly_lattice_coordinates,ph_poly_lattice_indexes,ph_solvent_loc);
            
            int loc_energy_before = 0;
            int loc_energy_after = 0;
            
            
            
           
            
            
    //        stdoutlog = fopen(stdoutlogname, "a");
    //        fprintf(stdoutlog, "Left totalenergyx2 increm: %lli\n   ep: %i   e:%i   ep-e: %i\n",increm%1,ep,e,ep-e);
    //        fprintf(stdoutlog, "local_energy_before: %i   local_energy_after: %i   loc_en_b-loc_en_a: %i\n",loc_energy_before,loc_energy_after,loc_energy_before-loc_energy_after);
    //        fclose(stdoutlog);
            
                
            if(MoveProposal==1)
                           {
                               if(poly_solvent_mov==0)
                               {
                                   //lattice_polymer_reset(poly_lattice_indexes , ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                   lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                   poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                               }
                               if(poly_solvent_mov==1)
                               {
                                   solvent_reset(solvent_loc,ph_solvent_loc);
                                   
                   //                stdoutlog = fopen(stdoutlogname, "a");
                   //                fprintf(stdoutlog, "solvent_loc[i]: %i ph_solvent_loc[i]: %i\n",solvent_loc[reset_indexes[0]+1],ph_solvent_loc[reset_indexes[0]+1]);
                   //                fclose(stdoutlog);
                               }
                           }
            
//            if(loc_energy_before-loc_energy_after != ep-e)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "\nMade it to local_energy\n");
//                fprintf(stdoutlog, "chosen_mov:%i \t poly_solv:%i\n",chosen_mov,poly_solvent_mov);
//                fprintf(stdoutlog, "reset_indexes[0]+1:%i \t reset_indexes[1]:%i\n",reset_indexes[0]+1,reset_indexes[1]);
//                fclose(stdoutlog);
//
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "2Difference  %i\n",cccounter);
//                fclose(stdoutlog);
//                 eye();
//                MPI_Abort(MPI_COMM_WORLD, 1);
//            }
    
    if(MoveProposal==1)
               {
                   loc_energy_before  = loc_en_before;
               loc_energy_after = local_energy(ph_poly_lattice_coordinates,ph_poly_lattice_indexes,ph_solvent_loc);
               }
               else
               {
                   loc_energy_before =0;
               }
    
    return(loc_energy_after - loc_energy_before);
}


int mov_crankshaft()
{
    int arm = rand() % (numbercores * pertrusionpercore);// choose the arm the crankshaft will be performed on
    int poly_index_p1 = arm * lengthpercore + (rand() % (lengthpercore - 2));
//    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
//    general_coord(point1, poly_lattice_indexes[poly_index_p1]);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov crankshaft()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);

    
    
    int poly_index_p2 = arm * lengthpercore + (rand() % lengthpercore);
    while (poly_index_p2 <= poly_index_p1 + 1)
    {
        poly_index_p2 = arm * lengthpercore + (rand() % lengthpercore);
    }
//    int point2[3]; // choose the end point on the arm to perform the movement (must be farther along the arm and 2 indexes from point 1
//    general_coord(point2, poly_lattice_indexes[poly_index_p2]);
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "2nd point defined\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
//    fclose(stdoutlog);
    
    int gap = poly_index_p2 - poly_index_p1 - 1; // defines the length of the gap between the monomers
    
    reset_indexes[0] = poly_index_p1;
    reset_indexes[1] = poly_index_p2;
    
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1, poly_index_p2); // clears all of the lattice points between the two axis points

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i  latticepoint_S[poly_lattice_indexes[poly_index_p1+1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p2]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1]],latticepoint_S[poly_lattice_indexes[poly_index_p1+1]],latticepoint_S[poly_lattice_indexes[poly_index_p2]]);
//    fclose(stdoutlog);
    
    double rotated_vector[3]; // holds the rotated point values
    int int_rotated_vector[3]; // holds the rounnded interger point values

    double starting_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    double end_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p2],(double)poly_lattice_coordinates[3 * poly_index_p2 + 1],(double)poly_lattice_coordinates[3 * poly_index_p2 + 2] };
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
//    fclose(stdoutlog);
    
    int assignment_possible = 0; // checks to ensure assignment is possible

    
    for (int i = 0; i < gap; i++)
    {
        int initial_point[3]={ poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i)],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 1],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 2] };


//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
//        fclose(stdoutlog);
        
        rotation_wrt_vector(starting_point, end_point, initial_point, rotated_vector);

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "rotated={%f,%f,%f}\n",rotated_vector[0],rotated_vector[1],rotated_vector[2]);
//        fclose(stdoutlog);
//
        for (int iii = 0; iii < 3; iii++) // rounds the outputted values from the rotated vector
        {
            if (rotated_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(rotated_vector[iii] + .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] + .5);
            }
            if (rotated_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(rotated_vector[iii] - .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] - .5);
            }
        }
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "int_rotated_vector={%i,%i,%i}\n",int_rotated_vector[0],int_rotated_vector[1],int_rotated_vector[2]);
//        fclose(stdoutlog);
        
        int new_index = general_index(int_rotated_vector);
                
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//        fclose(stdoutlog);
        
        assignment_possible = latticepoint_assign(new_index,1); // checks too see if assignment of the new lattice point wwas possible
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//        fclose(stdoutlog);
        
        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p1 + 1 + i] = new_index; // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+2] = index_z[new_index];
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }
        
        if (new_index >= L1dim_S_xy*L1dim_S_xy*plane_z_max || new_index < L1dim_S_xy*L1dim_S_xy*(plane_z_min+1))
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit1\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p1 + 1 + j)],poly_lattice_indexes[(poly_index_p1 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < no_solvent_sites; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\n", solvent_loc[i],ph_solvent_loc[i]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1 + 1 + i], ph_poly_lattice_indexes[poly_index_p1 + i]) >= int_max_dist && assignment_possible == 1) //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit2\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2], ph_poly_lattice_indexes[poly_index_p2 - 1]) >= int_max_dist)  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "early exit3\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 1;
    }

    return 1;
}

int mov_pivot()
{
    int arm = rand() % (numbercores * pertrusionpercore);// choose the arm the crankshaft will be performed on
    int poly_index_p1 = arm * lengthpercore + (rand() % (lengthpercore - 2));
//    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
//    general_coord(point1, poly_lattice_indexes[poly_index_p1]);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov pivot()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);

    int poly_index_p2 = poly_index_p1 + 1;
//    int point2[3]; // initializes end point that is one away from previous to create the pivot axis
//    general_coord(point2, poly_lattice_indexes[poly_index_p2]);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "2nd point defined\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
//    fclose(stdoutlog);
    
    int poly_index_arm_end = arm * lengthpercore + lengthpercore - 1;

    int gap = poly_index_arm_end - poly_index_p2;
    
    //int gap = poly_index_p2 - poly_index_p1 - 1; // defines the length of the gap between the monomers
    reset_indexes[0] = poly_index_p2-1;
    reset_indexes[1] = arm * lengthpercore + lengthpercore;
    
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p2, arm * lengthpercore + lengthpercore); // clears all of the lattice points between the two axis points
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i\n\n",latticepoint_S[poly_lattice_indexes[arm * lengthpercore + lengthpercore-1]],latticepoint_S[poly_lattice_indexes[arm * lengthpercore + lengthpercore-1]]);
//    fclose(stdoutlog);
    
        
    double rotated_vector[3]; // holds the rotated point values
    int int_rotated_vector[3]; // holds the rounnded interger point values

    double starting_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 1],(double)poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    double end_point[3] = { (double)poly_lattice_coordinates[3 * poly_index_p2],(double)poly_lattice_coordinates[3 * poly_index_p2 + 1],(double)poly_lattice_coordinates[3 * poly_index_p2 + 2] };

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "sp={%f,%f,%f} ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2],end_point[0],end_point[1],end_point[2]);
//    fclose(stdoutlog);
    
    int assignment_possible = 0; // checks to ensure assignment is possible

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "initial check\n\n");
//    fclose(stdoutlog);
    
//    for (int j = 0; j < gap; j++)
//    {
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "%i: ph_poly:%i  %i poly:%i  %i\n",poly_index_p2 + 1 + j,ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)]],poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[poly_lattice_indexes[(poly_index_p2 + 1 + j)]]);
//        fclose(stdoutlog);
//    }

    for (int i = 0; i < gap; i++)
    {
        int initial_point[3]={ poly_lattice_coordinates[3 * (poly_index_p2 + 1 + i)],poly_lattice_coordinates[3 * (poly_index_p2 + 1 + i) + 1],poly_lattice_coordinates[3 * (poly_index_p2 + 1 + i) + 2] };

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "%i: init={%i,%i,%i}\n",poly_index_p2 + 1 + i,initial_point[0],initial_point[1],initial_point[2]);
//        fclose(stdoutlog);
        
        rotation_wrt_vector(starting_point, end_point, initial_point, rotated_vector);

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "rotated={%f,%f,%f}\n",rotated_vector[0],rotated_vector[1],rotated_vector[2]);
//        fclose(stdoutlog);
        
        for (int iii = 0; iii < 3; iii++) // rounds the outputted values from the rotated vector
        {
            if (rotated_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p2 + 1 + i) + iii] = (int)(rotated_vector[iii] + .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] + .5);
            }
            if (rotated_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p2 + 1 + i) + iii] = (int)(rotated_vector[iii] - .5);
                int_rotated_vector[iii] = (int)(rotated_vector[iii] - .5);
            }

        }
        
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "int_rotated_vector={%i,%i,%i}\n",int_rotated_vector[0],int_rotated_vector[1],int_rotated_vector[2]);
//        fclose(stdoutlog);
        
        int new_index = general_index(int_rotated_vector);
       
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//        fclose(stdoutlog);
        
        assignment_possible = latticepoint_assign(new_index,1); // checks too see if assignment of the new lattice point wwas possible
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//        fclose(stdoutlog);
        
        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p2 + 1 + i] = new_index; // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i)*3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i)*3+1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p2 + 1 + i)*3+2] = index_z[new_index];
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            
            
            return 0;
        }
        
        if (new_index >= L1dim_S_xy*L1dim_S_xy*plane_z_max || new_index < L1dim_S_xy*L1dim_S_xy*(plane_z_min+1))
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            return 0;
        }

        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "%i: Distance Check %i  %i   %i\n",poly_index_p2 + 1 + i, reduced_lattice_optimization(ph_poly_lattice_coordinates[poly_index_p2 + 1 + i],ph_poly_lattice_coordinates[ poly_index_p2 + i]),ph_poly_lattice_coordinates[poly_index_p2 + 1 + i],ph_poly_lattice_coordinates[poly_index_p2 + i]);
//        fclose(stdoutlog);
        
        
        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2 + 1 + i],ph_poly_lattice_indexes[ poly_index_p2 + i]) >= int_max_dist && assignment_possible == 1) //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
            

//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  %i poly:%i  %i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)]],poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[poly_lattice_indexes[(poly_index_p2 + 1 + j)]]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            
            return 0;
        }

    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[arm * lengthpercore + lengthpercore - 2],ph_poly_lattice_indexes[ arm * lengthpercore + lengthpercore - 1]) >= int_max_dist)  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "early exit 2\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  %i poly:%i  %i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)]],poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[poly_lattice_indexes[(poly_index_p2 + 1 + j)]]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "pivot success\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  %i poly:%i  %i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)]],poly_lattice_indexes[(poly_index_p2 + 1 + j)],latticepoint_S[poly_lattice_indexes[(poly_index_p2 + 1 + j)]]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        
        return 1;
    }

    return 1;

}

int mov_inversion()
{
    return 0;
}

int mov_planar_reflection()
{
    int arm = rand() % (numbercores * pertrusionpercore);// choose the arm the crankshaft will be performed on
    int poly_index_p1 = arm * lengthpercore + (rand() % (lengthpercore - 2));
//    int point1[3]; // choose the initial point to perform the movement (must be 2 less than length of arm)
//    general_coord(point1, poly_lattice_indexes[poly_index_p1]);
//
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov planar reflection()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%i point1={%i,%i,%i} poly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,point1[0],point1[1],point1[2],poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);
    
    
    int poly_index_p2 = arm * lengthpercore + (rand() % lengthpercore);
    while (poly_index_p2 <= poly_index_p1 + 1)
    {
        poly_index_p2 = arm * lengthpercore + (rand() % lengthpercore);
    }
//    int point2[3]; // choose the end point on the arm to perform the movement (must be farther along the arm and 2 indexes from point 1
//    general_coord(point2, poly_lattice_indexes[poly_index_p2]);
//
//
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "2nd point defined\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p2=%i point2={%i,%i,%i} poly_lattice_indexes[poly_index_p2] = %i\n\n",arm,poly_index_p2,point2[0],point2[1],point2[2],poly_lattice_indexes[poly_index_p2]);
//    fclose(stdoutlog);
    
    
    int gap = poly_index_p2 - poly_index_p1 - 1; // defines the length of the gap between the monomers
    
    reset_indexes[0] = poly_index_p1;
    reset_indexes[1] = poly_index_p2;
    
    loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1, poly_index_p2); // clears all of the lattice points between the two axis points

    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear() gap: %i\n",gap);
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i poly_lattice_indexes[poly_index_p1+1] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1]],latticepoint_S[poly_lattice_indexes[poly_index_p1+1]] ,latticepoint_S[poly_lattice_indexes[poly_index_p2]]);
//    fclose(stdoutlog);
    
    
    int axis = rand() % 3;

    double reflected_vector[3]; // holds the rotated point values
    int int_reflected_vector[3]; // holds the rounnded interger point values

    double starting_point[3] = { (double)ph_poly_lattice_coordinates[3 * poly_index_p1],(double)ph_poly_lattice_coordinates[3 * poly_index_p1 + 1],(double)ph_poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    double end_point[3] = { (double)ph_poly_lattice_coordinates[3 * poly_index_p2], (double)ph_poly_lattice_coordinates[3 * poly_index_p2 + 1], (double)ph_poly_lattice_coordinates[3 * poly_index_p2 + 2] };
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "sp={%f,%f,%f} axis = %i ep={%f,%f,%f}\n",starting_point[0],starting_point[1],starting_point[2], axis,end_point[0],end_point[1],end_point[2]);
//    fclose(stdoutlog);

    int assignment_possible = 0; // checks to ensure assignment is possible

    for (int i = 0; i < gap; i++)
    {
        int initial_point[3]={ poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i)],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 1],poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + 2] };

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "init={%i,%i,%i}\n",initial_point[0],initial_point[1],initial_point[2]);
//        fclose(stdoutlog);
        
        reflection_wrt_plane(starting_point, end_point, initial_point, reflected_vector, axis);

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "reflected={%f,%f,%f}\n",reflected_vector[0],reflected_vector[1],reflected_vector[2]);
//        fclose(stdoutlog);
        
        for (int iii = 0; iii < 3; iii++) // rounds the outputted values from the rotated vector
        {
            if (reflected_vector[iii] >= 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(reflected_vector[iii] + .5);
                int_reflected_vector[iii] = (int)(reflected_vector[iii] + .5);
            }
            if (reflected_vector[iii] < 0)
            {
                ph_poly_lattice_coordinates[3 * (poly_index_p1 + 1 + i) + iii] = (int)(reflected_vector[iii] - .5);
                int_reflected_vector[iii] = (int)(reflected_vector[iii] - .5);
            }
        }
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "int_reflected_vector={%i,%i,%i}\n",int_reflected_vector[0],int_reflected_vector[1],int_reflected_vector[2]);

        //fclose(stdoutlog);
        
        int new_index = general_index(int_reflected_vector);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//        fclose(stdoutlog);
        
        assignment_possible = latticepoint_assign(new_index,1); // checks too see if assignment of the new lattice point wwas possible
//
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//        fclose(stdoutlog);
        
        if (assignment_possible == 1)
        {
            ph_poly_lattice_indexes[poly_index_p1 + 1 + i] = new_index; // generates an index from the rotated vector
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3] = index_x[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+1] = index_y[new_index];
            ph_poly_lattice_coordinates[(poly_index_p1 + 1 + i)*3+2] = index_z[new_index];
        }
        else
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit\n\n");
//            fclose(stdoutlog);

//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            return 0;
        }
        
        if (new_index >= L1dim_S_xy*L1dim_S_xy*plane_z_max || new_index < L1dim_S_xy*L1dim_S_xy*(plane_z_min+1))
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit0\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            return 0;
        }
        

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog,"ph_poly_lattice_indexes[poly_index_p1 + 1 + i]:%i   ph_poly_lattice_indexes[poly_index_p1 + i]:%i   max_dist: %f\n\n",poly_index_p1 + 1 + i,poly_index_p1 + i,max_dist);
//        fclose(stdoutlog);
        
        if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1 + 1 + i],ph_poly_lattice_indexes[ poly_index_p1 + i]) >= int_max_dist && assignment_possible == 1) //checks new distance between created monomer and previous monomer for every change
        {
            //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
            poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "early exit1\n\n");
//            fclose(stdoutlog);
//
//            for (int j = 0; j < gap; j++)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//                fclose(stdoutlog);
//            }
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//            {
//                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//            }
//            fclose(stdoutlog);
            return 0;
        }

    }

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "max_dist: %f\n\n",max_dist);
//    fclose(stdoutlog);
    
    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p2], ph_poly_lattice_indexes[poly_index_p2 - 1]) >= int_max_dist)  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "early exit2\n\n");
//        fclose(stdoutlog);
//
//        for (int j = 0; j < gap; j++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ph_poly:%i  poly:%i\n",ph_poly_lattice_indexes[(poly_index_p2 + 1 + j)],poly_lattice_indexes[(poly_index_p2 + 1 + j)]);
//            fclose(stdoutlog);
//        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
//        {
//            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
//        }
//        fclose(stdoutlog);
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        //poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "success\n\n");
//        fclose(stdoutlog);

        return 1;
    }

    return 1;
}

int mov_translation()
{
    int arm = rand() % (numbercores * pertrusionpercore);// choose the arm the translation mov will be performed on
    int poly_index_p1 = arm * lengthpercore + (rand() % (lengthpercore - 1)) + 1; // cant eqaul 0 otherwise the core would move
    //int point1[3]; // choose the initial point to perform the movement
    
    //general_coord(point1, poly_lattice_indexes[poly_index_p1]);
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to mov translation()\n");
//    fprintf(stdoutlog, "arm = %i poly_index_p1=%ipoly_lattice_indexes[poly_index_p1] = %i\n\n",arm,poly_index_p1,poly_lattice_indexes[poly_index_p1]);
//    fclose(stdoutlog);
    
    reset_indexes[0] = poly_index_p1 - 1;
    reset_indexes[1] = poly_index_p1 + 1;
    
    
    
    poly_latticepoint_clear(poly_lattice_indexes, poly_index_p1 - 1, poly_index_p1 + 1); // clears all of the lattice points between the two axis points
    
   
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called poly_clear()\n");
//    fprintf(stdoutlog, "latticepoint_S[poly_lattice_indexes[poly_index_p1-1]] = %i latticepoint_S[poly_lattice_indexes[poly_index_p1]] = %i\n\n",latticepoint_S[poly_lattice_indexes[poly_index_p1-1]],latticepoint_S[poly_lattice_indexes[poly_index_p1]]);
//    fclose(stdoutlog);

    int int_translated_vector[3] = { poly_lattice_coordinates[3 * poly_index_p1],poly_lattice_coordinates[3 * poly_index_p1 + 1],poly_lattice_coordinates[3 * poly_index_p1 + 2] };
    //double end_point[3] = { poly_lattice_coordinates[3 * poly_index_p2],poly_lattice_coordinates[3 * poly_index_p2 + 1],poly_lattice_coordinates[3 * poly_index_p2 + 2] };
    
    //int int_translated_vector[3] = { starting_point[0],starting_point[1],starting_point[2] }; // holds the translated interger point values

    
    
    int assignment_possible = 0; // checks to ensure assignment is possible
    
    int random_comp = (rand() % 3);
    int value;
    if((rand() % 2)==0)
    {
        value = -1;
    }
    else{
        value =  1;
    }
    
    int multiply = (rand() % 2)+1;
    
    int_translated_vector[random_comp] += value*multiply;

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "tv={%i,%i,%i}\n",int_translated_vector[0],int_translated_vector[1],int_translated_vector[2]);
//    fclose(stdoutlog);
    
    int new_index = general_index_op(int_translated_vector);
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called general_index() new_index=%i\n",new_index);
//    fclose(stdoutlog);
    
    assignment_possible = latticepoint_assign2(new_index,1); // checks too see if assignment of the new lattice point wwas possible
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Called assignment possible() assignment possible=%i\n\n",assignment_possible);
//    fclose(stdoutlog);
    
    if (assignment_possible == 1)
    {
         
        
        ph_poly_lattice_indexes[poly_index_p1] = new_index; // generates an index from the rotated vector
        ph_poly_lattice_coordinates[poly_index_p1*3] = index_x[new_index];
        ph_poly_lattice_coordinates[poly_index_p1*3+1] = index_y[new_index];
        ph_poly_lattice_coordinates[poly_index_p1*3+2] = index_z[new_index];
    }
    else
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }
    
    if (new_index >= lower_bounds || new_index < upper_bounds)
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }
    
    //assignment_possible = latticepoint_assign(new_index,1); // checks too see if assignment of the new lattice point wwas possible

    
    
    if ( (poly_index_p1+1)%lengthpercore != 0 && (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1],ph_poly_lattice_indexes[poly_index_p1+1])>= int_max_dist )) //checks new distance between created monomer and previous monomer for every change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }

    if (reduced_lattice_optimization(ph_poly_lattice_indexes[poly_index_p1],ph_poly_lattice_indexes[poly_index_p1-1]) >= int_max_dist  )  //checks new distance between last created monomer and considered end monomer for last change
    {
        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
        return 0;
    }
    else  //checks new distance between last created monomer and considered end monomer for last change
    { // if a move is successful do not immediately reset non ph_poly_lattice in case of movement rejection
        //lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, poly_index_p1, poly_index_p2);
        loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
        
        return 1;
    }

    return 1;
}




int reflection_wrt_plane(double(&s_pos)[3], double(&e_pos)[3], int(&init_p)[3], double(&refl_p)[3], int axis) // output index from any x , y , z components
{
    if (axis == 0) // axis_end_pos[0]+1
    {
        if ((e_pos[1] - s_pos[1])==0 && (e_pos[2] - s_pos[2])==0) {
            refl_p[0] = init_p[0];
            refl_p[1] = init_p[1];
            refl_p[2] = init_p[2];
            return 0;
        }
        
        
        refl_p[0] = init_p[0]; //Good
        
        refl_p[1]=((e_pos[1]*e_pos[1])*init_p[1]-(e_pos[2]*e_pos[2]) *(init_p[1]-2*s_pos[1])+init_p[1]* (s_pos[1]*s_pos[1])+2* e_pos[2] *init_p[1] *s_pos[2]+2 *init_p[2]* s_pos[1] *s_pos[2]-init_p[1]* (s_pos[2]*s_pos[2])-2*e_pos[2] *s_pos[1] *(init_p[2]+s_pos[2])+2 *e_pos[1]* (e_pos[2]* init_p[2]-init_p[1]* s_pos[1]-(e_pos[2]+init_p[2])* s_pos[2]+(s_pos[2]*s_pos[2])))/((e_pos[1]-s_pos[1])*(e_pos[1]-s_pos[1])+(e_pos[2]-s_pos[2])*(e_pos[2]-s_pos[2]));

        
        refl_p[2] = init_p[2]+(2 *(e_pos[1]-s_pos[1]) *(e_pos[2] *(init_p[1]-s_pos[1])+init_p[2]* s_pos[1]-init_p[1] *s_pos[2]+e_pos[1]* (-init_p[2]+s_pos[2])))/(((e_pos[1]-s_pos[1])*(e_pos[1]-s_pos[1]))+((e_pos[2]-s_pos[2])*(e_pos[2]-s_pos[2])));
    }
    if (axis == 1)
    {
        
        if ((e_pos[0] - s_pos[0])==0 && (e_pos[2] - s_pos[2])==0) {
            refl_p[0] = init_p[0];
            refl_p[1] = init_p[1];
            refl_p[2] = init_p[2];
            return 0;
        }
        
        refl_p[0] = ((e_pos[0]*e_pos[0])*
                     init_p[0] - (e_pos[2]*e_pos[2])* (init_p[0] - 2 *s_pos[0]) +
                    init_p[0] *(s_pos[0]*s_pos[0]) + 2 *e_pos[2] *init_p[0] *s_pos[2] +
                    2 *init_p[2]* s_pos[0] *s_pos[2] - init_p[0] *(s_pos[2]*s_pos[2]) -
                    2 *e_pos[2] *s_pos[0]* (init_p[2] + s_pos[2]) +
                    2 *e_pos[0] *(e_pos[2] *init_p[2] -
                       init_p[0] *s_pos[0] - (e_pos[2] + init_p[2]) *
                        s_pos[2] + (s_pos[2]*s_pos[2])));

        refl_p[0] /= (((e_pos[0] - s_pos[0])*(e_pos[0] - s_pos[0])) + ((e_pos[2] -
                                                                        s_pos[2])*(e_pos[2] - s_pos[2])));

        refl_p[1] = init_p[1];

        refl_p[2] = init_p[2]+(2*(e_pos[0]-s_pos[0])*(e_pos[2]*(init_p[0]-s_pos[0])+init_p[2]* s_pos[0]-init_p[0]*s_pos[2]+e_pos[0]* (-init_p[2]+s_pos[2])))/((e_pos[0]-s_pos[0])*(e_pos[0]-s_pos[0])+(e_pos[2]-s_pos[2])*(e_pos[2]-s_pos[2]));
        
    }
    if (axis == 2)
    {
        
        if ((e_pos[0] - s_pos[0])==0 && (e_pos[1] - s_pos[1])==0) {
            refl_p[0] = init_p[0];
            refl_p[1] = init_p[1];
            refl_p[2] = init_p[2];
            return 0;
        }
        
        refl_p[0] = ((e_pos[0]*e_pos[0])*
                     init_p[0] - (e_pos[1]*e_pos[1])*(init_p[0] - 2 *s_pos[0]) +
                    init_p[0] *(s_pos[0]*s_pos[0]) + 2 *e_pos[1]* init_p[0]* s_pos[1] +
                    2 *init_p[1] *s_pos[0]* s_pos[1] - init_p[0]* (s_pos[1]*s_pos[1]) -
                    2* e_pos[1]* s_pos[0]* (init_p[1] + s_pos[1]) +
                    2 *e_pos[0]*(e_pos[1] *init_p[1] -
                       init_p[0]* s_pos[0] - (e_pos[1] + init_p[1])*
                        s_pos[1] + (s_pos[1]*s_pos[1])));

        refl_p[0] /= (((e_pos[0] - s_pos[0])*(e_pos[0] - s_pos[0])) + ((e_pos[1] -
                                                                        s_pos[1])*(e_pos[1] - s_pos[1])));

        refl_p[1] = init_p[1]+(2*(e_pos[0]-s_pos[0])*(e_pos[1]*(init_p[0]-s_pos[0])+init_p[1]* s_pos[0]-init_p[0]*s_pos[1]+e_pos[0]* (-init_p[1]+s_pos[1])))/((e_pos[0]-s_pos[0])*(e_pos[0]-s_pos[0])+(e_pos[1]-s_pos[1])*(e_pos[1]-s_pos[1]));

        refl_p[2] = init_p[2];
    }
    //return ( arr[0]+arr[1]*L1dim_S+arr[2]*L1dim_S*L1dim_S)
    return 0;
}


int rotation_wrt_vector(double(&s_pos)[3], double(&e_pos)[3], int(&init_p)[3], double(&rot_p)[3]) // output index from any x , y , z components
{
    double angle = (((rand() % 23 + 1) * 15) * 3.14159265) / 180.0;

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "angle = %f\n",angle);
//    fclose(stdoutlog);
    
    double ph =
        sqrt((s_pos[0] - e_pos[0]) * (s_pos[0] - e_pos[0]) + (s_pos[1] -
            e_pos[1]) * (s_pos[1] - e_pos[1]) + (s_pos[2] -
                e_pos[2]) * (s_pos[2] - e_pos[2]));

    double phh = (s_pos[0] - e_pos[0]) * (s_pos[0] -
        e_pos[0]) + (s_pos[1] - e_pos[1]) * (s_pos[1] -
            e_pos[1]) + (s_pos[2] - e_pos[2]) * (s_pos[2] - e_pos[2]);
    double sin_sub = sin(angle);
    double cos_sub = cos(angle);


    rot_p[0] = ((s_pos[0] * s_pos[0]) * init_p[0] + (e_pos[0] * e_pos[0]) * init_p[0] +
        cos_sub *
        init_p[0] * ((s_pos[1] - e_pos[1]) * (s_pos[1] -
            e_pos[1]) + (s_pos[2] - e_pos[2]) * (s_pos[2] - e_pos[2])) -
        sin_sub * e_pos[1] * s_pos[2] * ph +
        sin_sub * init_p[1] * s_pos[2] * ph +
        sin_sub * s_pos[1] * ph * e_pos[2] -
        sin_sub * init_p[1] * ph * e_pos[2] +
        s_pos[0] * (-2 * e_pos[0] *
            init_p[0] + (-1 +
                cos_sub) * ((s_pos[1] - e_pos[1]) * (e_pos[1] -
                    init_p[1]) + (s_pos[2] - e_pos[2]) * (e_pos[2] -
                        init_p[2]))) - (-1 + cos_sub) *
        e_pos[0] * ((s_pos[1] - e_pos[1]) * (s_pos[1] -
            init_p[1]) + (s_pos[2] - e_pos[2]) * (s_pos[2] - init_p[2])) -
        sin_sub * s_pos[1] * ph * init_p[2] +
        sin_sub * e_pos[1] * ph * init_p[2]) / (phh);

    rot_p[1] = init_p[1] + (sin_sub * e_pos[0] * s_pos[2]) /
        ph - (sin_sub * init_p[0] * s_pos[2]) /
        ph - (sin_sub * s_pos[0] * e_pos[2]) /
        ph + (sin_sub * init_p[0] * e_pos[2]) / ph + (sin_sub * s_pos[0] * init_p[2]) / ph - ((sin_sub * e_pos[0] * init_p[2]) / ph) + ((-1 + cos_sub) * (e_pos[0] * init_p[0] * (s_pos[1] - e_pos[1]) +
            s_pos[0] * (init_p[0] * (-s_pos[1] + e_pos[1]) +
                e_pos[0] * (s_pos[1] + e_pos[1] - 2 * init_p[1])) +
            e_pos[0] * e_pos[0] * (-s_pos[1] + init_p[1]) +
            s_pos[0] *
            s_pos[0] * (-e_pos[1] + init_p[1]) - (s_pos[2] -
                e_pos[2]) * (init_p[1] * (-s_pos[2] + e_pos[2]) +
                    e_pos[1] * (s_pos[2] - init_p[2]) +
                    s_pos[1] * (-e_pos[2] + init_p[2])))) / phh;

    rot_p[2] = -(sin_sub * e_pos[0] * s_pos[1]) / ph + (sin_sub * init_p[0] * s_pos[1]) /
        ph + (sin_sub * s_pos[0] * e_pos[1]) /
        ph - (sin_sub * init_p[0] * e_pos[1]) /
        ph - (sin_sub * s_pos[0] * init_p[1]) /
        ph + (sin_sub * e_pos[0] * init_p[1]) +
        e_pos[2] + ((-1 + cos_sub) * (s_pos[2] -
            e_pos[2]) * ((s_pos[0] - e_pos[0]) * (e_pos[0] -
                init_p[0]) + (s_pos[1] - e_pos[1]) * (e_pos[1] -
                    init_p[1]) + (s_pos[2] - e_pos[2]) * (e_pos[2] - init_p[2])) /
            phh) + cos_sub * (-e_pos[2] + init_p[2]);
    
    return 0;
}


int solvent_translation() // attempts to move a solvent molecule to a different location
{
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to solvent translation  %i\n",rand()%numberspins_S );
//    fclose(stdoutlog);
    
    int random_index = rand()%numberspins_S; // chooses a random insex
    //int random_index = (rand()%solvent_max) + solvent_offset; // chooses a random insex
    
    int random_solvent_site = rand()%no_solvent_sites; // chooses random_solvent molecule
    
    
    
    while((random_index >= lower_bounds || random_index < upper_bounds)) //prevents the new index from being in certain planes
    {
        random_index = rand()%numberspins_S;

//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "random_index:%i  %i  %i\n",random_index,random_index<((L1dim_S_xy*L1dim_S_xy)),random_index>(((L1dim_S_xy)*(L1dim_S_xy)*(L1dim_S_z-1))));
//        fclose(stdoutlog);
    }
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "random_index:%i   random_solvent_site:%i   solvent_loc[random_solvent_site]: %i \n\n",random_index,random_solvent_site,solvent_loc[random_solvent_site]);
//    fclose(stdoutlog);
    reset_indexes[0] =random_solvent_site-1;
    reset_indexes[1] =random_solvent_site+1;
    
    
    
    poly_latticepoint_clear(solvent_loc[random_solvent_site]); // clears the solvent in the previous localtion
    
    
    
    int assignment_possible = latticepoint_assign2(random_index,2); // checks to see if assignment was possible
            
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "assignment possible %i\n",assignment_possible);
//    fclose(stdoutlog);
            
    if(assignment_possible) //if assignment is possible
    {
        
        loc_en_before = local_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
        
        reset_solvent_index = random_solvent_site; // reset solvent index tracked for wl
        ph_solvent_loc[random_solvent_site] = random_index; // new index is assigned to solvent_loc
        reset_indexes[0] =random_solvent_site-1;
        reset_indexes[1] =random_solvent_site+1;
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "random_index:%i   random_solvent_site:%i   solvent_loc[random_solvent_site]: %i \n\n",random_index,random_solvent_site,ph_solvent_loc[random_solvent_site]);
//        fclose(stdoutlog);
        
        return 1;
    }
    else{
        reset_indexes[0] =random_solvent_site-1;
        reset_indexes[1] =random_solvent_site+1;
        solvent_reset(ph_solvent_loc,solvent_loc);
        
    } return 0;
    
    return 0;
}
