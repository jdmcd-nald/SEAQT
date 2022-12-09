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


void solvent_reset(int* p_index,int* c_index)
{
    for (int i = reset_indexes[0]+1; i < reset_indexes[1]; i++) // clears the lattice polymer in the range set
    {
        latticepoint_S[p_index[i]] = 0;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L * p_index[i] + iii]] = 0;
        }
    }

    for (int i = reset_indexes[0]+1; i < reset_indexes[1]; i++) // clears the lattice polymer in the range set
    {
        latticepoint_S[c_index[i]] = 2;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L * c_index[i] + iii]] = 2;
        }

        p_index[i] = c_index[i];
    }
}


int anint(int a, int b)
{
    float c = (float)(a) / ((float)b);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to anint()  c:%f\n",c);
    //fclose(stdoutlog);
    
    if(c>=0)
    {
        return ((int)(c+.5));
    }
    if(c<0)
    {
        return ((int)(c-.5));
    }
    if (c <-0.5 && c>-1.5)
    {
        return -1;
    }
    if (c <= 0.5 && c >= -0.5)
    {
        return 0;
    }
    if (c < 1.5 && c > 0.5)
    {
        return 1;
    }
    return 0;
}

int anint(float c)
{
    
    if(c>=0)
    {
        return ((int)(c+.5));
    }
    if(c<0)
    {
        return ((int)(c-.5));
    }
    
    if (c <-0.5 && c>-1.5)
    {
        return -1;
    }
    if (c <= 0.5 && c >= -0.5)
    {
        return 0;
    }
    if (c < 1.5 && c > 0.5)
    {
        return 0;
    }
    return 0;
}

void init_poly_cores()
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to init poly_core()\n");
    fclose(stdoutlog);
    //holds the indexes of the chains

    poly_lattice_indexes[0] = 149400-L1dim_S_xy*L1dim_S_xy*50;//core definition

    poly_lattice_indexes[1 * lengthpercore] = 149978-L1dim_S_xy*L1dim_S_xy*50; //core definition

    poly_lattice_indexes[2 * lengthpercore] = 149238-L1dim_S_xy*L1dim_S_xy*50;//core definition
    
    poly_lattice_indexes[3 * lengthpercore] =149816-L1dim_S_xy*L1dim_S_xy*50;//core definition
    
    poly_lattice_indexes[4 * lengthpercore] = 149416-L1dim_S_xy*L1dim_S_xy*50;//core definition

    poly_lattice_indexes[5 * lengthpercore] = 149994-L1dim_S_xy*L1dim_S_xy*50; //core definition
    
    poly_lattice_indexes[6 * lengthpercore] = 149594-L1dim_S_xy*L1dim_S_xy*50;//core definition
    
    poly_lattice_indexes[7 * lengthpercore] = 150172-L1dim_S_xy*L1dim_S_xy*50;//core definition

    
    
    for(int j=0;j<numbercores;j++)
    {
        for(int i=1;i<lengthpercore;i++)
        {
            poly_lattice_indexes[(j*lengthpercore)+i] = poly_lattice_indexes[(j*lengthpercore)+i-1] - L1dim_S_xy*L1dim_S_xy*2;
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "%i\n", poly_lattice_indexes[(j*lengthpercore)+i] );
//            fclose(stdoutlog);
            
        }
    }
    
    if(myid<minmaxid)
    {
    if (fopen("Minimum_Config_Poly.txt","r") == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
        fclose(stdoutlog);
    }
    else // rewrite of prior file reading system to be more portable
    {
        FILE* fptr;
        int value;
        int valuee;
        
        fptr = fopen("Minimum_Config_Poly.txt", "r");

        while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
        {
            poly_lattice_indexes[value] = valuee;
        }

        fclose(fptr);
    }
    }
    
    if(myid>=minmaxid)
    {
    if (fopen("Maximum_Config_Poly.txt","r") == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
        fclose(stdoutlog);
    }
    else // rewrite of prior file reading system to be more portable
    {
        FILE* fptr;
        int value;
        int valuee;
        
        fptr = fopen("Maximum_Config_Poly.txt", "r");

        while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
        {
            poly_lattice_indexes[value] = valuee;
        }

        fclose(fptr);
    }
    }
    
    // writes the chain occupation to lattice point
    int poly_index = -1;
    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
    {
        latticepoint_S[poly_lattice_indexes[i]] = 1;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]] = 1;
        }
    }

 optomize();

//    int new_index;
//    for(int i=0;i<numbercores*pertrusionpercore;i++)
//    {
//        poly_latticepoint_clear(poly_lattice_indexes, i*lengthpercore*pertrusionpercore, (1+i)*lengthpercore*pertrusionpercore); // clears all of the lattice points between the two axis points
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "base %i\n",poly_lattice_indexes[i*pertrusionpercore*lengthpercore]);
//        fclose(stdoutlog);
//        for(int j=1;j<lengthpercore;j++)
//        {
//            int counter = 0;
//            int diff1= 46-73+1;
//            int randstart = rand()%(diff1);
//            int assignable=0;
//            while(counter<diff1 && assignable!=1)
//            {
//                new_index = local_en_indexes_op[poly_lattice_indexes[i*lengthpercore*pertrusionpercore+j-1]*120+(randstart+counter)%diff1+74];

//                if(new_index<(L1dim_S_xy*L1dim_S_xy*(L1dim_L_z-3)))
//                {
//                    assignable = latticepoint_assign2(new_index,1);
//                }
//                if(assignable)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "assigned 1 %i  %i\n",new_index,assignable);
//                    fclose(stdoutlog);
//                }

//                counter++;
//            }

//            int diff2= 74-93+1;
//            counter = 0;
//            randstart = rand()%(diff2);
//            while(counter<diff2 && assignable!=1)
//            {
//                new_index = local_en_indexes_op[poly_lattice_indexes[i*lengthpercore*pertrusionpercore+j-1]*120+(randstart+counter)%diff2+74];

//                if(new_index>(L1dim_S_xy*L1dim_S_xy*10))
//                {
//                    if(new_index<(L1dim_S_xy*L1dim_S_xy*(L1dim_L_z-3)))
//                    {
//                           assignable = latticepoint_assign2(new_index,1);
//                    }
//                }
//                else
//                {
//                    assignable = 0;
//                }

//                if(assignable)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "assigned 2 %i %i\n",new_index,assignable);
//                    fclose(stdoutlog);
//                }
//                counter++;
//            }

//            int diff3= 94-119+1;
//            randstart = rand()%(diff3);
//            while(counter<diff3 && assignable!=1)
//            {
//                new_index = local_en_indexes_op[poly_lattice_indexes[i*lengthpercore*pertrusionpercore+j-1]*120+(randstart+counter)%diff3+94];
//                if(new_index>
//                (L1dim_S_xy*L1dim_S_xy*10))
//                {
//                    if(new_index<
//                    (L1dim_S_xy*L1dim_S_xy*(L1dim_L_z-3)))
//                    {
//                    assignable = latticepoint_assign2(new_index,1);
//                     }
//                }
//                else{
//                    assignable = 0;
//                }

//                if(assignable) {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "assigned 3 %i  %i\n",new_index,assignable);
//                    fclose(stdoutlog);
//                }
//                counter++;
//            }

//            int diff4= 0-45+1;
//            counter = 0;
//            randstart = rand()%(diff4);
//            while(counter<diff4 && assignable!=1)
//            {
//                new_index = local_en_indexes_op[poly_lattice_indexes[i*lengthpercore*pertrusionpercore+j-1]*120+(randstart+counter)%diff4+0];
//                assignable = latticepoint_assign2(new_index,1);
//                if(assignable) {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "assigned 4 %i  %i\n",new_index,assignable);
//                    fclose(stdoutlog);
//                }
//                counter++;
//            }

//            if(assignable)
//            {
//                latticepoint_assign(new_index,1);
//                poly_lattice_indexes[(i*lengthpercore)+j]=new_index;
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "guess assign %i,%i,%i\n",i,j,new_index);
//                fclose(stdoutlog);
//            }
//            else{


//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "couldnt assign %i,%i,%i\n",i,j,new_index);
//                fclose(stdoutlog);

//                while(assignable==0)
//                {
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "trying again %i  %i\n",i,j);
//                    fclose(stdoutlog);

//                    poly_latticepoint_clear(poly_lattice_indexes, i*lengthpercore*pertrusionpercore, (1+i)*lengthpercore*pertrusionpercore); // clears all of the lattice points between the two axis points
//                    stdoutlog = fopen(stdoutlogname, "a");
//                    fprintf(stdoutlog, "ta base %i\n",poly_lattice_indexes[i*pertrusionpercore*lengthpercore]);
//                    fclose(stdoutlog);
//                    for(int k=1;k<=j;k++)
//                    {
//                        int counter = 0;
//                        int randstart = rand()%(120);
//                        assignable=0;
//                        while(counter<120 && assignable!=1)
//                        {
//                            new_index = local_en_indexes_op[poly_lattice_indexes[i*lengthpercore*pertrusionpercore+k-1]*120+(randstart+counter)%120+0];
//                            if(new_index< (L1dim_S_xy*L1dim_S_xy*(L1dim_L_z-3)))
//                            {
//                                if(new_index>L1dim_S_xy*L1dim_S_xy*6)
//                                {
//                                    assignable = latticepoint_assign2(new_index,1);
//                                }
//                            }
//                            if(assignable) {
//                                stdoutlog = fopen(stdoutlogname, "a");
//                                fprintf(stdoutlog, "assigned ta %i  %i\n",new_index,assignable);
//                                fclose(stdoutlog);
//                            }
//                            counter++;
//                        }

//                        if(assignable)
//                        {
//                            latticepoint_assign(new_index,1);
//                            poly_lattice_indexes[(i*lengthpercore)+k]=new_index;

//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "TA guess assign %i,%i,%i\n",i,k,new_index);
//                            fclose(stdoutlog);
//                        }
//                        else{
//                         stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "restart ...\n");
//                            fclose(stdoutlog);

//                            break;
//                        }
//                    }
//                }




//                //MPI_Abort(MPI_COMM_WORLD,1);
//            }
//        }
//    }


    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice point init poly_core()\n");
    fclose(stdoutlog);
    
    // writes the lattice polymer data
    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
    {
        int chain_position = i % lengthpercore;

        lattice_polymer[3 * poly_lattice_indexes[i]] = poly_lattice_indexes[i];

        if (chain_position == 0) // core definition needs to be removed or ignored if there are multiple petrusions from a single core
        {
            lattice_polymer[3 * poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 * poly_lattice_indexes[i] + 2] = -100;
        }
        else if (chain_position == (lengthpercore - 1)) // head of the poly chain
        {
            lattice_polymer[3 * poly_lattice_indexes[i] + 1] = -200;
            lattice_polymer[3 * poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
        else { // aqny intermediate monomer
            lattice_polymer[3 * poly_lattice_indexes[i] + 1] = poly_lattice_indexes[i + 1];
            lattice_polymer[3 * poly_lattice_indexes[i] + 2] = poly_lattice_indexes[i - 1];
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Populated lattice polymer init poly_core()\n");
    fclose(stdoutlog);
    
    // fills poly_lattice_coordinates
    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
    {
        poly_coord_decompose(poly_lattice_indexes[i], i);
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "plc:(%i,%i,%i)\n",poly_lattice_coordinates[3*i],poly_lattice_coordinates[3*i+1],poly_lattice_coordinates[3*i+2]);
//        fclose(stdoutlog);
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "assigned coord init poly_core()\n");
    fclose(stdoutlog);
    
    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
    {
        ph_poly_lattice_indexes[i] = poly_lattice_indexes[i];
        ph_poly_lattice_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
        ph_poly_lattice_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
        ph_poly_lattice_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "created place holders init poly_core()\n");
    fclose(stdoutlog);
    
    init_en_and_dist_array();
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Called distance array init poly_core()\n");
    fclose(stdoutlog);
    
    init_solvent(no_solvent_sites);
    
    optomize();
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Called init solvent()\n");
    fclose(stdoutlog);
    
    for(int i=0;i<hist_size;i++)
    {
        rog[i]=0;
        rog_buf[i]=0;
        rog_z[i]=0;
        rog_z_buf[i]=0;
        rog_xy[i]=0;
        rog_xy_buf[i]=0;
        visitdp[i]=0;
        visits[i]=0;
        visits_buf[i]=0;
        endtoend[i]=0;
        endtoend_buf[i]=0;
    }
    
}

void poly_coord_decompose(int reference_index, int polymer_position)
{
    int x = reference_index % L1dim_S_xy; // column
    int y = (reference_index % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    int z = (reference_index / (L1dim_S_xy * L1dim_S_xy)); // depth

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "\tri:%i  z:%i\n",reference_index,z);
//    fclose(stdoutlog);
    
    poly_lattice_coordinates[3 * polymer_position] = x;
    poly_lattice_coordinates[3 * polymer_position + 1] = y;
    poly_lattice_coordinates[3 * polymer_position + 2] = z;
}

double point_distance(int reference_index_1, int reference_index_2)
{
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "Made it to point_distance() reference_index_1:%i reference_index_2:%i\n",reference_index_1,reference_index_2);
//    fclose(stdoutlog);
    
    float x_1 = reference_index_1 % L1dim_S_xy; // column
    float y_1 = (reference_index_1 % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    float z_1 = (reference_index_1 / (L1dim_S_xy * L1dim_S_xy)); // depth

    float x_2 = reference_index_2 % L1dim_S_xy; // column
    float y_2 = (reference_index_2 % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    float z_2 = (reference_index_2 / (L1dim_S_xy * L1dim_S_xy));// depth

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "after dz:%f \t dy: %f \t dz: %f \n",(x_1 - x_2),(y_1 - y_2),(z_1 - z_2));
    
    int dx = abs((x_1 - x_2) - anint((x_1 - x_2), L1dim_S_xy) * L1dim_S_xy);
    int dy = abs((y_1 - y_2) - anint((y_1 - y_2), L1dim_S_xy) * L1dim_S_xy);
    int dz = abs((z_1 - z_2) - anint((z_1 - z_2), L1dim_S_z) * L1dim_S_z);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,sqrt(dx * dx + dy * dy + dz * dz));
//    fclose(stdoutlog);
    
    return(abs(sqrt(dx * dx + dy * dy + dz * dz)));
}

int int_point_distance(int reference_index_1, int reference_index_2)
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to point_distance() reference_index_1:%i reference_index_2:%i\n",reference_index_1,reference_index_2);
    //fclose(stdoutlog);
    
    int x_1 = reference_index_1 % L1dim_S_xy; // column
    int y_1 = (reference_index_1 % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    int z_1 = (reference_index_1 / (L1dim_S_xy * L1dim_S_xy)); // depth

    int x_2 = reference_index_2 % L1dim_S_xy; // column
    int y_2 = (reference_index_2 % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    int z_2 = (reference_index_2 / (L1dim_S_xy * L1dim_S_xy));// depth

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i \n",(x_1 - x_2),(y_1 - y_2),(z_1 - z_2));
    
    int dx = abs((x_1 - x_2) - anint((x_1 - x_2), L1dim_S_xy) * L1dim_S_xy);
    int dy = abs((y_1 - y_2) - anint((y_1 - y_2), L1dim_S_xy) * L1dim_S_xy);
    int dz = abs((z_1 - z_2) - anint((z_1 - z_2), L1dim_S_z) * L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,sqrt(dx * dx + dy * dy + dz * dz));
    //fclose(stdoutlog);
    
    return(dx * dx + dy * dy + dz * dz);
}

double point_distance(int reference_index_1, int reference_index_2,int(&arr)[3])
{
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "\t\t\tMade it to point_distance() reference_index_1:%i reference_index_2:%i\n",reference_index_1,reference_index_2);
//    fclose(stdoutlog);
    
    int x_1 = reference_index_1 % L1dim_S_xy; // column
    int y_1 = (reference_index_1 % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    int z_1 = (reference_index_1 / (L1dim_S_xy * L1dim_S_xy)); // depth

    int x_2 = reference_index_2 % L1dim_S_xy; // column
    int y_2 = (reference_index_2 % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    int z_2 = (reference_index_2 / (L1dim_S_xy * L1dim_S_xy));// depth

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "\t\t\tafter dz:%i \t dy: %i \t dz: %i\n",(x_1 - x_2),(y_1 - y_2),(z_1 - z_2));
//    fclose(stdoutlog);
    
    int dx = abs((x_1 - x_2) - anint((x_1 - x_2), L1dim_S_xy) * L1dim_S_xy);
    int dy = abs((y_1 - y_2) - anint((y_1 - y_2), L1dim_S_xy) * L1dim_S_xy);
    int dz = abs((z_1 - z_2) - anint((z_1 - z_2), L1dim_S_z) * L1dim_S_z);

//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "\t\t\tafter dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,sqrt(dx * dx + dy * dy + dz * dz));
//    fclose(stdoutlog);
    
    arr[0]=0;
    arr[1]=0;
    arr[2]=0;
    
    if(dx<4&&dy<4&&dz<4)
    {
        arr[0]=abs(dx);
        arr[1]=abs(dy);
        arr[2]=abs(dz);
    }
    
    return 0;
}

double poly_point_distance(int* poly_lattice_coordinates,int polymer_position_1, int polymer_position_2)
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to poly_point_distance()\n");
    //fclose(stdoutlog);
    
    int dx = poly_lattice_coordinates[3 * polymer_position_1] - poly_lattice_coordinates[3 * polymer_position_2];
    
    int dy = poly_lattice_coordinates[3 * polymer_position_1 + 1] - poly_lattice_coordinates[3 * polymer_position_2 + 1];
    
    int dz = poly_lattice_coordinates[3 * polymer_position_1 + 2] - poly_lattice_coordinates[3 * polymer_position_2 + 2];
    
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "befor dz:%i \t dy: %i \t dz: %i\n",dx,dy,dz);
    //fclose(stdoutlog);
    
    dx = abs(dx - anint(dx, L1dim_S_xy) * L1dim_S_xy);
    dy = abs(dy - anint(dy, L1dim_S_xy) * L1dim_S_xy);
    dz = abs(dz - anint(dz, L1dim_S_z) * L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,dist_array[dx][dy][dz]);
    //fclose(stdoutlog);
    
    if (dx < 4 && dy < 4 && dz < 4)
    {
        return dist_array[dx][dy][dz];
    }

    return dist_array[0][0][0];
}

double poly_point_distance(int* poly_lattice_coordinates,int polymer_position_1, int polymer_position_2,int(&arr)[3])
{
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "Made it to poly_point_distance()\n");
    //fclose(stdoutlog);
    
    int dx = poly_lattice_coordinates[3 * polymer_position_1] - poly_lattice_coordinates[3 * polymer_position_2];
    
    int dy = poly_lattice_coordinates[3 * polymer_position_1 + 1] - poly_lattice_coordinates[3 * polymer_position_2 + 1];

    int dz = poly_lattice_coordinates[3 * polymer_position_1 + 2] - poly_lattice_coordinates[3 * polymer_position_2 + 2];
    
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "before dz:%i \t dy: %i \t dz: %i \n",dx,dy,dz);
    //fclose(stdoutlog);
    
    dx = abs(dx - anint(dx, L1dim_S_xy) * L1dim_S_xy);
    dy = abs(dy - anint(dy, L1dim_S_xy) * L1dim_S_xy);
    dz = abs(dz - anint(dz, L1dim_S_z) * L1dim_S_z);

    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "after dz:%i \t dy: %i \t dz: %i    dist_check:%f\n",dx,dy,dz,dist_array[dx][dy][dz]);
    //fclose(stdoutlog);
    
    if (dx < 4 && dy < 4 && dz < 4)
    {
        arr[0]=dx;
        arr[1]=dy;
        arr[2]=dz;
        return dist_array[dx][dy][dz];
    }

    return dist_array[0][0][0];
}

void init_en_and_dist_array()
{
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "Made it to en array and distance()\n");
    fclose(stdoutlog);
    
    for (int i = 0; i < 11; i++)
    {
    int_en_array_pp[i]=0;
        int_en_array_ps[i]=0;
        int_en_array_ss[i]=0;
        int_en_array[i]=0;
        int_en_array_bond[i]=0;
    }
    
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
            
                en_array_pp[i][j][k] =0;
                en_array_ps[i][j][k] =0;
                en_array_ss[i][j][k] =0;
                en_array[i][j][k] =0;
                en_array_bond[i][j][k] =0;
                en_array_wall[i][j][k] =0;
                
                if (i * i + j * j + k * k == 4)
                {
                    en_array[i][j][k] =(int) 303.974;
                    en_array_bond[i][j][k] =(int) 1.5*2028.01;
                    en_array_wall[i][j][k] =(int) -95.7;
                    
                    en_array[i][j][k] =(int) -10;
                    en_array_bond[i][j][k] =1.5*235.;
                    en_array_wall[i][j][k] =(int) -9;
                    
                    int_en_array[i * i + j * j + k * k] = en_array[i][j][k];
                    int_en_array_bond[i * i + j * j + k * k] = en_array_bond[i][j][k];
                    
                                              en_array_pp[i][j][k] =-10;
                                             en_array_ps[i][j][k] =-13;
                                             en_array_ss[i][j][k] =-9;
                                             
                                                 int_en_array_pp[i * i + j * j + k * k]=en_array_pp[i][j][k];
                                              int_en_array_ps[i * i + j * j + k * k]=en_array_ps[i][j][k];
                                             int_en_array_ss[i * i + j * j + k * k]= en_array_ss[i][j][k];
                    
                   // en_array[i][j][k] =(int) 300;
                    //en_array_bond[i][j][k] =(int) 2000;
                    //en_array_wall[i][j][k] =(int) -100;
                }
                if (i * i + j * j + k * k == 5)
                {
                    en_array[i][j][k] = (int)16.3753;
                    en_array_bond[i][j][k] =(int) 2355.75;
                    en_array_wall[i][j][k] =(int) -75.9574;
                    
                    en_array[i][j][k] = (int)-8;
                    en_array_bond[i][j][k] = 1.5*315.;
                    en_array_wall[i][j][k] =(int) -8;
                    
                    
                    int_en_array[i * i + j * j + k * k] = en_array[i][j][k];
                    int_en_array_bond[i * i + j * j + k * k] = en_array_bond[i][j][k];
                    
                                         en_array_pp[i][j][k] =-8;
                                              en_array_ps[i][j][k] =-10;
                                              en_array_ss[i][j][k] =-7;
                                              
                                                  int_en_array_pp[i * i + j * j + k * k]=en_array_pp[i][j][k];
                                               int_en_array_ps[i * i + j * j + k * k]=en_array_ps[i][j][k];
                                              int_en_array_ss[i * i + j * j + k * k]= en_array_ss[i][j][k];
                    
                    //en_array[i][j][k] = (int) 20;
                    //en_array_bond[i][j][k] =(int) 2400;
                    //en_array_wall[i][j][k] =(int) -100;
                }
                if (i * i + j * j + k * k == 6)
                {
                    en_array[i][j][k] = (int) 3.50123;
                    en_array_bond[i][j][k] =(int) 3095.98;
                    en_array_wall[i][j][k] =(int) -60.2979;
                    
                    en_array[i][j][k] = (int) -5;
                    en_array_bond[i][j][k] = 1.5*410.;
                    en_array_wall[i][j][k] =(int) -6;
                    
                    int_en_array[i * i + j * j + k * k] = en_array[i][j][k];
                    int_en_array_bond[i * i + j * j + k * k] = en_array_bond[i][j][k];
                    
                                 en_array_pp[i][j][k] =-5;
                               en_array_ps[i][j][k] =-6;
                               en_array_ss[i][j][k] =-5;
                               
                                   int_en_array_pp[i * i + j * j + k * k]=en_array_pp[i][j][k];
                                int_en_array_ps[i * i + j * j + k * k]=en_array_ps[i][j][k];
                               int_en_array_ss[i * i + j * j + k * k]= en_array_ss[i][j][k];
                    
                    //en_array[i][j][k] = (int)10;
                    //en_array_bond[i][j][k] =(int) 3100;
                    //en_array_wall[i][j][k] =(int) -10;
                }
                if (i * i + j * j + k * k == 8)
                {
                    en_array[i][j][k] = (int) 43.1734;
                    en_array_bond[i][j][k] = (int)5475.03;
                    en_array_wall[i][j][k] =(int) - 40.462;
                    
                    en_array[i][j][k] = (int) -3;
                    en_array_bond[i][j][k] =  1.5*676.;
                    en_array_wall[i][j][k] =(int) -4;
                    
                    int_en_array[i * i + j * j + k * k] = en_array[i][j][k];
                    int_en_array_bond[i * i + j * j + k * k] = en_array_bond[i][j][k];
                    
                          en_array_pp[i][j][k] =-2;
                                     en_array_ps[i][j][k] =-3;
                                     en_array_ss[i][j][k] =-2;
                                     
                                         int_en_array_pp[i * i + j * j + k * k]=en_array_pp[i][j][k];
                                      int_en_array_ps[i * i + j * j + k * k]=en_array_ps[i][j][k];
                                     int_en_array_ss[i * i + j * j + k * k]= en_array_ss[i][j][k];
                    
                    //en_array[i][j][k] = (int) 50;
                    //en_array_bond[i][j][k] = (int)5500;
                    //en_array_wall[i][j][k] =(int) - 50;
                }
                if (i * i + j * j + k * k == 9)
                {
                    en_array[i][j][k] = (int) 57.6302;
                    en_array_bond[i][j][k] =(int) 7828.85;
                    en_array_wall[i][j][k] =(int) - 34.1454;
                    
                    en_array[i][j][k] = (int)-2;
                    en_array_bond[i][j][k] = 1.5*887.;
                    en_array_wall[i][j][k] =(int) -3;
                    
                    int_en_array[i * i + j * j + k * k] = en_array[i][j][k];
                    int_en_array_bond[i * i + j * j + k * k] = en_array_bond[i][j][k];
                    
                        en_array_pp[i][j][k] =-2;
                                     en_array_ps[i][j][k] =-2;
                                     en_array_ss[i][j][k] =-2;
                                     
                                         int_en_array_pp[i * i + j * j + k * k]=en_array_pp[i][j][k];
                                      int_en_array_ps[i * i + j * j + k * k]=en_array_ps[i][j][k];
                                     int_en_array_ss[i * i + j * j + k * k]= en_array_ss[i][j][k];
                    
                    //en_array[i][j][k] = (int)50;
                    //en_array_bond[i][j][k] =(int) 7800;
                    //en_array_wall[i][j][k] =(int) - 50;
                }

            
            if (i * i + j * j + k * k == 10)
                {
                    en_array[i][j][k] = (int) 57.6302;
                    en_array_bond[i][j][k] =(int) 7828.85;
                    en_array_wall[i][j][k] =(int) - 34.1454;
                    
                    en_array[i][j][k] = (int)-2;
                    en_array_bond[i][j][k] = 1.5*1247.;
                    en_array_wall[i][j][k] =(int) -3;
                    
                    int_en_array[i * i + j * j + k * k] = en_array[i][j][k];
                    int_en_array_bond[i * i + j * j + k * k] = en_array_bond[i][j][k];
                    
                        en_array_pp[i][j][k] =-1;
                                     en_array_ps[i][j][k] =-2;
                                     en_array_ss[i][j][k] =-1;
                                     
                                         int_en_array_pp[i * i + j * j + k * k]=en_array_pp[i][j][k];
                                      int_en_array_ps[i * i + j * j + k * k]=en_array_ps[i][j][k];
                                     int_en_array_ss[i * i + j * j + k * k]= en_array_ss[i][j][k];
                    
                    //en_array[i][j][k] = (int)50;
                    //en_array_bond[i][j][k] =(int) 7800;
                    //en_array_wall[i][j][k] =(int) - 50;
                }

            }
        }
    }

    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "en array created()\n");
    fclose(stdoutlog);
    
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                dist_array[i][j][k] = INT_MAX;
                dist_array[0][0][0] = INT_MAX;
                if (i * i + j * j + k * k < 11 && i * i + j * j + k * k >=4)
                {
                    dist_array[i][j][k] = sqrt(i*i + j*j + k*k);
                }

            }
        }
    }
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "dist array created()\n");
    fclose(stdoutlog);
}


int general_coord(int(&arr)[3], int reference_index) // output the x y and z compenent for any index
{
    float x = reference_index % L1dim_S_xy; // column
    float y = (reference_index % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
    float z = (reference_index / (L1dim_S_xy * L1dim_S_xy)); // depth

    arr[0] = x;
    arr[1] = y;
    arr[2] = z;

    return 0;
}

int general_index(int(&arr)[3]) // output index from any x , y , z components
{

    while (arr[0] < 0) {arr[0] = arr[0] + L1dim_S_xy;}
    while (arr[0] >= L1dim_S_xy) {arr[0] = arr[0] - L1dim_S_xy;}

    while (arr[1] < 0) {arr[1] = arr[1] + L1dim_S_xy;}
    while (arr[1] >= L1dim_S_xy) {arr[1] = arr[1] - L1dim_S_xy;}
    
    while (arr[2] < 0) {arr[2] = arr[2] + L1dim_S_z;}
    while (arr[2] >= L1dim_S_z) {arr[2] = arr[2] - L1dim_S_z;}

    return (arr[0] + arr[1] * L1dim_S_xy + arr[2] * L1dim_S_xy * L1dim_S_xy);
}

int general_index_op(int(&arr)[3]) // output index from any x , y , z components
{
    arr[0] =index_optomize_x[arr[0]+neg_xy];
    arr[1] =index_optomize_y[arr[1]+neg_xy];
    arr[2] =index_optomize_z[arr[2]+neg_z];
    return(calc_index_optomize_x[arr[0]] + calc_index_optomize_y[arr[1]] +calc_index_optomize_z[arr[2]]);
}

int poly_latticepoint_clear(int* c_chain, int start, int end) // set a range of monomers between two to zero for unimpeded movement
{
    for (int i = start + 1; i < end; i++)
    {
        
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "Made it to poly_clear()  start: %i  end: %i\n",start,end);
//        fprintf(stdoutlog, "c_chain: %i  latticepoint_S[c_chain[i]]: %i\t",c_chain[i],latticepoint_S[c_chain[i]]);
//        fprintf(stdoutlog, "neighbor_L[numberneighbors_L * c_chain[i] + iii]:");
//        fclose(stdoutlog);
        
        
        latticepoint_S[c_chain[i]] = 0;
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "c_chain: %i  %i\t",c_chain[i],latticepoint_L[neighbor_L[numberneighbors_L * c_chain[i] + iii]]);
//            fclose(stdoutlog);
            
            
            latticepoint_L[neighbor_L[numberneighbors_L * c_chain[i] + iii]] = 0;
        }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "\n");
//        fclose(stdoutlog);
    }
    
    return 0;
}
                    
int poly_latticepoint_clear(int index) // set a range of monomers between two to zero for unimpeded movement
{
    latticepoint_S[index] = 0;

    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        latticepoint_L[neighbor_L[numberneighbors_L * index + iii]] = 0;
    }
    
    return 0;
}

int latticepoint_assign(int index , int spin)
{
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "checking assign index: %i    %i\n",index,latticepoint_S[index]);
//    fclose(stdoutlog);
    
    if (latticepoint_S[index] == 0)
    {
        latticepoint_S[index] = spin;
    }
    else
    {
        return 0;
    }
    
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "checking assign index: %i    %i\n",index,latticepoint_S[index]);
//    fclose(stdoutlog);
    
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        
        if (latticepoint_L[neighbor_L[numberneighbors_L * index + iii]] == 0)
        {
            int x= latticepoint_L[neighbor_L[numberneighbors_L * index + iii]];
            
            latticepoint_L[neighbor_L[numberneighbors_L * index + iii]]= spin;
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "ni:%i\t%i\t%i\n",neighbor_L[numberneighbors_L * index + iii],x,latticepoint_L[neighbor_L[numberneighbors_L * index + iii]]);
//            fclose(stdoutlog);
            
        }
        else // resets the points switched to occupaton in the case where occupation is not allowed. Just performs cleanup itself so no other function has to
        {
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "checking assign index: %i  iii: %i\n",index,iii);
//            fprintf(stdoutlog, "b latticepoint_S[index]: %i\n",latticepoint_S[index]);
//            fclose(stdoutlog);
            
            latticepoint_S[index] = 0;
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "a latticepoint_S[index]: %i\n\n",latticepoint_S[index]);
//            fclose(stdoutlog);
            
            for (int j = 0; j < iii; j++)
            {
                int x= latticepoint_L[neighbor_L[numberneighbors_L * index + j]];
                
                latticepoint_L[neighbor_L[numberneighbors_L * index + j]] = 0;
                
                
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "%i\t%i\n",x,latticepoint_L[neighbor_L[numberneighbors_L * index + j]]);
//                fclose(stdoutlog);
                
            }
            
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "\n");
//            fclose(stdoutlog);
            
            return 0;
        }
    }

    return 1;
}

int latticepoint_assign2(int index , int spin)
{
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "checking assign index: %i    %i\n",index,latticepoint_S[index]);
//    fclose(stdoutlog);
    
    if (latticepoint_S[index] != 0)
    {
        return 0;
    }
    
    
//    stdoutlog = fopen(stdoutlogname, "a");
//    fprintf(stdoutlog, "checking assign index: %i    %i\n",index,latticepoint_S[index]);
//    fclose(stdoutlog);
    
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        
        if (latticepoint_L[neighbor_L[numberneighbors_L * index + iii]] != 0)
        {
            return 0;
        }
    }

    return 1;
}



void lattice_polymer_reset(int* p_chain, int* c_chain, int start, int end) //retrofit of reorganize function from previous code saved some time rewriting parts of the function
{
    for (int i = start + 1; i < end; i++) // clears the lattice polymer in the range set
    {
        lattice_polymer[3 * p_chain[i]] = p_chain[i];
        lattice_polymer[3 * p_chain[i] + 1] = -1; // next lattice point on chain (if not applicable -1 assigned)
        lattice_polymer[3 * p_chain[i] + 2] = -1;
    }
    //stdoutlog = fopen(stdoutlogname, "a");
    //fprintf(stdoutlog, "qq\n");
    //fclose(stdoutlog);
    for (int i = start + 1; i < end; i++)
    {

        lattice_polymer[3 * p_chain[i]] = p_chain[i];

        if (i % lengthpercore == 0) //core definition
        {
            lattice_polymer[3 * c_chain[i] + 1] = c_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
            lattice_polymer[3 * c_chain[i] + 2] = -100; //previous lattice point on chain (if applicable -1 assigned)
        }
        if (i % lengthpercore != 0 && i % (lengthpercore) != lengthpercore - 1) // intermediate monomer
        {
            lattice_polymer[3 * c_chain[i] + 1] = c_chain[i + 1]; // next lattice point on chain (if not applicable -1 assigned)
            lattice_polymer[3 * c_chain[i] + 2] = c_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
        }
        if (i % (lengthpercore) == lengthpercore - 1) // length of an arm
        {
            lattice_polymer[3 * c_chain[i] + 1] = -200; // next lattice point on chain (if not applicable -1 assigned)
            lattice_polymer[3 * c_chain[i] + 2] = c_chain[i - 1]; //previous lattice point on chain (if applicable -1 assigned)
        }
    }
}

int lattice_poly_index_reset(int* p_chain, int* c_chain, int start, int end) // reorganizes the latticepoints more efficient than previous scrubbing over entire lattice
{
    for (int i = start + 1; i < end; i++) // clears the lattice polymer in the range set
    {
        latticepoint_S[p_chain[i]] = 0;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L * p_chain[i] + iii]] = 0;
        }
    }

    for (int i = start + 1; i < end; i++) // clears the lattice polymer in the range set
    {
        latticepoint_S[c_chain[i]] = 1;

        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            latticepoint_L[neighbor_L[numberneighbors_L * c_chain[i] + iii]] = 1;
        }

        p_chain[i] = c_chain[i];
    }



    return 0;
}

int poly_coord_reset(int* p_chain, int* c_chain, int start, int end)
{
    for (int i = start + 1; i < end; i++) // clears the lattice polymer in the range set
    {
        p_chain[i * 3 + 0] = c_chain[i * 3 + 0];
        p_chain[i * 3 + 1] = c_chain[i * 3 + 1];
        p_chain[i * 3 + 2] = c_chain[i * 3 + 2];
    }
    return 0;
}


int init_solvent_sub(int index)
{
    for (int iii = 0; iii < numberneighbors_L; iii++)
    {
        if (latticepoint_L[neighbor_L[numberneighbors_L * index + iii]] == 0) // checks to ensure lattice is assignable
        {
            latticepoint_L[neighbor_L[numberneighbors_L * index + iii]]= 2;
        }
        else // resets the points switched to occupaton in the case where occupation is not allowed. Just performs cleanup itself so no other function has to
        {
            latticepoint_S[index] = 0;
            for (int j = 0; j < iii; j++)
            {
                latticepoint_L[neighbor_L[numberneighbors_L * index + j]] = 0;
            }
            return 0;
        }
    }
    return 1;
}

int init_solvent(int number) // randomly instance solvent atoms
{
    

    if(myid<minmaxid)
    {
    if (fopen("Minimum_Config_Solvent.txt","r") == NULL)
    {
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt %i \n", myid,numberspins_S);
        fclose(stdoutlog);
        
        int counter=0;// number of solvent placed
        for(int index =numberspins_S-1-L1dim_S_xy-L1dim_S_xy*3;index>=0;index--) // sequences through lattice
        {
        
        
            if(counter<number && index>=((upper_bounds)) && index<(L1dim_S_xy*L1dim_S_xy*L1dim_S_z-L1dim_S_xy*L1dim_S_xy))
            {// while all solvent has not been placed and while the solvent is within a given range
                
                if (latticepoint_S[index] == 0) // checks for if assignment is possible
                {
                
                    latticepoint_S[index] = 2; // assigns solvent spin
                    
                    int neighbor_check = 0; // checks if neighbor assignment is possible
                    neighbor_check =init_solvent_sub(index); // calls sub function to reduce size of function
                    
                    if(neighbor_check == 1) // if assignment was possible
                    {
                        solvent_loc[counter] = index; // assigns the spin to solvent_loc and increments
                        ph_solvent_loc[counter] = index;
                        wl_solvent_loc[counter] = index;
                        counter++;
                    }
                }
            }
            
        }
        return 0;
    }
    else // rewrite of prior file reading system to be more portable
    {
        FILE* fptr;
        int value;
        int valuee;
        
        fptr = fopen("Minimum_Config_Solvent.txt", "r");

        while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
        {
            solvent_loc[value] = valuee;
            ph_solvent_loc[value] = valuee;
            wl_solvent_loc[value] = valuee;
           
        }

        fclose(fptr);
    }
    
        
    reset_indexes[0]=-1;
    reset_indexes[1]=no_solvent_sites;
        
//         stdoutlog = fopen(stdoutlogname, "a");
//         for (int i = 0; i < no_solvent_sites; i++)
//         {
//             fprintf(stdoutlog, "0 i: %i\t%4i \t %i\n", i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
//         }
//         fclose(stdoutlog);
    
   
    solvent_reset(solvent_loc,ph_solvent_loc);
//        stdoutlog = fopen(stdoutlogname, "a");
//        for (int i = 0; i < no_solvent_sites; i++)
//        {
//            fprintf(stdoutlog, "1 i: %i\t%4i \t %i\n", i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
//        }
//        fclose(stdoutlog);
    solvent_reset(wl_solvent_loc,solvent_loc);
//        for (int i = 0; i < no_solvent_sites; i++)
//        {
//            fprintf(stdoutlog, "2 i: %i\t%4i \t %i\n", i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
//        }
//        fclose(stdoutlog);
return 0;
    }
    
    if(myid>=minmaxid)
       {
       if (fopen("Maximum_Config_Solvent.txt","r") == NULL)
       {
           stdoutlog = fopen(stdoutlogname, "a");
           fprintf(stdoutlog, "\nProc %i: Can't find file Access.txt \n", myid);
           fclose(stdoutlog);
           
           int counter=0;// number of solvent placed
           for(int index =numberspins_S-1-L1dim_S_xy-L1dim_S_xy*3;index>=0;index--) // sequences through lattice
        {
               if(counter<number && index>=((upper_bounds)) && index<(L1dim_S_xy*L1dim_S_xy*L1dim_S_z-L1dim_S_xy*L1dim_S_xy))
               {// while all solvent has not been placed and while the solvent is within a given range
                   
                   if (latticepoint_S[index] == 0) // checks for if assignment is possible
                   {
                       latticepoint_S[index] = 2; // assigns solvent spin
                       
                       int neighbor_check = 0; // checks if neighbor assignment is possible
                       neighbor_check =init_solvent_sub(index); // calls sub function to reduce size of function
                       
                       if(neighbor_check == 1) // if assignment was possible
                       {
                           solvent_loc[counter] = index; // assigns the spin to solvent_loc and increments
                           ph_solvent_loc[counter] = index;
                           wl_solvent_loc[counter] = index;
                           counter++;
                       }
                   }
               }
               
           }
           return 0;
       }
       else // rewrite of prior file reading system to be more portable
       {
           FILE* fptr;
           int value;
           int valuee;
           
           fptr = fopen("Maximum_Config_Solvent.txt", "r");

           while (fscanf(fptr, "%i\t%i\n", &value, &valuee) > 0)
           {
               solvent_loc[value] = valuee;
               ph_solvent_loc[value] = valuee;
               wl_solvent_loc[value] = valuee;
              
           }

           fclose(fptr);
       }
       
           reset_indexes[0]=-1;
        reset_indexes[1]=no_solvent_sites;
            
//             stdoutlog = fopen(stdoutlogname, "a");
//             for (int i = 0; i < no_solvent_sites; i++)
//             {
//                 fprintf(stdoutlog, "0 i: %i\t%4i \t %i\n", i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
//             }
//             fclose(stdoutlog);
        
       
        solvent_reset(solvent_loc,ph_solvent_loc);
//            stdoutlog = fopen(stdoutlogname, "a");
//            for (int i = 0; i < no_solvent_sites; i++)
//            {
//                fprintf(stdoutlog, "1 i: %i\t%4i \t %i\n", i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
//            }
//            fclose(stdoutlog);
        solvent_reset(wl_solvent_loc,solvent_loc);
//            for (int i = 0; i < no_solvent_sites; i++)
//            {
//                fprintf(stdoutlog, "2 i: %i\t%4i \t %i\n", i,solvent_loc[i],latticepoint_S[solvent_loc[i]]);
//            }
//            fclose(stdoutlog);
return 0;
       }
    return 0;
}

int reduced_lattice_optimization(int index_1,int index_2)
{
    const int reduced_optomize_value =  (L1dim_S_xy*L1dim_S_xy*3);
    int shifted_index_1 = lattice_dist_optomize_2[index_1];
    int index_1_difference = index_1 - shifted_index_1;
    int shifted_index_2 = index_2 - index_1_difference;
    
    if(shifted_index_2>=0 && shifted_index_2<7*L1dim_S_xy*L1dim_S_xy)
    {
        return lattice_dist_optomize_2_lattice[(shifted_index_1-(L1dim_S_xy*L1dim_S_xy*3))*(L1dim_S_xy * L1dim_S_xy*7)+shifted_index_2];
    }
        
    return INT_MAX;
}


void init_optomize_energy_index()
{
    int index = 1029;
        int c_x;
        int c_y;
        int c_z;
        int counter_2=0;
        c_x = index_x[index];
        c_y = index_y[index];
        c_z = index_z[index];
        
        int hold_indexes[120]={0};
        
    
        for(int a =-3;a<=3;a++)
        {
            for(int b = -3;b<=3;b++)
            {
                for(int c = -3 ; c<=3;c++)
                {
                    int checked_index = xyz_index[(a+c_z+5)*(L1dim_S_xy+10)*(L1dim_S_xy+10)+(b+c_y+5)*(L1dim_S_xy+10)+(c+c_x+5)];
                    
                    if(point_distance(checked_index,index)<max_en_dist && point_distance(checked_index,index)>=2)
                    {
                        local_en_indexes_op_x[counter_2]=c;
                        local_en_indexes_op_y[counter_2]=b;
                        local_en_indexes_op_z[counter_2]=a;
                        hold_indexes[counter_2]=checked_index;
                        counter_2++;
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "%i: {%i,%i,%i} {%i,%i,%i}\n\t",checked_index,c,b,a,local_en_indexes_op_x[counter_2-1],local_en_indexes_op_y[counter_2-1],local_en_indexes_op_z[counter_2-1]);
//                        fclose(stdoutlog);
//
//                        for(int i=0;i<numberneighbors_L;i++)
//                        {
//                            stdoutlog = fopen(stdoutlogname, "a");
//                            fprintf(stdoutlog, "%i  ",neighbor_L[checked_index*numberneighbors_L+i]);
//                            fclose(stdoutlog);
//                        }
//                        stdoutlog = fopen(stdoutlogname, "a");
//                        fprintf(stdoutlog, "\n");
//                        fclose(stdoutlog);
                    }
                }
            }
        }
        
    for(int q =0;q<120;q++)
    {
        for(int r=1;r<numberneighbors_L;r++)
        {
            what_to_exclude[q*numberneighbors_L+r]=120;
            for(int g = 0;g<120;g++)
            {
                if(hold_indexes[g] == neighbor_L[hold_indexes[q]*numberneighbors_L+r])
                {
                    what_to_exclude[q*numberneighbors_L+r]=g;
                }
            }
        }
    }
    
//    for(int q =0;q<96;q++)
//    {
//        stdoutlog = fopen(stdoutlogname, "a");
//                               fprintf(stdoutlog, "%i:  %i\n\t",q,hold_indexes[q]);
//                               fclose(stdoutlog);
//        for(int r=1;r<numberneighbors_L;r++)
//        {
//            stdoutlog = fopen(stdoutlogname, "a");
//            fprintf(stdoutlog, "%i  ",what_to_exclude[q*numberneighbors_L+r]);
//            fclose(stdoutlog);
//        }
//
//        stdoutlog = fopen(stdoutlogname, "a");
//                               fprintf(stdoutlog, "\n");
//                               fclose(stdoutlog);
//    }
//        stdoutlog = fopen(stdoutlogname, "a");
//        fprintf(stdoutlog, "\n\n\n%i\n",counter_2);
//        fclose(stdoutlog);
//
   //MPI_Abort(MPI_COMM_WORLD,1);
}

void optomize_energy_index()
{
    
    for(int y = 0;y<numberspins_S;y++)
    {
        int index = y;
        int c_x;
        int c_y;
        int c_z;
        int counter_2=0;
        c_x = index_x[index];
        c_y = index_y[index];
        c_z = index_z[index];
    
        for(int a =-3;a<=3;a++)
        {
            for(int b = -3;b<=3;b++)
            {
                for(int c = -3 ; c<=3;c++)
                {
                    int checked_index = xyz_index[(a+c_z+5)*(L1dim_S_xy+10)*(L1dim_S_xy+10)+(b+c_y+5)*(L1dim_S_xy+10)+(c+c_x+5)];
                    
                    if(point_distance(checked_index,index)<max_en_dist && point_distance(checked_index,index)>=2)
                    {
                        local_en_indexes_op[index*120+counter_2]=checked_index;
                        local_en_dist[counter_2]=int_point_distance(checked_index,index);
                        counter_2++;
                    }
                }
            }
        }
    }
}

void optomize()
{
    //lattice_dist_optomize = (int*)malloc(L1dim_S_xy * L1dim_S_xy * L1dim_S_z * L1dim_S_xy * L1dim_S_xy * L1dim_S_z * sizeof(int));
    
    index_optomize_x = (int*)malloc((L1dim_S_xy*12) * sizeof(int));
    index_optomize_y = (int*)malloc((L1dim_S_xy*12) * sizeof(int));
    index_optomize_z = (int*)malloc((L1dim_S_z*6) * sizeof(int));
    
    calc_index_optomize_x = (int*)malloc((L1dim_S_xy*12) * sizeof(int));
    calc_index_optomize_y = (int*)malloc((L1dim_S_xy*12) * sizeof(int));
    calc_index_optomize_z = (int*)malloc((L1dim_S_z*6) * sizeof(int));
    
    polymer_optomize = (int*)malloc(L1dim_S_xy * L1dim_S_xy * L1dim_S_z * sizeof(int));
    
//    for(int i=0;i<numberspins_S;i++)
//    {
//        polymer_optomize[i] = i*numberspins_S;
//        for(int j=0;j<numberspins_S;j++)
//        {
//            lattice_dist_optomize[i*numberspins_S +j]=int_point_distance(i,j);
//        }
//    }
    
    lattice_dist_optomize_2 = (int*)malloc(L1dim_S_xy * L1dim_S_xy * L1dim_S_z * sizeof(int));
    for(int j=0;j<numberspins_S;j++)
    {
        float x = j % L1dim_S_xy; // column
        float y = (j % (L1dim_S_xy * L1dim_S_xy)) / L1dim_S_xy; // row
        
        lattice_dist_optomize_2[j]=(x+y*L1dim_S_xy)+L1dim_S_xy*L1dim_S_xy*3;
    }
    
    lattice_dist_optomize_2_lattice = (int*)malloc((L1dim_S_xy*L1dim_S_xy)*(L1dim_S_xy * L1dim_S_xy*7) * sizeof(int));
    for(int i=0+(L1dim_S_xy*L1dim_S_xy*3);i<(L1dim_S_xy*L1dim_S_xy*4);i++)
    {
        for(int j=0;j<(L1dim_S_xy * L1dim_S_xy*7);j++)
        {
            lattice_dist_optomize_2_lattice[(i-(L1dim_S_xy*L1dim_S_xy*3))*(L1dim_S_xy * L1dim_S_xy*7) + j ]=int_point_distance(i,j);
        }
    }
    
    for(int i=0;i<L1dim_S_xy*12;i++)
    {
        index_optomize_x[i] = i%L1dim_S_xy;
        index_optomize_y[i] = i%L1dim_S_xy;
        
        calc_index_optomize_x[i] = (i%L1dim_S_xy);
        calc_index_optomize_y[i] = (i%L1dim_S_xy)*L1dim_S_xy;
    }
    for(int i=0;i<L1dim_S_z*6;i++)
    {
        index_optomize_z[i] = i%L1dim_S_z;
        calc_index_optomize_z[i] = (i%L1dim_S_z)*L1dim_S_xy*L1dim_S_xy;
    }
    
    xyz_index = (int*)malloc((L1dim_S_z+10)*(L1dim_S_xy+10)*(L1dim_S_xy+10) * sizeof(int));
    int xyz_array[3] = {0};
    int x = -5;
    int y = -5;
    int z = -5;
    for(int i =0;i<(L1dim_S_xy+10);i++)
    {
        y=-5;
        for(int j = 0; j <L1dim_S_xy+10;j++)
        {
            z=-5;
            for(int k = 0;k<(L1dim_S_z+10);k++)
            {
                
                xyz_array[0] = i-5;
                xyz_array[1] = j-5;
                xyz_array[2] = k-5;
                
                xyz_index[k*(L1dim_S_xy+10)*(L1dim_S_xy+10)+j*(L1dim_S_xy+10)+i] = general_index(xyz_array);
                z++;
            }
            y++;
        }
        x++;
    }
    
    index_x = (int*)malloc(numberspins_S * sizeof(int));
    index_y = (int*)malloc(numberspins_S * sizeof(int));
    index_z = (int*)malloc(numberspins_S * sizeof(int));
    int ph_array[3]={0};
    
    for(int i =0;i<numberspins_S;i++)
    {
        general_coord(ph_array,i);
        index_x[i]=ph_array[0];
        index_y[i]=ph_array[1];
        index_z[i]=ph_array[2];
    }
    
    init_optomize_energy_index();
    optomize_energy_index();
}


void eye()
{
    
//    for(int i=0;i<numberspins_S;i++)
//    {
//        if(latticepoint_S[i]==2)
//        {
//            int yes=0;
//            for(int ii=0;ii<no_solvent_sites;ii++)
//            {
//                if(i==solvent_loc[ii])
//                {yes=1;}
//
//
//            }
//            if(yes==0)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "solvent no\n");
//                fclose(stdoutlog);
//
//                MPI_Abort(MPI_COMM_WORLD,1);
//            }
//        }
//
//        if(latticepoint_S[i]==1)
//        {
//            int yes=0;
//            for(int ii=0;ii<lengthpercore*pertrusionpercore*numbercores;ii++)
//            {
//                if(i==poly_lattice_indexes[ii])
//                {yes=1;}
//            }
//            if(yes==0)
//            {
//                stdoutlog = fopen(stdoutlogname, "a");
//                fprintf(stdoutlog, "poly no\n");
//                fclose(stdoutlog);
//
//                MPI_Abort(MPI_COMM_WORLD,1);
//            }
//        }
//    }
    
    for(int i=0;i<lengthpercore*pertrusionpercore*numbercores;i++)
    {

        if (latticepoint_S[poly_lattice_indexes[i]] != 1)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "veye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i",i,latticepoint_S[poly_lattice_indexes[i]]);
            fclose(stdoutlog);
            
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        for (int iii = 0; iii < numberneighbors_L; iii++)
        {
            
            if (latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]] != 1)
            {
                int x= latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]];
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "neye iii:%i  i:%i\t %i\t %i\t %i\n",iii, i,neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii],x,latticepoint_L[neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii]]);
                fclose(stdoutlog);

                MPI_Abort(MPI_COMM_WORLD,1);
            }
        }
        
        for(int j=0;j<lengthpercore*pertrusionpercore*numbercores;j++)
        {
            if((reduced_lattice_optimization(poly_lattice_indexes[i],poly_lattice_indexes[j])<4 || poly_lattice_indexes[i]==poly_lattice_indexes[j]) && i!=j)
            {
                float x = point_distance(poly_lattice_indexes[i],poly_lattice_indexes[j]);
                int y = int_point_distance(poly_lattice_indexes[i],poly_lattice_indexes[j]);
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],poly_lattice_indexes[j]),poly_lattice_indexes[i],poly_lattice_indexes[j],x,y);
                fclose(stdoutlog);
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],poly_lattice_indexes[j]),poly_lattice_indexes[i],poly_lattice_indexes[j],x,y);
                fclose(stdoutlog);
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            
        }
        
        for (int j = 0; j < no_solvent_sites; j++) // calculatw=es energy for solvent interactions
        {
            if((reduced_lattice_optimization( poly_lattice_indexes[i],solvent_loc[j])<4 || poly_lattice_indexes[i]==solvent_loc[j]))
            {
                float x = point_distance(poly_lattice_indexes[i],solvent_loc[j]);
                int y = int_point_distance(poly_lattice_indexes[i],solvent_loc[j]);
                
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "sol i:%i    j:%i   dist:%i  p_l_i:%i    p_l_i:%i    x:%f   y:%i\n",i,j,reduced_lattice_optimization( poly_lattice_indexes[i],solvent_loc[j]),poly_lattice_indexes[i],solvent_loc[j],x,y);
                fclose(stdoutlog);
      
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            if (latticepoint_S[solvent_loc[j]] != 2)
            {
                stdoutlog = fopen(stdoutlogname, "a");
                fprintf(stdoutlog, "sveye i:%i latticepoint_S[poly_lattice_indexes[i]]: %i",i,latticepoint_S[solvent_loc[i]]);
                fclose(stdoutlog);
                
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            for (int iii = 0; iii < numberneighbors_L; iii++)
            {
                
                if (latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[j] + iii]] != 2)
                {
                    int x= latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[j] + iii]];
                    
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "sneye iii:%i  i:%i\t %i\t %i\t %i\n",iii, i,neighbor_L[numberneighbors_L * poly_lattice_indexes[i] + iii],x,latticepoint_L[neighbor_L[numberneighbors_L * solvent_loc[i] + iii]]);
                    fclose(stdoutlog);

                    MPI_Abort(MPI_COMM_WORLD,1);
                }
            }
            
        }
    }
    
}
