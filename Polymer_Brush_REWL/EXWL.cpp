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
#include "Mov.h"
#include "EXWL.h"


int histflat(int imin, int imax, double ratio)
{
    // checks flatness of histograms for _all_ walkers in energy window
    int myflat, otherflat;
    int flatproc, ioffset;
    int flatprocr, ioffsetr;
    int flatness_crit = 1;

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
            if ((lngE[x] > 0 || HE[x] > 0) && HE[x] <= minval)
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

            if (lngE[x] > 0 || HE[x] > 0)
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


void accessiblelevels()
{
    //stdoutlog = fopen(stdoutlogname, "a");

    
    
    if (myid == 0) // 'root' in energy window, receive individual g(E) and send merged g(E)
    {
        for (int i = Eminindex; i <= Eglobalwidth; i++)
        {
            if (lngE[i] > 0.005)
            {
                real_lngE[i] = 1;
            }
        }

        MPI_Recv(&maximum_en_config, 1, MPI_INT, numprocs-1, 71, MPI_COMM_WORLD, &status);
        MPI_Recv(&max_printed, 1,MPI_INT, numprocs-1, 72, MPI_COMM_WORLD, &status);
        
        for (int ii = 1; ii < numprocs; ii++)
        {
            MPI_Recv(&lngE_buf[0], hist_size, MPI_DOUBLE, ii, 76, MPI_COMM_WORLD, &status); // get other dens. of states
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "\naccess Proc %i: Received lngE from Proc. %i\n", myid, ii);
            fclose(stdoutlog);
            for (int i = 0; i <= hist_size; i++)
            {
                if (lngE_buf[i] > 0.005)
                {
                    real_lngE[i] = 1;
                }

            }
            for (int i = 0; i < hist_size; i++)
            {
                lngE_buf[i] = 0.0;
            }
        }

        //fclose(file);
        
        int w_l = 0;
        
        sprintf(filename, "Access_Levels.txt");
        if ((file = fopen(filename, "w")) == NULL)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Can not open file %s. Exit.\n", filename);
            fclose(stdoutlog);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            for (int i = Eminindex; i < Eglobalwidth; i++)
            {
                if (real_lngE[i] > 0.005)
                {
                    if(w_l==0)
                    {
                        if(i!=minimum_en_config)
                        {
                            printed=0;
                        }
                        
                        minimum_en_config=i;
                        w_l=1;
                    }
                    if(w_l==1)
                    {
                        if(i>maximum_en_config)
                        {
                            max_printed=0;
                        }
                        
                        maximum_en_config=i;
                    }
                        
                     fprintf(file, "%i\t%f\n", (i + (Eglobalmin)), real_lngE[i]);
                }
            }
            
        }
        fclose(file);
        
        MPI_Send(&maximum_en_config,1, MPI_INT, numprocs-1, 73, MPI_COMM_WORLD);
        MPI_Send(&max_printed,1, MPI_INT, numprocs-1, 74, MPI_COMM_WORLD);
        
    }
    else // send individual g(E) and receive merged g(E)
    {
        if(myid==numprocs-1)
        {
            MPI_Send(&maximum_en_config,1, MPI_INT, 0, 71, MPI_COMM_WORLD);
            MPI_Send(&max_printed,1, MPI_INT, 0, 72, MPI_COMM_WORLD);
        }
        
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "ufig");
        fclose(stdoutlog);
        MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, 0, 76, MPI_COMM_WORLD);
        stdoutlog = fopen(stdoutlogname, "a");
        fprintf(stdoutlog, "Proc %i: Sent lngE to Proc. %i\n", myid, 0);
        fclose(stdoutlog);
        
        if(myid==numprocs-1)
        {
            MPI_Recv(&maximum_en_config, 1, MPI_INT, 0, 73, MPI_COMM_WORLD, &status);
            MPI_Recv(&max_printed, 1, MPI_INT, 0, 74, MPI_COMM_WORLD, &status);
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "Recieved max details %i %i\n", max_printed,maximum_en_config);
            fclose(stdoutlog);
        }
    }
}


void recombine(double countd)
{
    
        sprintf(filename, "DP_%i.txt", myid);
                if ((file = fopen(filename, "w")) == NULL)
                {

                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                else
                {
                    for (int i = 0; i <= Eglobalwidth; i++)
                    {
                        if (lngE[i] > 0.5 && visitdp[i]>0)
                        {
                            fprintf(file, "%i:\t%lli\t", i,visitdp[i]);
                            for(int j=0;j<L1dim_S_z;j++)
                            {
                                fprintf(file, "%f\t ",((double)density_profile[i*L1dim_S_z+j])/((double)(pertrusionpercore*numbercores*lengthpercore))/((double)visitdp[i]));
                            }
                            fprintf(file, "\n");
                        }
                    }
                }
                fclose(file);
        //
                sprintf(filename, "DPS_%i.txt", myid);
                if ((file = fopen(filename, "w")) == NULL)
                {

                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                else
                {
                    for (int i = 0; i <= Eglobalwidth; i++)
                    {
                        if (lngE[i] > 0.5 && visitdp[i]>0){
                            fprintf(file, "%i:\t%lli\t", i,visitdp[i]);
                            for(int j=0;j<L1dim_S_z;j++)
                            {
                                fprintf(file, "%f\t ",((double)density_profile_solvent[i*L1dim_S_z+ j])/((double)visitdp[i])/((double)no_solvent_sites));
                            }
                            fprintf(file, "\n");
                        }
                    }
                }
                fclose(file);
        
        double init_dos = 19962000;

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
            sprintf(filename, "TestL%iq%i.HE.proc%04i.iter0", L1dim_S_xy, q, myid);
            if ((file = fopen(filename, "w")) == NULL)
            {

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


            if ((file = fopen(filename, "w")) == NULL)
            {

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

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                for (int i = Eminindex; i <= Eglobalwidth; i++)
                {
                    if (real_lngE[i] > 0.5) fprintf(file, "%i\t%f\n", (i + (Eglobalmin)), real_lngE[i]);
                }
                fclose(file);
            }
sprintf(filename, "DOS.txt");
            if ((file = fopen(filename, "w")) == NULL)
            {

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                for (int i = Eminindex; i <= Eglobalwidth; i++)
                {
                    if (real_lngE[i] > 0.5) fprintf(file, "%i\t%f\n", (i + (Eglobalmin)), real_lngE[i]);
                }
                fclose(file);
            }

            
            for (int ii = 1; ii < numprocs; ii++)
            {
                MPI_Recv(&tortuosity_buf[0], hist_size, MPI_DOUBLE, ii, 83, MPI_COMM_WORLD, &status);
                MPI_Recv(&rog_buf[0], hist_size, MPI_DOUBLE, ii, 84, MPI_COMM_WORLD, &status);
                MPI_Recv(&visits_buf[0], hist_size, MPI_LONG_LONG, ii, 85, MPI_COMM_WORLD, &status);

    //            MPI_Recv(&density_profile_buf[0], hist_size*L1dim_S_z, MPI_INT, ii, 97, MPI_COMM_WORLD, &status);
    //            MPI_Recv(&density_profile_solvent_buf[0], hist_size*L1dim_S_z, MPI_INT, ii, 98, MPI_COMM_WORLD, &status);
                
                MPI_Recv(&indiv_rog_buf[0], hist_size*numbercores, MPI_DOUBLE, ii, 99, MPI_COMM_WORLD, &status);
                MPI_Recv(&density_profile_buf[0], hist_size*L1dim_S_z, MPI_DOUBLE, ii, 100, MPI_COMM_WORLD, &status);
                MPI_Recv(&density_profile_solvent_buf[0], hist_size*L1dim_S_z, MPI_DOUBLE, ii, 101, MPI_COMM_WORLD, &status);

                MPI_Recv(&rog_z_buf[0], hist_size, MPI_DOUBLE, ii, 102, MPI_COMM_WORLD, &status);
                MPI_Recv(&rog_xy_buf[0], hist_size, MPI_DOUBLE, ii, 103, MPI_COMM_WORLD, &status);
                
                for (int i = Eminindex; i <= Eglobalwidth; i++)
                {
                    if (real_lngE[i] > 0.5)
                    {
                        tortuosity[i] += tortuosity_buf[i];
                        rog[i] += rog_buf[i];
                        rog_z[i] += rog_z_buf[i];
                        rog_xy[i] += rog_xy_buf[i];
                        visits[i] += visits_buf[i];
                        
    //                    for(int j=0;j<L1dim_S_z;j++)
    //                    {
    //                        density_profile[i*L1dim_S_z+j]+=density_profile_buf[i*L1dim_S_z+j];
    //                    }
    //                    for(int j=0;j<L1dim_S_z;j++)
    //                    {
    //                        density_profile_solvent[i*L1dim_S_z+ j]+=density_profile_solvent_buf[i*L1dim_S_z+ j];
    //                    }
                        for(int j=0;j<numbercores;j++)
                        {
                            indiv_rog[i*numbercores+j]+=indiv_rog_buf[i*numbercores+j];
                        }
                        
                        for(int j=0;j<L1dim_S_z;j++)
                        {
                            density_profile[i*L1dim_S_z+j]+=density_profile_buf[i*L1dim_S_z+j];
                        }
                        for(int j=0;j<L1dim_S_z;j++)
                        {
                            density_profile_solvent[i*L1dim_S_z+j]+=density_profile_solvent_buf[i*L1dim_S_z+j];
                        }
                        
                    }
                }
            }

            sprintf(filename, "SysParam_%e.txt", countd);
            if ((file = fopen(filename, "w")) == NULL)
            {

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                for (int i = Eminindex; i <= Eglobalwidth; i++)
                {
                    if (real_lngE[i] > 0.5) fprintf(file, "%i\t%1.13f\t%1.13f\t%lli\t%1.13f\t%1.13f\n", i, (float)rog[i] / ((double)visits[i]), (float)tortuosity[i] / ((double)visits[i]), visits[i],(float)rog_z[i] / ((double)visits[i]),(float)rog_xy[i] / ((double)visits[i]));
                }
            }
            fclose(file);
            
            sprintf(filename, "SysParam.txt");
            if ((file = fopen(filename, "w")) == NULL)
            {

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                for (int i = Eminindex; i <= Eglobalwidth; i++)
                {
                    if (real_lngE[i] > 0.5)  fprintf(file, "%i\t%1.13f\t%1.13f\t%lli\t%1.13f\t%1.13f\n", i, (float)rog[i] / ((double)visits[i]), (float)tortuosity[i] / ((double)visits[i]), visits[i],(float)rog_z[i] / ((double)visits[i]),(float)rog_xy[i] / ((double)visits[i]));
                }
            }
            fclose(file);
            
            sprintf(filename, "Indiv_Rog_%e.txt", countd);
            if ((file = fopen(filename, "w")) == NULL)
            {

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                    for (int i = Eminindex; i <= Eglobalwidth; i++)
                    {
                        if (real_lngE[i] > 0.5)
                        {
                        fprintf(file, "%i:\t", i);
                        for(int j=0;j<numbercores;j++)
                        {
                            fprintf(file, "%f\t ",indiv_rog[i*numbercores+j]/visits[i]);
                        }
                        fprintf(file, "\n");
                        }
                    }
                
            }
            fclose(file);
            sprintf(filename, "DP.txt");
            if ((file = fopen(filename, "w")) == NULL)
            {

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                    for (int i = Eminindex; i <= Eglobalwidth; i++)
                    {
                        if (real_lngE[i] > 0.5)
                        {
                        fprintf(file, "%i:\t%lli\t", i,visits[i]);
                        for(int j=0;j<L1dim_S_z;j++)
                        {
                            fprintf(file, "%f\t ",((double)density_profile[i*L1dim_S_z+j])/((double)(pertrusionpercore*numbercores*lengthpercore))/((double)visits[i]));
                        }
                        fprintf(file, "\n");
                        }
                    }
                
            }
            fclose(file);
            sprintf(filename, "DPS.txt");
            if ((file = fopen(filename, "w")) == NULL)
            {

                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                    for (int i = Eminindex; i <= Eglobalwidth; i++)
                    {
                        if (real_lngE[i] > 0.5)
                        {
                        fprintf(file, "%i:\t%lli\t", i,visits[i]);
                        for(int j=0;j<L1dim_S_z;j++)
                        {
                            fprintf(file, "%f\t ",((double)density_profile_solvent[i*L1dim_S_z+j])/((double)(pertrusionpercore*numbercores*lengthpercore))/((double)visits[i]));
                        }
                        fprintf(file, "\n");
                        }
                    }
                
            }
            
            fclose(file);
            
                    for (int ii = 1; ii < numprocs; ii++)
                    {
                        MPI_Recv(&tortuosity_buf[0], hist_size, MPI_DOUBLE, ii, 830, MPI_COMM_WORLD, &status);
                        MPI_Recv(&rog_buf[0], hist_size, MPI_DOUBLE, ii, 840, MPI_COMM_WORLD, &status);
                        MPI_Recv(&visits_buf[0], hist_size, MPI_LONG_LONG, ii, 850, MPI_COMM_WORLD, &status);

            //            MPI_Recv(&density_profile_buf[0], hist_size*L1dim_S_z, MPI_INT, ii, 97, MPI_COMM_WORLD, &status);
            //            MPI_Recv(&density_profile_solvent_buf[0], hist_size*L1dim_S_z, MPI_INT, ii, 98, MPI_COMM_WORLD, &status);
                        
                        MPI_Recv(&indiv_rog_buf[0], hist_size*numbercores, MPI_DOUBLE, ii, 990, MPI_COMM_WORLD, &status);
                        
                        MPI_Recv(&density_profile_buf[0], hist_size*L1dim_S_z, MPI_DOUBLE, ii, 1000, MPI_COMM_WORLD, &status);
                        MPI_Recv(&density_profile_solvent_buf[0], hist_size*L1dim_S_z, MPI_DOUBLE, ii, 1010, MPI_COMM_WORLD, &status);

                        MPI_Recv(&rog_z_buf[0], hist_size, MPI_DOUBLE, ii, 1020, MPI_COMM_WORLD, &status);
                        MPI_Recv(&rog_xy_buf[0], hist_size, MPI_DOUBLE, ii, 1030, MPI_COMM_WORLD, &status);
                        
                        for (int i = Eminindex; i <= Eglobalwidth; i++)
                        {
                            if (real_lngE[i] > 0.5)
                            {
                                tortuosity[i] -= tortuosity_buf[i];
                                rog[i] -= rog_buf[i];
                                rog_xy[i] -= rog_xy_buf[i];
                                rog_z[i] -= rog_z_buf[i];
                                visits[i] -= visits_buf[i];
                                
            //                    for(int j=0;j<L1dim_S_z;j++)
            //                    {
            //                        density_profile[i*L1dim_S_z+j]+=density_profile_buf[i*L1dim_S_z+j];
            //                    }
            //                    for(int j=0;j<L1dim_S_z;j++)
            //                    {
            //                        density_profile_solvent[i*L1dim_S_z+ j]+=density_profile_solvent_buf[i*L1dim_S_z+ j];
            //                    }
                                for(int j=0;j<numbercores;j++)
                                {
                                    indiv_rog[i*numbercores+j]-=indiv_rog_buf[i*numbercores+j];
                                }
                                for(int j=0;j<L1dim_S_z;j++)
                                {
                                    density_profile[i*L1dim_S_z+j]-=density_profile_buf[i*L1dim_S_z+j];
                                }
                                for(int j=0;j<L1dim_S_z;j++)
                                {
                                    density_profile_solvent[i*L1dim_S_z+j]-=density_profile_solvent_buf[i*L1dim_S_z+j];
                                }
                            }
                        }
                    }
            
        }
        else // send individual g(E) and receive merged g(E)
        {
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "5\n");
            fclose(stdoutlog);
            MPI_Send(&lngE[0], hist_size, MPI_DOUBLE, 0, 76, MPI_COMM_WORLD);

            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "6\n");
            fclose(stdoutlog);
            
            MPI_Send(&tortuosity[0], hist_size, MPI_DOUBLE, 0, 83, MPI_COMM_WORLD);
            MPI_Send(&rog[0], hist_size, MPI_DOUBLE, 0, 84, MPI_COMM_WORLD);
            MPI_Send(&visits[0], hist_size, MPI_LONG_LONG, 0, 85, MPI_COMM_WORLD);
            
    //        MPI_Send(&density_profile[0], hist_size*L1dim_S_z, MPI_INT, 0, 97, MPI_COMM_WORLD);
    //        MPI_Send(&density_profile_solvent[0], hist_size*L1dim_S_z, MPI_INT, 0, 98, MPI_COMM_WORLD);
            MPI_Send(&indiv_rog[0], hist_size*numbercores, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
            
            MPI_Send(&density_profile[0], hist_size*L1dim_S_z, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
            MPI_Send(&density_profile_solvent[0], hist_size*L1dim_S_z, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD);
            
            MPI_Send(&rog_z[0], hist_size, MPI_DOUBLE, 0, 102, MPI_COMM_WORLD);
            MPI_Send(&rog_xy[0], hist_size, MPI_DOUBLE, 0, 103, MPI_COMM_WORLD);
            
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "7\n");
            fclose(stdoutlog);
            
           
            
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "8\n");
            fclose(stdoutlog);

                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "61\n");
                    fclose(stdoutlog);
                    
                    MPI_Send(&tortuosity[0], hist_size, MPI_DOUBLE, 0, 830, MPI_COMM_WORLD);
                    MPI_Send(&rog[0], hist_size, MPI_DOUBLE, 0, 840, MPI_COMM_WORLD);
                    MPI_Send(&visits[0], hist_size, MPI_LONG_LONG, 0, 850, MPI_COMM_WORLD);
                    
            //        MPI_Send(&density_profile[0], hist_size*L1dim_S_z, MPI_INT, 0, 97, MPI_COMM_WORLD);
            //        MPI_Send(&density_profile_solvent[0], hist_size*L1dim_S_z, MPI_INT, 0, 98, MPI_COMM_WORLD);
                    MPI_Send(&indiv_rog[0], hist_size*numbercores, MPI_DOUBLE, 0, 990, MPI_COMM_WORLD);
                    MPI_Send(&density_profile[0], hist_size*L1dim_S_z, MPI_DOUBLE, 0, 1000, MPI_COMM_WORLD);
                    MPI_Send(&density_profile_solvent[0], hist_size*L1dim_S_z, MPI_DOUBLE, 0, 1010, MPI_COMM_WORLD);
                MPI_Send(&rog_z[0], hist_size, MPI_DOUBLE, 0, 1020, MPI_COMM_WORLD);
            MPI_Send(&rog_xy[0], hist_size, MPI_DOUBLE, 0, 1030, MPI_COMM_WORLD);
                    
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "71\n");
                    fclose(stdoutlog);
                    
                   
                    
                    stdoutlog = fopen(stdoutlogname, "a");
                    fprintf(stdoutlog, "81\n");
                    fclose(stdoutlog);
        
        }
    
}


void pseudowl() // A Fake WL Function used to explore the energy landscape before the main WL operation to better enesure all lnge values of the system have been explored
{
    int maxe = Eglobalmin; //holds the index of the maximum energy found by the system used to reset the for loop counter if a higher energy is found

    int qhold = q; // Holds the value of q from before (could be better coded but as far as I know would require significant retooling of some functions)
    q = 2; // minimum q needed for system to explore all of the possible values of the energy landscape (shape of lng(e) is not important in this test

     
    double lnf = 1.0;

    double wsk, dice; // wsk: probability
    int wiggle;
    int wiggletwo;

    // terminal modification factor
    double check_flatness_every = 1000;   // in number of sweeps changed to a large number to encourage visitation of all levels (reduces time in the long run)
    //double check_flatness_every = 10;   // in number of sweeps changed to a large number to encourage visitation of all levels (reduces time in the long run)
    int backup;
    int backuptwo;

    int eold, energie;
    energie = total_energy(poly_lattice_coordinates,poly_lattice_indexes,solvent_loc);
    energie -= Eglobalmin;
    eold = energie;

    int swtch;
    int found = 0;

    long long iqq=0;
    
    for (int i = 0; i < hist_size; i++) HE[i] = 0; //init H(E)

    for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
    {
        wl_pseudo_chain[i] = poly_lattice_indexes[i];
        wl_pseudo_chain_coordinates[3 * i] = poly_lattice_coordinates[3 * i];
        wl_pseudo_chain_coordinates[3 * i + 1] = poly_lattice_coordinates[3 * i + 1];
        wl_pseudo_chain_coordinates[3 * i + 2] = poly_lattice_coordinates[3 * i + 2];
    }
    for (int i = 0; i < no_solvent_sites; i++)
    {
        wl_solvent_loc[i] = solvent_loc[i];
    }

    int within = 0;
    
    if(((energie < Emaxindex) && (energie > Eminindex)))
    {
        within=1;
    }
    
    stdoutlog = fopen(stdoutlogname, "a");
    fprintf(stdoutlog, "iqq: %lli \n within check: %i  %i %i  %i \n",iqq,within,energie,Eminindex,Emaxindex);
    fclose(stdoutlog);
    
    for (int k = 0; k < check_flatness_every; k++)
    {
        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++) // this does 1 MC sweep
        {

            energie = eold + propose_update(); // calculate energy of updated configuration
            
//            if(MoveProposal==1)
//            {
//                if(poly_solvent_mov==0)
//                {
//                    //lattice_polymer_reset(poly_lattice_indexes , ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//                    lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//                    poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//                }
//                if(poly_solvent_mov==1)
//                {
//                    solvent_reset(solvent_loc,ph_solvent_loc);
//
////                    stdoutlog = fopen(stdoutlogname, "a");
////                    fprintf(stdoutlog, "solvent_loc[i]: %i ph_solvent_loc[i]: %i\n",solvent_loc[reset_indexes[0]+1],ph_solvent_loc[reset_indexes[0]+1]);
////                    fclose(stdoutlog);
//                }
//            }
            
            // reject because you got out of bounds of energy window

            if (MoveProposal == 0 || (within==1 && ((energie > Emaxindex) || (energie < Eminindex))) || (within==0 && myid<minmaxid && (energie > Emaxindex)) || (within==0 && myid>=minmaxid && (energie < Eminindex))) // boundary check
            {
                if(MoveProposal==1)
                                    {
                                    if(poly_solvent_mov==0)
                                    {
                                        //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                                        lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                                        poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
                                        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                                        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                                    }
                                    if(poly_solvent_mov==1)
                                    {
                                        solvent_reset(solvent_loc,wl_solvent_loc);
                //                        stdoutlog = fopen(stdoutlogname, "a");
                //                        for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
                //                        {
                //                            fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
                //                        }
                //                        for (int i = 0; i < no_solvent_sites; i++)
                //                        {
                //                            fprintf(stdoutlog, "%4i\t\t%4i\t\t%i\n", solvent_loc[i],ph_solvent_loc[i],wl_solvent_loc[i]);
                //                        }
                //
                //                        fclose(stdoutlog);
                                        solvent_reset(ph_solvent_loc,solvent_loc);

                                    }
                                    }
                
                if(within)
                {
                energie = eold;
                lngE[energie-energie%10] = 1; // the real lngE of the system is set to a value of 1 if visitation occured at a specific energy (can be subtracted later) to force visitation in the real wl process
                }
            }
            else // calculate acceptance propbability
            {
                if(within)
                {
                lngE[energie-energie%10] = 1;
                }
                
                if (energie > maxe) // checks to see if new energy is higher than max e
                {
                    maxe = energie;
                    k = 0; // reset for loop iterative
                }

                // the regular wl process is kept using a place holder pseudolngE to force the system to not move to move to higher energies
                if(within)
                               {
                dice = (1.0 * rand() / (RAND_MAX + 1.0)); // roll the dice
                wsk = exp(pseudolngE[eold-eold%10] - pseudolngE[energie-energie%10]); // calculate acceptance probability
                               }
                else{
                    wsk=1;
                    dice=0;
                }
                
                if (dice > wsk)  // reject
                {
                    if(poly_solvent_mov==0)
                    {
                        //lattice_polymer_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(poly_lattice_indexes, wl_pseudo_chain, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(poly_lattice_coordinates, wl_pseudo_chain_coordinates, reset_indexes[0], reset_indexes[1]);
                        //lattice_polymer_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(ph_poly_lattice_indexes, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(ph_poly_lattice_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                    }
                    if(poly_solvent_mov==1)
                    {
                        solvent_reset(solvent_loc,wl_solvent_loc);
                        solvent_reset(ph_solvent_loc,solvent_loc);

                    }

                    energie = eold;

                }
                else {
                    eold = energie; // accept
                    
                    
                    if(((energie < Emaxindex) && (energie > Eminindex)))
                    {
                        within=1;
                    }
                    
                   if(poly_solvent_mov==0)
                    {
                        //lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
                        poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
                    }
                    if(poly_solvent_mov==1)
                    {
                        solvent_reset(wl_solvent_loc,solvent_loc);

                    }
//                    lattice_polymer_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//                    lattice_poly_index_reset(poly_lattice_indexes, ph_poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//                    poly_coord_reset(poly_lattice_coordinates, ph_poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//                    solvent_reset(solvent_loc,ph_solvent_loc);
//
//
//                    lattice_polymer_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//                    lattice_poly_index_reset(wl_pseudo_chain, poly_lattice_indexes, reset_indexes[0], reset_indexes[1]);
//                    poly_coord_reset(wl_pseudo_chain_coordinates, poly_lattice_coordinates, reset_indexes[0], reset_indexes[1]);
//                    solvent_reset(wl_solvent_loc,wl_solvent_loc);
                    
                    //sysparam(energie-energie%10);
                }

                //if ((energie > eold))
                //{
                //    latticepoint[wiggle] = backup;
                //    latticepoint[wiggletwo] = backuptwo;
                //    energie = eold;
                //}
                //if ((energie < eold))
                //{
                //    k = 0;
                //    eold = energie;
                //}
            }
            if(within)
                           {
            pseudolngE[energie-energie%10] += lnf;
            lngE[energie-energie%10] = 1;
                           }


        }
        if( iqq%200000==0)
        {
            eye();
            sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim_S, q, myid, -1);
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
                    if(pseudolngE[i]>.5)
                    {
                        fprintf(file, "%i\t%i\t%e\t%e\t%e\t%f\t%f\n", i, (i + (Eglobalmin)), 0.0, lngE[i], pseudolngE[i], rog[i] / visits[i], tortuosity[i] / visits[i]);
                    }
                }
                
            }
            fclose(file);
            
            stdoutlog = fopen(stdoutlogname, "a");
            fprintf(stdoutlog, "iqq: %lli \n within check: %i  %i %i  %i \n",iqq,within,energie,Eminindex,Emaxindex);
            for (int i = 0; i < numberspins_S; i++)
            {
                fprintf(stdoutlog, "%4i", latticepoint_S[i]);
                if ((i + 1) % L1dim_S_xy == 0)
                    fprintf(stdoutlog, "\n");
                if ((i + 1) % (L1dim_S_xy*L1dim_S_xy) == 0)
                    fprintf(stdoutlog, "\n");

            }
            fprintf(stdoutlog, "\n");
            fprintf(stdoutlog, "\n");
            fclose(stdoutlog);
            
            stdoutlog = fopen(stdoutlogname, "a");
            for (int i = 0; i < numberspins_L; i++)
            {
                fprintf(stdoutlog, "%4i", latticepoint_L[i]);
                if ((i + 1) % L1dim_L_xy == 0)
                    fprintf(stdoutlog, "\n");
                if ((i + 1) % (L1dim_L_xy*L1dim_L_xy) == 0)
                    fprintf(stdoutlog, "\n");

            }
            fclose(stdoutlog);
            
            stdoutlog = fopen(stdoutlogname, "a");
            for (int i = 0; i < numbercores * lengthpercore * pertrusionpercore; i++)
            {
                fprintf(stdoutlog, "%4i\t\t%4i\t%4i\t%4i\t%4i\t\t%4i\t%4i\t%4i\n", poly_lattice_indexes[i], poly_lattice_coordinates[3 * i], poly_lattice_coordinates[3 * i + 1], poly_lattice_coordinates[3 * i + 2],ph_poly_lattice_indexes[i], ph_poly_lattice_coordinates[3 * i], ph_poly_lattice_coordinates[3 * i + 1], ph_poly_lattice_coordinates[3 * i + 2]);
            }
            fprintf(stdoutlog, "\n");
            fclose(stdoutlog);
        }
        iqq++;
    }

    q = qhold; // when the pseudowl process is done q is reset to the system value

    // The occucpied placeholder lnge and the values of the real lnge are out putted into a prelim file
    sprintf(filename, "Prelim.L%iq%i.HE.proc%04i.iter%i", L1dim_S, q, myid, 0);
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
            fprintf(file, "%i\t%i\t%e\t%e\t%e\t%f\t%f\n", i, (i + (Eglobalmin)), 0.0, lngE[i], pseudolngE[i], rog[i] / visits[i], tortuosity[i] / visits[i]);
        }
        fclose(file);
    }
}
