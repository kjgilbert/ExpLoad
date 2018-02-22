/* 
 * File:   main.cpp
 * Author: stephan peischl
 *
 * Created on February 27, 2012, 5:37 PM
 * Modified by Kim Gilbert, starting April 2017
 */

#include <cstdlib>
#include <iostream> // For input/output
#include <fstream>
#include <string>
#include <sstream>  // kjg adding to test input
#include <vector>   // kjg adding to test input

#include "range_expansion.h"
//#include "rng.h"
#include "rng2.h"

#ifdef  _GCC_
#include <sys/time.h>
#else
#include <time.h>
#endif

#include <unistd.h>

using namespace std;


double fitnessConstant; // a global variable for scaling fitness


int main(int argc, char* argv[]) {
 
  
 ///*  // comment this out for debug, and the part in parenthesis of main above  
    if ( argc != 2 )
    {                   // argc should be length 2 for correct execution
                        //spot 0 is the program name when running it in command line
                        //spot 1 is the name of the input file that is going to be read
       cout << "\n You must provide the input parameter file. \n";
       return 0;
    }
    
    
    string fn = argv[1]; //filename for parameter inputs
    cout << fn; 
//*/

      
//    fstream file;
//    file.open(fn);
    //your logic here
//    file.close(); 
    
    
/* find the path of the executable:
    // Get the last position of '/'
    std::string aux(argv[0]);
    // get '/' or '\\' depending on unix/mac or windows.
#if defined(_WIN32) || defined(WIN32)
    int pos = aux.rfind('\\');
#else
    int pos = aux.rfind('/');
#endif
    // Get the path and the name
    std::string path = aux.substr(0,pos+1);
    std::string name = aux.substr(pos+1);
    // show results
    std::cout << "Path: " << path << std::endl;
    std::cout << "Name: " << name << std::endl;
*/
    
    
// find the path of the paramfile:
// I think this one is the more useful because then I can put the outputs wherever the paramfile is stored
   // Get the last position of '/'
///*   // comment this out for debug
    std::string aux(argv[1]);       

    // get '/' or '\\' depending on unix/mac or windows.
#if defined(_WIN32) || defined(WIN32)
    int pos = aux.rfind('\\');
#else
    int pos = aux.rfind('/');
#endif

    // Get the path and the name
    std::string path = aux.substr(0,pos+1);
    std::string name = aux.substr(pos+1);
    // show results
//*/    
    
    
    
    
    
    
    
    
    
    
    
    
// make input file setting to either read in random seed or generate a random one
   
   //Loro_09_04_13 Random number initialization
   long randSeed;
   #ifdef  _GCC_
   srand(1);
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   randSeed = ((tv.tv_sec ^ tv.tv_usec) ^ getpid()) % 1000000; //Implies only  one million possible seeds!

   #else
   int random();    // used to be an error saying random not declared in this scope... check -KJG
   long curRand = rand() % 1000;
   long curTime = time(NULL);
   long curClock = time(0);     // I added this - check with Stephan on this line, what is curClock KJG
   randSeed = (long) (1.0 * curTime + curClock * curRand) % (200000 - curRand);
   #endif

   //Could be good to have the possibility to use an external seed for random number generator
   long manualSeed=0;
   bool seedFound=false;
   long curSeed;
   if (seedFound) curSeed=abs(manualSeed);
   else curSeed = randSeed;
   
   
   // read in inputs from file
   // otherwise give error that no input file found
   
   
   
   
   
   
   
   
   
   
   //__________________________________________________________________________
   
   
   initializeRan3(curSeed);

       
    // parameters 
    
    int anc_pop_size;               // size of ancestral population - exists outside landscape, it's the pop that burns in for the burnin time defined below
    int capacity;                   // carrying capacity of a deme
    int generations;                // number of generations that the expansion lasts
    int burnin_time;                // number of burn in phase of the ancestral population
    double wf;                      // ignore - it tells you what deme is most recently colonized
    double s;                       // selection coefficient
    double m;                       // migration rate
    double mu;                      // genome-wide mut rate U
    double phi;                     // the proportion of occurring mutations that are DELETERIOUS - the rest will have the opposite seln coefficient (i.e. be beneficial)
    int snapshot;                   // number of generations between two snapshots of the whole metapopulation - it outputs data every 'snapshot' generations
    int m1,m2;                      // size of the 2D grid x = m2, y = m1 if you expand across the x-axis
    int replicates;                 // number of replicates of the simulation
    int expansion_start;            // burnin time ON the grid - to mig-seln balance AFTER the ancestral pop burnin
    int selectionMode;              // SoftSelection = 0, HardSelection = 1
    int expansionMode;              // 0 = linear expansion, 1 = radial expansion, 2 linear - starting from both ends of the habitat, 
    int expansionModeKim;           // 0 = normal range expansion, 1 = open front shift with receding trailing edge, 2 = controlled shift at expanding and receding edge
    int starting_demes;             // used to calculate number of initially colonized demes as landscape height times this number
    int niche_width;           // the width for a 1-D shift of the demes undergoing the shift, aka the niche
    int theta;                  // the patch will shift right by one deme every theta number of generations
    int loci;                       // number of loci to simulate
//    int loci = 1000;    // 3 first and last deme are colonized (at opposite edges), 4 = four corners of the habitat are colonized
          
    // counters and other parameters for running
    
    int rep;                        // this is the counter for going through reps, see if runs after deleting this - we need the variable, but don't need it to have value 0
    int i,j,k;                      // other counters
    int tot_demes = m1*m2;          // total number of demes in the world
    int initial_colonized;          // number of initially colonized demes (location of demes is determined via mode)
    
//__________________________________________________________________________

    
    
   
    
    
// read in the input file from the command line argument

    //___________________________________________________________________
    // Stephan's old code:
    //      ifstream infile;
    //      // this is the default input file, looks in working directory for this file
    //      infile.open ("input_parameters.txt", ifstream::in);                            
    //___________________________________________________________________
//    char wd[250];           // the working directory defined in the input file

    
    // comment all this out for debug with hardcoded params
///* 
    double par;
    vector<double> params;          // character array to hold all the inputs from paramfile
    
    ifstream infile ( argv[1] );
    if ( !infile.is_open() )
    {
      cout << "Could not open file, exiting simulation.\n";
      return 0;
    }
    else
    {
        char contents;
        while ( infile >> par )
            params.push_back(par);
    }
    

    if(params.size() > 18)
    {
        m1 = params[0];
        m2 = params[1];
        starting_demes = params[2];
        niche_width = params[3];
        capacity = params[4];
        anc_pop_size = params[5];
        
        burnin_time = params[6];
        expansion_start = params[7];
        theta = params[8];
        generations = params[9];
        snapshot = params[10];
        replicates = params[11];

        expansionMode = params[12];
        expansionModeKim = params[13];
        selectionMode = params[14];
        mu = params[15];
        m = params[16];
        s = params[17];
        phi = params[18];

        tot_demes = m1*m2;
        initial_colonized = starting_demes*m1;  
    }
    else    // it's never going to get to this else statement the way I have it coded now, but it doesn't seem to count as only less than a number for param size?
    {
        cout << "\n MISSING REQUIRED NUMBER OF PARAMETERS, EXITING SIMULATION.\n";
    }
//*/ 

/*   uncomment this to debug and have hardcoded parameters
    m1 = 5;
    m2 = 40;
    starting_demes = 5;
    niche_width = 5;
    capacity = 20;
    anc_pop_size = 5000;
    
    burnin_time = 10;
    expansion_start = 400;
    theta = 5;
    generations = 10;
    snapshot = 10;
    replicates = 1;

    expansionMode = 0;
    expansionModeKim = 1;   // 0 = full expand, 1 = open front, 2 = controlled
    selectionMode = 1;      // 0 = soft, 1 = hard
    mu = 0.1;
    m = 0.05;
    s = -0.005;
    phi = 0.9;

    tot_demes = m1*m2;
    initial_colonized = starting_demes*m1;  
*/
    
   loci = 2000;     // right now number of loci has to be hard coded in 
    
/*
 * This code works, but I have yet to figure out how to parse it when it contains both text and integers, so for now sticking to Stephan's original implementation for input files

 *  // We assume argv[1] is a filename to open
    ifstream the_file ( argv[1] );
    // Always check to see if file opening succeeded
    if ( !the_file.is_open() )
      cout << "Could not open file\n";
    else {
      char contents;
      // the_file.get ( contents ) returns false if the end of the file
      //  is reached or an error occurs
      while ( the_file.get ( contents ) )
        cout<< contents;
    }
    // the_file is closed implicitly here
 */

    
    
    
    
    
// set up all the log and output filenames and files
    char filename[150]; 
    char filename2[150]; 
    char filename3[150];  
    char filename4[150];  
    char filename5[150];  
    
    // root dir for output files and  the starting name - fix this to come from param file
    // this doesn't work because it needs to know the size ahead of time :(   std::string full_path = path + "out_test_";
    const char base[] = "Outputs/out_";
   	// root dir for log files and the starting name - fix this to come from param file
    const char filename_log[] = "Outputs/out_log";


 
                                                           
    ofstream outputfile,outputfile2,outputfile3,outputfile4,outputfile5,logfile;                        // streams to outputfiles
    logfile.open(filename_log);
    logfile << "Random number generator initialized with seed " << curSeed << "\n";    
    logfile << "Simulating an expansion on a " << m1 << "x"<<m2<<" grid. \n";
    logfile << "Selection is ";
    if(selectionMode == 0)
        logfile << "soft.\n";
    else
        logfile << "hard.\n";
    
    logfile << "\nParameters: \n   Carrying capacity: "     << capacity << 
                             "\n   Migration rate: "        << m << 
                             "\n   Selection coefficient: " << s << 
                             "\n   Mutation rate: "         << mu << 
                             "\n   Burnin time: "           << burnin_time <<
                             "\n   Expansion start: "       << expansion_start;
    logfile << "\n Ancestral population size: " << anc_pop_size << "\n";
    logfile << "\n Expansion Mode: " << expansionMode << "\n";
    logfile << "\n Number of replicates: " << replicates<< "\n";
    
    logfile << endl;
    
    cout << "\nSimulating an expansion on a " << m1 << "x"<<m2<<" grid. \n";
    cout << "Selection is ";
    if(selectionMode == 0)
        cout << "soft.\n";
    else
        cout << "hard.\n";
    
    cout << "\nParameters: \n   Carrying capacity: "     << capacity << 
                             "\n   Migration rate: "        << m << 
                             "\n   Selection coefficient: " << s << 
                             "\n   Mutation rate: "         << mu << 
                             "\n   Burnin time: "           << burnin_time <<
                             "\n   Expansion start: "       << expansion_start;
    cout << "\n  Ancestral population size: " << anc_pop_size << "\n";
    cout << "\n  Expansion Mode: " << expansionMode;
    cout << "\n  Number of replicates: " << replicates;
    cout << "\n  Random number seed: " << curSeed;
    
    cout << endl;
 
   
    
    // set the global variable for fitness scaling to 1 before the burn-in starts
    // NO!! doesn't work here, it's been set in deme where reproduce SS burnin happens because the ancestral pop only ever reproduces with that function
    //    fitnessConstant = 1;
     
    // INITIALIZE THE SIMULATION (AND ITS LANDSCAPE)
    vector<double> outdata(tot_demes);  
    vector<double> tempdata1(tot_demes);  
    vector<double> tempdata2(tot_demes);  
    vector<double> trailing_edge(m1);              // the identity of the deme (for 1-D shifts) at the far left of the landscape, it will always be zero, as is defined below in the code
   
    World Grid2D(m1,m2,initial_colonized,anc_pop_size,burnin_time,capacity,expansionMode,mu,s,m,phi);   // initialize world: grid size (m1,m2), number of initially colonized demes, 
                                                                                                    // size of original population, burn in time of original population, capacity of demes, mode of intial colonization   

    //srand(time(NULL));  // add back in, maybe this was what made reps different
    
    // GO THROUGH REPS
    for (rep = 0;rep<replicates;rep++)                                  // loop that simulates replicates for the same set of parameters and initial conditions
    {
        cout << "Beginning replicate " << rep+1 << "/" << replicates << endl;
                
        Grid2D.setParams(capacity,mu,s,m);
// orig        sprintf(filename,"%s%d",base,rep);                              // these are the names of the outputs per rep, only need the one line the 2 below are for different cases, but the others could create separate file outputs for diff summ stats
// orig        sprintf(filename2,"%s%s%d",base,"_hom-wt_rep_",rep);
// orig        sprintf(filename3,"%s%s%d",base,"_het_rep_",rep);
        sprintf(filename,"%s%s%d%s%d%s%d",base,"wid",niche_width,"_speed",theta,"_rep_",rep);
        sprintf(filename2,"%s%s%d%s%d%s%d",base,"hom-mut_wid",niche_width,"_speed",theta,"_rep_",rep);
        sprintf(filename3,"%s%s%d%s%d%s%d",base,"het_wid",niche_width,"_speed",theta,"_rep_",rep); // this is just making a filename that is appended to and used later as filename 4, so don't need to open it at the moment
        sprintf(filename5,"%s%s%d%s%d%s%d",base,"popdens_wid",niche_width,"_speed",theta,"_rep_",rep);

 
        outputfile.open(filename);
        cout << " filename:" <<filename;
        outputfile2.open(filename2);
        outputfile4.open(filename5);
        
    

        for (k = 0;k<expansion_start;k++)                                      
        {                        
            Grid2D.migrate(initial_colonized);                                       // migration        
            Grid2D.reproduceBurnin(selectionMode, phi);                                // reproduction and selection     
        } 
        
    
        cout << "\n Burn-in finished, expansion into new territory starts.";
        
        // get mean fitness of the whole metapop at the end of the burning to scale fitness throughout the remainder of the simulation
        // and scale mean fitness by pop size in each deme
        tempdata1 = Grid2D.getMeanFit();
        tempdata2 = Grid2D.getDemeDensity();
        double fitnessSum;
        double popSizeSum;
        
        for(int it=0; it<tot_demes; it++)
        {
            outdata[it]=tempdata1[it] * tempdata2[it];
            fitnessSum += outdata[it];               // this sum will include the sum of all empty, size 0 demes, but won't matter as long as I only divide by the number of occupied demes during the burnin
            popSizeSum += tempdata2[it];
        }

        // get the average by dividing the sum by the number of occupied demes in the burnin
        // this resets the global variable for fitness scaling to the mean fitness of the metapop at the end of the burn-in
        fitnessConstant = fitnessSum/popSizeSum;
        

// GET THE MEAN HERE TO SCALE FITNESS
        
        
// ALSO BEGIN BENEFICIAL MUTATIONS HERE, and turn them off during the burn-in
        
        
        
// controlled width range shift in following section
    if(expansionModeKim == 2)
    {
        Grid2D.setCapacity(0);      // remove Migration-barrier completely (and any barrier that might've been drawn on the landscape) all are removed here
                                    // so here we set carrying capacity to 0 for all demes
                                    // then in the next line, we put the carrying capacity back up to what we want, but only in the "niche" or shifting range we want
        for (int row = 0; row < m1; row++)
        {
            int deme = row*m2;
            for (deme; deme<(niche_width+(row*m2)); deme++)
            {
                Grid2D.setDemeCapacity(deme,capacity);
            }
        }
    }
        
// open front range shift
    if(expansionModeKim == 1){
        Grid2D.setCapacity(capacity);       // set all demes up to size
    }
                                            //Grid2D.startExpansion((m1/2)*m2+(initial_colonized/m1)+1,0);  // old version for starting the expansion  
    
    if(expansionModeKim == 1 || expansionModeKim == 2)
    {
        trailing_edge = Grid2D.getEdgeDemes(0);    // the trailing edge always starts at the leftmost deme(s)
        
        for(i = 0; i< generations;i++)      // loop through all generations of the simulation
        {        
            if(i % snapshot == 0)           // write out all the data every snapshot
            {
                outdata = Grid2D.getMeanFit(); // get mean fitness of the whole population

                for (j = 0; j < tot_demes; j++) // write it to file
                {
                    outputfile << outdata[j] << " ";
                }
                outputfile << "\n";

                outdata = Grid2D.getGenotypeFrequencies(0, loci, 1); // get  heterozygotes

                sprintf(filename4, "%s%s%d", filename3, "_gen_", (i));
                
                outputfile3.open(filename4);

                for (j = 0; j < (tot_demes); j++) {
                    for (k = 0; k < loci; k++) {
                        outputfile3 << outdata[j * loci + k] << " ";
                    }
                    outputfile3 << "\n";
                }
                outputfile3.close();

//                outdata = Grid2D.getGenotypeFrequencies(0, loci, 0); // get  wt homozygotes
                outdata = Grid2D.getGenotypeFrequencies(0, loci, 2); // get  mutant homozygotes

                sprintf(filename4, "%s%s%d", filename2, "_gen_", (i));
                outputfile3.open(filename4);
                for (j = 0; j < (tot_demes); j++) {
                    for (k = 0; k < loci; k++) {
                        outputfile3 << outdata[j * loci + k] << " ";
                    }
                    outputfile3 << "\n";
                }
                outputfile3.close();
                
                
                /// GET POP DENSITY ACROSS DEMES
                outdata = Grid2D.getDemeDensity(); // get deme density across the whole population

                for (j = 0; j < tot_demes; j++) // write it to file
                {
                    outputfile4 << outdata[j] << " ";
                }
                outputfile4 << "\n";

            }


            if (i % theta == 0)             // move the 'patch' carrying capacity for a 1- or 2-D shift, it moves one deme/column of demes at a time
            {
                for(int j=0; j<m1; j++)     // go through all the demes along the trailing edge and reset those capacities then move the edge up by 1
                {
                    Grid2D.setDemeCapacity(trailing_edge[j],0);
                    if(expansionModeKim == 2) Grid2D.setDemeCapacity(trailing_edge[j]+niche_width,capacity);
                    if(trailing_edge[m1-1] < (((m2*m1)-1) - niche_width)){     // once the trailing edge reaches the far right of the landscape, stop so can look at recovery -- minus 1 because m2 is intuitive size, but C goes from 0 to that size - 1, so measure in C's language
                        trailing_edge[j] += 1;                        
                    }
                }
                
//                Grid2D.setDemeCapacity(trailing_edge,0);
//                if(expansionModeKim == 2) Grid2D.setDemeCapacity(trailing_edge+niche_width,capacity); // comment this out for open front expansion, include it for controlled expansion
//                if(trailing_edge < (m2 - niche_width)){     // once the trailing edge reaches the far right of the landscape, stop so can look at recovery -- CHECK IF I'M OFF BY 1 HERE AND NEED TO STOP IT SOONER?
//                    trailing_edge += 1;
//                }
            }
                    

            // then normal simulation routines proceed for the rest of this generation
            Grid2D.migrate(tot_demes); // migration        
            Grid2D.reproduce(selectionMode); // reproduction and selection     
 
        }


        outdata = Grid2D.getMeanFit();                                          // write all data at the end of sim to object
        
        for (j = 0;j<tot_demes;j++)                                             // put all the final data from that object into the output file
        { 
                outputfile << outdata[j] << " ";
                //finaloutputfile << outdata[j] << " ";
        }
        
        
        outdata = Grid2D.getDemeDensity(); // get deme density across the whole population at the end of the sim

        for (j = 0; j < tot_demes; j++) // write it to file
        {
             outputfile4 << outdata[j] << " ";
        }
        outputfile4 << "\n";
    }   // end expansionModeKim 1 or 2
        
        
 // this is the normal range expansion setup, replace all of above from the cout before looping through generations
    if(expansionModeKim == 0)
    { 
        Grid2D.setCapacity(capacity);      // remove Migration-barrier completely (and any barrier that might've been drawn on the landscape) all are removed here

        for(i = 0; i< (generations)/snapshot;i++)                                // loop through one set of generations to the first, 2nd, ... snapshot 
        {
 
                outdata = Grid2D.getMeanFit();                                  // get mean fitness of the whole population
     
                for (j = 0;j<tot_demes;j++)                                     // write it to file
                { 
                    outputfile << outdata[j] << " ";
                }
                outputfile << "\n";

                outdata = Grid2D.getGenotypeFrequencies(0, loci, 1); // get  heterozygotes

                sprintf(filename4, "%s%s%d", filename3, "_gen_", (i*snapshot));
                outputfile3.open(filename4);

                for (j = 0; j < (tot_demes); j++)
                {
                    for (k = 0; k < loci; k++)
                    {
                        outputfile3 << outdata[j * loci + k] << " ";
                    }
                    outputfile3 << "\n";
                }
                outputfile3.close();

//                outdata = Grid2D.getGenotypeFrequencies(0, loci, 0); // get  wt homozygotes
                outdata = Grid2D.getGenotypeFrequencies(0, loci, 2); // get  mutant homozygotes

                sprintf(filename4, "%s%s%d", filename2, "_gen_", (i*snapshot));
                outputfile3.open(filename4);
                for (j = 0; j < (tot_demes); j++)
                {
                    for (k = 0; k < loci; k++)
                    {
                        outputfile3 << outdata[j * loci + k] << " ";
                    }
                    outputfile3 << "\n";
                }
                outputfile3.close();

                /// GET POP DENSITY ACROSS DEMES
                outdata = Grid2D.getDemeDensity(); // get mean fitness of the whole population  // NEED TO CREATE THIS FILE, PUT IT SOMEWHERE IN WORLD.CPP

                for (j = 0; j < tot_demes; j++) // write it to file
                {
                    outputfile4 << outdata[j] << " ";
                }
                outputfile4 << "\n";


                for (k = 0;k<snapshot;k++)                                      // now go through the first set of gens before the next snapshot, etc
                {         
                        Grid2D.migrate(tot_demes);                              // migration        
                        Grid2D.reproduce(selectionMode);                        // reproduction and selection     
                }                  
        }
    }    // end expansionModeKim0    

        outputfile << "\n";
             
//        outdata = Grid2D.getGenotypeFrequencies(0,loci,0);              // get ancestral homozygotes
        outdata = Grid2D.getGenotypeFrequencies(0,loci,2);              // get mutant homozygotes

        sprintf(filename4,"%s%s%d",filename2,"_gen_",(i));
        outputfile3.open(filename4);
        for (j = 0;j<(tot_demes);j++)                                     
        { 
            for (k = 0;k< loci;k++)
            {
                outputfile3 << outdata[j*loci+k] << " ";
            }
            outputfile3 << "\n";
        }
        outputfile3.close();


        outdata = Grid2D.getGenotypeFrequencies(0,loci,1);              // get  heterozygotes

        sprintf(filename4,"%s%s%d",filename3,"_gen_",(i));
        outputfile3.open(filename4);
        for (j = 0;j<(tot_demes);j++)                                     
        { 
            for ( k = 0;k< loci;k++)
            {
                outputfile3 << outdata[j*loci+k] << " ";
            }
            outputfile3 << "\n";
        }
        outputfile3.close();

        outputfile.close();
        outputfile2.close();
        outputfile4.close();
     

        Grid2D.clear(m1,m2,initial_colonized,anc_pop_size,burnin_time,capacity,expansionMode,mu,s,m,phi); 		// clear for next rep             

        cout << "Finished replicate " << rep+1 << "/" << replicates << endl;
    }

    return 0;
}



