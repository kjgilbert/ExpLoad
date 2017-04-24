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





int main(int argc, char* argv[]) {
    
    if ( argc != 2 )
    {                   // argc should be length 2 for correct execution
                        //spot 0 is the program name when running it in command line
                        //spot 1 is the name of the input file that is going to be read
       cout << "\n You must provide the input parameter file. \n";
       return 0;
    }
    
    
    string fn = argv[1]; //filename for parameter inputs
    cout << fn; 
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
   int random();
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
    int snapshot;                   // number of generations between two snapshots of the whole metapopulation - it outputs data every 'snapshot' generations
    int m1,m2;                      // size of the 2D grid x = m2, y = m1 if you expand across the x-axis
    int replicates;                 // number of replicates of the simulation
    int expansion_start;            // burnin time ON the grid - to mig-seln balance AFTER the ancestral pop burnin
    int selectionMode;              // SoftSelection = 0, HardSelection = 1
    int expansionMode;              // 0 = linear expansion, 1 = radial expansion, 2 linear - starting from both ends of the habitat, 
    int starting_demes;             // used to calculate number of initially colonized demes as landscape height times this number
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
    

    if(params.size() < 16)
    {
        m1 = params[0];
        m2 = params[1];
        starting_demes = params[2];
        capacity = params[3];
        anc_pop_size = params[4];
        
        burnin_time = params[5];
        expansion_start = params[6];
        generations = params[7];
        snapshot = params[8];
        replicates = params[9];

        expansionMode = params[10];
        selectionMode = params[11];
        mu = params[12];
        m = params[13];
        s = params[14];
        loci = params[15];

        tot_demes = m1*m2;
        initial_colonized = starting_demes*m1;  
    }
    else if(params.size() > 16)
    {
           // then there must have been a random seed specified
        m1 = params[0];
        m2 = params[1];
        starting_demes = params[2];
        capacity = params[3];
        anc_pop_size = params[4];
        
        burnin_time = params[5];
        expansion_start = params[6];
        generations = params[7];
        snapshot = params[8];
        replicates = params[9];

        expansionMode = params[10];
        selectionMode = params[11];
        mu = params[12];
        m = params[13];
        s = params[14];
        loci = params[15];
        
        curSeed = params[16];

        tot_demes = m1*m2;
        initial_colonized = starting_demes*m1;  
    }
    else    // it's never going to get to this else statement the way I have it coded now, but it doesn't seem to count as only less than a number for param size?
    {
        cout << "\n MISSING REQUIRED NUMBER OF PARAMETERS, EXITING SIMULATION.\n";
    }

    
    
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
    
    // root dir for output files and  the starting name - fix this to come from param file
    // this doesn't work because it needs to know the size ahead of time :(   std::string full_path = path + "out_test_";
    const char base[] = "Outputs/out_";
   	// root dir for log files and the starting name - fix this to come from param file
    const char filename_log[] = "Outputs/out_log";


 
                                                           
    ofstream outputfile,outputfile2,outputfile3,outputfile4,logfile;                        // streams to outputfiles
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
    
    cout << "Simulating an expansion on a " << m1 << "x"<<m2<<" grid. \n";
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
    
    cout << endl;
 
    vector<double> outdata(tot_demes);  
    
    World Grid2D(m1,m2,initial_colonized,anc_pop_size,burnin_time,capacity,expansionMode,mu,s,m);   // initialize world: grid size (m1,m2), number of initially colonized demes, 
                                                                                                    // size of original population, burn in time of original population, capacity of demes, mode of intial colonization   
    
   // srand(time(NULL));		//might be an artifact -- CHECK - in case it overwrites the random number seed ???
    
    
    
    for (rep = 0;rep<replicates;rep++)                                  // loop that simulates replicates for the same set of parameters and initial conditions
    {
        cout << "Beginning replicate " << rep+1 << "/" << replicates << endl;
                
        Grid2D.setParams(capacity,mu,s,m);
        sprintf(filename,"%s%d",base,rep);                              // these are the names of the outputs per rep, only need the one line the 2 below are for different cases, but the others could create separate file outputs for diff summ stats
        //sprintf(filename2,"%s%s%d",base,"_hom_wt_",rep);
        //sprintf(filename3,"%s%s%d",base,"_het_",rep);

        
        outputfile.open(filename);
        cout << " filename:" <<filename;
//        outputfile2.open(filename2);
//        outputfile3.open(filename3);
        
         // migration barrier along expansion axis:
        
//        for(i = 10; i < m2; i++)
//        {
//                Grid2D.startExpansion((m1/2)*m2+i,0);
//        }
//       for (i = 0;i < m1; i++)
//        {
//             Grid2D.startExpansion((i)*m2+(initial_colonized/m1)+1,0);     //  Migration-barrier for burn in
//        }
     

        for (k = 0;k<expansion_start;k++)                                      
        {                        
            Grid2D.migrate(initial_colonized);                                       // migration        
            Grid2D.reproduce(selectionMode);                                // reproduction and selection     
        }  
    
        cout << "\n Burn-in finished, expansion into new territory starts.";

//        for (i = 0;i < m1; i++)
//        {
//             Grid2D.startExpansion((i)*m2+(initial_colonized/m1)+1,capacity);     // remove Migration-barrier after burn in
//        }
        
//        Grid2D.startExpansion((m1/2)*m2+(initial_colonized/m1)+1,0);     
        
        // Grid2D.startExpansion((m1/2)*m2+(initial_colonized/m1)+1,capacity);     // open 1 deme in Migration-barrier 
        
//        Grid2D.startExpansion((m1-1)*m2+(initial_colonized/m1)+1,capacity);             // open 1 deme in Migration-barrier, part of corridor ooA
//        Grid2D.startExpansion((m1-1)*m2+(initial_colonized/m1)+2,capacity);             // open 1 deme in Migration-barrier
//        Grid2D.startExpansion((m1-1)*m2+(initial_colonized/m1)+3,capacity);             // open 1 deme in Migration-barrier
//        
//        Grid2D.startExpansion((m1-2)*m2+(initial_colonized/m1)+1,capacity);             // open 1 deme in Migration-barrier, part of corridor ooA
//        Grid2D.startExpansion((m1-2)*m2+(initial_colonized/m1)+2,capacity);             // open 1 deme in Migration-barrier
//        Grid2D.startExpansion((m1-2)*m2+(initial_colonized/m1)+3,capacity);             // open 1 deme in Migration-barrier
//       
        
        Grid2D.setCapacity(capacity);                                           // remove Migration-barrier completely (and any barrier that might've been drawn on the landscape) all are removed here
      
        //Grid2D.startExpansion((m1/2)*m2+(initial_colonized/m1)+1,0);   

        for(i = 0; i< (generations)/snapshot;i++)                                // loop through one set of generations to the first, 2nd, ... snapshot 
        {            
                outdata = Grid2D.getMeanFit();                                  // get mean fitness of the whole population
    
                for (j = 0;j<tot_demes;j++)                                     // write it to file
                { 
                        outputfile << outdata[j] << " ";
                }
                outputfile << "\n";
                
//////                outdata = Grid2D.getGenotypeFrequencies(0,loci,1);              // get  heterozygotes
//////                
//////                sprintf(filename4,"%s%s%d",filename3,"_gen_",(i*snapshot));
//////                outputfile3.open(filename4);
//////                
//////                for (j = 0;j<(tot_demes);j++)                                     
//////                { 
//////                    for ( k = 0;k< loci;k++)
//////                    {
//////                        outputfile3 << outdata[j*loci+k] << " ";
//////                    }
//////                    outputfile3 << "\n";
//////                }
//////                outputfile3.close();
//////                
//////                outdata = Grid2D.getGenotypeFrequencies(0,loci,0);              // get  wt homozygotes
//////                
//////                sprintf(filename4,"%s%s%d",filename2,"_gen_",(i*snapshot));
//////                outputfile3.open(filename4);
//////                for (j = 0;j<(tot_demes);j++)                                     
//////                { 
//////                    for ( k = 0;k< loci;k++)
//////                    {
//////                        outputfile3 << outdata[j*loci+k] << " ";
//////                    }
//////                    outputfile3 << "\n";
//////                }
//////                outputfile3.close();
        
                
                for (k = 0;k<snapshot;k++)                                      // now go through the first set of gens before the next snapshot, etc
                {         
                        Grid2D.migrate(tot_demes);                              // migration        
                        Grid2D.reproduce(selectionMode);                        // reproduction and selection     
                }                  
        }
    
        outdata = Grid2D.getMeanFit();                                          // write data to output file
        
        for (j = 0;j<tot_demes;j++) 
        { 
                outputfile << outdata[j] << " ";
                //finaloutputfile << outdata[j] << " ";
        }
        
        outputfile << "\n";
             
//        outdata = Grid2D.getGenotypeFrequencies(0,loci,0);              // get ancestral homozygotes
//
//        sprintf(filename4,"%s%s%d",filename2,"_gen_",(i*snapshot));
//        outputfile3.open(filename4);
//        for (j = 0;j<(tot_demes);j++)                                     
//        { 
//            for ( k = 0;k< loci;k++)
//            {
//                outputfile3 << outdata[j*loci+k] << " ";
//            }
//            outputfile3 << "\n";
//        }
//        outputfile3.close();
//
//        outdata = Grid2D.getGenotypeFrequencies(0,loci,1);              // get  heterozygotes
//
//        sprintf(filename4,"%s%s%d",filename3,"_gen_",(i*snapshot));
//        outputfile3.open(filename4);
//        for (j = 0;j<(tot_demes);j++)                                     
//        { 
//            for ( k = 0;k< loci;k++)
//            {
//                outputfile3 << outdata[j*loci+k] << " ";
//            }
//            outputfile3 << "\n";
//        }
//        outputfile3.close();

        outputfile.close();
//        outputfile2.close();
//        outputfile3.close();
//        outputfile4.close();
     
        Grid2D.clear(m1,m2,initial_colonized,anc_pop_size,burnin_time,capacity,expansionMode,mu,s,m); 		// clear for next rep             
    }

    return 0;
}



