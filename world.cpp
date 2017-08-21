#include <iostream> // For input/output
#include <fstream> // For file input/output
#include <string>  // For strcpy
#include <time.h>  // For time
#include <stdlib.h>  // For toupper and tolower
#include <math.h>
#include <vector>
#include <list>
#include "range_expansion.h"
//#include <rng.h>
#include "rng2.h"
#include <algorithm>


using namespace std;

int World::number_demes = 700;
int World::colonized_demes = 1;
int World::m1 = 100;
int World::m2 = 100;
double World::phi;

const int bn_size = 1000;
const int bn_length = 10;

//inline int rand_n(int n)
//{
//        return rand()%n;
//}

World::World()
{
    
   int i;
   demes.resize(number_demes);
   migrants.resize(number_demes);
   
   for(i=0;i<colonized_demes;i++) {demes[i].colonize();}  
   
   
}
       
World::World(int length1,int length2,int initial_colonized,int initial_popsize,int burnin_time,int capacity,int mode,double mutation_rate,double s,double migration_rate,double mut_prop)
{
     
   int i,j;
   number_demes = length1*length2;
   m1 = length1;
   m2 = length2;
   phi = mut_prop;
   demes.resize(number_demes);
   migrants.resize(number_demes);
   Migrants propagule;
   vector<Individual>::iterator it;
   initial_population.resize(1);
  
   
   int width, start;
   
   initial_population[0].setParams(initial_popsize,mutation_rate,s,migration_rate);
      
   initial_population[0].colonize();
   initial_population[0].set_selection_dist(phi);
   demes[1].set_selection_dist(phi);
   
   
   for (i = 0; i < burnin_time ; i++)
   {
       initial_population[0].reproduceSSburnin(0, phi);
   }
   
   //bottleneck:
   cout << "\n Ancestral population evolved for " << burnin_time << " generations. \n";
   initial_population[0].printStat();      

   

 
   
   for(j=0;j<number_demes;j++)
   {    
       demes[j].setParams(capacity);
   }
   
     
   
   
   
   if (mode == 1)     // colonize a square in the center of grid for radial expansion
   {
        width = floor(pow(initial_colonized,0.5));
        start = (m1/2-floor(width/2))*m2;
   
        colonized_demes = initial_colonized;
   
        for(i=0;i<width;i++) 
        {
                for(j=0;j<width;j++)
                {
                    propagule = initial_population[0].sampleIndividuals(capacity);    
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                    
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {
                        demes[start+i*m2+m2/2+(j-width/2)].addMigrant((*it));
                    }
                }
        }
   }
   
   width = floor(initial_colonized/m1); // gives the number of columns being colonized at the left of a linear, 2D expansion

   if (mode == 0)     // colonize first initial_size demes at left edge of grid and set a barrier to gene flow to the right of colonized area 
   { 
       
        for(i=0;i<width;i++)   // iterate through columns being colonized
        {
            for(j=0;j<m1;j++)   // m1 is the y-axis, or width/height (not length) of the landscape, go through each row in the column being iterated
            {    
                propagule = initial_population[0].sampleIndividuals(capacity);      // get a propoagule - sample K inds    
                //demes[start+i*m2+m2/2+(j-width/2)].colonize();
               
                for(it = propagule.begin();it!=propagule.end();it++)                // iterate through all inds in propagule
                {    
                    demes[j*m2+i].addMigrant(*it);                                  // deme 0*300+0 = 0, 1*300+0 = 300, 2*300+0 = 600
                }
                
                demes[j*m2+width].setParams(0);                                   // make the barrier of carrying cap 0 at the edge of the starting demes
            }
        }  
        
   }
   
   if (mode == 6)     // colonize first initial_size demes at left edge  of grid, a barrier to gene flow to the right of colonized area, barrier 5 demes-wide (to mimic sahara)
   {    
        for(i=0;i<width;i++)
        {
                for(j=0;j<m1;j++)
                {    
                    propagule = initial_population[0].sampleIndividuals(capacity);      
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                   
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {    
                        demes[j*m2+i].addMigrant(*it);
                    }
                    
                     demes[j*m2+width+1].setParams(0);   
                     demes[j*m2+width+2].setParams(0); 
                     demes[j*m2+width+3].setParams(0); 
                }
        }  
   }
   
   if (mode == 5)     // all demes colonized at the beginning without any founder effects (empty demes are colonized by sampling from left neighbour deme with replacement)
   {
        for(i=0;i<m2;i++)
        {
                for(j=0;j<m1;j++)
                {    
                    if(i == 0)
                    {   
                        propagule = initial_population[0].sampleIndividuals(capacity);  
                    }
                    
                    if (i > 0)    
                    {   
                        propagule = demes[j*m2+i-1].sampleIndividuals(capacity);  
                    }
                       
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                   
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {    
                        demes[j*m2+i].addMigrant(*it);
                    }
                    
                    //demes[j*m2+width+1].setParams(0);
                }
        }  
   }
   
   if (mode == 2)     // colonize first initial_size demes at left and right edge  of grid 
   {
        for(i=0;i<width;i++)
        {
                for(j=0;j<m1;j++)
                {    
                    propagule = initial_population[0].sampleIndividuals(capacity);    
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                    
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {
                        demes[j*m2+i].addMigrant(*it);
                    }
                    
                    propagule = initial_population[0].sampleIndividuals(capacity);    
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                    
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {
                        demes[(j+1)*m2-1 -i].addMigrant(*it);
                    }
                }
        }  
   }
   
      
   if (mode == 3)     // colonize a single deme at opposing corners of the habitat
   {
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
            for(int i = 0;i < m1;i++)    
            {
                demes[i*m2].addMigrant((*it));
            }  
        }
        
        
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
            for(int i = 0;i < m1;i++)    
            {
                demes[(i+1)*m2-1].addMigrant((*it));
            }  
        }
        
        

   }
   
   if (mode == 4)     // colonize a single deme at opposing corners of the habitat
   {
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
                demes[0].addMigrant((*it));
        }
        
        
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
                demes[m1*m2-1].addMigrant(*it);
        }
        
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
                demes[m1-1].addMigrant(*it);
        }
        
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
                demes[m1*m2-1-m1].addMigrant(*it);
        }      
   }
   
   
   for(i=0;i<number_demes;i++) {demes[i].setID(i); }  
   //demes[0].colonize();
   //demes[number_demes-1].colonize();
   //demes[25].colonize();
}

                    
World::~World()
{
}

void World::setDemeCapacity(int deme_location,int capacity)
{
    demes[deme_location].setParams(capacity);
}


void World::clear(int length1,int length2,int initial_colonized,int initial_popsize,int burnin_time,int capacity,int mode,double mutation_rate,double s,double migration_rate,double mut_prop)
{
    int i,j;
    
    phi = mut_prop;
    
    demes.clear();
    migrants.clear();
    initial_population.clear();
    
    demes.resize(number_demes);
    migrants.resize(number_demes);
    initial_population.resize(1);
   
   int width, start;
   Migrants propagule;
   vector<Individual>::iterator it;
   

   cout << "\n Reinitializing ancestral population... ";
   initial_population[0].setParams(initial_popsize,mutation_rate,s,migration_rate);
      
   initial_population[0].colonize();
   initial_population[0].set_selection_dist(phi);
   demes[1].set_selection_dist(phi);
   
   for (i = 0; i < burnin_time ; i++)
   {
       initial_population[0].reproduceSSburnin(0, phi);    // make the ancestral burnin time also have no beneficial mutations (and it's always soft selection)
   }
   
   //bottleneck:
   cout << "\n Ancestral population evolved for " << burnin_time << " generations. \n";
   initial_population[0].printStat();      

   initial_population[0].setParams(bn_size);
   
   for(i = 0; i < bn_length;i++)
       initial_population[0].reproduceSS(wavefrontID);
   
   cout << "\n Ancestral went through bottleneck: " << bn_size <<  " individuals for " << bn_length << " generations. \n";
   
   initial_population[0].printStat();      

   
   for(j=0;j<number_demes;j++)
   {    
       demes[j].setParams(capacity);
   }
   
     
   
   
   
   if (mode == 1)     // colonize a square in the center of grid for radial expansion
   {
        width = floor(pow(initial_colonized,0.5));
        start = (m1/2-floor(width/2))*m2;
   
        colonized_demes = initial_colonized;
   
        for(i=0;i<width;i++) 
        {
                for(j=0;j<width;j++)
                {
                    propagule = initial_population[0].sampleIndividuals(capacity);    
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                    
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {
                        demes[start+i*m2+m2/2+(j-width/2)].addMigrant((*it));
                    }
                }
        }
   }
   
   width = floor(initial_colonized/m1);

   

   if (mode == 0)     // colonize first initial_size demes at left edge  of grid, a barrier to gene flow to the right of colonized area 
   {    
        for(i=0;i<width;i++)
        {
                for(j=0;j<m1;j++)
                {    
                    propagule = initial_population[0].sampleIndividuals(capacity);      
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                   
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {    
                        demes[j*m2+i].addMigrant(*it);
                    }
                     demes[j*m2+width].setParams(0);   
                }
        }  
   }
   
   
   if (mode == 2)     // colonize first initial_size demes at left and right edge  of grid 
   {
        for(i=0;i<width;i++)
        {
                for(j=0;j<m1;j++)
                {    
                    propagule = initial_population[0].sampleIndividuals(capacity);    
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                    
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {
                        demes[j*m2+i].addMigrant(*it);
                    }
                    
                    propagule = initial_population[0].sampleIndividuals(capacity);    
                    //demes[start+i*m2+m2/2+(j-width/2)].colonize();
                    
                    for(it = propagule.begin();it!=propagule.end();it++)
                    {
                        demes[(j+1)*m2-1 -i].addMigrant(*it);
                    }
                }
        }  
   }
   
      
   if (mode == 3)     // colonize a single deme at opposing corners of the habitat
   {
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
            for(int i = 0;i < m1;i++)    
            {
                demes[i*m2].addMigrant((*it));
            }  
        }
        
        
        propagule = initial_population[0].sampleIndividuals(capacity);    
        for(it = propagule.begin();it!=propagule.end();it++)
        {
            for(int i = 0;i < m1;i++)    
            {
                demes[(i+1)*m2-1].addMigrant((*it));
            }  
        }
   }
   for(i=0;i<number_demes;i++) {demes[i].setID(i); }  
}


void World::reproduceBurnin(int mode, double phi)
{
    vector<Deme>::iterator it;
    //updateWaveFront();
    if (mode == 0)
    {
        reproduceSSburnin(phi); 
    }
    else
    {   
        reproduceHSburnin(phi);
    }
   
    
}

void World::reproduce(int mode)
{
    vector<Deme>::iterator it;
    //updateWaveFront();
    if (mode == 0)
    {
        reproduceSS(); 
    }
    else
    {   
        reproduceHS1();
    }
   
    
}

void World::reproduceSS()
{
    vector<Deme>::iterator it;
    //updateWaveFront();

    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->reproduceSS(wavefrontID); 
    }
}


void World::reproduceSSAM()
{
    vector<Deme>::iterator it;
    //updateWaveFront();
    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->reproduceSSAM(wavefrontID); 
    }
}


void World::reproduceHS1()
{
    vector<Deme>::iterator it;
    double mean_fit;
    //updateWaveFront();
    for(it = demes.begin();it!=demes.end();it++)  
    {
        mean_fit = it->getMeanFit();
        it->reproduceHS1(mean_fit,wavefrontID); 
    }
}

void World::reproduceSSburnin(double phi)
{
    vector<Deme>::iterator it;
    //updateWaveFront();

    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->reproduceSSburnin(wavefrontID, phi); 
    }
}
void World::reproduceHSburnin(double phi)
{
    vector<Deme>::iterator it;
    double mean_fit;
    //updateWaveFront();
    for(it = demes.begin();it!=demes.end();it++)  
    {
        mean_fit = it->getMeanFit();
        it->reproduceHSburnin(mean_fit,wavefrontID, phi); 
    }
}


void World::migrate(int range)     
{ 
    vector<Deme>::iterator it;
    Migrants::iterator m_it;
    Individual ind;
   
    int i=0;
    int j;
    int width;
    int mig_distance,destination;   // number of demes an individual migrates and location to which an individual migrates
   
    width = floor(range/m1);
   //get a vector of migrants, migrants[i] contains the emigrants of deme i -- doesn't seem to let me iterate specific i's out of order
   it = demes.begin();
   for(i = 0; i<(m1*m2); i++)
   {
       migrants[i] = it->getMigrants();
       it++;
   }

    //distribute the emigrants according to migration pattern
    for(i=0;i<width;i++)
    {
        for(j=0;j<m1;j++)
        {
            for(m_it = migrants[(j*m2+i)].begin();m_it!=migrants[(j*m2+i)].end();m_it++)  
            {
                destination = (j*m2+i)+randint(0,1)*2-1;              // nearest neighbor migration: migrate to deme i - 1 or i + 1, each with prob 1/2
              
                //destination = (j*m2+i) + (randint(0,1)*2-1) * (1+(int)(randexp(1/5)));                 // Poisson distribution  
              
                // 2D grid 
                mig_distance = randint(0,1);
              
                mig_distance = (mig_distance * 2 -1);
              
                if(randreal(0,1)<0.5)
                {
                    mig_distance =  mig_distance*m2;
                }
              
                destination = (j*m2+i) + mig_distance;
              
                // reflecting boundaries     
                if (destination<0) {destination=(j*m2+i);}                                      // don't go negative, i.e. below deme 0 or off the top edge                        
                if (destination>=(m1*m2)) {destination=(j*m2+i);}                               // don't go beyond last deme on landscape
                if ((destination>=(j*m2+width))&&(mig_distance<m2)) {destination = (j*m2+i);}   // don't migrate beyond the burn-in barrier or off the bottom edge
                if ((((j*m2+i)%m2)==0)&&(mig_distance==-1)) {destination = (j*m2+i);}           // don't migrate from the left edge core to the right edge
                if ((((j*m2+i)%m2)==(m2-1))&&(mig_distance==1)) {destination = i;}              // don't migrate from the right edge back into the core

                //cout << "final dest " << destination << endl;
                ind = *m_it;
                demes[destination].addMigrant(*m_it);
            }
        }
    }
}


void World::select()
{
    vector<Deme>::iterator it;

    for(it = demes.begin();it!=demes.end();it++)  {it->select(); }
}
   

void World::print()
{
    vector<Deme>::iterator it;

    
    for(it = demes.begin();it!=demes.end();it++)  {it->print(); }
}


void World::printStat()
{
    vector<Deme>::iterator it;

    
    for(it = demes.begin();it!=demes.end();it++)  {it->printStat(); }
}


vector<double> World::getMeanFit()                                              // returns vector with mean fitness of all demes
{
    vector<Deme>::iterator it;
    vector<double> data;
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        data.push_back(it->getMeanFit()); 
    }
    
    return(data);
}

vector<double> World::getVarFit(vector<double> mean_fit)                                               // retuns vector with variance in fitness of all demes
{
    vector<Deme>::iterator it;
    vector<double> data;
    int i = 0;
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        data.push_back(it->getVarFit(mean_fit[i++])); 
    }
    
    return(data);
}

// KJG adding
///*
vector<double> World::getDemeDensity()                                              // returns vector with pop density of all demes
{
    vector<Deme>::iterator it;
    vector<double> data;
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        data.push_back(it->getDemeDensity());   // go to deme.cpp for this function
    }
    
    return(data);
}
//*/
// KJG adding
///*
vector<double> World::getEdgeDemes(int edgeColumn)         // returns vector with deme IDs of all demes in column x across m2's length of landscape
{
    vector<double> edgeDemes(m1);   // make the vector as long as the width of the landscape
    int landscapeWidth = m1;    // the "height" or "width" of the landscape - expansion happens over m2's axis
    int landscapeLength = m2;   // the length of the landscape over which expansion happens
        
    for(int i=0; i<landscapeWidth; i++)
    {
        edgeDemes[i] = edgeColumn + (i * landscapeLength);
    }    
    
    return(edgeDemes);
}
//*/

vector<double> World::getHeterozygosity(int loci_begin,int loci_end)                                               // returns vector with variance in fitness of all demes
{
    vector<Deme>::iterator it;
    vector<int>::iterator int_it;
    vector<double> data;
    vector<int> al1,al2,al3,al4; 
    int l1,l2;
    
    asc_loci.clear();
    //asc_loci = demes[m1/2*m2+1].getAscLoci();
    
    
    //ascertainment variant 1
    //al1 = demes[m1/2*m2+6].getAscLoci(loci); 
    //al2 = demes[(m1/2-1)*m2+6].getAscLoci(loci);
    //al3 = demes[(m1/2+1)*m2+6].getAscLoci(loci);
    
    // ascertainment variant 2
    //al1 = demes[m1/2*m2+2].getAscLoci(loci); 
    //al2 = demes[m1/2*m2+10].getAscLoci(loci); 
    //al3 = demes[m1/2*m2+50].getAscLoci(loci); 
    
    // ascertainment variant 3
    //al1 = demes[m1/2*m2+6].getAscLoci(loci); 
    //al2 = demes[m1/2*m2+25].getAscLoci(loci); 
    //al3 = demes[m1/2*m2+50].getAscLoci(loci);
    
    // ascertainment variant 4
    //al1 = demes[(m1/2+2)*m2+25].getAscLoci(loci); 
    //al2 = demes[(m1/2+2)*m2+50].getAscLoci(loci); 
    ////al3 = demes[(m1/2+2)*m2+50].getAscLoci(loci);
    
     // ascertainment variant 5
//    al1 = demes[(m1/2+2)*m2+25].getAscLoci(loci); 
//    al2 = demes[(m1/2+2)*m2+50].getAscLoci(loci); 
//    al3 = demes[(m1/2+2)*m2+75].getAscLoci(loci);
    
//    set_intersection(al1.begin(),al1.end(), al2.begin(),al2.end(), std::back_inserter(asc_loci)); 
//        
//    al1.clear();
//    
//    set_intersection(asc_loci.begin(),asc_loci.end(), al3.begin(),al3.end(), std::back_inserter(al1)); 
//    asc_loci = al1;
    
//    ascertainment variant 6
     asc_loci = demes[(m1/2)*m2+10].getAscLociSample(loci_begin,loci_end,20);
     
     
    //ascertainment variant 7
    //al1 = demes[m1/2*m2+10].getAscLoci(loci); 
    //al2 = demes[m1/2*m2+2].getAscLoci(loci);
    //set_union(al1.begin(),al1.end(), al2.begin(),al2.end(), std::back_inserter(asc_loci)); 
    
    
    // ascertainment variant 8
//    al1 = demes[m1/2*m2+10].getAscLoci(loci); 
//    al2 = demes[m1/2*m2+2].getAscLoci(loci);
//    
//    l1 = al1.size();
//    l2 = al2.size();
//    
//    
//    
//    set_union(al1.begin(),al1.end(), al2.begin(),al2.end(), std::back_inserter(asc_loci)); 
        

    for(it = demes.begin();it!=demes.end();it++)  
    {
        data.push_back(it->getHeterozygosity(asc_loci,loci_begin,loci_end)); 
    }
    
    return(data);
}

vector<double> World::getAlleleFrequenciesWF(int loci_begin,int loci_end,int age)                                        
{
    vector<Deme>::iterator it;
    vector<double> data;
    vector<double> frequencies;
    int wf_size = 0;

    int temp_loci = loci_end-loci_begin;
    data.resize(temp_loci);
    fill_n(data.begin(),temp_loci,0);
    frequencies.resize(temp_loci);


    for(it = demes.begin();it!=demes.end();it++)  
    {
        if (it->getAge() < age && it->getAge() > 2)
        {
                frequencies = (it->getFrequencies(loci_begin,loci_end)); 
                for(int i=0;i<temp_loci;i++) 
                { 
                    data[i] = data[i] + frequencies[i];
                }
                wf_size++;
        }
    }
    
         
    for(int i=0;i<temp_loci;i++) 
    { 
        data[i]=data[i]/wf_size;
    }

    return(data);
}

vector<double> World::getAlleleFrequencies(int loci_begin,int loci_end)                                        
{
    vector<Deme>::iterator it;
    vector<double> data;
    vector<double> frequencies;
    int deme_counter = 0;
    

    int temp_loci = loci_end-loci_begin;
    data.resize(temp_loci*m1*m2);
    frequencies.resize(temp_loci);


    for(it = demes.begin();it!=demes.end();it++)  
    {
        frequencies = (it->getFrequencies(loci_begin,loci_end)); 
                
        for(int i=0;i<temp_loci;i++) 
        { 
            data[(deme_counter*temp_loci)+i] = frequencies[i];      
        }
        deme_counter++;        
    }

    return(data);
}

vector<double> World::getGenotypeFrequencies(int loci_begin,int loci_end,int genotype)                                        
{
    vector<Deme>::iterator it;
    vector<double> data;
    vector<double> frequencies;
    int deme_counter = 0;
    

    int temp_loci = loci_end-loci_begin;
    data.resize(temp_loci*m1*m2);
    
    frequencies.resize(temp_loci);

    for(it = demes.begin();it!=demes.end();it++)  
    {
        frequencies = (it->getGenotypeFrequencies(loci_begin,loci_end,genotype)); 
                
        for(int i=0;i<temp_loci;i++) 
        { 
            data[(deme_counter*temp_loci)+i] = frequencies[i];
        }
        deme_counter++;        
    }

    return(data);
}

vector<double> World::getInversionFrequency()                                        
{
    vector<Deme>::iterator it;
    vector<double> data;
    double frequency;
    int deme_counter = 0;
    
    data.resize(m1*m2);

    for(it = demes.begin();it!=demes.end();it++)  
    {
        frequency = it->getInversionFrequency(); 
        data[deme_counter] = frequency;
        deme_counter++;        
    }
   
    return(data);
}


void World::setParams(int K,double mu,double s,double m)
{
    demes[0].setParams(K,mu,s,m);
}


void World::setCapacity(int K)
{
    vector<Deme>::iterator it;
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->setParams(K);
    }
}

/*void World::updateWaveFront()
{
    int i = 0;
    int wf=0;
    
    for(i=0;i<number_demes;i++)
    {
        if(demes[i].getSize()>0)
        {wf++;}
    }
    
    //while (demes[i].colonized())
    //{
        
     //  i++;
    
    //}
    
    wavefrontID = wf-1;
    //cout<< '\n' << wavefrontID;
    
    //updateDistance();
    
}

*/



/*vector<Count> World::getStatMut()
{
    vector<Deme>::iterator it;
    vector<Count> data;
    
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        data.push_back(it->getStatMut()); 
    }
    
    return(data);
}*/


/*void World::ResetMutationOrigin()
{
    vector<Deme>::iterator it;
    
    
    for(it = demes.begin();it!=demes.end();it++)  
    {
        it->ResetMutationOrigin(); 
    } 
    
}*/

double World::sample_wfID(int location)
{
    return(demes[location].sample_wfID(10));
}

bool World::isColonized(int deme)
{
    return(demes[deme].colonized());
}

