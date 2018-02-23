#include <iostream> // For input/output
#include <fstream> // For file input/output
#include <string>  // For strcpy
#include <time.h>  // For time
#include <stdlib.h>  // For toupper and tolower
#include <math.h>
#include <vector>
#include <list>
#include "range_expansion.h"
//#include "rng.h"
#include "rng2.h"

using namespace std;

int Individual::loci = 2000;
double Individual::rrate = 0.5; 
vector<int> Individual::used_loci;
vector<float> Individual::s_coeff;
//double phi = 0.9;                   // right now it is hard-coded that every mutation is deleterious
//unsigned long long Individual::counter = 0;

//inline double rand_unif(double x0, double x1)
//{
//return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);
//}

inline int max(int a, int b) { return (a < b) ? b : a; }
inline int min(int a, int b) { return (a < b) ? a : b; }


Individual::Individual()
{
    haplotypes.resize(2);
    //mutations.resize(2);
    
    used_loci.resize(loci);
    
    for(int i = 0; i < loci; i++)
    {
        used_loci[i] = i;
    }
    
    wf_ID = 1;
    
    ancestors = 0;

    /*mutations_b.resize(2);
    mb_front.resize(2);
    mutations_d.resize(2);
    md_front.resize(2);*/
    
    haplotypes[0].resize(loci);
    haplotypes[1].resize(loci);
    
    //mutations[0].resize(loci);
    //mutations[1].resize(loci);
    
    /*mutations_b[0].resize(loci);
    mb_front[0].resize(loci);
    mutations_d[0].resize(loci);
    md_front[0].resize(loci);
 
    mutations_b[1].resize(loci);
    mb_front[1].resize(loci);
    mutations_d[1].resize(loci);
    md_front[1].resize(loci); */
    
    fill_n(haplotypes[0].begin(),loci,0);
    fill_n(haplotypes[1].begin(),loci,0);

    
   /* fill_n(mutations_d[0].begin(),loci,0);
    fill_n(md_front[0].begin(),loci,0);
    fill_n(mutations_b[0].begin(),loci,0);
    fill_n(mb_front[0].begin(),loci,0);
    
    fill_n(mutations_d[1].begin(),loci,0);
    fill_n(md_front[1].begin(),loci,0);
    fill_n(mutations_b[1].begin(),loci,0);
    fill_n(mb_front[1].begin(),loci,0);*/
}


Individual::~Individual()
{
    
}


heritableUnit Individual::getNewGamete(double mu,double s,bool front)           // this function adds mutations to the genome. mutations and backmutations occur at the same rate. (this will have an effect on the DFE, need to investigate this)
{
    Loci hap_new;
    double rec_rate = rrate; 
        
    heritableUnit gam_new;
    
   int i;
    
   int nmutations;
   int site = 0;

   
//   if(haplotypes[0][loci] ==  haplotypes[1][loci])   // last locus is inversion marker
//       rec_rate = rrate; 
   

   site = randint(0,1);
   hap_new = haplotypes[site];
   
   if(rec_rate > 0)
   {
       for (i = 0; i < (loci);i++)   // recombination
       {
           if(randreal(0,1) < rec_rate)
           {
               site = (site+1)%2;
           }
           hap_new[i]=haplotypes[site][i]; 
       } 
   }
  
//   if(randreal(0,1) < 0.01)
//   {
//       hap_new[loci] = !hap_new[loci];
//   }
   
    nmutations = 0;
    if (mu > 0) {nmutations = randpois(mu); }       // draw poisson distributed number of mutations around mean mu
                                                    // this is the number of mutations per gamete, so mu is indeed U, the genome-wide mutation rate per generation
    
    for(i = 0; i < nmutations; i++)
    {   
        site = randint(0,loci-1); 
        if(site < loci/2){
            hap_new[site] = 1; // when 2000 loci (last 1000 neutral), selected loci have no back mutation but neutral do
        }else{
            hap_new[site] = 1 - hap_new[site];  // when 2000 loci (last 1000 neutral), neutral muts have back-mut
            //!hap_new[site];    // doesn't seem to properly work for back mutations?t
       }
        //hap_new[site] = 1;//!hap_new[site]; // this is where back-mutation happens, to get rid of it, set it as 1 so it doens't back mutate
    } 
    
        gam_new.haplotype = hap_new;
   
        return(gam_new); 
}


heritableUnit Individual::getNewGameteBurnin(double mu,double s,double phi)           // this function adds mutations to the genome. mutations and backmutations occur at the same rate. (this will have an effect on the DFE, need to investigate this)
{
    Loci hap_new;
    double rec_rate = rrate; 
        
    heritableUnit gam_new;
    
   int i;
    
   int nmutations;
   int site = 0;

   
//   if(haplotypes[0][loci] ==  haplotypes[1][loci])   // last locus is inversion marker
//       rec_rate = rrate; 
   

   site = randint(0,1);
   hap_new = haplotypes[site];
   
   if(rec_rate > 0)
   {
       for (i = 0; i < (loci);i++)   // recombination
       {
           if(randreal(0,1) < rec_rate)
           {
               site = (site+1)%2;
           }
           hap_new[i]=haplotypes[site][i]; 
       } 
   }
   
    nmutations = 0;
    //double scaledMu = mu/(phi*4);   // IMPORTANT - ONLY FOR 2000 LOCI VERSION - because in those inputs I halved phi AND doubled U, so here need to quadruple it to make same expected U just for 900 delet muts of 0.1111  (U/phi = 0.1/0.9 = 0.1111)
    double scaledMu = mu/((phi*2)+1);     // VERSION WITH NEUTRALS IN THE BURNIN PHASE - I want U=0.1 everywhere for 1000 loci, so U=0.2 for 2000 loci, so for 1900 loci, U=0.19 (U/phi = 0.2/0.9 = 0.1111) -- these rescalings may actually be wrong? because I want to downscale mu when I don't have as much of the genome that is mutating. But I don't htink it's an issue in general since it's only the burnin phase, and upping the mutation rate then helps to get to equil faster
    if (scaledMu > 0) {nmutations = randpois(scaledMu); }       // draw poisson distributed number of mutations around mean mu that is scaled down by the number of loci that are deleterious
                                                    // this is the number of mutations per gamete, so mu is indeed U, the genome-wide mutation rate per generation
    
    for(i = 0; i < nmutations; i++)
    {   
        // site = randint(0,(phi*loci)-1); // original for only counting delets during the burnin
        site = randint(0,((phi*loci)+(loci/2))-1);
        if(site > (phi*loci-1)) site = site + 100;
        
        if(site < loci/2){
            hap_new[site] = 1; // when 2000 loci (last 1000 neutral), selected loci have no back mutation but neutral do
        }else{
            hap_new[site] = 1 - hap_new[site];  // when 2000 loci (last 1000 neutral), neutral muts have back-mut
            //!hap_new[site];    // doesn't seem to properly work for back mutations?
        }
        //hap_new[site] = 1;//!hap_new[site]; // this is where back-mutation happens, to get rid of it, set it as 1 so it doens't back mutate
    } 
    
        gam_new.haplotype = hap_new;
   
        return(gam_new); 
}

heritableUnit Individual::getNewGameteMM2(double mu1,double mu2,double s)
{
    Loci hap_new;
    vector<int>::iterator it;
    heritableUnit gam_new;
    
    int i;
   
    int site = 0;
    site = randint(0,1);
    
    hap_new = haplotypes[0];
      
    for (i = 0; i < loci;i++)   // recombination
    {
        if(randint(0,1)<rrate)   
        {
            site = (site+1)%2;
        }
        
        hap_new[i]=haplotypes[site][i]; 
        
//         if (hap_new[i]==0)
//        {
//                if(randreal(0,1)<mu1) 
//                { 
//                        hap_new[i]=1;
//                }
//        }
    }
    
    for (it = used_loci.begin(); it < used_loci.end();)
    {
        if(randreal(0,1)<mu1) 
        { 
            hap_new[*it]=1;
            it = used_loci.erase(it); 
        }
        else    {it++; }
    }
//    
//        if (hap_new[i]==0)
//        {
//                if(randreal(0,1)<mu1) 
//                { 
//                        hap_new[i]=1;
//                }
//        }
//        
//        if (hap_new[i]==1)                                                    //back mutation
//        {
//            if(randreal(0,1)<(mu2))
//                { 
//                        hap_new[i]=0;
//                }
//        }
    
    gam_new.haplotype = hap_new;

    return(gam_new); 
}



void Individual::setGenotype(heritableUnit g1,heritableUnit g2)
{
    haplotypes[0] = g1.haplotype;
    haplotypes[1] = g2.haplotype;
    
//    mutations[0] = g1.muts;
//    mutations[1] = g2.muts;
    
    /*mutations_d[0] = g1.m_d;
    md_front[0] = g1.md_front;
    mutations_b[0] = g1.m_b;
    mb_front[0] = g1.mb_front;
    
    mutations_d[1] = g2.m_d;
    md_front[1] = g2.md_front;
    mutations_b[1] = g2.m_b;
    mb_front[1] = g2.mb_front;*/
}
 
 
 
double Individual::getFitness(double s)                                         // this is no longer used
{

}

double Individual::getRelativeFitness(double s)
{
    double w = 1;
    int i;
    float h = 0.3;                  // dominance is set here
//    float h;
    extern double fitnessConstant;  // this is the global variable for scaling fitness, defined in main.cpp
    
    // calculate fitness from genotype
    
    for (i=0;i<(loci/2);i++)                                                    // add up the effects of deleterious mutations        
    {
       
       if ((haplotypes[0][i] && !haplotypes[1][i])||(!haplotypes[0][i] && haplotypes[1][i]))                              // completely recessive mutations
       { 
           //h = 1/((1/0.5)-(2500*s_coeff[i])); // make h depend on the mutational effect size, params approx from Huber et al 2017 would be 0.99 and 20,000 but I changed to 0.5 because then neutral things are additive, and 2500 because then s-0.0005 gives h=0.307 and s=-0.005 gives h=0.068
           //if(s_coeff[i] > 0) h = 0.5;        // make all beneficials additive
           w *= (1 + h*s_coeff[i]);  
       }
       else if (haplotypes[0][i] && haplotypes[1][i])                              // completely recessive mutations
       { 
           w *= (1 + s_coeff[i]);  
       }
////                                                                                // mutliplicative effects  
//        if(haplotypes[0][i]) w *= (1-s);
//        if(haplotypes[1][i]) w *= (1-s); <- 
       
        
        
//         w *= (1-( (haplotypes[1][i]+haplotypes[0][i])* s_coeff[i]/2 ) );                       // additive effects
    }  
    


    return(w / fitnessConstant); //use the global variable to scale fitness (so that we can compare across h cases)
}



double Individual::getMaxFitness(double s)
{

}

void Individual::print()
{
    cout << "\n h1:";
    for(int i=0;i<haplotypes[0].size();i++) { cout << haplotypes[0][i] << " ";}
    cout << "\n h2:";
    for(int i=0;i<haplotypes[1].size();i++) { cout << haplotypes[1][i] << " ";}
}

void Individual::setParams(int number_loci)
{
    loci = number_loci;         // HERE KJG
    
    haplotypes[0].resize(loci);
    haplotypes[1].resize(loci);
    
 
    fill_n(haplotypes[0].begin(),loci,0);      // initialize haplotype 
    fill_n(haplotypes[1].begin(),loci,0);  
}


 




void Individual::setAncestors(int a)
{
    ancestors = a;
}
        

void Individual::setWFID(double id)
{
    wf_ID = id;
}


int Individual::getAncestors()
{
    return(ancestors);
}

double Individual::getWFID()
{
    return(wf_ID);
}


vector<double> Individual::getSumAlleles(int loci_begin,int loci_end)
{
    vector<double> p;
    
    
    p.resize(loci_end);
        
    fill_n(p.begin(),loci_end,0);
    
    for (int i = loci_begin;i<loci_end;i++)
    {
        p[i] += haplotypes[0][i]+haplotypes[1][i];
    }
    return(p);
}


vector<double> Individual::getSumGenotypes(int loci_begin,int loci_end,int genotype)
{
    vector<double> p;
    
    
    p.resize(loci_end-loci_begin);
        
    fill_n(p.begin(),loci_end-loci_begin,0);
  
    for (int i = loci_begin;i<loci_end;i++)
    {
        if(haplotypes[0][i]+haplotypes[1][i]==genotype)
            p[i] += 1;
    }

    return(p);
}

double Individual::getInversionCount()                                     // not used anymore
{
    return(haplotypes[0][loci]+haplotypes[1][loci]);
}

void Individual::normalizeFitness(double mean_fit)                              // not used anymore
{

}

unsigned long Individual::getNumberMutations()
{
 
}

void Individual::set_selection_dist(double s,double mut_prop)   // here we set the distribution of effect sizes for loci (bens and dels)
{       
    Individual::s_coeff.resize(loci);
    int i, j;
    
    for (i = 0;i<int(loci*mut_prop);i++)       // these are the deleterious mutations
    {
        s_coeff[i] = s;                                               // for constant fitness effects
        // s_coeff[i] = s * (-log( 1 - ((float)i/(loci*mut_prop)) )); // delet muts are now a negative exponential
    }
    
    j=0; // because we'll iterate i through the remaining loci that are bens, but j from 0 to number ben loci to get the correct quantile
    for (i = int(loci*mut_prop); i<loci; i++) // these are the beneficials
    {
        s_coeff[i] = -s;                                                            // for constant fitness effects
        //s_coeff[i] = -(s * (-log( 1 - ((float)j/(loci - (loci*mut_prop))) )));    // make beneficials reverse exponentially distributed from the negative s
        j=j+1;
    } 
 
}

//    THIS IS THE OLD CODE FOR MAKING THE DISTRIBUTION OF MUTATION EFFECTS EXPONENTIAL
//    for (i = 1;i<loci;i++)
//    {
//        s_coeff[i-1] = s;//*(-log(1-((float)i-1)/(loci)));
//    }
//    
//    s_coeff[loci-1] = s;//s*(-log(1-((float)loci-1)/loci));