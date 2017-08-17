/* 
 * File:   range_expansion.h
 * Author: stephan peischl
 *
 * Created on March 28, 2012, 5:38 PM
 */

#ifndef _RANGE_EXPANSION_H_
#define _RANGE_EXPANSION_H_


#include <vector>
#include <list>
#include <map>
#include "rng2.h"

using std::vector;
using std::list;
using std::map;

struct heritableUnit;
class World;
class Deme;
class Individual;


typedef map<int,float> MutationBag;

typedef vector<bool> Loci; 

typedef vector<MutationBag> MutationHaploList;

typedef vector<int> Distance;

typedef vector<double> Count;

typedef vector<Individual> Migrants; 



struct heritableUnit
{
    
 
        Loci haplotype;
        //Distance d;
        //MutationHaploList muts;
                
       /* Count m_d;          // total number of deleterious mutations
        Count md_front;    //  number of mutations that originated at the front
        
        Count m_b;         // same for beneficial ones
        Count mb_front;*/
         
        //double ass_mat;
        //double mig;
            
};



class World
{
 private:
         vector<Deme> demes; 
         vector<Deme> initial_population;
         vector<Migrants> migrants;
         static int number_demes;
         static int m1;
         static int m2;
         static int colonized_demes;
         int wavefrontID;
         int deltaWF;    
         vector<int> asc_loci;
         static double phi;
         
 public:
        World();
        World(int length1,int length2,int initial_colonized,int initial_popsize,int burnin_time,int capacity,int mode,double mu, double s, double m,double mut_prop);                
        ~World();
        bool isColonized(int deme);
        void clear(int length1,int length2,int initial_colonized,int initial_popsize,int burnin_time,int capacity,int mode,double mutation_rate,double s,double migration_rate,double mut_prop);
        void reproduceBurnin(int mode, double phi);
        void reproduceSSburnin(double phi);                     //reproduction plus soft selection
        void reproduceHSburnin(double phi);                    // reproduction plus hard selection v1 - growth rate proportional to mean fit
        void reproduce(int mode);
        void reproduceSS();                     //reproduction plus soft selection
        void reproduceSSAM(); 
        void reproduceHS1();                    // reproduction plus hard selection v1 - growth rate proportional to mean fit
        void select();
        void migrate(int range);
        void print();
        void printStat();
        void changeEnvironment(int demeID,double s0, double s1);
        vector<double> getMeanFit();
        vector<double> getDemeDensity(); // added KJG
        vector<double> getEdgeDemes(int edgeColumn); // added KJG
        vector<double> getHeterozygosity(int loci_begin,int loci_end);
        vector<double> getAlleleFrequenciesWF(int loci_begin,int loci_end, int age);     // returns vector with allele frequencies at the wave front (demes that are "age" generations old)
        vector<double> getAlleleFrequencies(int loci_begin,int loci_end);
        vector<double> getGenotypeFrequencies(int loci_begin,int loci_end,int genotype);
        vector<double> getInversionFrequency();
        vector<double> getVarFit(vector<double> mean_fit);   
        vector<double> getStatDist();
        //vector<Count> getStatMut();
        void setParams(int K,double mu, double s,double m);
        void setCapacity(int K);
        //void updateWaveFront();
        void updateDistance();
        void setDemeCapacity(int bn_location, int capacity);
        //void ResetMutationOrigin();
        double sample_wfID(int location);
      
};


class Deme
{
 private:
         list<Individual> this_generation; 
         list<Individual> next_generation;
         static double m;
         int capacity; 
         static double s;                                                                 
         static double mutation_rate;
         double max_fit;
         int ID;
         int age;
         
 public:
        Deme();
        ~Deme();
        void initialize();
        void colonize();     
        void set_selection_dist(double mut_prop);
        void reproduce(int wf); 
        void reproduceSS(int wf); 
        void reproduceSSAM(int wf); 
        void reproduceHS1(double mean_fit,int wf);
        void reproduceSSburnin(int wf,double phi); 
        void reproduceHSburnin(double mean_fit,int wf,double phi);
        void select();
        void migrate();
        void print();
        Migrants getMigrants();
        Migrants sampleIndividuals(int samplesize);
        void addMigrant(Individual);
        void printStat();
        double getMeanFit();
        double getDemeDensity(); // added KJG
        double getEdgeDemes(int edgeColumn); // added KJG
        double getHeterozygosity(vector<int> a_loci, int loci_begin,int loci_end);
        vector<double> getFrequencies(int loci_begin,int loci_end);
        vector<double> getGenotypeFrequencies(int loci_begin,int loci_end,int genotype);
        double getInversionFrequency();
        double getVarFit(double mean_fit);
        double getStatDist();
        int getAge();
        void setParams(int K,double mig, double s,double mu);
        void setParams(int K);
        void setID(int i);
        bool colonized();
        int getSize();
        //Count getStatMut();
        //void ResetMutationOrigin();
        double sample_wfID(int max_age);
        void normalizeFitness();
        vector<int> getAscLoci(int loci_begin,int loci_end);
        vector<int> getAscLociSample(int loci_begin,int loci_end,int n);
        
        
};



class Individual
{
 private:
         vector<Loci> haplotypes;
         //vector<MutationHaploList> mutations;
         
         //FitnessComponents haplotypes[2];
         //MutationHaploList mutations[2];
         
         //vector<Count> mutations_d;
         //vector<Count> md_front;
         //vector<Count> mutations_b;
         //vector<Count> mb_front;
         static int loci;
         static double rrate;                                                   // recombination rate between two consecutive loci, (no rec =) 0 <= rrate <= 0.5 (= free rec) 
         static vector<int> used_loci;
         static vector<float> s_coeff;
         double wf_ID;
         int ancestors;
         
 public:
        Individual();
        ~Individual();
        void set_selection_dist(double s,double mut_prop);
        heritableUnit getNewGamete(double mu,double s,bool front);         // infinite sites within recombining regions
        heritableUnit getNewGameteBurnin(double mu,double s,double phi);   // ignores beneficial mutations
        heritableUnit getNewGameteMM2(double mu1,double mu2,double s);     // mutation model with 2 alleles per locus
        void setGenotype(heritableUnit g1,heritableUnit g2);
        double getFitness(double s);
        double getRelativeFitness(double s);
        double getMaxFitness(double s);
        double getMeanDistance();
        void print();
        void setParams(int loci);
        void setAncestors(int a);
        void setWFID(double id);
        int getAncestors();
        vector<double> getSumAlleles(int loci_begin,int loci_end);
        vector<double> getSumGenotypes(int loci_begin,int loci_end,int genotypes);
        double getSquareSumAlleles();
        double getWFID();
        void normalizeFitness(double co);
        unsigned long getNumberMutations();
        double getInversionCount();
        
        //Count getMutationCount();
        //void ResetMutationOrigin();
        
        

        
};

//#ifdef	__cplusplus
//extern "C" {
//#endif
//
//#ifdef	__cplusplus
//}
//#endif

#endif	/* _RANGE_EXPANSION_H_ */

