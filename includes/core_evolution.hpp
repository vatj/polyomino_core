#include "core_phenotype.hpp"
#include <random>

namespace simulation_params {
  extern uint16_t population_size;
  extern double fitness_factor;
}

extern thread_local std::mt19937 RNG_Engine;

struct FitnessPhenotypeTable : PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
         
  void RelabelPhenotypes(std::vector<Phenotype_ID >& pids,std::map<Phenotype_ID, std::set<interaction_pair> >& p_ints);
  std::map<Phenotype_ID,uint16_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids);
  double GenotypeFitness(std::map<Phenotype_ID,uint16_t> ID_counter);
  double SingleFitness(Phenotype_ID pid,uint16_t commonness);

};

std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);
