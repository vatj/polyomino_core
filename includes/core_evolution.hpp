#pragma once

#include "core_assembly.hpp"
#include "core_phenotype.hpp"

namespace simulation_params {
  extern double fitness_factor; 
}

struct FitnessPhenotypeTable : PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0,0}}};
         
  void RelabelPhenotypes(std::vector<Phenotype_ID >& pids,std::map<Phenotype_ID, std::set<interaction_pair> >& p_ints);
  double GenotypeFitness(std::map<Phenotype_ID,uint16_t> ID_counter);
  double SingleFitness(Phenotype_ID pid,uint16_t commonness);

};

std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses);




