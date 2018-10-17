#include "core_evolution.hpp"
#include <algorithm>

std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  std::vector<uint16_t> selected_indices(fitnesses.size());
  std::uniform_real_distribution<double> random_interval(0,fitnesses.back());
  for(auto& sel_val : selected_indices) 
    sel_val=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(RNG_Engine))-fitnesses.begin());
  std::sort(selected_indices.begin(),selected_indices.end());
  return selected_indices;
}

void FitnessPhenotypeTable::RelabelPhenotypes(std::vector<Phenotype_ID >& pids,std::map<Phenotype_ID, std::set<interaction_pair> >& p_ints)  {
  for(auto x_iter=new_phenotype_xfer.begin();x_iter!=new_phenotype_xfer.end();x_iter+=3) {
    p_ints[std::make_pair(*x_iter,*(x_iter+2))].insert(p_ints[std::make_pair(*x_iter,*(x_iter+1))].begin(),p_ints[std::make_pair(*x_iter,*(x_iter+1))].end());
    phenotype_fitnesses[*x_iter].emplace_back(std::gamma_distribution<double>(*(x_iter)*2,.5*std::pow(*x_iter,-.5))(RNG_Engine));
  }
  PhenotypeTable::RelabelPhenotypes(pids);
}

double FitnessPhenotypeTable::GenotypeFitness(std::map<Phenotype_ID,uint16_t> ID_counter) {
  double fitness=0;
  for(auto kv : ID_counter)
    if(kv.second>=ceil(model_params::UND_threshold*model_params::phenotype_builds))
      fitness+=phenotype_fitnesses[kv.first.first][kv.first.second] * std::pow(static_cast<double>(kv.second)/model_params::phenotype_builds,simulation_params::fitness_factor);     
  return fitness;
}

double FitnessPhenotypeTable::SingleFitness(Phenotype_ID pid,uint16_t commonness)  {
  return phenotype_fitnesses[pid.first][pid.second] * std::pow(static_cast<double>(commonness)/model_params::phenotype_builds,simulation_params::fitness_factor);     
}
