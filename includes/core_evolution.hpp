#pragma once
#include <functional>
#include <algorithm>

namespace simulation_params {
  extern double fitness_factor; 
}

struct FitnessPhenotypeTable : PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
  std::function<double(uint8_t)> fit_func;
  FitnessPhenotypeTable(void) {fit_func=[](double s) {return std::gamma_distribution<double>(s*2,.5*std::pow(s,-.5))(RNG_Engine);};};
         
  inline void RelabelPhenotypes(std::vector<Phenotype_ID >& pids,std::map<Phenotype_ID, std::set<interaction_pair> >& p_ints);
  inline double GenotypeFitness(std::map<Phenotype_ID,uint16_t> ID_counter);
  inline double SingleFitness(Phenotype_ID pid,uint16_t commonness);
  inline void ExtendFitness();

};

inline std::vector<uint16_t> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::vector<uint16_t> selected_indices(fitnesses.size());
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  if(fitnesses.back()==0.) {
    std::iota(selected_indices.begin(),selected_indices.end(),0);
  }
  else{
    std::uniform_real_distribution<double> random_interval(0,fitnesses.back());
    for(auto& sel_val : selected_indices) 
      sel_val=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(RNG_Engine))-fitnesses.begin());
    std::sort(selected_indices.begin(),selected_indices.end());
  }
  return selected_indices;
}

void FitnessPhenotypeTable::ExtendFitness() {
  for(auto x_iter=new_phenotype_xfer.begin();x_iter!=new_phenotype_xfer.end();) {
    phenotype_fitnesses[*x_iter].emplace_back(fit_func(*x_iter));
    if(*(x_iter+1)==*(x_iter+2))
      x_iter=new_phenotype_xfer.erase(x_iter,x_iter+3);
    else
      x_iter+=3;
  }
}
 
void FitnessPhenotypeTable::RelabelPhenotypes(std::vector<Phenotype_ID >& pids,std::map<Phenotype_ID, std::set<interaction_pair> >& p_ints)  {
  for(auto x_iter=new_phenotype_xfer.begin();x_iter!=new_phenotype_xfer.end();x_iter+=3)
    p_ints[std::make_pair(*x_iter,*(x_iter+2))].insert(p_ints[std::make_pair(*x_iter,*(x_iter+1))].begin(),p_ints[std::make_pair(*x_iter,*(x_iter+1))].end());
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



