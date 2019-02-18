#pragma once
#include <functional>
#include <algorithm>

struct FitnessPhenotypeTable : PhenotypeTable {
  inline static double fitness_factor=1;

  std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};
  std::function<double(uint8_t)> fit_func;
  FitnessPhenotypeTable(void) {fit_func=[](double s) {return std::gamma_distribution<double>(s*2,.5*std::pow(s,-.5))(RNG_Engine);};};

  inline void UpdateFitnesses() {
	for(const auto& kv : known_phenotypes)
	while(kv.second.size()>phenotype_fitnesses[kv.first].size())
	  phenotype_fitnesses[kv.first].emplace_back(fit_func(kv.first));
  }

  inline double GenotypeFitness(std::map<Phenotype_ID,uint16_t> ID_counter) {
    double fitness=0;
    for(auto& kv : ID_counter)
        fitness+=phenotype_fitnesses[kv.first.first][kv.first.second] * std::pow(static_cast<double>(kv.second)/phenotype_builds,fitness_factor);
    return fitness;
  }
  
  inline double SingleFitness(Phenotype_ID pid,uint16_t commonness) {
    return phenotype_fitnesses[pid.first][pid.second] * std::pow(static_cast<double>(commonness)/phenotype_builds,fitness_factor);     
  }
  
  inline void LoadTable(std::string f_name) {
    PhenotypeTable::LoadTable(f_name);
    for(auto& kv : known_phenotypes)
      phenotype_fitnesses[kv.first].insert(phenotype_fitnesses[kv.first].end(),kv.second.size(),0); 
  }  
    
};

//fitness proportional selection, or equal selection if net zero fitness
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
