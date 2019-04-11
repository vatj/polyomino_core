#pragma once
#include <functional>
#include <algorithm>

using PopulationSize = uint16_t;

//simple extension of PhenotypeTable to include fitnesses accessed via pid
struct FitnessPhenotypeTable : PhenotypeTable {
  inline static double fitness_factor=1;

  std::unordered_map<uint8_t,std::vector<double> > phenotype_fitnesses{{0,{0}}};

  //user defined fitness function, by default some random variable based on size
  std::function<double(uint8_t)> fit_func=[](double s) {return std::gamma_distribution<double>(s*2,.5*std::pow(s,-.5))(RNG_Engine);};


  //add as many new fitnesses as newly discovered phenotypes
  inline void UpdateFitnesses() {
    for(const auto& kv : known_phenotypes)
      while(kv.second.size()>phenotype_fitnesses[kv.first].size())
        phenotype_fitnesses[kv.first].emplace_back(fit_func(kv.first));
  }

  //sum over all phenotypes according to weighted fitness
  inline double GenotypeFitness(std::map<Phenotype_ID,uint16_t> ID_counter) {
    double fitness=0;
    for(auto& kv : ID_counter)
      fitness+=phenotype_fitnesses[kv.first.first][kv.first.second] * std::pow(static_cast<double>(kv.second)/phenotype_builds,fitness_factor);
    return fitness;
  }

  //fitness based on a single given pid and its commonness
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
inline std::vector<PopulationSize> RouletteWheelSelection(std::vector<double>& fitnesses) {
  std::vector<PopulationSize> selected_indices(fitnesses.size());
  std::partial_sum(fitnesses.begin(), fitnesses.end(), fitnesses.begin());
  
  //no fitness in population, select from all equally
  if(fitnesses.back()==0.) {
    std::uniform_int_distribution<PopulationSize> dist{0, static_cast<PopulationSize>(fitnesses.size()-1)};
    std::generate(selected_indices.begin(),selected_indices.end(), [&dist](){return dist(RNG_Engine);});
  }
  //select proportional to weighted fitness intervals
  else{ 
    std::uniform_real_distribution<double> random_interval(0,fitnesses.back());
    for(auto& sel_val : selected_indices) 
      sel_val=static_cast<uint16_t>(std::lower_bound(fitnesses.begin(),fitnesses.end(),random_interval(RNG_Engine))-fitnesses.begin());
    std::sort(selected_indices.begin(),selected_indices.end());
  }
  
  return selected_indices;
}
