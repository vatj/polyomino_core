#include "core_genotype.hpp"
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <cmath>


/*! free vs one-sided polyominoes and tile vs orientation determinism */
constexpr uint8_t FREE_POLYOMINO = true ? 2 : 1, DETERMINISM_LEVEL=3;
/*! 
  (true) free polyominoes are not chirally distinct
  (false) one-sided polyominoes are

determinism levels as follows:
    shape       : 1
    tile        : 2
    orientation : 3
*/

/*! phenotype definitions */
struct Phenotype {
  uint8_t dx=1,dy=1;
  std::vector<uint8_t> tiling{1};
};

typedef std::pair<uint8_t,uint16_t> Phenotype_ID;
constexpr Phenotype_ID NULL_pid{0,0};

/*! print phenotype to stdout */
void PrintShape(Phenotype& phen);

/*! phenotype rotations */
void ClockwiseRotation(Phenotype& phen);
void ClockwisePiRotation(Phenotype& phen);
void ChiralFlip(Phenotype& phen);

/*! phenotype comparisons*/
bool ComparePolyominoes(Phenotype& phen1, const Phenotype& phen2);
void MinimizePhenRep(std::vector<uint8_t>& tiling);
void GetMinPhenRepresentation(Phenotype& phen);

namespace model_params
{
  extern bool FIXED_TABLE;
  extern double UND_threshold;
  extern uint16_t phenotype_builds;
}

struct PhenotypeTable {
  std::unordered_map<uint8_t,std::vector<Phenotype> > known_phenotypes,undiscovered_phenotypes;
  std::vector<size_t> new_phenotype_xfer;
  std::unordered_map<uint8_t, std::vector<uint16_t> > undiscovered_phenotype_counts;
  
  Phenotype_ID GetPhenotypeID(Phenotype& phen) {
    uint8_t phenotype_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;});
    
    /*! compare against existing table entries*/
    for(uint16_t phenotype_index=0; phenotype_index != known_phenotypes[phenotype_size].size();++phenotype_index)
      if(ComparePolyominoes(phen,known_phenotypes[phenotype_size][phenotype_index])) 
	return std::make_pair(phenotype_size,phenotype_index);
    if(model_params::FIXED_TABLE)
      return NULL_pid;
    /*! compare against temporary table entries*/
    uint16_t new_phenotype_index=0;
    for(Phenotype phen_p : undiscovered_phenotypes[phenotype_size]) {
      if(ComparePolyominoes(phen,phen_p)) {
	if(++undiscovered_phenotype_counts[phenotype_size][new_phenotype_index]>=std::ceil(model_params::UND_threshold*model_params::phenotype_builds)) {
	  known_phenotypes[phenotype_size].emplace_back(phen);
          new_phenotype_xfer.insert(new_phenotype_xfer.end(),{phenotype_size,known_phenotypes[phenotype_size].size()-1+new_phenotype_index+model_params::phenotype_builds,known_phenotypes[phenotype_size].size()-1});            
	  return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()-1);
	}
	else
	  return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()+new_phenotype_index+model_params::phenotype_builds);
      }
      ++new_phenotype_index;
    }

    /*! brand new phenotype, get minimal representation and add to temporary table (or existing if limit was 1) */
    GetMinPhenRepresentation(phen);
    if(static_cast<uint16_t>(std::ceil(model_params::UND_threshold*model_params::phenotype_builds))<=1) {
      known_phenotypes[phenotype_size].emplace_back(phen);
      return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()-1);
    }
    
    undiscovered_phenotypes[phenotype_size].emplace_back(phen);
    undiscovered_phenotype_counts[phenotype_size].emplace_back(1);
    return std::make_pair(phenotype_size,known_phenotypes[phenotype_size].size()+new_phenotype_index+model_params::phenotype_builds);
  }


  void RelabelPhenotypes(std::vector<Phenotype_ID >& pids) {
    /*! relabels stored in tuple (size, swap_from, swap_to) */
    for(auto x_iter=new_phenotype_xfer.begin(); x_iter!=new_phenotype_xfer.end();x_iter+=3)
      std::replace(pids.begin(),pids.end(),std::make_pair(static_cast<uint8_t>(*x_iter),static_cast<uint16_t>(*(x_iter+1))),std::make_pair(static_cast<uint8_t>(*x_iter),static_cast<uint16_t>(*(x_iter+2))));
	
    undiscovered_phenotypes.clear();
    undiscovered_phenotype_counts.clear();
    new_phenotype_xfer.clear();
  }

  /* Count each ID frequency */
  std::map<Phenotype_ID,uint16_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids, bool& rare_phenotypes) {
    std::map<Phenotype_ID, uint16_t> ID_counter;
    for(std::vector<Phenotype_ID >::const_iterator ID_iter = pids.begin(); ID_iter!=pids.end(); ++ID_iter) {
      if(ID_iter->second < known_phenotypes[ID_iter->first].size())
	++ID_counter[std::make_pair(ID_iter->first,ID_iter->second)];
      else
	rare_phenotypes=true;
    }
    return ID_counter;
  }
        
  void PrintTable(std::ofstream& fout) {
    for(auto known_phens : known_phenotypes) {
      uint16_t n_phen=0;
      for(Phenotype known : known_phens.second) {
	fout<<+known_phens.first<<" "<<+n_phen++<<" "<<+known.dx<<" "<<+known.dy<<" ";
	for(uint8_t tile : known.tiling)
	  fout<<+tile<<" ";
	fout<<"\n";
      }
    }
  }

};
