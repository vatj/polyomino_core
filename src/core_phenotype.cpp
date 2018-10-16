#include "core_phenotype.hpp"
#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

void PrintShape(Phenotype& phen) {
  for(uint8_t y=0;y<phen.dy;++y) {
    for(uint8_t x=0;x<phen.dx;++x)
      std::cout<<+phen.tiling[y*phen.dx+x]<<" ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

bool ComparePolyominoes(Phenotype& phen1, const Phenotype& phen2) {
  if(phen1.dy > phen1.dx)
    ClockwiseRotation(phen1);
  MinimizePhenRep(phen1.tiling);

  /*! test same size of phenotypes*/
  if(phen1.tiling.size()!=phen2.tiling.size() || std::count(phen1.tiling.begin(),phen1.tiling.end(),0)!=std::count(phen2.tiling.begin(),phen2.tiling.end(),0))
    return false; 
  
  /*! square phenotypes*/
  if(phen1.dx==phen2.dx && phen1.dy==phen2.dy && phen1.dx==phen2.dy) {
    for(uint8_t flip=0; flip<FREE_POLYOMINO;++flip) {
      if(phen1.tiling==phen2.tiling) 
        return true;
      for(int rotation=0;rotation<3;++rotation) {
        ClockwiseRotation(phen1);
        MinimizePhenRep(phen1.tiling);
        if(phen1.tiling==phen2.tiling) 
          return true;
      }
      if(flip==(FREE_POLYOMINO-1))
        return false; //square phenotype, but never matching
      ChiralFlip(phen1);
      MinimizePhenRep(phen1.tiling);
    }
  }
  /*! rectangular phenotypes*/
  if(phen1.dx==phen2.dx && phen1.dy==phen2.dy) {
    for(uint8_t flip=0; flip<FREE_POLYOMINO;++flip) {
      if(phen1.tiling==phen2.tiling)
        return true;
      ClockwisePiRotation(phen1);
      MinimizePhenRep(phen1.tiling);
      if(phen1.tiling==phen2.tiling)
        return true;
      if(flip==(FREE_POLYOMINO-1))
        return false; //rectangular phenotype, but never matching
      ChiralFlip(phen1);
      MinimizePhenRep(phen1.tiling);
    }
  }
  return false; //catch all
}

void ClockwiseRotation(Phenotype& phen) {
  std::vector<uint8_t> swapper;
  swapper.reserve(phen.tiling.size());
  for(uint8_t column=0;column<phen.dx;++column) 
    for(uint8_t row=phen.dy;row!=0;--row)  
      swapper.emplace_back(phen.tiling[(row-1)*phen.dx+column]);                       
  std::swap(phen.dx,phen.dy);
  phen.tiling=swapper;
}

void ClockwisePiRotation(Phenotype& phen) {
  std::reverse(phen.tiling.begin(),phen.tiling.end());
}

void ChiralFlip(Phenotype& phen) {
  for(uint8_t row=0;row<phen.dy;++row)
    std::reverse(phen.tiling.begin()+row*phen.dx,phen.tiling.begin()+(row+1)*phen.dx);
  if(DETERMINISM_LEVEL==3)
    for(uint8_t& element : phen.tiling)
      if(element && element%2==0)
	element+=-(element-1)%4+((element-1)%4+2)%4;
}

void MinimizePhenRep(std::vector<uint8_t>& tiling) {
  if(tiling.size()==1) {
    tiling={1};
    return;
  }
  if(DETERMINISM_LEVEL==1) {
    for(uint8_t& t : tiling)
      t= t ? 1:0;
    return;
  }
  for(uint8_t& t:tiling)
    t+=128*(t!=0);
  uint8_t swap_count=1;  
  for(std::vector<uint8_t>::iterator t_iter=std::find_if(tiling.begin(),tiling.end(),[](const int s) { return s>0; });t_iter!=tiling.end();) {
    const uint8_t static_swap=*t_iter;
    for(uint8_t cyclic = 0; cyclic< (DETERMINISM_LEVEL==3?4:1); ++cyclic) {
      const uint8_t pre_swap=DETERMINISM_LEVEL==3 ? (static_swap-(static_swap-1)%4)+((static_swap-1)%4+cyclic)%4 : static_swap;
      std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
      std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
      std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
      ++swap_count;
    }  
    t_iter=std::find_if(tiling.begin(),tiling.end(),[swap_count](const int s) { return (s>0 && (DETERMINISM_LEVEL==3?s-(s-1)%4:s) >= swap_count); });
  }
}

void GetMinPhenRepresentation(Phenotype& phen) {
  std::vector< std::vector<uint8_t> > min_tilings;
  if(phen.dy > phen.dx)
    ClockwiseRotation(phen);
  
  for(uint8_t rot=0;rot<FREE_POLYOMINO;++rot) {
    MinimizePhenRep(phen.tiling);
    min_tilings.emplace_back(phen.tiling);
      
    if(phen.dy != phen.dx) {
      ClockwisePiRotation(phen);
      MinimizePhenRep(phen.tiling);
      min_tilings.emplace_back(phen.tiling);
    }
    else {
      for(uint8_t rot=0;rot<3;++rot) {
        ClockwiseRotation(phen);
        MinimizePhenRep(phen.tiling);
        min_tilings.emplace_back(phen.tiling);
      }
    }
    if(rot==(FREE_POLYOMINO-1))
      break;
    ChiralFlip(phen);
  }
  std::nth_element(min_tilings.begin(),min_tilings.begin(),min_tilings.end());
  phen.tiling=min_tilings.front();
}

Phenotype GetPhenotypeFromGrid(std::vector<int8_t>& placed_tiles) {
  std::vector<int8_t> x_locs, y_locs,tile_vals;
  x_locs.reserve(placed_tiles.size()/3);y_locs.reserve(placed_tiles.size()/3);tile_vals.reserve(placed_tiles.size()/3);
  
  for(std::vector<int8_t>::iterator check_iter = placed_tiles.begin();check_iter!=placed_tiles.end();check_iter+=3) {
    x_locs.emplace_back(*check_iter);
    y_locs.emplace_back(*(check_iter+1));
    tile_vals.emplace_back(*(check_iter+2));
  }
  std::vector<int8_t>::iterator x_left,x_right,y_top,y_bottom;
  std::tie(x_left,x_right)=std::minmax_element(x_locs.begin(),x_locs.end());
  std::tie(y_bottom,y_top)=std::minmax_element(y_locs.begin(),y_locs.end());
  uint8_t dx=*x_right-*x_left+1,dy=*y_top-*y_bottom+1;
  std::vector<uint8_t> spatial_grid(dx*dy);
  
  for(uint16_t tileIndex=0;tileIndex<x_locs.size();++tileIndex) {
    uint8_t tile_detail=0;
    switch(DETERMINISM_LEVEL) {
    case 1:
      tile_detail=tile_vals[tileIndex] > 0 ? 1 : 0;
      break;
    case 2:
      tile_detail=tile_vals[tileIndex] > 0 ? (tile_vals[tileIndex]-1)/4+1 : 0;
      break;
    case 3:
      tile_detail=tile_vals[tileIndex];
    }
    spatial_grid[(*y_top-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*x_left)]=tile_detail;
  }
  return Phenotype{dx,dy,spatial_grid};
}


Phenotype_ID PhenotypeTable::GetPhenotypeID(Phenotype& phen) {
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

void PhenotypeTable::RelabelPhenotypes(std::vector<Phenotype_ID >& pids) {
  /*! relabels stored in tuple (size, swap_from, swap_to) */
  for(auto x_iter=new_phenotype_xfer.begin(); x_iter!=new_phenotype_xfer.end();x_iter+=3)
    std::replace(pids.begin(),pids.end(),std::make_pair(static_cast<uint8_t>(*x_iter),static_cast<uint16_t>(*(x_iter+1))),std::make_pair(static_cast<uint8_t>(*x_iter),static_cast<uint16_t>(*(x_iter+2))));
	
  undiscovered_phenotypes.clear();
  undiscovered_phenotype_counts.clear();
  new_phenotype_xfer.clear();
}
std::map<Phenotype_ID,uint16_t> PhenotypeTable::PhenotypeFrequencies(std::vector<Phenotype_ID >& pids, bool& rare_phenotypes) {
  std::map<Phenotype_ID, uint16_t> ID_counter;
  for(std::vector<Phenotype_ID >::const_iterator ID_iter = pids.begin(); ID_iter!=pids.end(); ++ID_iter) {
    if(ID_iter->second < known_phenotypes[ID_iter->first].size())
      ++ID_counter[std::make_pair(ID_iter->first,ID_iter->second)];
    else
      rare_phenotypes=true;
  }
  return ID_counter;
}

void PhenotypeTable::PrintTable(std::ofstream& fout) {
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


