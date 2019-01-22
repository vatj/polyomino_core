#pragma once
#include <cstdint>
#include <fstream>
#include <iostream>

#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>



/*! free vs one-sided polyominoes and tile vs orientation determinism */
/*! 
  (true) free polyominoes are not chirally distinct
  (false) one-sided polyominoes are

determinism levels as follows:
    shape       : 1
    tile        : 2
    orientation : 3
*/


/*! phenotype definitions */
using Phenotype_ID = std::pair<uint8_t,uint16_t>;
constexpr Phenotype_ID NULL_pid{0,0}, UNBOUND_pid{255,0};

struct Phenotype {
  uint8_t dx,dy;
  std::vector<uint8_t> tiling;

  Phenotype(uint8_t tdx=1, uint8_t tdy=1, std::vector<uint8_t> tt={1}) {dx=tdx;dy=tdy;tiling=tt;}
  
  bool operator==(const Phenotype& rhs) {return this->dx==rhs.dx && this->dy==rhs.dy && this->tiling==rhs.tiling;}
  
  inline static bool FREE_POLYOMINO=true;
  inline static uint8_t DETERMINISM_LEVEL=3;
  inline friend std::ostream& operator<<(std::ostream& out, Phenotype& phen);
  
};

inline std::ostream & operator<<(std::ostream& os,Phenotype& phen) {
  os << "Phenotype dx: " << +phen.dx<<", dy:  " << +phen.dy << "\n";
  for(uint8_t y=0;y<phen.dy;++y) {
    for(uint8_t x=0;x<phen.dx;++x)
      os << +phen.tiling[y*phen.dx+x] << " ";
    os << "\n";
  }
  return os;
}

inline std::ostream & operator<<(std::ostream& os,Phenotype_ID& pid) {
  os <<"pID: "<< +pid.first<<", " << +pid.second;
  return os;
}


inline void ClockwiseRotation(Phenotype& phen) {
  std::vector<uint8_t> swapper;
  swapper.reserve(phen.tiling.size());
  for(uint8_t column=0;column<phen.dx;++column) 
    for(uint8_t row=phen.dy;row!=0;--row)  
      swapper.emplace_back(phen.tiling[(row-1)*phen.dx+column]);                       
  std::swap(phen.dx,phen.dy);
  phen.tiling=swapper;
}

inline void ChiralFlip(Phenotype& phen) {
  for(uint8_t row=0;row<phen.dy;++row)
    std::reverse(phen.tiling.begin()+row*phen.dx,phen.tiling.begin()+(row+1)*phen.dx);
  if(Phenotype::DETERMINISM_LEVEL==3)
    for(uint8_t& element : phen.tiling)
      if(element && element%2==0)
	element+=-(element-1)%4+((element-1)%4+2)%4;
}

inline void MinimizePhenRep(std::vector<uint8_t>& tiling) {
  if(tiling.size()==1 || Phenotype::DETERMINISM_LEVEL==1) {
    for(uint8_t& t : tiling)
      t= t ? 1:0;
    return;
  }
  for(uint8_t& t:tiling)
    t+=128*(t!=0);
  uint8_t swap_count=1;  
  for(std::vector<uint8_t>::iterator t_iter=std::find_if(tiling.begin(),tiling.end(),[](const int s) { return s>0; });t_iter!=tiling.end();) {
    const uint8_t static_swap=*t_iter;
    for(uint8_t cyclic = 0; cyclic< (Phenotype::DETERMINISM_LEVEL==3?4:1); ++cyclic) {
      const uint8_t pre_swap=Phenotype::DETERMINISM_LEVEL==3 ? (static_swap-(static_swap-1)%4)+((static_swap-1)%4+cyclic)%4 : static_swap;
      std::replace(tiling.begin(),tiling.end(),swap_count,uint8_t(255));
      std::replace(tiling.begin(),tiling.end(),pre_swap,swap_count);
      std::replace(tiling.begin(),tiling.end(),uint8_t(255),pre_swap);
      ++swap_count;
    }  
    t_iter=std::find_if(tiling.begin(),tiling.end(),[swap_count](const int s) { return (s>0 && (Phenotype::DETERMINISM_LEVEL==3?s-(s-1)%4:s) >= swap_count); });
  }
}

inline void GetMinPhenRepresentation(Phenotype& phen) {
  std::vector< std::vector<uint8_t> > min_tilings;
  if(phen.dy > phen.dx)
    ClockwiseRotation(phen);
  
  for(uint8_t flip=0;flip<=Phenotype::FREE_POLYOMINO;++flip) {
    MinimizePhenRep(phen.tiling);
    min_tilings.emplace_back(phen.tiling);
      
    if(phen.dy != phen.dx) {
      std::reverse(phen.tiling.begin(),phen.tiling.end());
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
    if(!Phenotype::FREE_POLYOMINO)
      break;
    ChiralFlip(phen);
  }
  phen.tiling=*std::min_element(min_tilings.begin(),min_tilings.end());
}

 
inline Phenotype GetPhenotypeFromGrid(std::vector<int8_t>& placed_tiles) {
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
    switch(Phenotype::DETERMINISM_LEVEL) {
    case 1:
      tile_detail=tile_vals[tileIndex] > 0 ? 1 : 0;
      break;
    case 2:
      tile_detail=tile_vals[tileIndex] > 0 ? (tile_vals[tileIndex]-1)/4+1 : 0;
      break;
    case 3:
	[[fallthrough]]
    default:
      tile_detail=tile_vals[tileIndex];
    }
    spatial_grid[(*y_top-y_locs[tileIndex])*dx + (x_locs[tileIndex]-*x_left)]=tile_detail;
  }
  Phenotype phen{dx,dy,spatial_grid};
  GetMinPhenRepresentation(phen);
  return phen;
}

struct PhenotypeTable {
  bool FIXED_TABLE=false;
  inline static double UND_threshold=0;
  inline static uint16_t phenotype_builds=10;

  std::unordered_map<uint8_t,std::vector<Phenotype> > known_phenotypes;
  
  inline Phenotype_ID GetPhenotypeID(Phenotype& phen) {
    uint8_t phenotype_size=std::count_if(phen.tiling.begin(),phen.tiling.end(),[](const int c){return c != 0;});
    
    //compare against existing table entries
    for(size_t phenotype_index=0; phenotype_index != known_phenotypes[phenotype_size].size();++phenotype_index)
      if(phen==known_phenotypes[phenotype_size][phenotype_index]) 
        return Phenotype_ID{phenotype_size,phenotype_index};
    if(FIXED_TABLE)
      return NULL_pid;
  
    //compare against temporary table entries
    uint16_t new_phenotype_index=0;
    while(new_phenotype_index<undiscovered_phenotypes[phenotype_size].size()) { 
      if(phen==undiscovered_phenotypes[phenotype_size][new_phenotype_index])
        goto found_phen;
      ++new_phenotype_index;
    }

    undiscovered_phenotypes[phenotype_size].emplace_back(phen);
    undiscovered_phenotype_counts[phenotype_size].emplace_back(0);

  found_phen:
    ++undiscovered_phenotype_counts[phenotype_size][new_phenotype_index];
    return Phenotype_ID{phenotype_size,known_phenotypes[phenotype_size].size()+new_phenotype_index+phenotype_builds};
  }
  
  inline void RelabelPhenotypes(std::vector<Phenotype_ID >& pids) {
    const uint16_t thresh_val=std::ceil(UND_threshold*phenotype_builds);
    for(auto& kv : undiscovered_phenotype_counts) {
      const size_t table_size=known_phenotypes[kv.first].size(); 
      for(size_t nth=0; nth<kv.second.size(); ++nth)
        if(kv.second[nth] >= thresh_val) {
          std::replace(pids.begin(),pids.end(),Phenotype_ID{kv.first,table_size+phenotype_builds+nth},Phenotype_ID{kv.first,known_phenotypes[kv.first].size()});
          known_phenotypes[kv.first].emplace_back(undiscovered_phenotypes[kv.first][nth]);
        }
    }
    undiscovered_phenotypes.clear();
    undiscovered_phenotype_counts.clear();
  }
  
  inline std::map<Phenotype_ID,uint16_t> PhenotypeFrequencies(std::vector<Phenotype_ID >& pids,const Phenotype_ID RARE_pid=NULL_pid) {
    const uint16_t thresh_val=std::ceil(UND_threshold*phenotype_builds);
    std::map<Phenotype_ID, uint16_t> ID_counter;
    for(auto& pid : pids)
      ++ID_counter[pid];

    const size_t all_size= ID_counter.size();
    for(auto map_it = ID_counter.begin(); map_it != ID_counter.end();) {
      if(map_it->second < thresh_val) 
        map_it=ID_counter.erase(map_it);
      else
        ++map_it;
    }
    if(ID_counter.size()!=all_size)
      ++ID_counter[RARE_pid];

    return ID_counter;    
  }

  inline void LoadTable(std::string f_name) {
    std::ifstream fin(f_name);
    std::string temp;
    while (std::getline(fin, temp)) {
      std::istringstream buffer(temp);
      std::vector<uint8_t> vec{std::istream_iterator<double>(buffer), std::istream_iterator<double>()};
      
      uint8_t phenotype_size=std::count_if(vec.begin()+2,vec.end(),[](const int c){return c != 0;});
      known_phenotypes[phenotype_size].emplace_back(vec[0],vec[1],std::vector<uint8_t>{vec.begin()+2,vec.end()});
    }
  }
    
  inline void PrintTable(std::string f_name) {
    std::ofstream fout(f_name);
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
  
protected:
  std::unordered_map<uint8_t,std::vector<Phenotype>> undiscovered_phenotypes;
  std::unordered_map<uint8_t, std::vector<uint16_t>> undiscovered_phenotype_counts;

};
