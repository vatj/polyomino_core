#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include <tuple>
#include <map>
#include <set>

#include <array>
#include <random>

extern thread_local std::mt19937 RNG_Engine;

using Genotype = std::vector<int8_t>;
using interaction_pair = std::pair<uint8_t,uint8_t>;

struct PotentialTileSite {
  interaction_pair bonding_pair;
  std::array<int8_t,3> site_information;
  PotentialTileSite(interaction_pair bp,  int8_t x,int8_t y, int8_t f)
  {bonding_pair=bp;site_information={x,y,f};}
};
struct PotentialTileSites {
  std::vector<PotentialTileSite> sites;
  std::vector<double> strengths;
};

template<class Q>
class PolyominoAssembly {
public:
  
  template<typename T, typename A>
  static uint8_t CountActiveInterfaces(const std::vector<T,A>& genotype) {
    uint8_t N_interfaces=0;
    for(uint8_t b1=0;b1<genotype.size() ;++b1)
      for(uint8_t b2=b1;b2<genotype.size();++b2)
        if(Q::InteractionMatrix(genotype[b1],genotype[b2])) 
          ++N_interfaces;
    return N_interfaces;
  }

  template<typename T, typename A>
  static void StripNoncodingGenotype(std::vector<T,A>& genotype) {
    std::vector<uint8_t> coding{0},noncoding(genotype.size()/4-1);
    std::iota(noncoding.begin(), noncoding.end(), 1);
  
    for(uint8_t c_in=0;c_in<coding.size();++c_in)
      for(uint8_t nc_in=0;nc_in<noncoding.size();++nc_in) {
        for(uint8_t cface=0;cface<4;++cface)
          for(uint8_t ncface=0;ncface<4;++ncface)
            if(Q::InteractionMatrix(genotype[coding[c_in]*4+cface],genotype[noncoding[nc_in]*4+ncface])){
              coding.emplace_back(noncoding[nc_in]);
              noncoding.erase(noncoding.begin()+nc_in--);
              goto newtile;
            }
      newtile: ;
      }

    for(uint8_t rm=0;rm<noncoding.size();++rm)
      genotype.erase(genotype.begin()+(noncoding[rm]-rm)*4,genotype.begin()+(1+noncoding[rm]-rm)*4);
  }

  template <typename T, typename A>
  static std::vector<std::pair<interaction_pair,double> > GetEdgePairs(const std::vector<T,A>& genotype) {
    std::vector<std::pair<interaction_pair,double> > edge_pairs;
    for(uint8_t b1=0;b1<genotype.size();++b1)
      for(uint8_t b2=b1;b2<genotype.size();++b2) 
        if(double B_ij=Q::InteractionMatrix(genotype[b1],genotype[b2]))
          edge_pairs.emplace_back(std::minmax(b1,b2),B_ij);
    
    return edge_pairs;
  }

	static inline std::vector<int8_t> AssemblePolyomino(const std::vector<std::pair<interaction_pair,double> > edges,const int8_t seed,const size_t UNBOUND_LIMIT, std::set<interaction_pair>& interacting_indices) {
  //uint16_t accumulated_time=0;
     
  std::vector<int8_t> placed_tiles{0,0,seed},growing_perimeter;
  std::vector<double> strengths_cdf;
  PotentialTileSites perimeter_sites;
    
  ExtendPerimeter(edges,seed,0,0,placed_tiles,perimeter_sites);

  while(!perimeter_sites.strengths.empty()) {
    /*! select new site proportional to binding strength */
    strengths_cdf.resize(perimeter_sites.strengths.size());
    std::partial_sum(perimeter_sites.strengths.begin(), perimeter_sites.strengths.end(), strengths_cdf.begin());
    std::uniform_real_distribution<double> random_interval(0,strengths_cdf.back());
    size_t selected_choice=static_cast<size_t>(std::lower_bound(strengths_cdf.begin(),strengths_cdf.end(),random_interval(RNG_Engine))-strengths_cdf.begin());
      
    auto chosen_site=perimeter_sites.sites[selected_choice];
    /*
    if(true) {
      accumulated_time+=1+std::geometric_distribution<uint16_t>(perimeter_sites.strengths[selected_choice])(RNG_Engine);
      if(accumulated_time>2)
	break;
    }
    */
    /*! place new tile */
    placed_tiles.insert(placed_tiles.end(),chosen_site.site_information.begin(),chosen_site.site_information.end());
    interacting_indices.insert(chosen_site.bonding_pair);
    if(placed_tiles.size()>UNBOUND_LIMIT)
      return {};
      
    /*! remove all further options in same tile location */
    auto [f_x, f_y, f_t] = chosen_site.site_information;
    for(size_t cut_index=0;cut_index<perimeter_sites.strengths.size();) {
      if(f_x==perimeter_sites.sites[cut_index].site_information[0] && f_y==perimeter_sites.sites[cut_index].site_information[1]) {
	perimeter_sites.strengths.erase(perimeter_sites.strengths.begin()+cut_index);
	perimeter_sites.sites.erase(perimeter_sites.sites.begin()+cut_index);
      }
      else
	++cut_index;
    }
    ExtendPerimeter(edges,f_t,f_x,f_y,placed_tiles,perimeter_sites);
  }
  return placed_tiles;
}
 
static inline void ExtendPerimeter(const std::vector<std::pair<interaction_pair,double> >& edges,uint8_t tile_detail, int8_t x,int8_t y, std::vector<int8_t>& placed_tiles,PotentialTileSites& perimeter_sites) {
  int8_t dx=0,dy=0,tile=(tile_detail-1)/4,theta=(tile_detail-1)%4;
  for(uint8_t f=0;f<4;++f) {
    switch(f) {
    case 0:dx=0;dy=1;break;
    case 1:dx=1;dy=0;break;
    case 2:dx=0;dy=-1;break;
    case 3:dx=-1;dy=0;break;
    }
    /*! site already occupied, move on */
    for(std::vector<int8_t>::reverse_iterator tile_info=placed_tiles.rbegin();tile_info!=placed_tiles.rend();tile_info+=3) 
      if((x+dx)==*(tile_info+2) && (y+dy)==*(tile_info+1)) 
	goto nextloop;

    {//empty scope to prevent cross-initilisation errors from goto
      uint8_t g_index=static_cast<uint8_t>(tile*4+(f-theta+4)%4);
      std::vector<std::pair<interaction_pair,double> >::const_iterator iter = edges.begin();
      while ((iter = std::find_if(iter, edges.end(),[&g_index](const auto& edge){ return (edge.first.first == g_index || edge.first.second == g_index);})) != edges.end()) {
	int8_t base = iter->first.first==g_index ? iter->first.second : iter->first.first;
	perimeter_sites.sites.emplace_back(iter->first,static_cast<int8_t>(x+dx),static_cast<int8_t>(y+dy),static_cast<int8_t>(base-base%4+((f+2)%4-base%4+4)%4+1));
	perimeter_sites.strengths.emplace_back(iter->second);
	++iter;
      }
    }
    /*! find all above threshold bindings and add to new sites */
  nextloop: ; //continue if site already occupied
  }   
}
  
};


/*
void Clean_Genome(Genotype& genome,int secondNonInteracting,bool Remove_Duplicates);
void Minimize_Tile_Set(Genotype& genome);
std::map<uint8_t,uint8_t> DuplicateGenes(Genotype& genome);
bool Disjointed_Check(Genotype& genome);
void Search_Next_Tile(Genotype& genome, std::vector<uint8_t>& Unvisited, std::vector<uint8_t>& Connected_Components,uint8_t tile);
*/



