#include "core_assembly.hpp"

std::vector<int8_t> AssemblePolyomino(const std::vector<std::pair<interaction_pair,double> > edges,const int8_t seed,const size_t UNBOUND_LIMIT, std::set<interaction_pair>& interacting_indices) {
  uint16_t accumulated_time=0;
     
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
 
void ExtendPerimeter(const std::vector<std::pair<interaction_pair,double> >& edges,uint8_t tile_detail, int8_t x,int8_t y, std::vector<int8_t>& placed_tiles,PotentialTileSites& perimeter_sites) {
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
