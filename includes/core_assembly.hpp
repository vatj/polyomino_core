#pragma once
#include "core_phenotype.hpp"
#include "core_genotype.hpp"

#include <array>
#include <random>

extern thread_local std::mt19937 RNG_Engine;

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

std::vector<int8_t> AssemblePolyomino(const std::vector<std::pair<interaction_pair,double> > edges,const int8_t seed,const size_t UNBOUND_LIMIT, std::set<interaction_pair>& interacting_indices);
void ExtendPerimeter(const std::vector<std::pair<interaction_pair,double> >& edges,uint8_t tile_detail, int8_t x,int8_t y, std::vector<int8_t>& placed_tiles,PotentialTileSites& perimeter_sites);
