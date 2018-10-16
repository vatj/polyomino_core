#include <cstdint>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include <tuple>
#include <map>

using Genotype = std::vector<uint8_t>;
using interaction_pair = std::pair<uint8_t,uint8_t>;

bool InteractionMatrix(const uint8_t face_1,const uint8_t face_2);
double BindingStrength(const uint8_t face_1,const uint8_t face_2);

template <typename T>
uint8_t CountActiveInterfaces(const T& genotype) {
  uint8_t N_interfaces=0;
  for(uint8_t b1=0;b1<genotype.size() ;++b1)
    for(uint8_t b2=b1;b2<genotype.size();++b2)
      if(InteractionMatrix(genotype[b1],genotype[b2]))
        ++N_interfaces;
  return N_interfaces;
}
template <typename T>
void StripNoncodingGenotype(T& genotype) {
  std::vector<uint8_t> coding{0},noncoding(genotype.size()/4-1);
  std::iota(noncoding.begin(), noncoding.end(), 1);
  
  for(uint8_t c_in=0;c_in<coding.size();++c_in)
    for(uint8_t nc_in=0;nc_in<noncoding.size();++nc_in) {
      for(uint8_t cface=0;cface<4;++cface)
        for(uint8_t ncface=0;ncface<4;++ncface)
          if(InteractionMatrix(genotype[coding[c_in]*4+cface],genotype[noncoding[nc_in]*4+ncface])){
            coding.emplace_back(noncoding[nc_in]);
            noncoding.erase(noncoding.begin()+nc_in--);
            goto newtile;
          }
    newtile: ;
    }

  for(uint8_t rm=0;rm<noncoding.size();++rm)
    genotype.erase(genotype.begin()+(noncoding[rm]-rm)*4,genotype.begin()+(1+noncoding[rm]-rm)*4);
}
template <typename T>
std::vector<std::pair<interaction_pair,double> > GetEdgePairs(const T& genotype) {
  std::vector<std::pair<interaction_pair,double> > edge_pairs;
  for(uint8_t b1=0;b1<genotype.size();++b1)
    for(uint8_t b2=b1;b2<genotype.size();++b2)
      if(InteractionMatrix(genotype[b1],genotype[b2]))
        edge_pairs.emplace_back(std::minmax(b1,b2),BindingStrength(genotype[b1],genotype[b2]));

  return edge_pairs;
}

uint8_t Interaction_Matrix(uint8_t input_face);
void Clean_Genome(Genotype& genome,int secondNonInteracting,bool Remove_Duplicates);
void Minimize_Tile_Set(Genotype& genome);
std::map<uint8_t,uint8_t> DuplicateGenes(Genotype& genome);
bool Disjointed_Check(Genotype& genome);
void Search_Next_Tile(Genotype& genome, std::vector<uint8_t>& Unvisited, std::vector<uint8_t>& Connected_Components,uint8_t tile);
