#include "core_phenotype.hpp"

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
