/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _mgvec_h
#define _mgvec_h


#include "gdvec.h"
#include "omp.h"
#include "structure_df.h"

using namespace std;


template <class T> 
class mgvec{
public:
  // constructor 
  // value_type : "imageF" for reference image, "imageG" for image to be correlated 
  // "node" for nodal type value
  mgvec(Stack *U, const int maxlevel, const int d); 
/*
  /// copy contructor
  mgvec(const mgvec <T> & in):nblevels(in.nblevels), level_max(in.level_max), 
  level_min(in.level_min), Levels(in.Levels) {};
  */
  // get the grid on level l 
  gdvec<T>& getLevel(const int l)  {return *Levels[l];};
  // const get the grid on level l 
  const gdvec<T>& getLevel(const int l) const {return *Levels[l];};
  int nbl() const {return nblevels;};
  // deconstructor 
  ~mgvec();
private:
  int nblevels;
  std::vector<gdvec<T> * >  Levels;
};
using mgvecdouble = mgvec <double>; 
template <class T> 
mgvec<T>::mgvec(Stack *U, const int _maxlevel, const int d):nblevels(_maxlevel+1),
Levels(nblevels)
{    
      Level  *L;
      for (int l = 0; l < nblevels; ++l)
      {
        L=U->Ll+l;
        vector <size_t> dim = {(size_t)L->nbcolox,(size_t)L->nbcoloy,(size_t)L->nbcoloz};
        Levels[l] = new gdvec<T>(d,dim);        
     }
}

template <class T> 
mgvec<T>::~mgvec()
{
    #pragma omp parallel for
    for (int l = 0; l < nblevels; ++l)
      delete Levels[l];
}


#endif
