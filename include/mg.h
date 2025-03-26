/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _MG_H
#define _MG_H
#include "grid.h"
#include "omp.h"
#include "structure_df.h"

using namespace std;


// mg
template <class T> 
class mg{
public:
  // constructor 
  // value_type : "imageF" for reference image, "imageG" for image to be correlated 
  // "node" for nodal type value
  mg(Stack *U, const int maxlevel); 
/*
  /// copy contructor
  mg(const mg <T> & in):nblevels(in.nblevels), level_max(in.level_max), 
  level_min(in.level_min), Levels(in.Levels) {};
  */
  // get the grid on level l 
  grid<T>& getLevel(const int l)  {return *Levels[l];};
  const grid<T>& getLevel(const int l) const {return *Levels[l];};
  int nbl() const {return nblevels;};
  // deconstructor 
  ~mg();
private:
  int nblevels;
  std::vector<grid<T> * >  Levels;
};

using mgdouble = mg <double>;
using mgshort = mg <short>;

template <class T> 
//For mg
mg<T>::mg(Stack *U, const int _maxlevel): nblevels(_maxlevel+1),Levels(nblevels)
{    
      Level  *L;
      for (int l = 0; l < nblevels; ++l)
      {
        L=U->Ll+l;
        vector <size_t> dim = {(size_t)L->nbcolox,(size_t)L->nbcoloy,(size_t)L->nbcoloz};
        Levels[l] = new grid<T>(dim);
       }
}

template <class T> 
mg<T>::~mg()
{
    #pragma omp parallel for
    for (int l = 0; l < nblevels; ++l)
      delete Levels[l];
}

#endif
