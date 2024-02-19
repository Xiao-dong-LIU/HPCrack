/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _fd_h
#define _fd_h
#include "grid.h"
#include "structure_df.h"

//------------------- Compute the new fd dure to the change of H
void fd_h(Stack* U, griddouble & fd, const griddouble & H, const int l);


#endif // _fd_h