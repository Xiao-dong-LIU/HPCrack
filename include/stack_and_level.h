/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _STACK_AND_LEVEL_H
#define _STACK_AND_LEVEL_H

#include <cstdlib>  // For free()

typedef struct
{
	double hx;             /// mesh size x
	double hy;             /// mesh size y
	double hz;			   /// mesh size z
	long    ii;             /// number of nodes x
	long    jj;             /// number of nodes y
	long    kk;             /// number of nodes z
	int    nbcolox;         /// number of nodes x in each processor
	int    nbcoloy;         /// number of nodes y in each processor
	int    nbcoloz;         /// number of nodes z in each processor
} Level;

typedef struct
{
	int    maxlevel;        /// NB of levels
	int    nx0;             /// number of nodes coarsest grid
	int    ny0;             /// number of nodes coarsest grid
	int    nz0;             /// number of nodes coarsest grid
	double xa,xb;           /// begin,end of computational domain
	double ya,yb;           /// begin,end of computational domain
	double za,zb;           /// begin,end of computational domain
	double wu;              /// work unit counter
	Level  *Ll ;            /// array of grid levels
} Stack;

// Inline function to free memory of Stack
inline void finalize(Stack* U) {
    /// Free memory at the end of the program
    free(U->Ll);
}


#endif
