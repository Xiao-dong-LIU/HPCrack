/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include "initialisation.h"
#include "omp.h"
#include <cmath>
#include "iostream"

using namespace std;

/// ------------- initialize values in datastructure
void initialize(Stack *U, const Parameters & para, const MPI_Setting &M, const int maxlevel)
{
  double xa = para.X_start[0];
  double ya = para.X_start[1];
  double za = para.X_start[2];
  double xb = xa + para.IDM.coefxb;
  double yb = ya + para.IDM.coefyb;
  double zb = za + para.IDM.coefzb;
  double hx, hy, hz;
  Level *L;
  int ii, jj, kk;

  U->maxlevel = maxlevel;
  U->xa = xa;
  U->xb = xb;
  U->ya = ya;
  U->yb = yb;
  U->za = za;
  U->zb = zb;
  U->wu = 0.0;

  U->Ll = (Level *)calloc(maxlevel+1, sizeof(Level));

  hx = (xb - xa) / para.element_nb[0];
  hy = (yb - ya) / para.element_nb[1];
  hz = (zb - za) / para.element_nb[2];
  ii = para.element_nb[0];
  jj = para.element_nb[1];
  kk = para.element_nb[2];

  int npx = M.dims[0], npy = M.dims[1], npz = M.dims[2];
  int nbcolox = (ii + 1) / npx + 2;
  int nbcoloy = (jj + 1) / npy + 2;
  int nbcoloz = (kk + 1) / npz + 2;
  if (npx == 1)
    nbcolox = ii + 1;
  else
  {
    if (M.coordinates[0] == 0)
      nbcolox = (ii + 1) / npx + 1;
    if (M.coordinates[0] == npx - 1)
      nbcolox = (ii + 1) / npx + (ii + 1) % npx + 1;
  }
  if (npy == 1)
    nbcoloy = jj + 1;
  else
  {
    if (M.coordinates[1] == 0)
      nbcoloy = (jj + 1) / npy + 1;
    if (M.coordinates[1] == npy - 1)
      nbcoloy = (jj + 1) / npy + (jj + 1) % npy + 1;
  }
  if (npz == 1)
    nbcoloz = kk + 1;
  else
  {
    if (M.coordinates[2] == 0)
      nbcoloz = (kk + 1) / npz + 1;
    if (M.coordinates[2] == npz - 1)
      nbcoloz = (kk + 1) / npz + (kk + 1) % npz + 1;
  }
  for (int l = 0; l <= maxlevel; l++)
  {
    L = U->Ll + l;
    L->hx = hx;
    L->hy = hy;
    L->hz = hz;
    L->ii = ii;
    L->jj = jj;
    L->kk = kk;
    if (M.coordinates[0] == 0 && M.coordinates[1] == 0 && M.coordinates[2] == 0 && l == maxlevel)
      cout << "ii=" << ii << "\t jj=" << jj << "\t kk=" << kk << endl;
    L->nbcolox = nbcolox;
    L->nbcoloy = nbcoloy;
    L->nbcoloz = nbcoloz;
    hx *= 0.5;
    hy *= 0.5;
    hz *= 0.5;
    ii *= 2;
    jj *= 2;
    kk *= 2;
    if (npx == 1)
      nbcolox = 2 * (nbcolox - 1) + 1;
    else
    {
      if (M.coordinates[0] == 0)
        nbcolox = 2 * (nbcolox - 1) + 1;
      else if (M.coordinates[0] == npx - 1)
        nbcolox = 2 * (nbcolox - 1);
      else
        nbcolox = 2 * (nbcolox - 2) + 2;
    }
    if (npy == 1)
      nbcoloy = 2 * (nbcoloy - 1) + 1;
    else
    {
      if (M.coordinates[1] == 0)
        nbcoloy = 2 * (nbcoloy - 1) + 1;
      else if (M.coordinates[1] == npy - 1)
        nbcoloy = 2 * (nbcoloy - 1);
      else
        nbcoloy = 2 * (nbcoloy - 2) + 2;
    }
    if (npz == 1)
      nbcoloz = 2 * (nbcoloz - 1) + 1;
    else
    {
      if (M.coordinates[2] == 0)
        nbcoloz = 2 * (nbcoloz - 1) + 1;
      else if (M.coordinates[2] == npz - 1)
        nbcoloz = 2 * (nbcoloz - 1);
      else
        nbcoloz = 2 * (nbcoloz - 2) + 2;
    }
  }
}