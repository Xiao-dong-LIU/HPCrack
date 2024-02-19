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
void initialize(Stack *U, const inputdomain *IDM, const int maxlevel, vector <int> const& element_nb, 
vector <double> const& X_start, const MPI_Setting &M)
{
  double xa = X_start[0];
  double ya = X_start[1];
  double za = X_start[2];
  double xb = xa + IDM->coefxb;
  double yb = yb + IDM->coefyb;
  double zb = zb + IDM->coefzb;
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

  hx = (xb - xa) / element_nb[0];
  hy = (yb - ya) / element_nb[1];
  hz = (zb - za) / element_nb[2];
  ii = element_nb[0];
  jj = element_nb[1];
  kk = element_nb[2];

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
  //  cout << "l="<<l<<" x="<<L->ii<<" y="<<L->jj<<" z="<<L->kk<<endl;
  }
}
/*
///------------initialisation of u,f,d,H------------///

void init_ufdH(Stack *U, multigrid3d & u, multigrid3d & v, multigrid3d & w,
               multigrid3d & fu, multigrid3d & fv, multigrid3d & fw,
               multigrid3d & d, multigrid3d & fd, multigrid3d & H,
               int l,MPI_Setting & M)

/// initializes u,
/// initializes right hand side values on grid level l
{
Level *L;
double x,y,z;
L =U->Ll+l;
#pragma omp parallel for private(x,y,z)
for (int k=0;k<=L->nbcoloz-1;k++)
  {
      /// z in each processor
    if (M.coordinates[2]==0) z=U->za+k*L->hz;
    else if (M.coordinates[2]==(M.dims[2]-1))
    z = U->za+(L->kk-L->nbcoloz+1+k)*L->hz;
    else
    z = U->za+(M.coordinates[2]*(L->nbcoloz-2) + k-1)*L->hz;
    for (int j=0;j<=L->nbcoloy-1;j++)
      {
            /// z in each processor
        if (M.coordinates[1]==0) y=U->ya+j*L->hy;
        else if (M.coordinates[1]==(M.dims[1]-1))
        y = U->ya+(L->jj-L->nbcoloy+1+j)*L->hy;
        else
        y = U->ya+(M.coordinates[1]*(L->nbcoloy-2) + j-1)*L->hy;
        for (int i=0;i<=L->nbcolox-1;i++)
        {
              /// z in each processor
          if (M.coordinates[0]==0) x=U->xa+i*L->hx;
          else if (M.coordinates[0]==(M.dims[0]-1))
          x = U->xa+(L->ii-L->nbcolox+1+i)*L->hx;
          else
          x = U->xa+(M.coordinates[0]*(L->nbcolox-2) + i-1)*L->hx;
          /// Initialize

          u(l-1,i,j,k)=0.0*(U->xb-U->xa)*(x+U->xb);
          v(l-1,i,j,k)=0.0*(U->yb-U->ya)*(y+U->yb);
          w(l-1,i,j,k)=0.0*(U->zb-U->za)*(z+U->zb);
          fu(l-1,i,j,k)=0.0;
          fv(l-1,i,j,k)=0.0;
          fw(l-1,i,j,k)=0.0;
          d(l-1,i,j,k)=0.0;
          fd(l-1,i,j,k)=0.0;
          H(l-1,i,j,k)=0.0;
        /*
          double r=sqrt((x-0.5*(U->xa+U->xb))*(x-0.5*(U->xa+U->xb))+
              (y-0.5*(U->ya+U->yb))*(y-0.5*(U->ya+U->yb))+
              (z-0.5*(U->za+U->zb))*(z-0.5*(U->za+U->zb)));
          if (r<=0.25*(U->xb-U->xa))
          {
            d(l - 1, i, j, k) = 1.0;
          }

        }
      }
  }
}
*/