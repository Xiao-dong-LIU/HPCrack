/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _ROI_SETTING_H
#define _ROI_SETTING_H
//imput for image
typedef struct
{
 int coefxa;            ///begin of domain in X, start by 0
 int coefxb;			///end of domain in X, start by 0
 int coefya;			///begin of domain in Y, start by 0
 int coefyb;			///end of domain in Y, start by 0
 int coefza;			///begin of domain in Z, start by 0
 int coefzb;			///end of domain in Z, start by 0
 int voxelnx;            ///NB of voxels in image on X
 int voxelny;				///NB of voxels in image on Y
 int voxelnz;				///NB of voxels in image on Z
} inputdomain;

#endif 