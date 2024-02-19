/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include "mg_cycle_parameter.h"

void mg_parameter(MG & mgp, const int nu0, const int nu1, const int nu2, 
const int gamma, const int ncy)
{
    mgp.nu0=nu0;
    mgp.nu1=nu1;
    mgp.nu2=nu2;
    mgp.gamma=gamma;
    mgp.ncy=ncy;
}