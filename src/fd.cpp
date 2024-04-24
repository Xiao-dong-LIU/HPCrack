/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "fd.h"
#include "shape_function.h"

//------------------- Compute the new fd dure to the change of H
void fd_h(Stack* U, griddouble& fd, const griddouble& H, const int l)
{
	Level* L;
	double hx, hy, hz;
    L = U->Ll + l;
	hx = L->hx;
	hy = L->hy;
	hz = L->hz;
	grid2d N(8, 8);
    vector<size_t> DNsize = {3,8,8};
    grid <double> DN(DNsize);
	fd.Assign_zero();
	//   get shape function and its gradient
	Dshape_function(N, DN, hx, hy, hz);
	int Cx[8] = { 0,1,1,0,0,1,1,0 };
	int Cy[8] = { 0,0,1,1,0,0,1,1 };
	int Cz[8] = { 0,0,0,0,1,1,1,1 };
	double wg = 0.125 * hx * hy * hz;
	   // compute the right hand d
	for (int k = 0; k <= L->nbcoloz - 2; k++)   /// Nb. of Elements
		for (int j = 0; j <= L->nbcoloy - 2; j++)
			for (int i = 0; i <= L->nbcolox - 2; i++)
			
				for (int tf = 0; tf < 8; tf++)   ///  Test function
                {
                    double fd_ite=0.0;
                    for (int g = 0; g < 8; g++)    /// Gauss Points
					{
					fd_ite += wg *2.0 * H(i, j, k) * N(g, tf);
					}
                    fd(Cx[tf] + i, Cy[tf] + j, Cz[tf] + k)+=fd_ite;
                }
					
}