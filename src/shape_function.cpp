/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "shape_function.h"
#include <cmath>

using namespace std;

///    Shape function  and Derivative of shape function
void Dshape_function(grid2d& N, grid <double>& DN, const double hx, const double hy, const double hz)
{
	double x[8] = { -1., 1., 1.,-1.,-1., 1., 1.,-1. };
	double y[8] = { -1.,-1., 1., 1.,-1.,-1., 1., 1. };
	double z[8] = { -1.,-1.,-1.,-1., 1., 1., 1., 1. };
	double w[8];
	double x1, x2, y1, y2, z1, z2, detJ;
	for (int i = 0; i < 8; i++)
	{
		///   Gauss points
		x[i] = x[i] / sqrt(3);
		y[i] = y[i] / sqrt(3);
		z[i] = z[i] / sqrt(3);
		N(i, 0) = 0.125 * (1 - x[i]) * (1 - y[i]) * (1 - z[i]);
		N(i, 1) = 0.125 * (1 + x[i]) * (1 - y[i]) * (1 - z[i]);
		N(i, 2) = 0.125 * (1 + x[i]) * (1 + y[i]) * (1 - z[i]);
		N(i, 3) = 0.125 * (1 - x[i]) * (1 + y[i]) * (1 - z[i]);
		N(i, 4) = 0.125 * (1 - x[i]) * (1 - y[i]) * (1 + z[i]);
		N(i, 5) = 0.125 * (1 + x[i]) * (1 - y[i]) * (1 + z[i]);
		N(i, 6) = 0.125 * (1 + x[i]) * (1 + y[i]) * (1 + z[i]);
		N(i, 7) = 0.125 * (1 - x[i]) * (1 + y[i]) * (1 + z[i]);
		for (int c = 0; c < 3; c++)
		{
			///Derivation on each direction
			if (c == 0)
			{
				x1 = 2, x2 = 0;
				y1 = y[i], y2 = y[i];
				z1 = z[i], z2 = z[i];
				detJ = 2. / hx;
			}
			if (c == 1)
			{
				x1 = x[i], x2 = x[i];
				y1 = 2, y2 = 0;
				z1 = z[i], z2 = z[i];
				detJ = 2. / hy;
			}
			if (c == 2)
			{
				x1 = x[i], x2 = x[i];
				y1 = y[i], y2 = y[i];
				z1 = 2, z2 = 0;
				detJ = 2. / hz;
			}
			/// 8 Shape Functions
			w[0] = 0.125 * (1 - x1) * (1 - y1) * (1 - z1);
			w[1] = 0.125 * (1 + x2) * (1 - y1) * (1 - z1);
			w[2] = 0.125 * (1 + x2) * (1 + y2) * (1 - z1);
			w[3] = 0.125 * (1 - x1) * (1 + y2) * (1 - z1);
			w[4] = 0.125 * (1 - x1) * (1 - y1) * (1 + z2);
			w[5] = 0.125 * (1 + x2) * (1 - y1) * (1 + z2);
			w[6] = 0.125 * (1 + x2) * (1 + y2) * (1 + z2);
			w[7] = 0.125 * (1 - x1) * (1 + y2) * (1 + z2);
			for (int j = 0; j < 8; j++) {
				DN(c, i, j) = w[j] * detJ;
			}
		}
	}
}