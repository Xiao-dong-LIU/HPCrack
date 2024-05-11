/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#include "shape_function.h"
#include "epsilon.h"

//---------------- Compute epsilon at one gauss point
void epsilon_routine(grid2d& epsilon, const gdvecdouble & u, const griddouble & DN,  
	 const int i, const int j, const int k, const int g) 
{
	double dc = 0.0;   // displacement in direction c
	double dm = 0.0;   // displacement in direction m
	int Cx[8] = {0,1,1,0,0,1,1,0};
	int Cy[8] = {0,0,1,1,0,0,1,1};
	int Cz[8] = {0,0,0,0,1,1,1,1};
	for (int c=0; c<3; c++)  //row of epsilon
		for (int m = 0; m < 3; m++)  //column of epsilon
		{
			double epsilon_temp=0.0;
			for (int tf = 0; tf < 8; tf++) // sum on Shape function 
			{
				dc = u.getgrid(c)(Cx[tf] + i, Cy[tf] + j, Cz[tf] + k);
				dm = u.getgrid(m)(Cx[tf] + i, Cy[tf] + j, Cz[tf] + k);
				epsilon_temp += 0.5*(DN(c, g, tf) * dm + DN(m, g, tf) * dc);
			}
			epsilon(c, m)=epsilon_temp;
		}
}


//---------------- Compute the trace of epsilon at one gauss point
double trace_epsilon(const grid2d& epsilon)
{
	double tr_epsilon = 0.0;
	tr_epsilon = epsilon(0, 0) + epsilon(1, 1) + epsilon(2, 2);
	return tr_epsilon;
}

//---------------- Double dot product of the epsilon_dev at gauss integration point
double epsilon_dev(const grid2d& epsilon)
{
	double ddp_epsilon=0.0;
	double third_trace = 0.0;  // One third of trace of epsilon
	third_trace = (1.0/3.0)*trace_epsilon(epsilon);
	for (int c = 0; c < 3; c++)
		for (int m = 0; m < 3; m++)
		{
			if (c==m)
			ddp_epsilon += (epsilon(c, m) - third_trace)* (epsilon(c, m) - third_trace);
			else 
			ddp_epsilon += epsilon(c, m) * epsilon(c, m) ;
		}
	return ddp_epsilon;
}