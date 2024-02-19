/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include <stdio.h>
#include <string.h>
#include <cmath>
#include "vtk_serial.h"

using namespace std;



///------------VTK XML OUTPUT------------///
void head(FILE * fp, const vector<int> & pindex, const vector<double> & axis, 
const vector<double> & h, const int extension)
{
	
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=");

	//------- for vti and pvti
	if (extension == 0)
	fprintf(fp, "\"ImageData\" ");
	else if (extension == 1)
	fprintf(fp, "\"PImageData\" ");

	fprintf(fp, "version=\"0.1\" byte_order=\"LittleEndian\">\n");
	//------------ Extent
	if (extension == 0)
	fprintf(fp, "<ImageData ");
	else if (extension == 1)
	fprintf(fp, "<PImageData ");

	fprintf(fp, "WholeExtent = \"");
	for (int i=0; i<pindex.size(); i++)
	{
		fprintf(fp,"%d ", pindex[i]);
	}
	fprintf(fp, "\" ");
	//------------ Axis Origin
	if (extension == 1)
	fprintf(fp, "GhostLevel=\"0\" ");
	fprintf(fp, "Origin = \"");
	for (int i=0; i<axis.size(); i++)
	{
		fprintf(fp,"%f ", axis[i]);
	}
	fprintf(fp, "\" ");
	//------------ Spacing
	fprintf(fp, "Spacing = \"");
	for (int i=0; i<h.size(); i++)
	{
		fprintf(fp,"%f ", h[i]);
	}
	fprintf(fp, "\">\n");
	//------------ Piece Extent 
	if (extension == 0)
	{
		fprintf(fp, "<Piece Extent = \"");
		for (int i=0; i<pindex.size(); i++)
		{
			fprintf(fp,"%d ", pindex[i]);
		}
		fprintf(fp, "\">\n");
	}
	
	
}


void point_data_null (FILE * fp)
{
	
	fprintf(fp,"<PointData> ");
}

void point_data(FILE * fp, const vector<float> & dat, const string & field_name)
{
	fprintf(fp,"<PointData Vectors=\"%s\">\n",field_name.c_str());
	fprintf(fp,"<DataArray Name=\"%s\" type=\"Float32\" ", field_name.c_str());
	fprintf(fp,"NumberOfComponents=\"3\" format=\"binary\">\n");
	write_data(fp,dat);
}

void point_data_end(FILE *fp)
{
	fprintf(fp,"</PointData>\n");
}

void celledata_end(FILE * fp)
{
	fprintf(fp,"</CellData>\n");
}

void write_end(FILE * fp)
{
	fprintf(fp,"</Piece>\n");
	fprintf(fp,"</ImageData>\n");
	fprintf(fp,"</VTKFile>");
}


