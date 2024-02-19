/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/

#ifndef _vtk_serial_H
#define _vtk_serial_H

#include <string>
#include <vector>
#include "vtk_file.h"
#include "vtkBase64Utilities.h"
#include <typeinfo>
using namespace std;


// extension : 0 for vti, 1 for pvti
void head(FILE * fp, const vector<int> & pindex, const vector<double> & axis, 
const vector<double> & h, const int extension);
// if there is no pointdata
void point_data_null (FILE * fp);
// write the end of pointfata
void point_data_end (FILE * fp);
// write pointdata 
void point_data(FILE * fp, const vector<float> & dat, const string & field_name);
// wrtite binary data
template <class T> 
void write_data (FILE * fp,  const vector<T> & dat)
{
	  long t; 
    T dtype; 
    string c= typeid(dtype).name();
   // cout<<c<<endl;
    int inttype = c.compare("i");
    int floattype = c.compare("f");
    int shorttype = c.compare("s");
    if (floattype==0)
    {
      t = dat.size()*sizeof(float);
    }
    if (shorttype==0)
    {
      t = dat.size()*sizeof(short);
    }
    if (inttype==0)
    {
      t = dat.size()*sizeof(int);
    }
    long taille,taille_entete ;
    const unsigned char *ubuffer;
    unsigned char *buf_code,*buf_entete;
    ubuffer = (unsigned char *)(&dat[0]);
    buf_code=new unsigned char[t*2];
    taille = vtkBase64Utilities::Encode(ubuffer,t,buf_code,'\0');

    ubuffer = reinterpret_cast<unsigned char*>(&t);
    buf_entete=new unsigned char[8];
    taille_entete=vtkBase64Utilities::Encode(ubuffer,sizeof(int),buf_entete,'\0');

    fwrite((char*)buf_entete,taille_entete,1,fp);
    fwrite((char*)buf_code,taille,1,fp);

    delete [] buf_code ;
    delete [] buf_entete ;

	fprintf(fp,"\n</DataArray>\n");
}
// write celledata
template <class T> 
void celledata (FILE * fp, const vector<T> & dat, const string & field_name)
{
	T dtype; 
	string c= typeid(dtype).name();
	int inttype=c.compare("i");
	int floattype=c.compare("f");
  int shorttype = c.compare("s");
	fprintf(fp,"<CellData Scalar = \"%s\">\n",field_name.c_str());
	fprintf(fp,"<DataArray Name = \"%s\" type = ",field_name.c_str());
	if (floattype==0)
	{
		fprintf(fp,"\"Float32\" ");
	}
	if (inttype==0)
	{
		fprintf(fp,"\"Int32\" ");
	}
  if (shorttype==0)
	{
		fprintf(fp,"\"Int16\" ");
	}
	fprintf(fp,"NumberOfComponents = \"1\" format = \"binary\">\n");
  write_data(fp,dat);
}
// write the the end line of celldata
void celledata_end(FILE * fp);
// write the end of file
void write_end(FILE * fp);
// function to write the voxel type vti
template <class T> 
void vkt_total_voxel(const string & filename, const string & field_name, const vector<T> & dat,
const vector<int> & pindex, const vector<double> & axis, const vector<double> & h)
{
    FILE * fp;
    fp = open_file(filename);
	  head(fp,pindex,axis,h,0);
    point_data_null(fp);
    point_data_end(fp);
    celledata(fp,dat,field_name);
    celledata_end(fp);
    write_end(fp);
    close_file(fp);
}
// function to write the elementary type vti
template <class T> 
void vkt_total_displacement(const string & filename, const string & field_name, 
const vector<T> & datmask, const vector<float> & datu, const vector<int> & pindex, 
const vector<double> & axis, const vector<double> & h)
{
    FILE * fp;
    fp = open_file(filename);
	  head(fp,pindex,axis,h,0);
    point_data(fp,datu,field_name);
    point_data_end(fp);
    celledata(fp,datmask,"Mask");
    celledata_end(fp);
    write_end(fp);
    close_file(fp);
}

#endif // _vtk_serial_H
