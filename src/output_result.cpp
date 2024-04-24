/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include <stdio.h>
#include <string.h>
#include <cmath>
#include "vtkBase64Utilities.h"
#include "output_result.h"
#include "shape_function.h"
#include "imput_image.h"
#include "boundary_conditions.h"


using namespace std;

static FILE *fp = NULL;

///------------VTK XML OUTPUT------------///

void open_file(const char *filename, int l, const int myid, int t)
{
    char full_filename[1024];
    if (strstr(filename, ".vti") != NULL){
      strcpy(full_filename, filename);}
    else{
      sprintf(full_filename, "Output/%s_l_%d_t_%d_Proc_%d.vti", filename,l,t,myid);}
    fp = fopen(full_filename, "wb+");
}

void close_file(void)
{
    fclose(fp);
    fp = NULL;
}

void write_header(Stack *U, const MPI_Setting & M, const int l)
{
    Level *L;
	L =U->Ll+l;
	double x0, y0, z0;
	/// Get the number of points in each direction at each processor
	int nxs = (L->nbcolox-2)*M.coordinates[0];
	int nx = (L->nbcolox-2)*(M.coordinates[0]+1);
	if (M.dims[0]==1) {
		nxs=0;
		nx = L->ii;}
	else{
		if(M.coordinates[0]==0) {
			nxs=0;
			nx = L->nbcolox-1;}
		if (M.coordinates[0]==M.dims[0]-1) {
			nxs = L->ii-L->nbcolox+2;
			nx = L->ii;}}
	int nys = (L->nbcoloy-2)*M.coordinates[1];
	int ny = (L->nbcoloy-2)*(M.coordinates[1]+1);
	if (M.dims[1]==1) {
		nys =0;
		ny = L->jj;}
	else{
		if (M.coordinates[1]==0) {
			nys = 0;
			ny = L->nbcoloy-1;}
		if (M.coordinates[1]==M.dims[1]-1) {
			nys = L->jj-L->nbcoloy+2;
			ny = L->jj;}}
	int nzs = (L->nbcoloz-2)*M.coordinates[2];
	int nz = (L->nbcoloz-2)*(M.coordinates[2]+1);
	if (M.dims[2]==1) {
		nzs = 0;
		nz = L->kk;}
	else{
		if (M.coordinates[2]==0) {
			nzs = 0;
			nz = L->nbcoloz-1;}
		if (M.coordinates[2]==M.dims[2]-1) {
			nzs = L->kk-L->nbcoloz+2;
			nz = L->kk;}}
	/// Get the original of X at each processor
	if (M.coordinates[0]==0) x0=U->xa;
	else if (M.coordinates[0]==(M.dims[0]-1))
		x0 = U->xa+(L->ii-L->nbcolox+2)*L->hx;
	else
		x0 = U->xa+(M.coordinates[0]*(L->nbcolox-2))*L->hx;
	/// Get the original of Y at each processor
	if (M.coordinates[1]==0) y0=U->ya;
	else if (M.coordinates[1]==(M.dims[1]-1))
		y0 = U->ya+(L->jj-L->nbcoloy+2)*L->hy;
	else
		y0 = U->ya+(M.coordinates[1]*(L->nbcoloy-2))*L->hy;
	/// Get the original of Z at each processor
	if (M.coordinates[2]==0) z0=U->za;
	else if (M.coordinates[2]==(M.dims[2]-1))
		z0 = U->za+(L->kk-L->nbcoloz+2)*L->hz;
	else
		z0 = U->za+(M.coordinates[2]*(L->nbcoloz-2))*L->hz;
    /// Open file and write

    char str[1024];
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	sprintf(str, "<ImageData WholeExtent=\"%d %d %d %d %d %d \" \
				Origin=\"%f %f %f\" Spacing=\"%f %f %f\">\n",
				nxs,nx,nys,ny,nzs,nz,x0,y0,z0,L->hx,L->hy,L->hz);
	fprintf(fp,"%s",str);
	sprintf(str,"<Piece Extent=\"%d %d %d %d %d %d\" >\n",nxs,nx,nys,ny,nzs,nz);
	fprintf(fp,"%s",str);
}

//------------- point data for different fields
void write_pointdata_head(const string & output_option, const char * field_name)
{
	int V_comp=output_option.compare("V");		// output field of vector
	int S_comp=output_option.compare("S");		// output field of scalar
	int VS_comp=output_option.compare("VS");		// output vector and scalar 
	if(S_comp==0) 
	{
		fprintf(fp,"<PointData Scalars=\"%s\">\n",field_name);
		fprintf(fp,"<DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\" format=\"binary\">\n",field_name);
	}
	if(V_comp==0) 
	{
		fprintf(fp,"<PointData Vectors=\"%s\">\n",field_name);
		fprintf(fp,"<DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n",field_name);
	}
	if (VS_comp==0)
	{
		fprintf(fp,"<PointData Scalars=\"E\" Vectors=\"Diplacement\">\n");
		fprintf(fp,"<DataArray Name=\"E\" type=\"Float32\" NumberOfComponents=\"1\" format=\"binary\">\n");
	}
}

void write_scalar(Stack *U, const griddouble& d, const MPI_Setting & M, const int l)
 
{
	Level *L;
	L =U->Ll+l;
	int is =1,js=1,ks=1;
	if(M.dims[0]==1) is =0;
	else if (M.coordinates[0]==0) is =0;
	if(M.dims[1]==1) js =0;
	else if (M.coordinates[1]==0) js =0;
	if(M.dims[2]==1) ks =0;
	else if (M.coordinates[2]==0) ks =0;
	int ie =L->nbcolox-1,je=L->nbcoloy-1,ke=L->nbcoloz-1;
	vector <size_t> WUsize = {(size_t)(ie-is+1),(size_t)(je-js+1),(size_t)(ke-ks+1)};
	gridfloat WU(WUsize);
	for (int i = is; i<=ie; i++)
	  for (int j = js; j<=je; j++)
	    for (int k = ks; k<=ke; k++)
	    {
			int x = i-1, y = j-1, z = k-1;
			if (M.coordinates[0]==0) x =i;
			if (M.coordinates[1]==0) y =j;
			if (M.coordinates[2]==0) z =k;
			WU(x,y,z)=(float)d(i,j,k);
		}
    long t = (ie-is+1)*(je-js+1)*(ke-ks+1)*sizeof(float);
    long taille,taille_entete ;
    const unsigned char *ubuffer;
    unsigned char *buf_code,*buf_entete;
    ubuffer = (unsigned char *)(&WU(0,0,0));
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

void write_vector(Stack *U, const gdvecdouble & u, const MPI_Setting & M, 
const int l, double u_proportion)
{
	Level *L;
	L =U->Ll+l;

	int is =1,js=1,ks=1;
	if(M.dims[0]==1) is =0;
	else if (M.coordinates[0]==0) is =0;
	if(M.dims[1]==1) js =0;
	else if (M.coordinates[1]==0) js =0;
	if(M.dims[2]==1) ks =0;
	else if (M.coordinates[2]==0) ks =0;
	int ie =L->nbcolox-1,je=L->nbcoloy-1,ke=L->nbcoloz-1;
	vector <size_t> WUsize = {(size_t)(3*(ie-is+1)),(size_t)(je-js+1),(size_t)(ke-ks+1)};
	gridfloat WU(WUsize);
	for (int i = is; i<=ie; i++)
	  for (int j = js; j<=je; j++)
	    for (int k = ks; k<=ke; k++)
	    {
			int x = i-1, y = j-1, z = k-1;
			if (M.coordinates[0]==0) x =i;
			if (M.coordinates[1]==0) y =j;
			if (M.coordinates[2]==0) z =k;
			WU(3*x  ,y,z)=(float)(u_proportion*u.getgrid(0)(i,j,k));
			WU(3*x+1,y,z)=(float)(u_proportion*u.getgrid(1)(i,j,k));
			WU(3*x+2,y,z)=(float)(u_proportion*u.getgrid(2)(i,j,k));
		}
    long t = 3*(ie-is+1)*(je-js+1)*(ke-ks+1)*sizeof(float);
    long taille,taille_entete ;
    const unsigned char *ubuffer;
    unsigned char *buf_code,*buf_entete;
    ubuffer = (unsigned char *)(&WU(0,0,0));
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

void write_pointdata_end(void)
{
	fprintf(fp,"</PointData>\n");
}

void write_celltdata(void)
{

	fprintf(fp,"<CellData></CellData>\n");
}


void write_celltdata_H(Stack *U, const griddouble &H, const MPI_Setting & M, const int l, 
const char *field_name)
{

	fprintf(fp,"<PointData></PointData>\n");
	fprintf(fp,"<CellData Scalar = \"%s\">\n",field_name);
	fprintf(fp,"<DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\" format=\"binary\">\n",field_name);
	Level *L;
	L =U->Ll+l;
	int is =1,js=1,ks=1;
	if (M.coordinates[0]==0) is =0;
	if (M.coordinates[1]==0) js =0;
	if (M.coordinates[2]==0) ks =0;
	int ie =L->nbcolox-2,je=L->nbcoloy-2,ke=L->nbcoloz-2;
	vector <size_t> WUsize = {(size_t)(ie-is+1),(size_t)(je-js+1),(size_t)(ke-ks+1)};
	gridfloat WU(WUsize);
	for (int i = is; i<=ie; i++)
	  for (int j = js; j<=je; j++)
	    for (int k = ks; k<=ke; k++)
	    {
			int x = i-1, y = j-1, z = k-1;
			if (M.coordinates[0]==0) x =i;
			if (M.coordinates[1]==0) y =j;
			if (M.coordinates[2]==0) z =k;
			WU(x,y,z)=(float)H(i,j,k);
			//std::cout<<"mpi="<<M.coordinates[0]<<",Output:"<<"i="<<i<<",j="<<j<<",k="<<k<<",H="<<H(l-1,i,j,k)<<std::endl;
		}
    long t = (ie-is+1)*(je-js+1)*(ke-ks+1)*sizeof(float);
    long taille,taille_entete ;
    const unsigned char *ubuffer;
    unsigned char *buf_code,*buf_entete;
    ubuffer = (unsigned char *)(&WU(0,0,0));
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
	fprintf(fp,"</CellData>\n");
}

void write_end(void)
{
	fprintf(fp,"</Piece>\n");
	fprintf(fp,"</ImageData>\n");
	fprintf(fp,"</VTKFile>");
	close_file();
}

void write_pvtk(Stack *U, const char *filename, const MPI_Setting & M, const string output_option, 
const char *field_name, const int l, const int myid, const int nbprocs, const int t)
{
	Level *L;
	L =U->Ll+l;
	int *nxt,*nxts,*nyts,*nyt,*nzts,*nzt;
	nxt = new int [nbprocs];
	nyt = new int [nbprocs];
	nzt = new int [nbprocs];
	nxts = new int [nbprocs];
	nyts = new int [nbprocs];
	nzts = new int [nbprocs];
	/// Get the number of points in each direction at each processor
	int nxs = (L->nbcolox-2)*M.coordinates[0];
	int nx = (L->nbcolox-2)*(M.coordinates[0]+1);
	if (M.dims[0]==1) {
		nxs=0;
		nx = L->ii;}
	else{
		if(M.coordinates[0]==0) {
			nxs=0;
			nx = L->nbcolox-1;}
		if (M.coordinates[0]==M.dims[0]-1) {
			nxs = L->ii-L->nbcolox+2;
			nx = L->ii;}}
	int nys = (L->nbcoloy-2)*M.coordinates[1];
	int ny = (L->nbcoloy-2)*(M.coordinates[1]+1);
	if (M.dims[1]==1) {
		nys =0;
		ny = L->jj;}
	else{
		if (M.coordinates[1]==0) {
			nys = 0;
			ny = L->nbcoloy-1;}
		if (M.coordinates[1]==M.dims[1]-1) {
			nys = L->jj-L->nbcoloy+2;
			ny = L->jj;}}
	int nzs = (L->nbcoloz-2)*M.coordinates[2];
	int nz = (L->nbcoloz-2)*(M.coordinates[2]+1);
	if (M.dims[2]==1) {
		nzs = 0;
		nz = L->kk;}
	else{
		if (M.coordinates[2]==0) {
			nzs = 0;
			nz = L->nbcoloz-1;}
		if (M.coordinates[2]==M.dims[2]-1) {
			nzs = L->kk-L->nbcoloz+2;
			nz = L->kk;}}
	MPI_Gather(&nx,1,MPI_INT,nxt,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&ny,1,MPI_INT,nyt,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&nz,1,MPI_INT,nzt,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&nxs,1,MPI_INT,nxts,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&nys,1,MPI_INT,nyts,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&nzs,1,MPI_INT,nzts,1,MPI_INT,0,MPI_COMM_WORLD);

	if (myid == 0)
	{
	static FILE *fpt;
	/// Open file
	char full_filename[1024];
    if (strstr(filename, ".pvti") != NULL){
      strcpy(full_filename, filename);}
    else{
      sprintf(full_filename, "Output/%s_l_%d_t_%d.pvti", filename,l,t);}
    fpt = fopen(full_filename, "wb+");
    /// Write the head
    char str[1024];
    fprintf(fpt, "<?xml version=\"1.0\"?>\n");
    fprintf(fpt, "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	sprintf(str, "<PImageData WholeExtent=\"0 %ld 0 %ld 0 %ld \" \
				GhostLevel=\"0\" Origin=\"%f %f %f\" Spacing=\"%f %f %f\">\n",
				L->ii,L->jj,L->kk,U->xa,U->ya,U->za,L->hx,L->hy,L->hz);
	fprintf(fpt,"%s",str);
	/// Write Point data

	int V_comp=output_option.compare("V");		// output field of vector
	int S_comp=output_option.compare("S");		// output field of scalar
	int VS_comp=output_option.compare("VS");		// output vector and scalar
	int Cell_comp=output_option.compare("C");		// output cell data
	if (S_comp==0)
	{
		fprintf(fpt,"<PPointData Scalars=\"%s\">\n",field_name);
		fprintf(fpt,"<PDataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n",field_name);
	}
	if (V_comp==0)
	{
		fprintf(fpt,"<PPointData Vectors=\"%s\">\n",field_name);
		fprintf(fpt,"<PDataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"3\"/>\n",field_name);
	}
	if (VS_comp==0)
	{
		fprintf(fpt,"<PPointData Scalars=\"E\" Vectors=\"Diplacement\">\n");
		fprintf(fpt,"<PDataArray Name=\"E\" type=\"Float32\" NumberOfComponents=\"1\"/>\n");
		fprintf(fpt,"<PDataArray Name=\"Diplacement\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
	}
	/// Write Cell Data
	if (Cell_comp==0)
	{
		fprintf(fpt,"<PPointData></PPointData>\n");
		fprintf(fpt,"<PCellData Scalar = \"%s\">\n",field_name);
		fprintf(fpt,"<PDataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n",field_name);
		fprintf(fpt,"</PCellData>\n");
	}
	else
	{
		fprintf(fpt,"</PPointData>\n");
		fprintf(fpt,"<PCellData></PCellData>\n");
	}	
	for (int i =0; i<nbprocs; i++)
	{
		sprintf(str,"<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s_l_%d_t_%d_Proc_%d.vti\"/>\n",
				nxts[i],nxt[i],nyts[i],nyt[i],nzts[i],nzt[i],filename,l,t,i);
		fprintf(fpt,"%s",str);
	}
	/// Write end
	fprintf(fpt,"</PImageData>\n");
	fprintf(fpt,"</VTKFile>\n");
    fclose(fpt);
    fpt = NULL;
	}
	delete [] nxt;
	delete [] nyt;
	delete [] nzt;
	delete [] nxts;
	delete [] nyts;
	delete [] nzts;
}
//-------- output scalar field: d or E
void write_total_de(Stack *U, const griddouble & de, const MPI_Setting & M, 
const string & output_option, const char *filename, const char *field_name, 
const int l, const int myid, const int nbprocs, const int t)
{
	open_file(filename,l,myid,t);
	write_header(U,M,l);
	write_pointdata_head(output_option,field_name);
	write_scalar(U,de,M,l);
	write_pointdata_end();
	write_celltdata();
	write_end();
	write_pvtk(U,filename,M,output_option,field_name,l,myid,nbprocs,t);
}
//-------- output scalar field: H
void write_total_H(Stack *U, const griddouble &H, const MPI_Setting & M, const string & output_option, 
const char *filename, const char *field_name, const int l, const int myid, const int nbprocs, 
const int t)
{
	open_file(filename,l,myid,t);
	write_header(U,M,l);
	write_celltdata_H(U,H,M,l,field_name);
	write_end();
	write_pvtk(U,filename,M,output_option,field_name,l,myid,nbprocs,t);
}
//-------- output vector field: u
void write_total_u(Stack *U, const gdvecdouble &u, const MPI_Setting & M, const int l, 
const int myid, const int nbprocs, const char *filename, const char *field_name, const int t,
const double u_proportion)
{
	std::string output_option="V";
	open_file(filename,l,myid,t);
	write_header(U,M,l);
	write_pointdata_head(output_option,field_name);
	write_vector(U,u,M,l,u_proportion);
	write_pointdata_end();
	write_celltdata();
	write_end();
	write_pvtk(U,filename,M,output_option,field_name,l,myid,nbprocs,t);
}
//-------- output scalar and vector field: E and u
void write_total_Eu(Stack *U, const gdvecdouble &u, const griddouble &de, const int l, 
const MPI_Setting & M, const int myid, const int nbprocs, const char *filename,
const char * field_name, const int t,const double u_proportion)
{
	std::string output_option="VS";
	open_file(filename,l,myid,t);
	write_header(U,M,l);
	write_pointdata_head(output_option,field_name);
	write_vector(U,u,M,l,u_proportion);
	write_scalar(U,de,M,l);
	write_pointdata_end();
	write_celltdata();
	write_end();
	write_pvtk(U,filename,M,output_option,field_name,l,myid,nbprocs,t);
}
