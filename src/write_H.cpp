/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#include "write_H.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h> 

using namespace std;
// --------- write solution to binary data to avoid previous time step computing
void write_vector_3D(griddouble & H, const string & filename, const int myid, const int t)
{
  /// Open file
	char full_filename[1024];
  if (strstr(filename.c_str(), ".dat") != NULL)
  {
      strcpy(full_filename, filename.c_str());
  }
  else
  {
      sprintf(full_filename, "Time_step_data/%s_t_%d_Proc_%d.dat", filename.c_str(),t,myid);
  }

  ofstream wf(full_filename, ios::out | ios::binary);
   if(!wf) 
   {
      if (myid==0) cout << "Cannot open file!" << endl;
      exit(-1);
   }
    for (int k=0;k<H.axis()[2];k++)
      for (int j=0;j<H.axis()[1];j++)
        for (int i=0;i<H.axis()[0];i++)
       {
           wf.write(reinterpret_cast<char*> (&H(i,j,k)), sizeof(double));
          // if (myid==3) cout<<"i="<<i<<",j="<<j<<",k="<<k<<",value="<<H(l-1,i,j,k)<<endl;
       }
   wf.close();
   if(!wf.good()) {
     if (myid==0) cout << "Error occurred at writing !" << endl;
      exit(-1);
   }
   else 
   {
     if (myid==0) cout<<"Write completed!"<<endl;
   }
}
// --------- read solution from binary data stocked by the previous computation
void read_vector_3D(griddouble & H, const string & filename, const int myid, const int t)
{
  /// Open file
	char full_filename[1024];
    if (strstr(filename.c_str(), ".dat") != NULL){
      strcpy(full_filename, filename.c_str());}
    else{
      sprintf(full_filename, "Time_step_data/%s_t_%d_Proc_%d.dat", filename.c_str(),t,myid);}

  ifstream rf(full_filename, ios::out | ios::binary);
   if(!rf) 
   {
      if (myid==0) cout << "Cannot open file!" << endl;
      exit(-1);
   }
  for (int k=0;k<H.axis()[2];k++)
      for (int j=0;j<H.axis()[1];j++)
        for (int i=0;i<H.axis()[0];i++)
       {
           rf.read(reinterpret_cast<char*> (&H(i,j,k)), sizeof(double));
         //  if (myid==3) cout<<"i="<<i<<",j="<<j<<",k="<<k<<",value="<<H(l-1,i,j,k)<<endl;
       }
   rf.close();
   if(!rf.good()) {
      if (myid==0) cout << "Error occurred at reading!" << endl;
      exit(-1);
   }
   else 
   {
     if (myid==0) cout<<"Read completed!"<<endl;
   }
}
