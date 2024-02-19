/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _gdvec_h
#define _gdvec_h
#include "grid.h"


template <class T>
class gdvec {

    public :
    //constructor
    gdvec(const int d, const vector<size_t> & dims);
    /// copy contructor
    gdvec(const gdvec <T> & in);
    // get the gird for each direction
    grid<T>& getgrid (const int di) {return *gv[di];};
    // const get the gird for each direction
    const grid<T>& getgrid (const int di) const {return *gv[di];} ;
    // return the dimenssion of gridvec
    int dims() const {return d;};
    
    // operator =
    gdvec<T> & operator=(const gdvec<T>& A);
    // operator +=
    gdvec<T> & operator+=(const gdvec<T>& A);
    // operator -=
    gdvec<T> & operator-=(const gdvec<T>& A);
    // operator *=
    gdvec<T> & operator*=(const T & a);
    // operator /=
    gdvec<T> & operator/=(const gdvec<T>& A);

    //deconstructor
    ~gdvec();
    private : 
    int d;
    vector<grid<T> *> gv;
};

using gdvecdouble = gdvec <double>; 

template<class T>
gdvec<T>::gdvec(const int _d, const vector<size_t> & dim):
d(_d), gv(d)
{
    for (int di=0; di<d; di++)
    gv[di] = new grid<T>(dim);
}

template<class T>
gdvec<T>::gdvec(const gdvec<T> & in):d(in.d)
{
    for (int di=0; di<d; di++)
    {
         gv[di] = new grid<T>(in.getgrid(0).axis());
         gv[di] = in.gv[di];
    }
}

template <class T>
gdvec<T> & gdvec<T>::operator=(const gdvec<T>& A)
{   
    if (this == &A) return *this;
    else 
    {
        d = A.d;
        gv.resize(d);
        for (int di=0; di<d; di++)
            gv[di] = new grid <T> (*A.gv[di]) ;
         return *this;
    }
}

template <class T>
gdvec<T> & gdvec<T>::operator+=(const gdvec<T>& A)
{   
    for (int di=0; di<d; di++)
        (*gv[di]) += (*A.gv[di]);
    return *this;
}

template <class T>
gdvec<T> & gdvec<T>::operator-=(const gdvec<T>& A)
{   
    for (int di=0; di<d; di++)
        gv[di] = gv[di] - A.gv[di];
    return *this;
}

template <class T>
gdvec<T> & gdvec<T>::operator*=(const T & a)
{   
    for (int di=0; di<d; di++)
        gv[di] = a*gv[di] ;
    return *this;
}

template <class T>
gdvec<T> & gdvec<T>::operator/=(const gdvec<T>& A)
{   
    for (int di=0; di<d; di++)
        gv[di] = gv[di] / A.gv[di];
    return *this;
}

template<class T>
gdvec<T>::~gdvec()
{
    for (int di=0; di<d; di++)
        delete gv[di];
}

///------------ no-member functions
/*
template <class T>
gdvec<T> operator+(const gdvec <T> & A, const gdvec <T> & B)
{
    gdvec <T> C(A);
    C += B;
    return C;
}

template <class T>
gdvec<T> operator-(const gdvec <T> & A, const gdvec <T> & B)
{
    gdvec <T> C(A);
    C -= B;
    return C;
}

// operator * for number*class gdvec
template <class T>
gdvec<T> operator*(const T & a, const gdvec <T> & A)
{
    gdvec<T> B(A);
    return B*=a;
}

// operator * for class gdvec*number
template <class T>
gdvec<T> operator*(const gdvec <T> & A, const T & a)
{
    return a*A;
}

template <class T>
gdvec<T> operator/(const gdvec <T> & A, const gdvec <T> & B)
{
    gdvec <T> C(A);
    C /= B;
    return C;
}
*/
#endif
