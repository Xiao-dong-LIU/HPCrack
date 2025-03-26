/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _GRID_H
#define _GRID_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "xDatatype.h"
#include <math.h>
#include "omp.h"
using namespace std;

// 3d class 
template <class T> 
class grid{
private:
  const std::vector<size_t> dim;
  std::vector<T> val;
public:
  /// Construct a dense nx*ny*nz matrix
  grid(const vector<size_t> & dim);
  /// copy contructor
  grid(const grid <T> & in):dim(in.dim), val(in.val){};
  grid(grid<T>&& in) noexcept : dim(std::move(in.dim)), val(std::move(in.val)) {};
  // operator ()
  T & operator()(const size_t i, const size_t j, const size_t k);
  const T & operator()(const size_t i, const size_t j, const size_t k) const;
  // operator []
  T & operator[](const size_t i) { return val[i];};
  const T & operator[](const size_t i) const { return val[i];};
  // return the dimenssion of the table
  const vector<size_t> & axis() const { return dim;};
  // assign One
  grid<T> & Assign_one ();
   // assign zero
  grid<T> & Assign_zero ();
  // sqrt
  grid<T> & gridsqrt ();
  // return the size of grid
  std::size_t size();
  // operator =
  grid<T> & operator=(const grid<T>& A);
  // operator +=
  grid<T> & operator+=(const T& a);
  grid<T> & operator+=(const grid<T>& A);
  // operator -=
  grid<T> & operator-=(const T& a);
  grid<T> & operator-=(const grid<T>& A);
  // operator *=
  grid<T> & operator*=(const T&a );
  grid<T> & operator*=(const grid<T>& A );
  // operator /=
//  grid<T> & operator/=(const T& a );
  grid<T> & operator/=(const grid<T>& A );
};

//------------- Alias
using griddouble = grid <double>;
using gridshort = grid <short>;
using gridfloat = grid <float>;
//------------- member functions

// Constructor
template <class T> 
grid<T>::grid(const vector<size_t> & _dim):dim(_dim)
{
  size_t total_size = 1;
    for (size_t d : dim) {
        total_size *= d;
    }
    val = vector<T>(total_size, xDataType<T>::zero());  
}
template <class T> 
T & grid<T>::operator()(const size_t i, const size_t j, const size_t k) 
{
    return val[i+j*dim[0]+k*dim[0]*dim[1]];
}
template <class T>  
const T & grid<T>::operator()(const size_t i, const size_t j, const size_t k) const
{
    return val[i+j*dim[0]+k*dim[0]*dim[1]];
}

template <class T>
grid <T> & grid<T>::Assign_one ()
{
    #pragma omp parallel for
      for (size_t i=0; i<val.size(); i++)
        val[i] = xDataType<T>::one();
    return *this;
}
template <class T>
grid <T> & grid<T>::Assign_zero ()
{
    #pragma omp parallel for
      for (size_t i=0; i<val.size(); i++)
        val[i] = xDataType<T>::zero();
    return *this;
}

template <class T>
grid <T> & grid<T>::gridsqrt ()
{
   // std::transform(val.begin(), val.end(), val.begin(), (double(*)(double)) sqrt);

    #pragma omp parallel for
      for (size_t i=0; i<val.size(); i++)
        val[i] = static_cast<T>(std::sqrt(static_cast<double>(val[i])));
    return *this;
}

template <class T>  
std::size_t grid<T>::size() 
{
    return val.size();
}

template <class T> 
grid<T> &grid<T>::operator=(const grid<T>& A)
{
  if (this == &A) return *this;
  else
    throw std::runtime_error("Grid dimensions mismatch in = operator");
}


template <class T> 
grid<T> &grid<T>::operator+=(const T& a)
{
    #pragma omp parallel for
    for (size_t i=0;i<val.size();i++)
    {
                val[i]+=a;
    }
    return *this;
}

template <class T> 
grid<T> &grid<T>::operator+=(const grid<T>& A)
{
   if (A.dim==dim)
   {
       #pragma omp parallel for
        for (size_t i=0;i<val.size();i++)
            {
                val[i]+=A.val[i];
            }
        return *this;
   }
   else 
    throw std::runtime_error("Grid dimensions mismatch in += operator");
}


template <class T> 
grid<T> &grid<T>::operator-=(const T& a)
{
    #pragma omp parallel for
    for (size_t i=0;i<val.size();i++)
    {
                val[i]-=a;
    }
    return *this;
}

template <class T> 
grid<T> &grid<T>::operator-=(const grid<T>& A)
{
   if (A.dim==dim)
   {
       #pragma omp parallel for
        for (size_t i=0;i<val.size();i++)
            {
                val[i]-=A.val[i];
            }
        return *this;
   }
   else 
    throw std::runtime_error("Grid dimensions mismatch in -= operator");
}

template <class T>
grid<T> & grid<T>::operator*=(const T&a )
{
    #pragma omp parallel for
    for  (size_t i = 0; i < val.size(); ++i) 
        val[i] *= a;
    return *this;
}

template <class T>
grid<T> & grid<T>::operator*=(const grid<T>& A )
{
    #pragma omp parallel for
    for  (size_t i = 0; i < val.size(); ++i) 
        val[i] *= A.val[i];
    return *this;
}

template <class T>
grid<T> & grid<T>::operator/=(const grid<T>& A )
{
    #pragma omp parallel for
    for  (size_t i = 0; i < val.size(); ++i) 
    {
      if (A.val[i] == xDataType<T>::zero()) {
        throw std::runtime_error("Division by zero in grid /= operation");
    }
      val[i] =  val[i]/A.val[i];
    }
        
    return *this;
}
/*
template <class T>
grid<T> & grid<T>::operator/=(const T& a )
{
    #pragma omp parallel for
    for  (size_t i = 0; i < val.size(); ++i) 
        val[i] =  a/val[i];
    return *this;
}
*/
///------------ no-member functions
/// Print the matrix content, line by line in the stream.
/*!
    typical usage : std::cout << A;
    !*/
template <class T>
std::ostream& operator<<(std::ostream& os, const grid<T> & A)
{
    for (int j=0; j < A.axis()[1]; ++j)
     for (int i=0; i < A.axis()[0]; ++i)
     for (int k=0; k < A.axis()[2]; ++k)
     { 
        os <<"i="<<i<<",j="<<j<<",k="<<k<<": "<< A(i, j, k) <<endl;
     }
  return os;
}


// operator +
template <class T> 
grid<T> operator+(const grid<T>& A, const grid<T>& B) 
{
  if (B.axis()==A.axis())
  {
    grid<T> C(A);
    C+=B;
    return C;
  }
  else
    throw;
}


// operator -
template <class T> 
grid<T> operator-(const grid<T>& A, const grid<T>& B) 
{
  if (B.axis()==A.axis())
  {
    grid<T> C(A);
    C-=B;
    return C;
  }
  else
    throw;
}

// operator * for number*class grid
template <class T>
grid<T> operator*(const T &a, const grid <T> & A)
{
    grid<T> B(A);
    return B*=a;
}
// operator * for class grid*number
template <class T>
grid<T> operator*(const grid <T> & A, const T &a)
{
    return a*A;
}

// operator * for class grid*grid 
template <class T> 
grid<T> operator*(const grid<T>& A, const grid<T>& B) 
{
  if (B.axis()==A.axis())
  {
    grid<T> C(A);
    C*=B;
    return C;
  }
  else
    throw;
}


template <class T> 
grid<T> operator/(const grid<T>& A, const grid<T>& B) 
{
  if (B.axis()==A.axis())
  {
    grid<T> C(A);
    C/=B;
    return C;
   }
  else
    throw;
}

#endif