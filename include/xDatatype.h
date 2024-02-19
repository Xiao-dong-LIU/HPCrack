/*=========================================================================

Copyright (c) 2022 CNRS, Ecole Centrale de Nantes, Nantes Université

Author : Xiaodong LIU  xiaodong.liu@cnrs.fr

Institut de Recherche en Génie Civil et Mécanique (GeM) UMR6183

=========================================================================*/
#ifndef _XDATATYPE_H
#define _XDATATYPE_H
#include <iostream>
#include <numeric>
using namespace std;

template < typename T >
class xDataType{

public:
    static T zero() { throw; }
    static T one() { throw; }
    static std::string  stype() { throw; }

};

/// Specialization of xDatatype for short
template < >
class xDataType < short >
{
    public:
        static short zero() { return 0; }
        static short one() { return 1; }
        static std::string  stype() { return "int"; }
};

/// Specialization of xDatatype for int
template < >
class xDataType < int >
{
    public:
        static int zero() { return 0; }
        static int one() { return 1; }
        static std::string  stype() { return "int"; }
};

/// Specialization of xDatatype for float
template < >
class xDataType < float >
{
    public:
        static float zero() { return 0.f; }
        static float one() { return 1.f; }
        static std::string  stype() { return "float"; }
};

/// Specialization of xDatatype for double
template < >
class xDataType < double >
{
    public:
        static double zero() { return 0.; }
        static double one() { return 1.; }
        static std::string  stype() { return "double"; }

};



#endif