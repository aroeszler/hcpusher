//
//  vector_add.hpp
//  hcsecondtry.cpp
//
//  Created by Annegret Roeszler on 06.04.20.
//  Copyright © 2020 Annegret Roeszler. All rights reserved.
//

/**
* Header file implementing simple arithmetic vector operations vector
* for std::vector<double> type.
**/

#pragma once

#include <vector>
#include <cassert>
#include <sstream>

using Vec3D = std::vector<double>;

std::string vec_print (const Vec3D& v)
{
  std::stringstream ss;
  ss << "[";
  for( unsigned i = 0; i < v.size(); ++i )
  {
    if(i != 0)
      ss << ", ";
    ss << v[i];
  }
  ss << "]";
  return ss.str();
}

Vec3D vec_add (const Vec3D& v1, const Vec3D& v2)
{
  assert( v1.size() == v2.size() );
  Vec3D vec_sum( v1.size() );
  for ( unsigned i=0; i < v1.size(); ++i )
    {
      vec_sum[i] = v1[i] + v2[i];
    }
  return vec_sum;
}

Vec3D scalmultip (const Vec3D& v1, const double d1)
{
    assert( v1.size() == 3 );
    Vec3D vec_scal( v1.size() );
    for ( unsigned i=0; i < v1.size(); ++i )
      {
        vec_scal[i] = v1[i] * d1;
      }
    return vec_scal;
}

double innerprod (const Vec3D& v1, const Vec3D& v2)
{
    assert( v1.size() == v2.size() );
    double scal_prod = 0;
    for ( unsigned i=0; i < v1.size(); ++i )
      {
        scal_prod += v1[i] * v2[i];
      }
    return scal_prod;
    
}

Vec3D crossprod (const Vec3D& v1, const Vec3D& v2)
{
    assert( v1.size() == v2.size() );
    assert( v1.size() == 3 );
    Vec3D vec_prod( v1.size() );
    vec_prod[0] = v1[1] * v2[2] - v1[2] * v2[1];
    vec_prod[1] = v1[2] * v2[0] - v1[0] * v2[2];
    vec_prod[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return vec_prod;
}
