// ----------------------------------------------------------------------------
// Matrix3x3.hpp
//
//  Created on: 13 Feb 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: 3x3 Matrix (DO NOT DISTRIBUTE!)
//
// Copyright 2020 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _MATRIX3X3_HPP_
#define _MATRIX3X3_HPP_

#include <cassert>

#include "typedefs.hpp"
#include "Vector3.hpp"

template<typename T> class Vector3;

template<typename T>
class Matrix3x3 {
public:
  enum {ROW=3, COL=3};

  typedef T ValueT;

  explicit Matrix3x3(
    const T &p00=1, const T &p01=0, const T &p02=0,
    const T &p10=0, const T &p11=1, const T &p12=0,
    const T &p20=0, const T &p21=0, const T &p22=1) {
    v[0][0]=p00; v[0][1]=p01; v[0][2]=p02;
    v[1][0]=p10; v[1][1]=p11; v[1][2]=p12;
    v[2][0]=p20; v[2][1]=p21; v[2][2]=p22;
  }

  explicit Matrix3x3(const Vector3<T> &diag) {
    v[0][0]=diag.x; v[0][1]=0; v[0][2]=0;
    v[1][0]=0; v[1][1]=diag.y; v[1][2]=0;
    v[2][0]=0; v[2][1]=0; v[2][2]=diag.z;
  }

  Matrix3x3(const Vector3<T> &c0, const Vector3<T> &c1, const Vector3<T> &c2) {
    v[0][0]=c0.x; v[0][1]=c1.x; v[0][2]=c2.x;
    v[1][0]=c0.y; v[1][1]=c1.y; v[1][2]=c2.y;
    v[2][0]=c0.z; v[2][1]=c1.z; v[2][2]=c2.z;
  }

  // assignment operators
  Matrix3x3& operator+=(const Matrix3x3 &m) {
    v00 += m.v00; v01 += m.v01; v02 += m.v02;
    v10 += m.v10; v11 += m.v11; v12 += m.v12;
    v20 += m.v20; v21 += m.v21; v22 += m.v22;
    return *this;
  }
  Matrix3x3& operator-=(const Matrix3x3 &m) {
    v00 -= m.v00; v01 -= m.v01; v02 -= m.v02;
    v10 -= m.v10; v11 -= m.v11; v12 -= m.v12;
    v20 -= m.v20; v21 -= m.v21; v22 -= m.v22;
    return *this;
  }
  Matrix3x3& operator*=(const T s) {
    v00 *= s; v01 *= s; v02 *= s;
    v10 *= s; v11 *= s; v12 *= s;
    v20 *= s; v21 *= s; v22 *= s;
    return *this;
  }
  Matrix3x3& operator/=(const T s) {
    v00 /= s; v01 /= s; v02 /= s;
    v10 /= s; v11 /= s; v12 /= s;
    v20 /= s; v21 /= s; v22 /= s;
    return *this;
  }

  // binary operators
  Matrix3x3 operator+(const Matrix3x3 &m) const { return Matrix3x3(*this)+=m; }
  Matrix3x3 operator-(const Matrix3x3 &m) const { return Matrix3x3(*this)-=m; }
  Matrix3x3 operator*(const Matrix3x3 &m) const {
    return Matrix3x3(
      v00*m.v00 + v01*m.v10 + v02*m.v20,
      v00*m.v01 + v01*m.v11 + v02*m.v21,
      v00*m.v02 + v01*m.v12 + v02*m.v22,

      v10*m.v00 + v11*m.v10 + v12*m.v20,
      v10*m.v01 + v11*m.v11 + v12*m.v21,
      v10*m.v02 + v11*m.v12 + v12*m.v22,

      v20*m.v00 + v21*m.v10 + v22*m.v20,
      v20*m.v01 + v21*m.v11 + v22*m.v21,
      v20*m.v02 + v21*m.v12 + v22*m.v22);
  }
  Matrix3x3 operator*(const T s) const { return Matrix3x3(*this)*=s; }
  Vector3<T> operator*(const Vector3<T> &v) const {
    return Vector3<T>(
      v00*v.x+v01*v.y+v02*v.z,
      v10*v.x+v11*v.y+v12*v.z,
      v20*v.x+v21*v.y+v22*v.z);
  }

  bool operator==(const Matrix3x3 &m) const {
    return
      (v00==m.v00 && v01==m.v01 && v02==m.v02 &&
       v10==m.v10 && v11==m.v11 && v12==m.v12 &&
       v20==m.v20 && v21==m.v21 && v22==m.v22);
  }

  const T& operator()(const tIndex r, const tIndex c) const { return v[r][c]; }
  T& operator()(const tIndex r, const tIndex c) {
    return const_cast<T &>(const_cast<const Matrix3x3 &>(*this)(r, c));
  }

  static Matrix3x3 I() { return Matrix3x3(1,0,0, 0,1,0, 0,0,1); }
  static Matrix3x3 zero() { return Matrix3x3(0,0,0, 0,0,0, 0,0,0); }
    
  Matrix3x3 getCofactor(int p, int q, int n){
    int i = 0, j = 0;
    Matrix3x3 temp = zero();
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++){
        for (int col = 0; col < n; col++){
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q){
                temp(i, j++) = v[row][col];
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1){
                    j = 0;
                    i++;
                }
            }
        }
    }
    return temp;
  }
    
  Real determinant(int n) {
    Real res = 0; // Initialize result
    //  Base case : if matrix contains single element
    if (n == 1)
        return v[0][0];
    Matrix3x3 temp = zero(); // To store cofactors
    int sign = 1;  // To store sign multiplier
    // Iterate for each element of first row
    for (int f = 0; f < n; f++){
        // Getting Cofactor of A[0][f]
        temp = getCofactor(0, f, n);
        res += sign * v[0][f] * temp.determinant(n-1);
        // terms are to be added with alternate sign
        sign = -sign;
    }
    return res;
 }

 Matrix3x3 adjoint(){
    // temp is used to store cofactors
     int sign = 1;
     Matrix3x3 temp = zero();
     Matrix3x3 adj = zero();
     for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            // Get cofactor of A[i][j]
            temp = getCofactor(i, j, 3);
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj(j, i) = (sign)*(temp.determinant(2));
        }
    }
    return adj;
 }
  
 Matrix3x3 inverse() {
    Real det = determinant(3);
    if (det == 0){
        std::cout << "Singular matrix, can't find its inverse";
        throw;
    }
    // Find adjoint
    Matrix3x3 adj = adjoint();
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    adj /= det;
    return adj;
  }
    
    void printMat() {
        std::cout << v00 << " " << v01 << " " << v02 << std::endl;
        std::cout << v10 << " " << v11 << " " << v12 << std::endl;
        std::cout << v20 << " " << v21 << " " << v22 << std::endl;
    }
  

  union {
    struct { T v00, v01, v02, v10, v11, v12, v20, v21, v22; };
    T v[3][3];
    T v1[9];
  };
};

typedef Matrix3x3<Real> Mat3f;

template <typename T>
inline
Matrix3x3<T>
operator*(const T s, const Matrix3x3<T> &m) {
  return m*s;
}

#endif  /* _MATRIX3X3_HPP_ */
