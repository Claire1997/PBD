#ifndef _CALCULATIONS_HPP_
#define _CALCULATIONS_HPP_

#include "typedefs.hpp"
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include <vector>
#include <iostream>
#include <map>

using namespace std;

Real k_damping = 0.01;

namespace calculations {
    Real sumM(vector<Real> M) {
        Real s = 0;
        for (int i=0; i<M.size(); i++) {
            s += M[i];
        }
        return s;
    }
    
    void dampVelocities_simple(vector<Vec3f> &V) {
        for (int i=1; i<V.size(); i++) {
            V[i] *= k_damping;
        }
    }
    
    void dampVelocities(vector<Vec3f> X, vector<Vec3f> &V, vector<Real> M) {
        // calculate center of mass / center of velocity / overall mass
        Real m_sum = sumM(M);
        Vec3f x_cm = Vec3f(0, 0, 0);
        Vec3f v_cm = Vec3f(0, 0, 0);
        for (int i=0; i<X.size(); i++) {
            x_cm += X[i] * M[i];
            v_cm += V[i] * M[i];
        }
        x_cm /= m_sum;
        v_cm /= m_sum;
        cout << "x_cm" << x_cm << endl;
        cout << "v_cm" << v_cm << endl;
        
        // computing global linear velocity (L = sum ri x (mivi); I = sum ?????)
        Vec3f L = Vec3f(0, 0, 0);
        Mat3f I = Mat3f::zero();
        Vec3f r = Vec3f(0, 0, 0);
        vector<Vec3f> r_mat;
        
        for (int i=0; i<X.size(); i++) {
            r = X[i] - x_cm;
            r_mat.push_back(r);
            L += r.cross(V[i] * M[i]);
            Mat3f I_temp = Mat3f(0, -r[2], r[1], -r[2], 0, r[0], r[1], -r[0], 0);
            Mat3f I_temp_tr = Mat3f(0, r[2], -r[1], r[2], 0, -r[0], -r[1], r[0], 0);
            I += I_temp * I_temp_tr * M[i];
        }
        
        cout << L << endl;
        I.printMat();
        Vec3f omega = I.inverse() * L;
        
        for (int i=1; i<X.size(); i++) {
            Vec3f delta_v = v_cm + omega.cross(r_mat[i]) - V[i];
            V[i] += k_damping * delta_v;
        }
        
    }
}


namespace mapCalculations {
    bool mapContainsKeyValuePair(int i, int j, map<int, vector<int> > &map1) {
        // always want i to be smaller than j for map iterations [removes the issue of duplicating elements when checking
        if (j < i) { swap(i, j); }
        if (map1.find(i) == map1.end()) {
            return false;
        }
        vector<int> vectorAt = map1.at(i);
        return (find(vectorAt.begin(), vectorAt.end(), j) == vectorAt.end());
    }
    
    void addToMap(int i, int j, map<int, vector<int> > &map1) {
        if (j < i) { swap(i, j); }
        if ( mapContainsKeyValuePair(i, j, map1)) { throw; }
        vector<int> adding = (map1.find(i) == map1.end()) ? vector<int>() : map1.at(i);
        adding.push_back(j);
        map1[i] = adding;
    }
}

#endif  /* _CALCULATIONS_HPP_ */
