#ifndef _CALCULATIONS_HPP_
#define _CALCULATIONS_HPP_

#include "typedefs.hpp"
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include <vector>
#include <iostream>
#include <map>

using namespace std;

Real k_damping = 0.999;
Real EPSILON = 1e-6;

class Face {
public:
    vector<int> indices;
    vector<Face*> attachedFaces;
    int id;
    
    Face(vector<int> indicesInput, int inputId)
    : indices(indicesInput), attachedFaces(vector<Face*>()), id(inputId) {}
    ~Face() { indices.clear(); attachedFaces.clear(); }
    
    bool containsIndex(int i) {
        return (indices[0] == i || indices[1] == i || indices[2] == i );
    }
    
    bool containsIndices(int i, int j) {
        return containsIndex(i) && containsIndex(j);
    }
    
    bool containsOrderedIndices(int i, int j) {
        return ((indices[0] == i && indices[1] == j)
                || (indices[1] == i && indices[2] == j)
                || (indices[2] == i && indices[0] == j));
    }
    
    bool adjacentToFaceById(int faceId) {
        // returns true if my attached faces vector already contains this face
        for (int i = 0; i < (int) attachedFaces.size(); ++i) {
            if (attachedFaces[i]->id == faceId) {
                return true;
            }
        }
        return false;
    }
    
    int getAdjacentFaceIndexByIndices(int sharedI, int sharedJ) {
        for (int i = 0; i < (int) attachedFaces.size(); ++i) {
            if (attachedFaces[i]->containsIndices(sharedI, sharedJ)) {
                return i;
            }
        }
        // if (testing) { throw; }
        return -1;
    }
    
    bool adjacentToFace(Face* adjacentFace) {
        return adjacentToFaceById(adjacentFace->id);
    }
    
    bool shouldBeAdjacentToFace(Face* face) {
        // if not already adjacent and contains a shared edge with this face
        // adds face to vector of attachedFaces vector and returns true; false otherwise
        if (!adjacentToFace(face)) {
            if (containsIndices(face->indices[0], face->indices[1])
                || containsIndices(face->indices[1], face->indices[2])
                || containsIndices(face->indices[2], face->indices[0])) {
                attachedFaces.push_back(face);
                face->shouldBeAdjacentToFace(this);
                
                return true;
            }
        }
        return false;
    }
    
    bool removeAdjacentFaceById(int id) {
        if (!adjacentToFaceById(id)) {
            return false;
        }
        
        bool ret = false;
        vector<Face*> attachedNew = vector<Face*>();
        for (int i = 0; i < (int) attachedFaces.size() && !ret; ++i) {
            if (attachedFaces[i]->id != id) {
                attachedNew.push_back(attachedFaces[i]);
                ret = true;
            }
        }
        attachedFaces = attachedNew;
        return ret;
    }
    
    Face* getAdjacentById(int faceId) {
        if (!(adjacentToFaceById(faceId))) { throw; }
        
        for (int i = 0; i < (int) attachedFaces.size(); ++i) {
            if (attachedFaces[i]->id == faceId) {
                return (attachedFaces[i]);
            }
        }
        // if (testing) { throw; }
        return nullptr;
    }
    
    int getNonListedIndex(int pi_index, int pj_index) {
        // returns index of point making face's triangle that is not one of the listed inputs
        for (int i = 0; i < (int) indices.size(); ++i) {
            if (indices[i] != pi_index && indices[i] != pj_index) {
                return indices[i];
            }
        }
        // if (testing) { throw; }
        return -1;
    }
    
    void resetIndices(int index0, int index1, int index2) {
        indices[0] = index0;
        indices[1] = index1;
        indices[2] = index2;
    }
};

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
        // cout << "x_cm" << x_cm << endl;
        // cout << "v_cm" << v_cm << endl;
        
        // computing global linear velocity (L = sum ri x (mivi); I = sum ?????)
        Vec3f L = Vec3f(0, 0, 0);
        Mat3f I = Mat3f::zero();
        Vec3f r = Vec3f(0, 0, 0);
        vector<Vec3f> r_mat;
        
        for (int i=0; i<X.size(); i++) {
            r = X[i] - x_cm;
            r_mat.push_back(r);
            L += r.cross(V[i] * M[i]);
            Mat3f I_temp = Mat3f(0, -r[2], r[1], r[2], 0, -r[0], -r[1], r[0], 0);
            Mat3f I_temp_tr = Mat3f(0, r[2], -r[1], -r[2], 0, r[0], r[1], -r[0], 0);
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
    
    void mathForFaceBendingConstraintCreation(Face* f1, Face* f2, vector<int> facesIndices, vector<Vec3f> p, Real& phi) {
        // using same p ordering as in update for this constraint
        Vec3f n1 = (p[2] - p[0]).cross(p[3] - p[0]); n1.norm();
        Vec3f n2 = (p[3] - p[1]).cross(p[2] - p[1]); n2.norm();
        
        Real d = n1.dot(n2);
        // d = clamp(n1.dot(n2), -1, 1);
        if (d<-1) d = -1;
        if (d>1)  d = 1;
        if (isnan(d)) { throw; }
        
        phi = acos(d);
        if (abs(acos(d) - phi) < EPSILON) { return; }
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
