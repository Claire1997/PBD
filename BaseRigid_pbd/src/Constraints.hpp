#ifndef _CONSTRAINTS_HPP_
#define _CONSTRAINTS_HPP_

#include "typedefs.hpp"
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "Calculations.hpp"
#include <vector>
#include <iostream>
#include <map>

using namespace std;

class Constraint {
protected:
    bool remove;
    bool allowedToBreak;
public:
    Constraint(bool inputAllowedToBreak) : remove(false), allowedToBreak(inputAllowedToBreak) {}
    ~Constraint() {}
    
    virtual void update(vector<Vec3f> &p) {}
    bool broken() {
        return (allowedToBreak && remove);
    }
    Real compressionStiffness=1.0, stretchingStiffness=1.0;
    int constraintIterations = 10;
};

class StretchConstraint : public Constraint {
protected:
    /// parental items
    // bool remove;
    // bool allowedToBreak;
public:
    int pi_index;
    int pj_index;
    
    Real originalLength;
    Real wi;
    Real wj;
    
    StretchConstraint(Real inputLen, Real inputwi, Real inputwj, int inputpi_index, int inputpj_index, bool inputAllowedToBreak)
    : Constraint(inputAllowedToBreak), pi_index(inputpi_index), pj_index(inputpj_index), originalLength(inputLen), wi(inputwi), wj(inputwj) {}
    
    ~StretchConstraint() { Constraint::~Constraint(); }
    
    virtual void update(vector<Vec3f> &p){ //, map<int, vector<int>>& broken) {
        // check if any fixed points - then dont do constraints calculations
        if (wi == 0 && wj == 0) { return; }
        // if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(pi_index, pj_index, broken))) { remove = true; return; }
        Vec3f piMinpj = p[pi_index] - p[pj_index];
        Real len = piMinpj.length();
        if (len == 0) {
            cout << "prevent dividing by zero: length is 0 error in calculations.cpp" << pi_index << pj_index << endl;
            throw;
        }
        
        Real diff = len - originalLength;
        
        /*if (allowedToBreak && diff > 1.5 * originalLength) {
            mapCalculations::addToMap(pi_index, pj_index, broken);
        }*/
        
        Vec3f updateValue = diff * piMinpj / ((wi + wj) * len);
        // compression or stretching stiffness
        Real stiffness = (len < originalLength) ? compressionStiffness : stretchingStiffness;
        stiffness = 1 - pow((1 - stiffness), constraintIterations);
        if (pi_index==0) {
            p[pj_index] += (wi + wj) * updateValue * (stiffness);
        }else {
            p[pi_index] += -wi * updateValue * (stiffness);
            p[pj_index] += wj * updateValue * (stiffness);
        }
    }
};


class FaceBendingConstraint : public Constraint {
public:
    // notEdge, edge, edge, notEdge
    vector<int> facesIndices;
    // weights for each of these vertices
    vector<Real> w;
    // original angle between the faces
    Real phi;
    
    FaceBendingConstraint(vector<int> inputFacesIndices, vector<Real> inputW, Real inputPhi, bool inputAllowedToBreak)
    : Constraint(inputAllowedToBreak), facesIndices(inputFacesIndices), w(inputW), phi(inputPhi) {}
    
    ~FaceBendingConstraint() { Constraint::~Constraint(); }
    
    virtual void update(vector<Vec3f> &p, map<int, vector<int>>& broken) {
        if (w[0] == 0 || w[1] == 0 || w[2] == 0 || w[3] == 0) { return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[0], facesIndices[1], broken))) { remove = true; return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[1], facesIndices[2], broken))) { remove = true; return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[2], facesIndices[3], broken))) { remove = true; return; }
        
        V p0 = p.row(facesIndices[0]).transpose();
        V p1 = p.row(facesIndices[1]).transpose();
        V p2 = p.row(facesIndices[2]).transpose();
        V p3 = p.row(facesIndices[3]).transpose();
        
        V n1 = (p2 - p0).cross(p3 - p0); n1 /= n1.squaredNorm(); //n1.normalize();
        V n2 = (p3 - p1).cross(p2 - p1); n2 /= n2.squaredNorm(); //n2 *= T(-1); n2.normalize();
        
        V p32 = p3 - p2;
        T len_p32 = p32.norm();
        T len_p12 = p1.cross(p2).norm();
        T len_p13 = p1.cross(p3).norm();
        
        if (len_p12 < EPSILON  || len_p13 < EPSILON) { return; }
        
        DM q = DM(dim, 4);
        
        q.col(0) = len_p32 * n1;
        q.col(1) = len_p32 * n2;
        q.col(2) = (p0 - p3).dot(p32) / len_p32 * n1 + (p1 - p3).dot(p32) / len_p32 * n2;
        q.col(3) = (p2 - p0).dot(p32) / len_p32 * n1 + (p2 - p1).dot(p32) / len_p32 * n2;
        
        n1.normalize();
        n2.normalize();
        
        T d = clamp(n1.dot(n2), T(-1), T(1));
        if (testing && isnan(d)) { throw; }
        T testPhi = acos(d);
        if (testing && abs(testPhi - phi) < EPSILON) { return; }
        
        T weighting = 0;
        for (int i = 0; i < 4; ++i) {
            weighting += w[i] * q.col(i).squaredNorm();
        }
        if (weighting == T(0)) { return; }
        
        // note in below line: (1 - sim.stiffness) --> 1 - (1-k)^iterations bc more than one iter inside the
        // simulation loop [reduces error] (M. Muller et al. / Position Based Dynamics Section 3.3)
        T stiffness = sim.bendingStiffness;
        stiffness = 1 - pow((1 - stiffness), sim.constraintIterations);
        weighting = (testPhi - phi) / weighting * stiffness;
        if (n1.cross(n2).dot(p32) > 0) { weighting *= T(-1); }
        
        for (int i = 0; i < 4; ++i) {
            p.row(facesIndices[i]) +=  - w[i] * weighting * q.col(i).transpose();
            if (testing) {
                for (int test_i = 0; test_i < 3; ++test_i) {
                    if ( isnan(p.row(facesIndices[i])(test_i)) ) { throw; }
                }
            }
        }
    }
};
/*
class BalloonVolumeConstraint : public Constraint {
 public:
 BalloonVolumeConstraint(bool inputAllowedToBreak) : Constraint(inputAllowedToBreak) {}
 ~BalloonVolumeConstraint() { Constraint::~Constraint(); }
 
 virtual void update(vector<Vec3f>& p, map<int, vector<int>>& broken) {
 // if allowedToBreak and broken -- return
 }
};

class CollisionConstraint {
public:
    // index of the position needing to be updated
    int index;
    bool ground;
    bool sphere;
    V_horizontal qc;
    V_horizontal nc;
    T w;
    
    CollisionConstraint(int inputIndex, bool inputGround, bool inputSphere, T inputWeight)
    : index(inputIndex), ground(inputGround), sphere(inputSphere), w(inputWeight) {}
    
    ~CollisionConstraint() { }
    
    void update(vector<Vec3f> &p) {}
    
 };
 */


class Constraints {
public:
    int constraintIterations = 10;
    
    vector<Constraint*> stretch;
    // vector<CollisionConstraint*> collisions;
    
    Constraints() : stretch(vector<Constraint*>()){};
    
    ~Constraints() { stretch.clear(); }
    
    void update(vector<Vec3f>& p ) {  // , map<int, vector<int>>& broken) {
        for (Constraint* c : stretch) {
            c->update(p);
        }
    }
    
    void createStretchConstraints(vector<Vec3f> &positions, vector<Real> w, bool allowedToBreak) {
        // for string
        Real len = (positions[1] - positions[0]).length();
        cout << "len" << len << " " << w.size() << endl;
        for (int i=0; i<w.size()-1; i++){
            stretch.push_back(new StretchConstraint(len, w[i], w[i+1], i, i+1, allowedToBreak) );
        }
    }
    // void createFaceBendingConstraints(vector<Vec3f>& p, vector<Real> w, bool allowedToBreak) {}
    // void createBalloonVolumeConstraints(vector<Vec3f> &p, bool allowedToBreak) {}
    // void createCollisionConstraints(vector<Vec3f> p, vector<Real> w) {}
};


#endif  /* _CALCULATIONS_HPP_ */
