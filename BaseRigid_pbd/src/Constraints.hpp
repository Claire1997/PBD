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
};


#endif  /* _CALCULATIONS_HPP_ */
