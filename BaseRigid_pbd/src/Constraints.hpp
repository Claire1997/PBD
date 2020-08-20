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
    
    virtual void update(vector<Vec3f> &p, map<int, vector<int> >& broken) {}
    bool broken() {
        return (allowedToBreak && remove);
    }
    Real compressionStiffness=1.0, stretchingStiffness=1.0, bendingStiffness = 1.0;;
    int constraintIterations = 10;
};

/*
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
    
    virtual void update(vector<Vec3f> &p), map<int, vector<int>>& broken) {
        // check if any fixed points - then dont do constraints calculations
        if (wi == 0 && wj == 0) { return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(pi_index, pj_index, broken))) { remove = true; return; }
        Vec3f piMinpj = p[pi_index] - p[pj_index];
        Real len = piMinpj.length();
        if (len == 0) {
            cout << "prevent dividing by zero: length is 0 error in calculations.cpp" << pi_index << pj_index << endl;
            throw;
        }
        
        Real diff = len - originalLength;
        
        if (allowedToBreak && diff > 1.5 * originalLength) {
            mapCalculations::addToMap(pi_index, pj_index, broken);
        }
        
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
*/

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
    
    virtual void update(vector<Vec3f> &p, map<int, vector<int> >& broken){ //, map<int, vector<int>>& broken) {
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
    
    virtual void update(vector<Vec3f> &p, map<int, vector<int> >& broken) {
        if (w[0] == 0 || w[1] == 0 || w[2] == 0 || w[3] == 0) { return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[0], facesIndices[1], broken))) { remove = true; return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[1], facesIndices[2], broken))) { remove = true; return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[2], facesIndices[3], broken))) { remove = true; return; }
        
        Vec3f p0 = p[facesIndices[0]];
        Vec3f p1 = p[facesIndices[1]];
        Vec3f p2 = p[facesIndices[2]];
        Vec3f p3 = p[facesIndices[3]];
        
        Vec3f n1 = (p2 - p0).cross(p3 - p0); n1 /= n1.lengthSquare();
        Vec3f n2 = (p3 - p1).cross(p2 - p1); n2 /= n2.lengthSquare();
        
        Vec3f p32 = p3 - p2;
        Real len_p32 = p32.length();
        Real len_p12 = p1.cross(p2).length();
        Real len_p13 = p1.cross(p3).length();
        
        if (len_p12 < EPSILON  || len_p13 < EPSILON) { return; }
        
        vector<Vec3f> q;
        
        q.push_back(len_p32 * n1);
        q.push_back(len_p32 * n2); // q.col(1) = len_p32 * n2;
        q.push_back( (p0 - p3).dot(p32) / len_p32 * n1 + (p1 - p3).dot(p32) / len_p32 * n2);
        q.push_back( (p2 - p0).dot(p32) / len_p32 * n1 + (p2 - p1).dot(p32) / len_p32 * n2);
        
        n1.norm();
        n2.norm();
        
        Real d = n1.dot(n2);
        // d = clamp(n1.dot(n2), -1, 1);
        if (d<-1) d = -1;
        if (d>1)  d = 1;
        if (isnan(d)) { throw; }
        Real testPhi = acos(d);
        if (abs(testPhi - phi) < EPSILON) { return; }
        
        Real weighting = 0;
        for (int i = 0; i < 4; ++i) {
            weighting += w[i] * q[i].lengthSquare();
        }
        if (weighting == 0) { return; }
        
        // note in below line: (1 - sim.stiffness) --> 1 - (1-k)^iterations bc more than one iter inside the
        // simulation loop [reduces error] (M. Muller et al. / Position Based Dynamics Section 3.3)
        Real stiffness = bendingStiffness;
        stiffness = 1 - pow((1 - stiffness), constraintIterations);
        weighting = (testPhi - phi) / weighting * stiffness;
        if (n1.cross(n2).dot(p32) > 0) { weighting *= -1; }
        
        for (int i = 0; i < 4; ++i) {
            p[facesIndices[i]] +=  - w[i] * weighting * q[i];
           // if (testing) {
                for (int test_i = 0; test_i < 3; ++test_i) {
                    if ( isnan(p[facesIndices[i]][test_i]) ) { throw; }
                }
            // }
        }
    }
};

class BalloonVolumeConstraint : public Constraint {
 public:
 BalloonVolumeConstraint(bool inputAllowedToBreak) : Constraint(inputAllowedToBreak) {}
 ~BalloonVolumeConstraint() { Constraint::~Constraint(); }
 
 virtual void update(vector<Vec3f>& p, map<int, vector<int> >& broken) {
 // if allowedToBreak and broken -- return
 }
};

class CollisionConstraint {
public:
    // index of the position needing to be updated
    int index;
    bool ground;
    bool sphere;
    Vec3f qc;
    Vec3f nc;
    Real w;
    
    CollisionConstraint(int inputIndex, bool inputGround, bool inputSphere, Real inputWeight)
    : index(inputIndex), ground(inputGround), sphere(inputSphere), w(inputWeight) {}
    
    ~CollisionConstraint() { }
    
    void update(vector<Vec3f> &p) {
        // collisions are an inequality constraint
        
        // zero length string - if p is within an object OR if p is below ground
        
        // collision is inequality constraint with stiffness = 1;
        
        // static: compute the closest surface point qc and surface normal nc at this location
        // continuous: we test for each vertix i the ray xi -> p.
        //      if the ray enters an object, we calculate the constraint between
        //      p and the entry point qc and the surface normal at this point, nc.
        // (diff between static and continuous is where the qc will be places)
        
        Real Cp = 0;
        Vec3f diff = Vec3f(0, 0, 0);
        
        if (ground) {
            // GROUND
            nc = Vec3f(0, 1, 0);
            qc = Vec3f(p[index][0], 0, p[index][2]); // hardcoding qc for ground constraint
            diff = p[index] - qc;
            Cp = diff.dot(nc);
            
            if (Cp < 0) {
                p[index] += -Cp * 1 * nc; // - angle * k-collisionrestitution * in direction of normal //qc;
            }
        } else if (sphere) {
            // SPHERE
            Vec3f center = Vec3f(5, 6.75, -3.5);
            Vec3f centerVelocity = Vec3f(-14 / 300, 0, 0);
            Real radius = 3;
            nc = (p[index] - center); nc /= nc.length();
            qc = radius*nc + center; // position - center
            diff = p[index] - qc;
            Cp = diff.dot(nc); // dot product
            
            // -angle * k-collisionrestitution * in direction of normal
            if (Cp < 0) {
                p[index] += (-Cp * 1 * nc);
            }
        } // else { if (testing) { throw; }}
    }
    
 };
 


class Constraints {
public:
    int constraintIterations = 10;
    
    vector<Constraint*> stretchbend;
    vector<CollisionConstraint*> collisions;
    
    Constraints() : stretchbend(vector<Constraint*>()),
                    collisions(vector<CollisionConstraint*>()){};
    
    ~Constraints() { stretchbend.clear(); collisions.clear();}
    
    void update(vector<Vec3f>& p, map<int, vector<int> >& broken) {  // , map<int, vector<int>>& broken) {
        for (Constraint* c : stretchbend) {
            c->update(p, broken);
        }
        for (CollisionConstraint* c : collisions) {
            c->update(p);
        }
    }
    
    /*
     void createStretchConstraints(vector<Face*> &faces, vector<Vec3f> &positions, vector<Real> &w, bool allowedToBreak) {
        // for box
        vector<vector<int> > check(sim.numParticles, vector<int>(sim.numParticles));
        for(int i=0; i<sim.numParticles; i++) {
            for(int j=0; j<sim.numParticles; j++){
                check[i][j] = 0;
            }
        }
        
        for (int i = 0; i < (int) faces.size(); ++i) {
            int size = faces[i]->indices.size();
            for (int j = 0; j < size; ++j) {
                int p1 = faces[i]->indices[j];
                int p2 = (j == size - 1) ? faces[i]->indices[0] : faces[i]->indices[j+1];
                
                if (check[p1][p2] != 1 && p1 != p2) {
                    if (check[p2][p1] == 1) { throw; }
                    check[p1][p2] = 1;
                    check[p2][p1] = 1;
                    Real len = ( positions[p1] - positions[p2]).length();
                    stretchbend.push_back(new StretchConstraint(len, w[p1], w[p2], p1, p2, allowedToBreak) );
                }
            }
        }
    }*/
    
    void createStretchConstraints(vector<Face*> &faces, vector<Vec3f> &positions, vector<Real> &w, bool allowedToBreak) {
        // for box
        int index_list[][2] =
        {
            {0, 1},
            {2, 3},
            {4, 5},
            {6, 7},
            {0, 3},
            {1, 2},
            {4, 7},
            {5, 6},
            {0, 4},
            {1, 5},
            {7, 3},
            {2, 6}
        };
        
        for (int i = 0; i < 12; ++i) {
            int p1 = index_list[i][0], p2 = index_list[i][1];
            Real len = ( positions[p1] - positions[p2]).length();
            stretchbend.push_back(new StretchConstraint(len, w[p1], w[p2], p1, p2, allowedToBreak) );
        }
    }
    
    /*
    void createStretchConstraints(vector<Vec3f> &positions, vector<Real> &w, bool allowedToBreak) {
        // for string
        Real len = (positions[1] - positions[0]).length();
        cout << "len" << len << " " << w.size() << endl;
        for (int i=0; i<w.size()-1; i++){
            stretchbend.push_back(new StretchConstraint(len, w[i], w[i+1], i, i+1, allowedToBreak) );
        }
    }*/
    /*
    void createFaceBendingConstraints(vector<Face*>& faces, vector<Vec3f>& p, vector<Real> &w, bool allowedToBreak){
        vector<int> edge = vector<int>();
        vector<int> notEdge = vector<int>();
        vector<int> facesIndices = vector<int>();
        vector<Vec3f> positions = vector<Vec3f>();
        vector<Real> weights = vector<Real>();
        Real phi = 0;
        cout << "facesize " << faces.size() << endl;
        vector<vector<int> > check(faces.size(), vector<int>(faces.size()));
        for(int i=0; i<faces.size(); i++) {
            for(int j=0; j<faces.size(); j++){
                check[i][j] = 0;
            }
        }
        
        for (int i = 0; i < faces.size(); ++i) {
            for (int j = 0; j < faces.size(); ++j) {
                if (i != j && check[i][j] != 1 && faces[i]->adjacentToFaceById(j)) {
                    check[i][j] = 1;
                    check[j][i] = 1; // ensures i will always be less than j
                     //cout << "i" << i << "j" << j <<endl;
                    edge.clear();
                    notEdge.clear();
                    facesIndices.clear();
                    positions.clear();
                    weights.clear();
                    //cout << "finish" <<endl;
                    // add in the first two indices for the shared edge using face j indices
                    for (int k = 0; k < 4; ++k) { // dim = 3
                        if (find(faces[i]->indices.begin(), faces[i]->indices.end(), faces[j]->indices[k]) != faces[i]->indices.end()) {
                            edge.push_back(faces[j]->indices[k]);
                        }
                    }
                    // cout << "edge" << edge[0] << edge[1] <<endl;
                    // add in the indices not on the shared edge
                    notEdge.push_back(faces[i]->getNonListedIndex(edge[0], edge[1]));
                    notEdge.push_back(faces[j]->getNonListedIndex(edge[0], edge[1]));
                    //cout << "notedge" << notEdge[0] << notEdge[1] <<endl;
                    if (edge.size() != 2) { throw; }
                    if (notEdge.size() != 2) { throw; }
     
                    // index ordering for bending = nonEdge, edge, edge, nonEdge
                    
                    weights.push_back(w[notEdge[0]]);  weights.push_back(w[edge[0]]);  weights.push_back(w[edge[1]]);  weights.push_back(w[notEdge[1]]);
                    facesIndices.push_back(notEdge[0]);  facesIndices.push_back(edge[0]);  facesIndices.push_back(edge[1]);  facesIndices.push_back(notEdge[1]);
                    positions.push_back(p[facesIndices[0]]); positions.push_back(p[facesIndices[1]]); positions.push_back(p[facesIndices[2]]); positions.push_back(p[facesIndices[3]]);
                    phi = acos(0.0);
     
                    calculations::mathForFaceBendingConstraintCreation(faces[i], faces[j], facesIndices, positions, phi);
                    cout << facesIndices[0] << facesIndices[1] << facesIndices[2] << facesIndices[3]  << " " << weights[0] << " "<< phi << " " << allowedToBreak << endl;
                    stretchbend.push_back(new FaceBendingConstraint(facesIndices, weights, phi, allowedToBreak));
                }
            }
        }
        cout << "end bend" <<endl;
    }*/
    void createFaceBendingConstraints(vector<Face*>& faces, vector<Vec3f>& p, vector<Real> &w, bool allowedToBreak){
        // for box
        double pi = 2*acos(0.0);
        int fi[24][4] =
        {
            {0,2,1,3},
            {0,7,3,4},
            {4,6,5,7},
            {1,6,2,5},
            {0,5,1,4},
            {2,7,3,6},
            
            {0,1,2,5},
            {0,3,2,7},
            {0,4,7,5},
            {1,2,0,6},
            {1,5,6,0},
            {2,3,0,7},
            {2,6,1,7},
            {3,7,0,2},
            {4,5,6,0},
            {4,7,0,6},
            {5,6,4,1},
            {6,7,4,2},
        };
        vector<int> facesIndices = vector<int>();
        vector<Vec3f> positions = vector<Vec3f>();
        vector<Real> weights = vector<Real>();
        Real phi;
        for (int i=0; i<18; i++) {
            facesIndices.clear();
            positions.clear();
            weights.clear();

            facesIndices.push_back(fi[i][2]);facesIndices.push_back(fi[i][0]);facesIndices.push_back(fi[i][1]);facesIndices.push_back(fi[i][3]);
            for (int j=0; j<4; j++) {
                positions.push_back(p[facesIndices[j]]);
                weights.push_back(w[facesIndices[j]]);
            }
            if (i<6) phi = 0;
            else phi = pi/2;
            // cout << facesIndices[0] << facesIndices[1] << facesIndices[2] << facesIndices[3]  << " " << weights[0] << " "<< phi << " " << allowedToBreak << endl;
            stretchbend.push_back(new FaceBendingConstraint(facesIndices, weights, phi, allowedToBreak));
        }
    }
    
    // void createBalloonVolumeConstraints(vector<Vec3f> &p, bool allowedToBreak) {}
    // void createCollisionConstraints(vector<Vec3f> p, vector<Real> w) {}
    void createCollisionConstraints(vector<Vec3f> &p, vector<Real> &w) {
        collisions.clear();
        bool sphereCollision = false;
        bool groundCollision = false;
        for (int i = 0; i < sim.numParticles; ++i) {
            // GROUND
            Vec3f qc = Vec3f(0, 0, 0);
            Vec3f nc = Vec3f(0, 1, 0);
            Real Cp = 0;
            Vec3f diff = Vec3f(0, 0, 0);
            
            qc = Vec3f(p[i][0], 0, p[i][2]); // hardcoding qc bc just grounding constraint for now
            diff = p[i] - qc;
            Cp = diff.dot(nc);
            
            if (Cp < 0) {
                groundCollision = true;
            }
            
            // SPHERE
            //  nc = (p[i] - center); nc.normalize();
            //  qc = radius*nc + center; // position - center
            //  diff = p[i] - qc;
            //  Cp = diff.dot(nc); // dot product
            //
            //  if (Cp < 0) {
            //      sphereCollision = true;
            //  }
            
            // w != 0 is so dont create a constraint between all restricted points
            if (groundCollision || sphereCollision) {
                collisions.push_back(new CollisionConstraint(i, groundCollision, sphereCollision, w[i]));
            }
            
            groundCollision = false;
            sphereCollision = false;
        }
    }
    
    void updateVelocitiesOfCollisions(vector<Vec3f> &v) {
        for (CollisionConstraint* c : collisions) {
            Vec3f nc = c->nc;
            Vec3f velo = v[c->index];
            
            // Reflect velocity in the direction of the collision normal.
            // [2 * k-restitution * lambert's law * normal]
            Real restitution = 0.01;
            velo -= 2 * restitution * (velo.dot(nc)) * nc;
            Vec3f friction = -(velo - velo.dot(nc) * nc);
            
            // frictional update
            v[c->index] += friction;
        }
    }
    
    void updateRigidBodies(vector<Vec3f> &p, vector<Vec3f> &x) {
        //        for (CollisionConstraint* c : collisions) {
        //            vec3f deltaP = p[c->index] - x[c->index];
        //            Real deltaT = 1e-3;
        //            deltaP /= (c->w * deltaT);
        //            // accel = (p[c->index] - x[c->index])/ (c->w * deltaT);
        //            // v = v0 + at
        //            // centerVelocity = centerVelocity + accel * deltaT;
        //            // shortens to: deltaP * deltaT / c->w
        //            // all weights are 1 or zero so fine
        //            centerVelocity += (p[c->index] - x[c->index]) * deltaT;
        //            //cout<<x[c->index]<<" << x, p >>"<<p[c->index]<<" deltaP: "<<deltaP<<endl;
        //            cout<<centerVelocity<<endl;
        //        } // impulse for sphere position
    }
    /*
    void updateEdges(map<int, vector<int> >& broken, vector<Vec3f> *x, vector<Vec3f> *v, vector<Real> &w, vector<Vec3f> &p, vector<Face*>& faces) {
        
        if (!(broken.size() > 0)) { return; }
        
        std::map<int, vector<int> >::iterator it = broken.begin();
        vector<int> edgePoints;
        int key;
        
        map<int, tuple<Vec3f, Vec3f>> newParticles = map<int, tuple<vec3f, vec3f>>();
        vector<Constraint*> addingConstraints = vector<Constraint*>();
        
        int numPositions = sim.numParticles;
        
        /// loop over all the key values
        while( it != broken.end()) {
            key = it->first;
            edgePoints = it->second;
            
            /// loop over each connected edgePoint to each key value
            for (int i = 0; i < (int) edgePoints.size(); ++i) {
                /// create edges DE, DF, and CE, CF [using indices]
                int d = key; // to match the notation in my notes
                int c = edgePoints[i];
                int e = numPositions;
                int f = numPositions + 1;
                numPositions += 2;
                
                /// store new particles x and v components
                // find midpoint of removed edge
                Vec3f midpoint = (x[d] + x[c]) / 2;
                Vec3f aveVelocity = (v[d] + v[c]) / 2;
                // create points E and F at midpoint with aveVelocity
                newParticles[e] = make_tuple(midpoint, aveVelocity);
                newParticles[f] = make_tuple(midpoint, aveVelocity);
                
                /// create and find indices for computations
                int indexFace0 = -1;
                int indexFace1 = -1;
                for (int k = 0; k < (int) faces.size() && indexFace0 == -1; ++k) {
                    if (faces[k]->containsIndices(c, d)) {
                        // on one of the correct faces
                        indexFace0 = k;
                    }
                }
                if (indexFace0 == -1) { throw; }
                Face *face0 = faces[indexFace0];
                indexFace1 = face0->getAdjacentFaceIndexByIndices(c, d);
                if (indexFace1 == -1) { throw; }
                if (indexFace0 == indexFace1) { throw; }
                Face *face1 = faces[indexFace1];
                
                // must be triangular faces for our implementation
                if (face0->indices.size() != 3 || face1->indices.size() != 3) { throw; }
                
                /// fill in a and b variables
                // a is the non-face index on face 1; b is the non-face index on face 0
                int a = -1;
                int b = -1;
                for (int index = 0; index < 3; ++index) {
                    if (face1->indices[index] != d && face1->indices[index] != c) {
                        a = face1->indices[index];
                    }
                    if (face0->indices[index] != d && face0->indices[index] != c) {
                        b = face0->indices[index];
                    }
                }
                if (a == -1 || b == -1) { throw; }
                
                /// check swap for proper indexing based on my faces framework
                // for proper triangulation for our notation face 0 goes c->d, face 1 goes d->c
                if (faces[indexFace1]->containsOrderedIndices(c, d)) {
                    if (testing && faces[indexFace0]->containsOrderedIndices(c, d)) { throw; }
                    
                    // face0 and face1 were labeled backwards for our indexing scheme so swap them
                    swap(indexFace0, indexFace1);
                }

                /// update faces indices in right order - CEB[update 0], AFC[update 1], ADF[new 2], BED[new 3]
                // dont update face0 and face1 indices yet
                Face *newFace2 = new Face({a, d, f}, faces.size());
                faces.push_back(newFace2);
                Face *newFace3 = new Face({b, e, d}, faces.size());
                faces.push_back(newFace3);
                
                /// update new faces' attached faces components for the one prev shared by face0 or face1 but not
                /// involved in other calculations
                // 2 takes 1's old face across DA edge connection [face 4]
                // 3 takes 0's old face across BD edge connection [face 5]
                int indexFace4 = face1->getAdjacentFaceIndexByIndices(d, a);
                int indexFace5 = face0->getAdjacentFaceIndexByIndices(b, d);
                if (indexFace5 == -1 || indexFace4 == -1) { throw; }
                Face* face4 = faces[indexFace4];
                Face* face5 = faces[indexFace5];
                
                /// delete face0 and face1's old attached faces across BD and DA edges
                bool adjacentsRemovedCorrectly = true;
                // delete 0's old face across BD edge connection [face 5]
                // delete 5's shared face connection with 0
                // delete 1's old face across DA edge connection [face 4]
                // delete 4's shared face connection with 1
                adjacentsRemovedCorrectly &= face0->removeAdjacentFaceById(indexFace5);
                adjacentsRemovedCorrectly &= face5->removeAdjacentFaceById(indexFace0);
                adjacentsRemovedCorrectly &= face1->removeAdjacentFaceById(indexFace4);
                adjacentsRemovedCorrectly &= face4->removeAdjacentFaceById(indexFace1); // TODO: THIS PART MIGHT CAUSE ISSUES
                if (!adjacentsRemovedCorrectly) { throw; }
                
                /// update face0 and face1 have proper position index values
                // resetting for simplicity [dont need to search for changed index]
                face0->resetIndices(c, e, b);
                face1->resetIndices(a, f, c);
                // 0 and 1 SHOULD STILL BE adjacents
                if (!(face0->adjacentToFace(face1))) { throw; }
                
                /// make sure all other faces in this arrangement [including newly created ones]
                /// are properly adjacents
                bool adjacentsAddedCorrectly = true;
                // new adjacents = 2 3; 0 3; 2 1 [currently connecting these adjacent faces]
                adjacentsAddedCorrectly &= newFace2->shouldBeAdjacentToFace(newFace3);
                adjacentsAddedCorrectly &= face0->shouldBeAdjacentToFace(newFace3);
                adjacentsAddedCorrectly &= face1->shouldBeAdjacentToFace(newFace2);
                if (!adjacentsAddedCorrectly) { throw; }
                // TODO: ADD IN CASE FOR IF FACES G AND H DONT EXIST
                /// update mesh face bending constraints
                // remove bend between 05; 14; - done in loop outside of this iter loop
                // add bend between 0 3; 1 2; 3 5; 2 4;
                // face indices is in order of notEdge, edge, edge, notEdge
                Real phi_03 = 0; Real phi_12 = 0; Real phi_35 = 0; Real phi_24 = 0;
                int g = face4->getNonListedIndex(d, a);
                int h = face5->getNonListedIndex(b, d);
                
                // fill in needed weight values and point values for constraint creations
                // if point is created before this edge method then grab it, if not, then set to 1
                // num - 1 bc indexing
                // int numParticles = w.size()
                Real wa = (a > sim.numParticles - 1) ? 1 : w[a];
                Real wb = (b > sim.numParticles - 1) ? 1 : w[b];
                Real wc = (c > sim.numParticles - 1) ? 1 : w[c];
                Real wd = (d > sim.numParticles - 1) ? 1 : w[d];
                Real we = 1; //new edge for this iteration
                Real wf = 1; //new edge for this iteration
                Real wg = (g > sim.numParticles - 1) ? 1 : w[g];
                Real wh = (h > sim.numParticles - 1) ? 1 : w[h];
                Vec3f pa = (a > sim.numParticles - 1) ? get<0>(newParticles[a]) : p[a];
                Vec3f pb = (b > sim.numParticles - 1) ? get<0>(newParticles[b]) : p[b];
                Vec3f pc = (c > sim.numParticles - 1) ? get<0>(newParticles[c]) : p[c];
                Vec3f pd = (d > sim.numParticles - 1) ? get<0>(newParticles[d]) : p[d];
                Vec3f pe = get<0>(newParticles[e]); //new edge for this iteration
                Vec3f pf = get<0>(newParticles[f]); //new edge for this iteration
                Vec3f pg = (g > sim.numParticles - 1) ? get<0>(newParticles[g]) : p[g];
                Vec3f ph = (h > sim.numParticles - 1) ? get<0>(newParticles[h]) : p[h];
                
                // facei, facej, faceIndices, phi to be filled in
                calculations::mathForFaceBendingConstraintCreation(face0, newFace3, {c, e, b, d}, {pc, pe, pb, pd}, phi_03);
                calculations::mathForFaceBendingConstraintCreation(face1, newFace2, {d, a, f, c}, {pd, pa, pf, pc}, phi_12);
                calculations::mathForFaceBendingConstraintCreation(newFace3, face5, {e, b, d, h}, {pe, pb, pd, ph}, phi_35);
                calculations::mathForFaceBendingConstraintCreation(newFace2, face4, {f, d, a, g}, {pf, pd, pa, pg}, phi_24);
                
                // inputFacesIndices, inputW, inputPhi
                addingConstraints.push_back(new FaceBendingConstraint({c, e, b, d}, {wc, we, wb, wd}, phi_03, true));
                addingConstraints.push_back(new FaceBendingConstraint({d, a, f, c}, {wd, wa, wf, wc}, phi_12, true));
                addingConstraints.push_back(new FaceBendingConstraint({e, b, d, h}, {we, wb, wd, wh}, phi_35, true));
                addingConstraints.push_back(new FaceBendingConstraint({f, d, a, g}, {wf, wd, wa, wg}, phi_24, true));
                
                /// update mesh stretch constraints
                // remove stretch from d,c - done in loop outside of this iter loop
                // add stretch to be; ce; cf; af; de; df;
                stretchBendVolume.push_back(new StretchConstraint((pb - pe).norm(), wb, we, b, e, true) );
                stretchBendVolume.push_back(new StretchConstraint((pc - pe).norm(), wc, we, c, e, true) );
                stretchBendVolume.push_back(new StretchConstraint((pc - pf).norm(), wc, wf, c, f, true) );
                stretchBendVolume.push_back(new StretchConstraint((pa - pf).norm(), wa, wf, a, f, true) );
                stretchBendVolume.push_back(new StretchConstraint((pd - pe).norm(), wd, we, d, e, true) );
                stretchBendVolume.push_back(new StretchConstraint((pd - pf).norm(), wd, wf, d, f, true) );
                
                // TODO FOR LATER : UPDATE VOLUME CONSTRAINTS
                
            } // end: going through all attached vertices of this key
            
            /// step the iterator
            ++it;
        }
        
        /// remove all the stretchConstraints and faceBendingConstraints that should be for this step
        stretchbend.erase(remove_if(stretchbend.begin(), stretchbend.end(),
                                          [](Constraint* c) { return c->broken(); }), stretchbend.end()); // TODO: DOES THIS ACTUALLY DELETE THE CONSTRAINT* AND ALL INSTANCES?
        
        /// add in all the new stretchConstraints and faceBendingConstraints
        stretchbend.insert(end(stretchbend), begin(addingConstraints), end(addingConstraints));
        
        /// reclear the edges that were broken since they were fixed
        broken.clear();
        
        /// calculate updates for all inputs x, v, p, w, faces for new added positions [includes resizing]
        // copy over x into newX
        // copy over v into newV
        vector<vec3f> newX;
        vector<vec3f> newV;
        for (int i = 0; i < sim.numParticles; ++i) {
            newX.push_back(x[i]);
            newV.push_back(v[i]);
        }
        map<int, tuple<vec3f, vec3f>>::iterator newParticlesIter = newParticles.begin();
        while (newParticlesIter != newParticles.end()) {
            newX[newParticlesIter->first] = get<0>(newParticlesIter->second);
            newV[newParticlesIter->first] = get<1>(newParticlesIter->second);
        }
        sim.numParticles = numPositions;
        
        /// actually updating the inputs
        // x->resize(sim.numParticles, dim);
        *x = newX;
        //v->resize(sim.numParticles, dim);
        *v = newV;
        
        // p.resize(sim.numParticles, dim);
        // below line is the same as mesh.calculate weights
        // w.resize(sim.numParticles, sim.numParticles);
        for (int i = 0; i < (int) sim.staticParticles.size(); ++i) {
            w[sim.staticParticles[i]] = 0;
        }
    }
     */


};


#endif  /* _CALCULATIONS_HPP_ */
