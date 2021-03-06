// ----------------------------------------------------------------------------
// rigid.cpp
//
//  Created on: 13 Feb 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: main codes (DO NOT DISTRIBUTE!)
//
// Copyright 2020 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>
#include <map>

#ifdef __APPLE__
#include <OpenGL/gl.h> //OS x libs
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
// #include <GLUT/freeglut.h>
#else
#include <GL/glut.h>
#include <GL/freeglut.h>
#endif

#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "Calculations.hpp"
#include "Constraints.hpp"
#include "FileIO.hpp"

bool gPause = false;
bool gSaveFile = false;
int gGlutTime;                  // glut elapsed time
int gSavedCnt = 0;


struct BodyAttributes {
  BodyAttributes() :
    X(0, 0, 0), R(Mat3f::I()), P(0, 0, 0), L(0, 0, 0),
    V(0, 0, 0), omega(0, 0, 0), F(0, 0, 0), tau(0, 0, 0) {}

  Vec3f worldCoordOf(const tIndex vidx) const {
    // TODO:
      return (X + R * vdata0[vidx]);
  }

  virtual void renderGl() const = 0;

  Real M;                       // mass
  Mat3f I0, I0inv;              // inertia tensor and its inverse in body space

  // rigid body state
  Vec3f X;                      // position
  Mat3f R;                      // rotation
  Vec3f P;                      // linear momentum
  Vec3f L;                      // angular momentum

  // auxiliary quantities
  Mat3f Iinv;                   // inverse of inertia tensor
  Vec3f V;                      // linear velocity
  Vec3f omega;                  // angular velocity

  // force and torque
  Vec3f F;                      // force
  Vec3f tau;                    // torque

  // mesh's vertices in body space
  std::vector<Vec3f> vdata0;
};


static const GLint index_list[][2] =
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

class Box: public BodyAttributes {
public:
    int n_points = 8;
    Real width, height, depth; // x, y, z
    std::vector<Vec3f> V_POINT;
    std::vector<Real> M_POINT;
    std::vector<Vec3f> pos;
    std::vector<Vec3f> projected;
    std::vector<Face*> faces;
    Constraints allConstraints;
    map<int, vector<int> > brokenEdges;
   
    explicit Box(Real w, Real he, Real de):  width(w), height(he), depth(de) {
        // vertices data (8 vertices)
        sim.numParticles = n_points;
        pos.push_back(Vec3f(-0.5*w, -0.5*he, -0.5*de));
        pos.push_back(Vec3f( 0.5*w, -0.5*he, -0.5*de));
        pos.push_back(Vec3f( 0.5*w,  0.5*he, -0.5*de));
        pos.push_back(Vec3f(-0.5*w,  0.5*he, -0.5*de));
        pos.push_back(Vec3f(-0.5*w, -0.5*he,  0.5*de));
        pos.push_back(Vec3f( 0.5*w, -0.5*he,  0.5*de));
        pos.push_back(Vec3f( 0.5*w,  0.5*he,  0.5*de));
        pos.push_back(Vec3f(-0.5*w,  0.5*he,  0.5*de));
        for (int i=0;i<n_points;i++) {
            V_POINT.push_back(Vec3f(0, 0, 0));
            M_POINT.push_back(1.0);
            projected.push_back(Vec3f(0, 0, 0));
        }
        M_POINT[0] = 0.0; // fixed
        int f[12][3] = {
            {0,1,2}, {0,2,3},
            {0,3,7}, {0,7,4},
            {4,5,6}, {4,6,7},
            {1,2,6}, {1,6,5},
            {0,1,5}, {0,5,4},
            {2,3,7}, {2,7,6}
        };
        //通过数组a的地址初始化，注意地址是从0到3（左闭右开区间）
        for (int i=0; i<12; i++) {
            vector<int> b = vector<int>();
            for (int j=0; j<3; j++) {
                b.push_back(f[i][j]);
            }
            faces.push_back(new Face(b, faces.size()) );
        }
        calculateAdjacentFaces();
        // [removes the issue of checking based on duplicate edgings in the map - i,j and j,i]
        brokenEdges = map<int, vector<int> >();
        bool allowedToBreak = false;
        /// setup overall constraints
        allConstraints = Constraints();
        allConstraints.createStretchConstraints(faces, pos, M_POINT, allowedToBreak);
        cout<<"finished stretch constraints"<<endl;
        // allConstraints.createFaceBendingConstraints(faces, pos, M_POINT, allowedToBreak);
        cout<<"finished facebend constraints"<<endl;
        // writeToFile(n_points, allowedToBreak);
        writeToFile(pos, n_points);
    }
    
    void calculateAdjacentFaces(){
        vector<vector<int> > check(faces.size(), vector<int>(faces.size()));
        for(int i=0; i<faces.size(); i++) {
            for(int j=0; j<faces.size(); j++){
                check[i][j] = 0;
            }
        }
        for (int i = 0; i < faces.size(); ++i) {
            for (int j = 0; j < faces.size(); ++j) {
                if (i != j && check[i][j] != 1) {
                    check[i][j] = 1;
                    check[j][i] = 1;
                    if ( faces[i]->shouldBeAdjacentToFace(faces[j]) ) {
                        faces[j]->shouldBeAdjacentToFace(faces[i]);
                    }
                }
            }
        }
    }
                                  
    virtual void renderGl() const {
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();

        // Render the vertices
        int i,j;
        glBegin(GL_LINES);
        for(i=0; i<12; ++i){
            for(j=0; j<2; ++j){
                GLfloat point[] = {pos[index_list[i][j]][0],pos[index_list[i][j]][1],pos[index_list[i][j]][2]};
                //cout<< i << " " << j << " " << pos[index_list[i][j]][0]<<endl;
                glVertex3fv(point);
            }
        }
        glEnd();
    
        glPopMatrix();
    }
    
};
/*
class String: public BodyAttributes {
public:
    explicit String(int n): n_points(n) {
        int i;
        for (i=0;i<n_points;i++) {
            V_POINT.push_back(Vec3f(0, 0, 0));
            pos.push_back(Vec3f(i*0.2, 0, 0));
            M_POINT.push_back(1.0);
            projected.push_back(Vec3f(0, 0, 0));
        }
        // [removes the issue of checking based on duplicate edgings in the map - i,j and j,i]
        // map<int, vector<int>> brokenEdges = map<int, vector<int>>();
        bool allowedToBreak = false;
        /// setup overall constraints
        allConstraints = Constraints();
        allConstraints.createStretchConstraints(pos, M_POINT, allowedToBreak);
        cout<<"finished stretch constraints"<<endl;
        // writeToFile(n_points, allowedToBreak);
        writeToFile(pos, n_points);
    }
    
    virtual void renderGl() const {
        glPushMatrix();
        // Render the vertices
        glEnableClientState(GL_VERTEX_ARRAY);
        glPointSize(5);
        glVertexPointer(3, GL_FLOAT, sizeof(Vec3f), &pos[0]);
        glDrawArrays(GL_LINE_STRIP, 0, pos.size());
        glDisableClientState(GL_VERTEX_ARRAY);
        
        glPopMatrix();
    }
    
    int n_points;
    std::vector<Vec3f> V_POINT;
    std::vector<Real> M_POINT;
    std::vector<Vec3f> pos;
    std::vector<Vec3f> projected;
    Constraints allConstraints;
};
 
 */

class RigidSolver {
public:
  explicit RigidSolver(
    // String *body0,   // BodyAttributes *body0,
    Box *body0,
    const Real dt=0.01, const Vec3f g=Vec3f(0, -0.1, 0)) :
    body(body0), dt(dt), _g(g), _step(0), _sim_t(0) {
  }

  void initScene(Box *body0) { // BodyAttributes *body0
    body = body0;
    _step = 0;
    _sim_t = 0;
  }

  void step() {
      
      if (_step>10001) return;
    // if (_step % 10==0) gSaveFile = true;
    // else gSaveFile = false;
    std::cout << "t=" << _sim_t << " (dt=" << dt << ")" << std::endl;
      
    /// (5) for all vertices - update velocity by doing vi = vi + deltat*wi*fext(xi)
    // for string we only consider the gravity
      for (int i=1;i<body->n_points;i++) { // first point fixed
          body->V_POINT[i] += _g * dt;
          std::cout << "velocity " << body->V_POINT[i] << std::endl;
      }
      
    /// (6) dampVelocites(v1,...vN)
     
    if (_step>0) calculations::dampVelocities(body->pos, body->V_POINT, body->M_POINT);
    else calculations::dampVelocities_simple(body->V_POINT);
      
    for (int i=0;i<body->n_points;i++) {
        std::cout << "pos " << i << body->pos[i] << std::endl;
    }
    
    /// (7) for all verticies i find projected point assuming no collisions pi = xi + deltat vi
    for (int i=0;i<body->n_points;i++) { // the first point will not move
        body->projected[i] = body->pos[i] + body->V_POINT[i] * dt;
    }
    
    /// (8) for all velocities generate collision constraints (xi -> pi), here stretch constraint already created
    // allConstraints.createCollisionConstraints(p, w);

    /// (9) loop the solver the number of desired iterations projecting the constraints
    for (int i = 0; i < body->allConstraints.constraintIterations; i++) {
        body->allConstraints.update(body->projected, body->brokenEdges); //!!!!
    }
      
    /// (12) for all vertices update vi and xi for next overall loop simulation step
    for (int i=1;i<body->n_points;i++) { // first point
        body->V_POINT[i] = (body->projected[i] - body->pos[i]) / dt;
    }
    body->pos = body->projected;

    /// (16) velocities update (for friction etc - only on particles involved in this iterations collisions)
      // since collision constraints are the last to be updated in allConstraints update - dont need to reconstrain here
    // body->allConstraints.updateVelocitiesOfCollisions(body->V_POINT);
      
    // write to file
    // writeToFile(_sim_t, dt, body->n_points);
    writeToFile(body->pos, body->n_points);
     
    ++_step;
    _sim_t += dt;
      
  }

  Box *body; // BodyAttributes *body;
  Real dt;                      // simulation time and time step size

private:
  Mat3f computeR(Vec3f omega) {
      Mat3f res = Mat3f::zero();
      res(0, 0) = cos(omega[1]) * cos(omega[2]);
      res(0, 1) = -sin(omega[2]) * cos(omega[0]) + cos(omega[2]) * sin(omega[1]) * sin(omega[0]);
      res(0, 2) = sin(omega[0]) * sin(omega[2]) + cos(omega[2]) * sin(omega[1]) * cos(omega[0]);
      res(1, 0) = sin(omega[2]) * cos(omega[1]);
      res(1, 1) = cos(omega[2]) * cos(omega[0]) + sin(omega[2]) * sin(omega[1]) * sin(omega[0]);
      res(1, 2) = -sin(omega[0]) * cos(omega[2]) + sin(omega[2]) * sin(omega[1]) * cos(omega[0]);
      res(2, 0) = -sin(omega[1]);
      res(2, 1) = cos(omega[1]) * sin(omega[0]);
      res(2, 2) = cos(omega[1]) * cos(omega[0]);
    
      return res;
  }
  /*
  void computeForceAndTorque() {
    // TODO:
    // F
    body->F = Vec3f(10.0, 0, 0);
    // T
    Vec3f r = Vec3f(body->width/2, body->height/2, body->depth/2);
    body->tau = Vec3f(r[0]*body->F[1] - r[1] * body->F[0], r[1] * body->F[2] - r[2] * body->F[1], r[2] * body->F[0] - r[0] * body->F[2]);
  }
    */
    
  // simulation parameters
  Vec3f _g;                     // gravity
  tIndex _step;                 // step count
  Real _sim_t;                  // simulation time
};

Box abox(0.8, 1.2, 0.8);
// String astring(10);
RigidSolver solver(&abox);

// simple rendering
void render() {
  glClearColor(.4f, .4f, .4f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
    
  glLoadIdentity();
  gluLookAt(0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

  glColor3f(1.0, 1.0, 1.0);
  solver.body->renderGl();

  glutSwapBuffers();

  if(gSaveFile) {
    std::stringstream fpath;
    fpath << "s" << std::setw(4) << std::setfill('0') << gSavedCnt++ << ".tga";

    std::cout << "Saving file " << fpath.str() << " ... " << std::flush;
    const short int w = glutGet(GLUT_WINDOW_WIDTH);
    const short int h = glutGet(GLUT_WINDOW_HEIGHT);
    std::vector<int> buf(w*h*3, 0);
    glReadPixels(0, 0, w, h, GL_BGR, GL_UNSIGNED_BYTE, &(buf[0]));

    FILE *out = fopen(fpath.str().c_str(), "wb");
    short TGAhead[] = {0, 2, 0, 0, 0, 0, w, h, 24};
    fwrite(&TGAhead, sizeof(TGAhead), 1, out);
    fwrite(&(buf[0]), 3*w*h, 1, out);
    fclose(out);
    gSaveFile = false;

    std::cout << "Done" << std::endl;
  }
}

void update() {
  if(!gPause) solver.step();
  glutPostRedisplay();
}

void reshape(int w, int h) {
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(70.0, static_cast<GLfloat>(w)/h, 1.0, 11.0);
  glMatrixMode(GL_MODELVIEW);
}

void restart() {
  abox = Box(2.0, 1.0, 1.0);
  solver.initScene(solver.body);
}

void keyboardFunc(unsigned char key, int, int) {
  switch(key) {
  // case 'q': case 'Q': glutLeaveMainLoop(); break;
  case 'p': case 'P': gPause = !gPause; break;
  case 's': case 'S': gSaveFile = true; break;
  case 'r': case 'R': gPause = true; restart(); break;
  }
}

int main(int argc, char **argv) {
  glutInitWindowSize(1200, 900);
  glutInit(&argc, argv);
  glutInitDisplayString("samples stencil>=3 rgb double depth");
  glutCreateWindow("Rigid Body Solver");
  glutDisplayFunc(render);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboardFunc);
  glutIdleFunc(update);

  glutMainLoop();

  return EXIT_SUCCESS;
}
