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

class Box : public BodyAttributes {
public:
  explicit Box(
    const Real w=1.0, const Real h=1.0, const Real d=1.0,
    const Vec3f v0=Vec3f(0, 0, 0), const Vec3f omega0=Vec3f(0, 0, 0)) :
    width(w), height(h), depth(d) {

    // TODO:
    M = 100.0;
    I0 = Mat3f::I();
    I0(0, 0) = M * (height * height + depth * depth) / 12;
    I0(1, 1) = M * (width * width + depth * depth) / 12;
    I0(2, 2) = M * (height * height + width * width) / 12;
    I0inv = I0.inverse();
    Iinv = I0inv;

    // vertices data (8 vertices)
    vdata0.push_back(Vec3f(-0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h, -0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h, -0.5*d));

    vdata0.push_back(Vec3f(-0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h,  0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h,  0.5*d));
  }

  virtual void renderGl() const {
    glPushMatrix();

    float mat[] = {              // OpenGL is column major
      R(0,0), R(1,0), R(2,0), 0, // column 0
      R(0,1), R(1,1), R(2,1), 0, // column 1
      R(0,2), R(1,2), R(2,2), 0, // column 2
      X[0],   X[1],   X[2],   1  // column 3
    };
    glMultMatrixf(mat);         // M

    // Render the vertices
    glEnableClientState(GL_VERTEX_ARRAY);
    glPointSize(5);
    glVertexPointer(3, GL_FLOAT, sizeof(Vec3f), &vdata0[0]);
    glDrawArrays(GL_POINTS, 0, vdata0.size());
    glDisableClientState(GL_VERTEX_ARRAY);

    glScalef(width, height, depth); // MS
    glutWireCube(1.0);              // MSx

    glPopMatrix();
  }

  // rigid body property
  Real width, height, depth;
};

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

class RigidSolver {
public:
  explicit RigidSolver(
    String *body0,   // BodyAttributes *body0,
    const Real dt=0.01, const Vec3f g=Vec3f(0, -9.8, 0)) :
    body(body0), dt(dt), _g(g), _step(0), _sim_t(0) {
  }

  void initScene(String *body0) { // BodyAttributes *body0
    body = body0;
    _step = 0;
    _sim_t = 0;
  }

  void step() {
      if (_step>1001) return;
    // if (_step % 10==0) gSaveFile = true;
    // else gSaveFile = false;
    std::cout << "t=" << _sim_t << " (dt=" << dt << ")" << std::endl;
      
    /// (5) for all vertices - update velocity by doing vi = vi + deltat*wi*fext(xi)
    // for string we only consider the gravity
      for (int i=1;i<body->n_points;i++) {
          body->V_POINT[i] += _g * dt;
          std::cout << "velocity " << body->V_POINT[i] << std::endl;
      }
      
    /// (6) dampVelocites(v1,...vN)
    if (_step>0) calculations::dampVelocities(body->pos, body->V_POINT, body->M_POINT);
    else calculations::dampVelocities_simple(body->V_POINT);
      
    for (int i=1;i<body->n_points;i++) {
        // body->pos[i] += body->V_POINT[i] * dt;
        std::cout << "pos " << body->pos[i] << std::endl;
    }
    
    /// (7) for all verticies i find projected point assuming no collisions pi = xi + deltat vi
    for (int i=1;i<body->n_points;i++) { // the first point will not move
        body->projected[i] = body->pos[i] + body->V_POINT[i] * dt;
    }
    
    /// (8) for all velocities generate collision constraints (xi -> pi), here stretch constraint already created
    // allConstraints.createCollisionConstraints(p, w);

    /// (9) loop the solver the number of desired iterations projecting the constraints
    for (int i = 0; i < body->allConstraints.constraintIterations; i++) {
        body->allConstraints.update(body->projected);// , brokenEdges); //!!!!
    }
      
    /// (12) for all vertices update vi and xi for next overall loop simulation step
    for (int i=1;i<body->n_points;i++) {
        body->V_POINT[i] = (body->projected[i] - body->pos[i]) / dt;
    }
    body->pos = body->projected;

    /// (16) velocities update (for friction etc - only on particles involved in this iterations collisions)


    // write to file
    // writeToFile(_sim_t, dt, body->n_points);
    writeToFile(body->pos, body->n_points);
      
    // computeForceAndTorque();

    // TODO:
    /* // simple Box demo
    body->V += _g * dt + body->F / body->M * dt;
    
    body->X += body->V * dt;
    
    body->omega += body->Iinv * body->tau * dt;
    
    body->R = computeR(body->omega);
     */
      
    ++_step;
    _sim_t += dt;
  }

  String *body; // BodyAttributes *body;
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

Box abox(2.0, 1.0, 1.0);
String astring(10);
RigidSolver solver(&astring);

// simple rendering
void render() {
  glClearColor(.4f, .4f, .4f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glLoadIdentity();
  gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

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
