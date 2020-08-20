

import math
from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGL.GLU import *
from OpenGL.GLUT import *
#import OpenGL.GLUT as glut
import numpy as ny

gPause = False
gSaveFile = False
_step = 0
_sim_t = 0
dt = 0.01
# pos_curr
# pos_an

def ReSizeGLScene(Width, Height): 
    glViewport(0, 0, Width, Height)        
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(70.0, float(Width)/float(Height), 1.0, 11.0)
    glMatrixMode(GL_MODELVIEW)


def draw():
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    #Render the vertices
    glBegin(GL_LINES)
    # for
    GLfloat point[] = {pos[index_list[i][j]][0],pos[index_list[i][j]][1],pos[index_list[i][j]][2]}
    glVertex3fv(point)
    glEnd()
    glPopMatrix()



def render():
    glClearColor(.4f, .4f, .4f, 1.0)
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()
    gluLookAt(0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    glColor3f(1.0, 1.0, 1.0)
    draw()
    glutSwapBuffers()
    if(gSaveFile) {
        ##
    }


def step():
    if (_step>10001): return
    # pos_curr
    # pos_an
    _step += 1
    _sim_t += dt



def update():
  if (gPause==False): step()
  glutPostRedisplay()



def KeyboardFunc(key, x, y):
    if k in ('p', 'P'): 
        gPause = True
    if k in ('s', 'S'):
        gSaveFile = True


def main():
    global window
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
    glutInitWindowSize(1200, 900)
    glutInitWindowPosition(800,400)
    window = glutCreateWindow("Rigid Body Solver")
    glutDisplayFunc(render)
    glutReshapeFunc(ReSizeGLScene)
    glutKeyboardFunc(KeyboardFunc)
    glutIdleFunc(update)
    InitGL(640, 480)
    glutMainLoop()

main()