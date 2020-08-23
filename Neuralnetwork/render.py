import tensorflow as tf
import math
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import numpy as np
import test2 as te
import make_my_data as md
import time
from PIL import Image
from PIL import ImageOps
import sys

gPause = False
gSaveFile = False
_step = 0
n_particle = 10
n_pic = 0

(X, Y, pos_int) = md.import_data()

(n_step, n_ver) = X.shape
n_test = int(n_step/5)

X_TRAIN = X[0:n_step-n_test, :]
Y_TRAIN = Y[0:n_step-n_test, :]
X_TEST = X[n_step-n_test:n_step, :]
Y_TEST = Y[n_step-n_test:n_step, :]

pos_cur = pos_int



#model = te.createModel(X_TRAIN, Y_TRAIN, X_TEST, Y_TEST, n_ver, True)
model = te.createModel(X_TRAIN, Y_TRAIN, X_TEST, Y_TEST, n_ver, False)

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
    global pos_cur
    for i in range(n_particle-1):
        for j in range(2):
            point = [pos_cur[i*3+j*3], pos_cur[i*3+j*3+1], pos_cur[i*3+j*3+2]]
            glVertex3fv(point)
    glEnd()
    glPopMatrix()



def render():
    glClearColor(0.4, 0.4, 0.4, 1.0)
    glClear(GL_COLOR_BUFFER_BIT)
    glLoadIdentity()
    gluLookAt(0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    glColor3f(1.0, 1.0, 1.0)
    draw()
    glutSwapBuffers()
    global gSaveFile
    if gSaveFile:
        global n_pic
        print("save file ", n_pic)
        
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        data = glReadPixels(0, 0, 1200, 900, GL_RGBA, GL_UNSIGNED_BYTE)
        image = Image.frombytes("RGBA", (1200, 900), data)
        image = ImageOps.flip(image)
        image.save("glutout"+n_pic + ".png", 'PNG')
        n_pic = n_pic + 1
        gSaveFile = False
        return


def step():
    global _step, pos_int, pos_cur
    if (_step>10000): 
        return

    pos_cur = model.predict(pos_int.reshape((1, 30))).reshape((30,))
    pos_int = pos_cur
    _step = _step + 1
    time.sleep(0.01)


def update():
    global gPause
    if (gPause==False): step()
    glutPostRedisplay()


def KeyboardFunc(key, x, y):
    global gSaveFile, gPause
    if key=='p': 
        print("========== key p ===========")
        gPause = not gPause
    if key=='s':
        gSaveFile = True


def main():
    global window
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
    glutInitWindowSize(1200, 900)
    glutInitWindowPosition(800,400)
    window = glutCreateWindow("Rigid Body Solver")
    time.sleep(1)
    glutDisplayFunc(render)
    glutReshapeFunc(ReSizeGLScene)
    glutKeyboardFunc(KeyboardFunc)
    glutIdleFunc(update)

    glutMainLoop()

main()
