#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import random
import WaveletCavity
import TriangleCavity
import math
from PyQt4 import Qt, QtGui, QtCore
from PyQt4.QtOpenGL import *
from OpenGL.GL import *

class CavityViewerMain(QtGui.QMainWindow):
    def __init__(self, glWindow):
        QtGui.QMainWindow.__init__(self)

        self.resize(800,600)
        self.setCentralWidget(glWindow)
        self.cavityView = glWindow

        menu = self.menuBar()

        fileMenu = menu.addMenu('&File')
        screenShotItem = fileMenu.addAction('&Take a screenshot')
        self.connect(screenShotItem, QtCore.SIGNAL('triggered()'), self.screenshot)
        fileMenu.addSeparator()
        quitItem = fileMenu.addAction('&Quit')
        self.connect(quitItem, QtCore.SIGNAL('triggered()'), self.quit)

        viewMenu = menu.addMenu('&View')

        self.cavityView.cavity.addMenuActions(viewMenu)

        viewMenu.addSeparator()

        self.useLightingItem = viewMenu.addAction('Use &lighting')
        self.useLightingItem.setCheckable(True)
        self.useLightingItem.setChecked(self.cavityView.useLighting)
        self.connect(self.useLightingItem, QtCore.SIGNAL('triggered()'), self.changeLighting)    

        self.orthogonalProjectionItem = viewMenu.addAction('&Orthogonal projection')
        self.orthogonalProjectionItem.setCheckable(True)
        self.orthogonalProjectionItem.setChecked(self.cavityView.ortho)
        self.connect(self.orthogonalProjectionItem, QtCore.SIGNAL('triggered()'),
                     self.changeProjection)

        viewMenu.addSeparator()

        self.resetToOriginItem = viewMenu.addAction('&Reset to origin')
        self.connect(self.resetToOriginItem, QtCore.SIGNAL('triggered()'), self.resetPosition)

        return

    def quit(self):
        sys.exit()
        return

    def screenshot(self):
        # Supported image formats
        imageFormats = ''
        for format in Qt.QImageWriter.supportedImageFormats():
            imageFormats += '*.' + str(format) + ' '
        imageFormats = 'Images ('+imageFormats[0:-1]+')'

        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save screenshot as', 'screenshot.png', imageFormats)
        if filename.isEmpty():
            return

        image = self.cavityView.grabFrameBuffer(True)
        image.save(filename)
        return

    def changeLighting(self):
        self.cavityView.useLighting = self.useLightingItem.isChecked()
        self.cavityView.updateGL()
        return

    def changeProjection(self):
        self.cavityView.ortho = self.orthogonalProjectionItem.isChecked()
        self.cavityView.resizeGL(self.cavityView.width(), self.cavityView.height())
        self.cavityView.updateGL()
        return
    
    def resetPosition(self):
        self.cavityView.resetToOrigin()
        self.cavityView.updateGL()
        return



class CavityView(QGLWidget):
    def __init__(self, cavity, parent = None):
        QGLWidget.__init__(self, parent, None)

        self.setMouseTracking(False)
        self.leftButtonDown = False
        self.rightButtonDown = False
        self.middleButtonDown = False
        self.lastPos = None

        self.cavity = cavity

        self.resetToOrigin()

        self.useLighting = True

        # set to > 0
        self.tesselationCoarsity = 15
        self.useLists = True

        self.ortho = False

        return

    def resetToOrigin(self):
        self.factor = 3.0
        self.xpos = 0.0
        self.ypos = 0.0
        self.zpos = -10.0
        glLoadIdentity()
        self.rotation = glGetFloatv(GL_MODELVIEW_MATRIX)
        return


    ############### OpenGL stuff ###############


    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        if self.useLighting:
            pos = [0.0, 100.0, 100.0]
            glLightfv(GL_LIGHT0, GL_POSITION, pos)
            glEnable(GL_LIGHTING)
        else:
            glDisable(GL_LIGHTING)
 
        glTranslate(self.xpos, self.ypos, self.zpos)

        glMultMatrixf(self.rotation)
        
        if self.ortho:
            glScale(self.factor/10.0, self.factor/10.0, self.factor/10.0)
        else:
            glScale(self.factor, self.factor, self.factor)


        self.cavity.drawOpenGL(self.useLighting)

        glFlush()
        return    

    def resizeGL(self, width, height):
        self.h = height
        self.w = width

        glViewport(0, 0, width, height) 
        glMatrixMode(GL_PROJECTION)    
        glLoadIdentity();
        aspectRatio = 1.0*width / height

        side = 0.5

        if width <= height:
            self.xScale = 1.0/width
            self.yScale = aspectRatio/height
        else:
            self.xScale = aspectRatio/width
            self.yScale = 1.0/height

        if not self.ortho:
            if width <= height: 
                glFrustum(-side, side, -side/aspectRatio, side/aspectRatio, 1.0, 20.0)
            else:        
                glFrustum(-side*aspectRatio, side*aspectRatio, -side, side, 1.0, 20.0)
        else:
            if width <= height:
                glOrtho(-side, side, -side/aspectRatio, side/aspectRatio, 1.0, 100.0)
            else:
                glOrtho(-side*aspectRatio, side*aspectRatio, -side, side, 1.0, 100.0)
                
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        return


    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glDepthFunc(GL_LESS)
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LINE_SMOOTH)
        glShadeModel(GL_SMOOTH)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        glClearDepth(1.0)
        glEnable(GL_RESCALE_NORMAL)

        # Light definitions
        a = 0.2
        ambient = [a, 0.0, a, 1.0]
        d = 1.0
        diff = [d, d, d, 1.0]
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diff)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        
        self.cavity.initializeOpenGL()
                
        return


    ############### Event handlers ###############


    def mousePressEvent(self, event):
        if event.button() == 1:
            self.leftButtonDown = True
            self.rightButtonDown = False
            self.middleButtonDown = False
            self.lastPos = event.pos()
            # Copy the stored rotation matrix
            glLoadMatrixf(self.rotation)
            self.lastRotation = glGetFloatv(GL_MODELVIEW_MATRIX)
            return
        if event.button() == 2:
            self.leftButtonDown = False
            self.rightButtonDown = True
            self.middleButtonDown = False
            self.lastPos = event.pos()
            return
        if event.button() == 4:
            self.leftButtonDown = False
            self.rightButtonDown = False
            self.middleButtonDown = True
            self.lastPos = event.pos()
        return

    def mouseReleaseEvent(self,event):
        if event.button() == 1:
            self.leftButtonDown = False
        if event.button() == 2:
            self.rightButtonDown = False
        if event.button() == 4:
            self.middleButtonDown = False
        return

    def mouseMoveEvent(self, event):
        diff = event.pos() - self.lastPos

        if self.leftButtonDown:

            # Construct rotation matrix by taking the plane of mouse draw and 
            # applying a rotation with an angle that depends on the draw length
            
            diff = [diff.x(), diff.y()]
            diffLen = math.sqrt(sum([d*d for d in diff]))
            if diffLen < 1e-1:
                return

            diff = [d/diffLen for d in diff]
            ax =  diff[1]
            ay =  diff[0]

            glLoadIdentity()
            glRotate(0.25 * diffLen, ax, ay, 0.0)
            glMultMatrixf(self.lastRotation)

            self.rotation = glGetFloatv(GL_MODELVIEW_MATRIX)
            # don't update the last position

        if self.middleButtonDown:
            self.factor += 0.1 * -diff.y()
            if self.factor < 1.0:
                self.factor = 1.0
            self.lastPos = event.pos()
        

        if self.rightButtonDown:
            if self.ortho:
                self.xpos += self.xScale * diff.x()
                self.ypos -= self.yScale * diff.y()
            else:
                self.xpos += -1.0 * self.zpos * self.xScale * diff.x()
                self.ypos -= -1.0 * self.zpos * self.yScale * diff.y()
            self.lastPos = event.pos()

        self.updateGL()
        return

    def wheelEvent(self, event):
        # One wheel 'click' is 120.0 units in most cases...
        self.factor += 0.5 * event.delta()/120.0
        if self.factor < 1.0:
            self.factor = 1.0
        self.updateGL()
        return


# The main program
        

def getCavity(filename):
    try:
        infile = open(filename)
        infile.close()
    except:
        print "The file '" + filename + "' does not exist."
        sys.exit()

    found = False

    try:
        cavity = WaveletCavity.WaveletCavity(filename)
        found = True
    except:
        pass
    if found == True:
        return cavity

    # Wasn't that type, try something else...
    try:
        cavity = TriangleCavity.TriangleCavity(filename)
        found = True
    except:
        pass
    if found == True:
        return cavity

    print "Cavity file is of unknown type!"
    sys.exit()
    return


if __name__=='__main__':

    app = QtGui.QApplication(sys.argv)

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'molec_dyadic.dat'

    cavity = getCavity(filename)

    glWidget = CavityView(cavity)
    mainApp = CavityViewerMain(glWidget)
    mainApp.show()
    sys.exit(app.exec_())
