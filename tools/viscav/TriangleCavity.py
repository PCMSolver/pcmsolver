""" Cavity composed of triangles """

import Cavity
from OpenGL.GL import *
import random

class TriangleCavity(Cavity.Cavity):
    def __init__(self, filename):
        """ File format:
        number of triangles
        triangle 0: vertex0[x,y,z]
        triangle 0: vertex1[x,y,z]
        triangle 0: vertex2[x,y,z]
        triangle 1: vertex0[x,y,z]
        ... 

        Vertices must be arranged so that the normal (v1-v0) x (v2-v0) points outside.
        """

        infile = open(filename, 'r')

        numTriangles = int(infile.readline().split()[0])
        self.mesh = []

        for i in xrange(numTriangles):
            verts = []
            for j in xrange(3):
                line = ''
                while(line.strip() == ''):
                    line = infile.readline()
                line = line.split()
                verts.append([float(line[0]), float(line[1]), float(line[2])])
            self.mesh.append(verts)

        self.scaleMesh()

        return

    def scaleMesh(self):
        import math
        """ Translates the center of mesh to origin and scales the mesh inside a unit ball. """
        
        # Center the mesh first
        center = [0.0, 0.0, 0.0]
        numPoints = 0
        for triangle in self.mesh:
            for point in triangle:
                for i in xrange(3):
                    center[i] += point[i]
                numPoints += 1
        center = [c/numPoints for c in center]

        # Translate and compute maximum radius (squared)
        maxRad = 0.0
        for triangle in self.mesh:
            for point in triangle:
                for i in xrange(3):
                    point[i] -= center[i]
                rad = sum([p*p for p in point])
                if rad > maxRad:
                    maxRad = rad
        scaleFactor = 1.0/math.sqrt(maxRad)

        # Scale
        for triangle in self.mesh:
            for point in triangle:
                for i in xrange(3):
                    point[i] *= scaleFactor

        return

    def initializeOpenGL(self):

        self.meshDL = glGenLists(1)
        glLoadIdentity()
        try:
            glNewList(self.meshDL, GL_COMPILE)
        except:
            pass

        for triangle in self.mesh:
            #normal = [0.0, 0.0, 0.0]
            #v1 = [a-b for a,b in zip(triangle[0], triangle[1])]
            #v2 = [a-b for a,b in zip(triangle[0], triangle[2])]
            #normal[0] = 
            glBegin(GL_TRIANGLES)
            glColor3f(random.random()*0.3+0.3, random.random()*0.3+0.3, random.random()*0.3+0.3)
            for vert in triangle:
                glVertex3f(vert[0], vert[1], vert[2])
            glEnd()

        try:
            glEndList()
        except:
            pass

        return

    def drawOpenGL(self, useLighting):
        try:
            glCallList(self.meshDL)
        except:
            pass
        return

    def finalizeOpenGL(self):
        glFreeLists(self.meshDL, 1)
        return

    def addMenuActions(self, menu):
        """ Add possible cavity-specific menu options. """
        return
