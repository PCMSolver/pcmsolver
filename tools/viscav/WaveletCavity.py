""" Implements Cavity class for Helmut's cavity definitons. """

from OpenGL.GL import *

import Cavity
from QuadraticInterpolation import interpolate1D, gen2DInterpolationCoefficients, interpolate2D, interpolate2DNormal


class WaveletCavity(Cavity.Cavity):
    def __init__(self, filename):

        infile = open(filename,"r")

        self.patchLevel = int(infile.readline().split()[0])
        self.nPatches = int(infile.readline().split()[0])

        pl = 2**self.patchLevel + 1
        self.mesh = []
        patchNumber = 0
        patch = [[None for i in xrange(pl)] for j in xrange(pl)]
        
        # Error checking method
        def checkPatch(patch):
            for i in patch:
                for j in patch:
                    if j == None:
                        raise Exception('Broken file/bad format')
            return True

        # Back to file reading
        for line in infile:

            parts = line.split()

            if int(parts[0]) != patchNumber:

                # Check everything is ok
                checkPatch(patch)

                if int(parts[0]) != (patchNumber + 1):
                    raise Exception('Broken file')

                # Seems so, more or less
                self.mesh.append(patch)
                patch = [[None for i in xrange(pl)] for j in xrange(pl)]

                patchNumber += 1
            
            patch[int(parts[1])][int(parts[2])] = [float(parts[3]), float(parts[4]), float(parts[5])]

        # Handle the end of file

        checkPatch(patch)
        self.mesh.append(patch)
        if (patchNumber + 1) != self.nPatches:
            raise Exception('Broken file')
        infile.close()

        # Other initialization

        self.tesselationLevel = 5

        self.drawPatchLines = True
        self.drawSubPatchLines = True
        self.drawPatches = True

        self.scaleMesh()

        return


    def scaleMesh(self):
        import math
        """ Translates the center of mesh to origin and scales the mesh inside a unit ball. """
        
        # Center the mesh first
        center = [0.0, 0.0, 0.0]
        numPoints = 0
        for mesh in self.mesh:
            for line in mesh:
                for point in line:
                    for i in xrange(3):
                        center[i] += point[i]
                    numPoints += 1
        center = [c/numPoints for c in center]

        # Translate and compute maximum radius (squared)
        maxRad = 0.0
        for mesh in self.mesh:
            for line in mesh:
                for point in line:
                    for i in xrange(3):
                        point[i] -= center[i]
                    rad = sum([p*p for p in point])
                    if rad > maxRad:
                        maxRad = rad
        self.scaleFactor = 1.0/math.sqrt(maxRad)

        # Scale
        for mesh in self.mesh:
            for line in mesh:
                for point in line:
                    for i in xrange(3):
                        point[i] *= self.scaleFactor

        return


    def initializeOpenGL(self):
        """ Initialize display lists for patch lines, sub-patch lines and solid panels. """

        # Display list calls are wrapped in try-except because VirtualBox's 
        # OpenGL "translation" generates unnecessary errors.
        
        # Lines between patches

        self.patchLinesDL = glGenLists(1)
        glLoadIdentity()
        try:
            glNewList(self.patchLinesDL, GL_COMPILE)        
        except:
            pass

        glLineWidth(2.0)        
        for mesh in self.mesh:
            line1 = mesh[0]
            line2 = [m[0] for m in mesh]
            line3 = [m[-1] for m in mesh]
            line4 = mesh[-1]
            
            self.drawLine(line1)
            self.drawLine(line2)
            self.drawLine(line3)
            self.drawLine(line4)
            
        try:
            glEndList()
        except:
            pass

        # Lines between sub-patches

        self.subPatchLinesDL = glGenLists(1)
        glLoadIdentity()
        try:
            glNewList(self.subPatchLinesDL, GL_COMPILE)        
        except:
            pass
        
        glLineWidth(1.0)
        for mesh in self.mesh:
            for i in xrange(len(mesh)):
                self.drawLine(mesh[i])
                horizontalLine = [m[i] for m in mesh]
                self.drawLine(horizontalLine)

        try:
            glEndList()
        except:
            pass

        # Solid surface

        self.patchesDL = glGenLists(1)
        glLoadIdentity()
        try:
            glNewList(self.patchesDL, GL_COMPILE)        
        except:
            pass

        uvstep = 1.0/self.tesselationLevel
        for mesh in self.mesh:
            pi = 0
            while (pi+2) < len(mesh[0]):
                pj = 0
                while (pj+2) < len(mesh):
                    controlPoints = [mesh[pi][pj:pj+3], mesh[pi+1][pj:pj+3], mesh[pi+2][pj:pj+3]]
                    coefs = gen2DInterpolationCoefficients(controlPoints)

                    # Generate vertices

                    vertices = []
                    normals = []
                    u = 0.0
                    for i in xrange(self.tesselationLevel+1):
                        v = 0.0
                        tmpv = []
                        tmpn = []
                        for j in xrange(self.tesselationLevel+1):
                            vert = interpolate2D(u, v, coefs)
                            norm = interpolate2DNormal(u, v, coefs)

                            # Substract position against the normal a bit so that 
                            # patch lines etc. won't be drawn inside the surface.
                            vert = [ve - 1e-3*n*self.scaleFactor for ve,n in zip(vert, norm)]

                            tmpv.append(vert)
                            tmpn.append(norm)
                            v += uvstep

                        vertices.append(tmpv)
                        normals.append(tmpn)
                        u += uvstep

                    # Draw panels

                    for i in xrange(0, len(vertices) - 1):
                        glBegin(GL_TRIANGLE_STRIP)
                        for j in xrange(len(vertices[i])):
                            glNormal3f(normals[i+1][j][0], normals[i+1][j][1], normals[i+1][j][2])
                            glVertex3f(vertices[i+1][j][0], vertices[i+1][j][1], vertices[i+1][j][2])
                            glNormal3f(normals[i][j][0], normals[i][j][1], normals[i][j][2])
                            glVertex3f(vertices[i][j][0], vertices[i][j][1], vertices[i][j][2])
                        glEnd()

                    pj += 2
                pi += 2
        
        try:
            glEndList()
        except:
            pass

        return



    def drawOpenGL(self, useLighting):

        if useLighting:
            glDisable(GL_LIGHTING)

        if self.drawPatchLines:
            glColor3f(0.0, 0.5, 0.5)
            try:
                glCallList(self.patchLinesDL)
            except:
                pass

        if self.drawSubPatchLines:
            glColor3f(0.3, 0.3, 0.0)
            try:
                glCallList(self.subPatchLinesDL)                
            except:
                pass

        if useLighting:
            glEnable(GL_LIGHTING)

        if self.drawPatches:
            try:
                glColor3f(0.0, 0.0, 0.4)
                glCallList(self.patchesDL)                
            except:
                pass

        return



    def finalizeOpenGL(self):

        glDeleteLists(self.patchLinesDL, 1)
        glDeleteLists(self.subPatchLinesDL, 1)
        glDeleteLists(self.PatchesDL, 1)

        return


    def drawLine(self, line):
        """ Helper routine for line drawing """

        glBegin(GL_LINE_STRIP)
        pointIndex = 0
        sstep = 1.0/self.tesselationLevel
        while (pointIndex+2) < len(line):
            s = 0.0
            for i in xrange(self.tesselationLevel+1):
                p = interpolate1D(s, line[pointIndex], line[pointIndex+1], line[pointIndex+2])
                glVertex3f(p[0], p[1], p[2])
                s += sstep
            pointIndex += 2
        glEnd()
        return


    def addMenuActions(self, menu):
        from PyQt4 import Qt, QtGui, QtCore

        patchLineMenuItem = menu.addAction('Draw patch lines')
        patchLineMenuItem.setCheckable(True)
        patchLineMenuItem.setChecked(self.drawPatchLines)
        menu.connect(patchLineMenuItem, QtCore.SIGNAL('triggered()'), self.changeDrawPatchLines)

        subPatchLineMenuItem = menu.addAction('Draw sub-patch lines')
        subPatchLineMenuItem.setCheckable(True)
        subPatchLineMenuItem.setChecked(self.drawSubPatchLines)
        menu.connect(subPatchLineMenuItem, QtCore.SIGNAL('triggered()'), self.changeDrawSubPatchLines)

        patchMenuItem = menu.addAction('Draw patches')
        patchMenuItem.setCheckable(True)
        patchMenuItem.setChecked(self.drawPatches)
        menu.connect(patchMenuItem, QtCore.SIGNAL('triggered()'), self.changeDrawPatches)
        
        return

    def changeDrawPatchLines(self):
        self.drawPatchLines = not self.drawPatchLines
        return

    def changeDrawSubPatchLines(self):
        self.drawSubPatchLines = not self.drawSubPatchLines
        return

    def changeDrawPatches(self):
        self.drawPatches = not self.drawPatches
        return


if __name__ == '__main__':
    cavity = WaveletCavity('molec_dyadic.dat')
    print cavity.patchLevel, cavity.nPatches, cavity.mesh[-1][-1][-1]
