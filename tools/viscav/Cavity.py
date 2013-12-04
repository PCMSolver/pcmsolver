""" This file describes a cavity class that takes care of input reading and OpenGL drawing. """

class Cavity:
    def __init__(self, filename):
        """ Open the file 'filename' and read data. Throw an error if inconsistencies are found. """
        raise Exception('Not implemented!')
        return

    def initializeOpenGL(self):
        """ Convert the data so that it's easily drawn (e.g., create display lists). """
        """ The object should be located at origin with a 'normalized' size of 1.0. """
        raise Exception('Not implemented!')
        return

    def drawOpenGL(self, useLighting):
        """ Draw the object (orientation is set by the caller). """
        raise Exception('Not implemented!')
        return

    def finalizeOpenGL(self):
        """ Free all OpenGL specific data (display lists and so on). """
        raise Exception('Not implemented!')
        return

    def addMenuActions(self, menu):
        """ Add possible cavity-specific menu options. """
        return
