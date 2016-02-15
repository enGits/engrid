#!/usr/bin/python

import drnum
import math

R = 0.1
H = 0.1
n = 400
s = 1.0
h = 2*R/n

nof = 0.0
gf  = 100
Nm  = 10

x0 = (R + 5*h)
z0 = (H/2 + 5*s*h)
N  = 10

x = [-x0, x0]
y = x
z = [-z0, z0]

mesh = drnum.Mesh()
patch = mesh.createRectGrid(x, y, z)
mesh.setNeighbourOverlapFactor(nof)

patch[0][0][0].setResI(h)
patch[0][0][0].setResJ(h)
patch[0][0][0].setResK(s*h)

mesh.setGrading2(gf)
mesh.setMinDim(Nm)

mesh.update()
mesh.save('patches/standard.grid')
mesh.printInfo()

print
print str(H/h) + ' cells over the domain height'
print str(n) + ' cells over the diameter'
print 'H   = ' + str(H)
print 'H/2 = ' + str(H/2)
print 'H/R = ' + str(H/R)
