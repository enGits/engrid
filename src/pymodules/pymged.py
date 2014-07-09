#!/usr/bin/python
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2014 enGits GmbH                                      +
# +                                                                      +
# + enGrid is free software: you can redistribute it and/or modify       +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import commands
import sys
import math
import time

class DataBase:
   
  def __init__(self, file_name):
    self.file_name = file_name
    self.object_type = {}
    self.setResolution(10)
    self.setTolerance(0.005)
    self.update()
    
  def setResolution(self, N):
    self.resolution = N
    
  def setTolerance(self, t):
    self.tolerance = t
    
  def update(self):
    cmd = "mged -c " + self.file_name + " ls -l"
    lines = commands.getoutput(cmd).split("\n");
    for line in lines:
      words = line.split()
      if len(words) > 1:
        self.object_type[words[0]] = words[1]
        
  def clear(self):
    cmd = "rm -rf " + self.file_name
    commands.getoutput(cmd)
    
  def mged(self, *args):
    cmd = "mged -c " + self.file_name + " '"
    space = 0
    for arg in args:
      if space == 1:
        cmd += " "
      space = 1
      cmd += str(arg)
    cmd += "'"
    pcmd = cmd
    if len(cmd) > 80:
      pcmd = cmd[:80] + ' ...'
    print pcmd
    self.update()
    output = commands.getoutput(cmd)
    print output
    return output
  
  def getObjects(self):
    return self.object_type.keys()

  def getSubObjects(self, object_name):
    sub_objects = []
    lines = self.mged("l", object_name).split("\n")
    for line in lines:
      words = line.split()
      if len(words) >= 2:
        if words[0] == "u" or words[0] == "-" or words[0] == "+":
          sub_objects.append(words[1])
    return sub_objects
    
  def getSolids(self, object_name):
    solids = []
    if self.object_type[object_name] == "region" or self.object_type[object_name] == "comb":
      for obj in self.getSubObjects(object_name):
        solids += self.getSolids(obj)
    else:
      solids.append(object_name)
    return solids
        
  def writeStl(self, object_name, file_name="undefined"):
    if file_name == "undefined":
      file_name = object_name.split(".")[0] + ".stl"
    cmd = "g-stl -b -D " + str(self.tolerance) + " -n " + str(math.pi/self.resolution) + " -o " + file_name + " " + self.file_name + " " + object_name;
    #cmd = "g-stl -b -n " + str(math.pi/self.resolution) + " -o " + file_name + " " + self.file_name + " " + object_name;
    #cmd = "g-stl -b -r " + str(self.tolerance) + " -o " + file_name + " " + self.file_name + " " + object_name;
    #cmd = "g-stl -b -r " + str(self.resolution) + " -o " + file_name + " " + self.file_name + " " + object_name;
    print cmd
    return commands.getoutput(cmd)
    
  def exportStl(self, object_name):
    self.writeStl(object_name, self.file_name + ".stl")
    
  def exportToEngrid(self, object_name):
    time.sleep(5)
    solids = self.getSolids(object_name)
    dir_name = object_name.split(".")[0] + ".gegc"
    cmd = "rm -rf " + dir_name;
    print cmd
    commands.getoutput(cmd)
    cmd = "mkdir " + dir_name
    print cmd
    commands.getoutput(cmd)
    file_name = dir_name + "/volume.stl"
    self.writeStl(object_name, file_name)
    for solid in solids:
      file_name = dir_name + "/" + solid.split(".")[0] + ".s.stl"
      self.writeStl(solid, file_name)
      
  def createSphere(self, name, x, y, z, radius):
    self.mged("in", name, "sph", x, y, z, radius)
    
  def createCylinder(self, name, x1, y1, z1, x2, y2, z2, radius):
    self.mged("in", name, "rcc", x1, y1, z1, x2-x1, y2-y1, z2-z1, radius);
    
  def createCone(self, name, x1, y1, z1, x2, y2, z2, radius1, radius2):
    self.mged("in", name, "trc", x1, y1, z1, x2-x1, y2-y1, z2-z1, radius1, radius2);
    
  def createConeWithNose(self, name, x1, y1, z1, x2, y2, z2, R, r):
    nx = x1-x2
    ny = y1-y2
    nz = z1-z2
    H = math.sqrt(nx*nx + ny*ny + nz*nz)
    nx /= H
    ny /= H
    nz /= H
    xe = R;
    ye = r - H
    scal = 1
    a1 = 0
    a2 = 0.5*math.pi
    xt = 0
    yt = 0
    count = 0
    while math.fabs(scal) > 1e-6 and count < 1000:
      a = 0.5*(a1 + a2)
      xt = r*math.cos(a)
      yt = r*math.sin(a)
      vx = xt - xe
      vy = yt - ye
      L = math.sqrt(vx*vx + vy*vy)
      scal = (vx*xt + vy*yt)/(L*r)
      if scal > 0:
        a2 = a
      else:
        a1 = a
      count += 1
    xs = x1 - r*nx
    ys = y1 - r*ny
    zs = z1 - r*nz
    xc = xs + yt*nx
    yc = ys + yt*ny
    zc = zs + yt*nz
    
    #self.mged("in", name + "_part.s", "part", x1, y1, z1, xc-x1, yc-y1, zc-z1, R, xt)
    #self.createCylinder(name + "_rcc", x1, y1, z1, x1 - 1.1*R*nx, y1 - 1.1*R*ny, z1 - 1.1*R*nz, 1.1*R)
    
    self.createCone(name + "_cone", x2, y2, z2, xc, yc, zc, R, xt)
    h = 2*self.tolerance*r
    self.createSphere(name + "_sphere", xs, ys, zs, math.sqrt(r*r + 0.25*h*h))
    self.mged("r", name, "u", name + "_cone", "u", name + "_sphere")

  def createBox(self, name, x1, y1, z1, x2, y2, z2):
    self.mged("in", name, "rpp", x1, x2, y1, y2, z1, z2)
    
  def importBot(self, name, file_name, scale, Dx, Dy, Dz):
    f = open(file_name)
    line = f.readline()
    num_nodes = int(line.split()[0])
    num_faces = int(line.split()[1])
    cmd = ""
    for i in range(num_nodes):
      line = f.readline()
      x = scale*float(line.split()[0]) + Dx
      y = scale*float(line.split()[1]) + Dy
      z = scale*float(line.split()[2]) + Dz
      cmd += " " + str(x) + " " + str(y) + " " + str(z)
    for i in range(num_faces):
      words = f.readline().split()
      for word in words:
        cmd += " " + word
    self.mged("in", name, "bot", num_nodes, num_faces, 2, 2, cmd)
    
  def rotX(self, name, angle, x=0, y=0, z=0):
    self.mged("e", name, ";", "oed /", name, ";", "keypoint", x, y, z, ";", "rot", angle, 0, 0, ";", "accept;")
    
  def createWedge(self, name, x1, y1, z1, x2, y2, z2, x3, y3, z3, L):
    ux = x2 - x1
    uy = y2 - y1
    uz = z2 - z1
    vx = x3 - x1
    vy = y3 - y1
    vz = z3 - z1
    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx
    H = math.sqrt(nx*nx + ny*ny + nz*nz)
    nx *= L/H
    ny *= L/H
    nz *= L/H
    self.mged("in", name, "arb6", x1, y1, z1, x2, y2, z2, x2+nx, y2+ny, z2+nz, x1+nx, y1+ny, z1+nz, x3, y3, z3, x3+nx, y3+ny, z3+nz)
    
    