import numpy as np
import polyscope as ps
import igl
from TensorData import *
import json

g = np.array([0, -9.8, 0])
h=0.01
breadDensity=0.25
spreadDensity=1
spreadThickness=0
spreadFaces = 1

# RigidBody class

class RigidBody:
    def __init__(self, objectFile):
        #self.name = "inertial frame"

        self.v = np.zeros(3)
        self.x = np.zeros(3)
        if type(objectFile) == list:
            objectFile = objectFile[0]
        output = DoEverything(objectFile,breadDensity,spreadDensity,spreadThickness,spreadFaces)
        data = json.load(open(objectFile["file_root"]+".json"))
        self.x0 = np.array(objectFile['x0'], dtype=float)
        self.R0 = np.array(objectFile.get('R0',np.eye(3)), dtype=float)
        self.v0 = np.array(objectFile.get('v0',(0,0,0)), dtype=float)
        self.omega0 = np.array(objectFile.get('omega0',(0,0,0)), dtype=float)
        self.R = np.asarray(self.R0)
        self.omega = np.asarray(self.omega0)

        #self.x0 = output[0]
        self.Vertices = output[1] + output[0]
        self.Faces = output[2]
        self.Volume = output[3][0]
        self.J0 = output[4]
        self.cm = output[0]
        self.Jinv0 = np.linalg.inv(self.J0)
        self.mass = breadDensity*self.Volume+spreadDensity*spreadThickness
        self.mass_inv = 1/self.mass
        self.rad_crit = 0
        for v in self.Vertices:
            n = np.linalg.norm(self.x0-np.array(v))
            if (n)>self.rad_crit:
                self.rad_crit = n
        self.reset()
        self.E = 0
        self.run = True
        self.firstContact = True
        self.ps_mesh = ps.register_surface_mesh("input mesh",self.Vertices, self.Faces, smooth_shade=False)
        self.update_display()
    def reset(self):
        self.R = self.R0
        self.x = self.x0
        self.v = self.v0
        self.omega = self.omega0
        self.elapsed = 0
    def update_display(self):
        # Construct and set the homogeneous transformation matrix
        T = np.eye(4)
        T[0:3,0:3] = self.R
        T[0:3,3] = self.x
        self.ps_mesh.set_transform( T )

    def resolveContact(self):
        contactList = self.groundContact([])
        #print(len(contactList))
        if len(contactList) != 0:
            self.run = False
        else:
            self.run = True
            '''
            if (self.firstContact):
                self.E = 0.5*np.cross(self.omega,self.J0 @ self.omega) + self.mass*g*self.x[1]*np.array([0,1,0])
                self.firstContact = False
            N = len(contactList)/3
            contactList = np.array(contactList)
            vprior=self.v[1]
            self.v[1] = 0
            F = -vprior*self.mass/(h*N)*np.array([0,1,0])
            print(F)
            F_g = self.mass*g/N
            #print(F)
            damping = 0.1
            Tau_1 = np.sum(np.cross(contactList,F)) - N*(np.cross(self.x,F))
            Tau_2 = (-1)*np.sum(np.cross(contactList,F_g)) + N*(np.cross(self.x,F_g))
            Tau_net = damping*(Tau_1 + Tau_2)
            norm = np.linalg.norm(2*self.Jinv0@(self.E-self.mass*g*self.x[1]*np.array([0,1,0])))
            omega1 = self.omega + h*self.Jinv0 @ (Tau_net - np.cross(self.omega, self.J0 @ self.omega))
            if norm!=0:
                omega1 = omega1/norm
            oHat = np.array([[0,-omega1[2],omega1[1]],[omega1[2],0,-omega1[0]],[-omega1[1],omega1[0],0]])
            self.omega = omega1
            self.R = np.eye(3) + h*oHat + (h**2)*(oHat @ oHat)
            self.R = self.R/np.linalg.det(self.R)
            '''

    def groundContact(self,contactList):
        for v in self.Vertices:
            #print(v[1])

            if v[1]+self.x[1]<= 0:
                contactList.append(v)
        return contactList

    def updateVelocity(self):
        self.v = self.v + h*g

    def updatePosition(self):
        self.x = self.x + self.v*h
        if np.linalg.norm(self.omega)>10**(-10):
            o = self.omega/np.linalg.norm(self.omega)
        else:
            o = self.omega
        oHat = np.array([[0,-o[2],o[1]],[o[2],0,-o[0]],[-o[1],o[0],0]])
        #oHat = np.array([[0,-self.omega[2],self.omega[1]],[self.omega[2],0,-self.omega[0]],[-self.omega[1],self.omega[0],0]])
        self.R = self.R + h*(-np.matmul(oHat,self.R))
        #print(self.Vertices)
        toNormalize = self.R.T
        toNormalize[0] = toNormalize[0]/np.linalg.norm(toNormalize[0])
        toNormalize[1] = toNormalize[1]/np.linalg.norm(toNormalize[1])
        toNormalize[2] = toNormalize[2]/np.linalg.norm(toNormalize[2])
        self.R = toNormalize.T
        self.Vertices = (self.Vertices @ self.R.T)


# function update positions 
tElapsed=0


    

