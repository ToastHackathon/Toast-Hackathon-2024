# TODO: YOUR NAME AND STUDENT NUMBER HERE

import numpy as np

class Contact:
    # p is point of contact
    # n is normal
    # d is penetration depth
    def __init__(self, body1, body2, p, n, d):
        self.body1 = body1
        self.body2 = body2
        self.p = p
        self.n = n
        if (abs(n[0]) < abs(n[1] and abs(n[0]) < abs(n[2]))):
            tmp = np.array([1,0,0])
        elif (abs(n[1]) < abs(n[2])):
            tmp = np.array([0,1,0])
        else:
            tmp = np.array([0,0,1])
        self.t2 = np.cross(self.n, tmp)
        self.t1 = np.cross(self.t2, self.n)
        self.d = d
        self.lamb = np.zeros(3)
        self.massINV = np.zeros(1)
        self.Jacobian = np.zeros(1)
        self.AUDL = []

    def compute_jacobian(self):
        row1 = np.array([-self.n.T,-np.cross(self.body1.x-self.p,self.n).T,self.n.T,np.cross(self.body2.x-self.p,self.n).T])

        t1 = np.array([-self.n[1],self.n[0],0])
        if (np.linalg.norm(t1)!=0):
            t1 = t1/np.linalg.norm(t1)
        t2 = np.array([self.n[0]*self.n[2],self.n[1]*self.n[2],self.n[0]**2-self.n[1]**2])
        if (np.linalg.norm(t2)!=0):
            t2 = t2/np.linalg.norm(t2)

        row2 = np.array([-t1.T,-np.cross(self.body1.x-self.p,t1).T,t1.T,np.cross(self.body2.x-self.p,t1).T])
        row3 = np.array([-t2.T,-np.cross(self.body1.x-self.p,t2).T,t2.T,np.cross(self.body2.x-self.p,t2).T])
        return np.array([row1.flatten(),row2.flatten(),row3.flatten()])

    def compute_jacobian_derivative(self):
        row1 = np.array([-self.n.T,-np.cross(self.body1.v,self.n).T,self.n.T,np.cross(self.body2.v,self.n).T])

        t1 = np.array([-self.n[1],self.n[0],0])
        if (np.linalg.norm(t1)!=0):
            t1 = t1/np.linalg.norm(t1)
        t2 = np.array([self.n[0]*self.n[2],self.n[1]*self.n[2],self.n[0]**2-self.n[1]**2])
        if (np.linalg.norm(t2)!=0):
            t2 = t2/np.linalg.norm(t2)

        row2 = np.array([-t1.T,-np.cross(self.body1.v,t1).T,t1.T,np.cross(self.body2.v,t1).T])
        row3 = np.array([-t2.T,-np.cross(self.body1.v,t2).T,t2.T,np.cross(self.body2.v,t2).T])
        return np.array([row1.flatten(),row2.flatten(),row3.flatten()])
    def compute_inv_effective_mass(self):
        m1Inverse = self.body1.mass_inv
        m2Inverse = self.body2.mass_inv
        arr = [[m1Inverse,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
               [0,m1Inverse,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
               [0,0,m1Inverse,0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],

               [0.00, 0.00, 0.00,self.body1.Jinv0[0][0],self.body1.Jinv0[0][1],self.body1.Jinv0[0][2],0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
               [0.00, 0.00, 0.00,self.body1.Jinv0[1][0],self.body1.Jinv0[1][1],self.body1.Jinv0[1][2],0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
               [0.00, 0.00, 0.00,self.body1.Jinv0[2][0],self.body1.Jinv0[2][1],self.body1.Jinv0[2][2],0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
               #now starts the second body
               [0.00, 0.00, 0.00, 0.00, 0.00, 0.00,m2Inverse,0.00, 0.00, 0.00, 0.00, 0.00],
               [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,m2Inverse,0.00, 0.00, 0.00, 0.00],
               [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,m2Inverse,0.00, 0.00, 0.00],
               [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,self.body2.Jinv0[0][0],self.body2.Jinv0[0][1],self.body2.Jinv0[0][2]],
               [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,self.body2.Jinv0[1][0],self.body2.Jinv0[1][1],self.body2.Jinv0[1][2]],
               [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,self.body2.Jinv0[2][0],self.body2.Jinv0[2][1],self.body2.Jinv0[2][2]]]
        mInv = np.array(arr)
        return mInv