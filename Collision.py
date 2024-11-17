# TODO: YOUR NAME AND STUDENT NUMBER HERE

import igl
import numpy as np
import polyscope as ps
from PARKER_RIGID import *
from contact import *
import os

class Collision:

    def __init__(self):
        self.contacts = []
        parser = argparse.ArgumentParser()
        parser.add_argument("--file", type=str, default = "scenes/scene2.json")
        args = parser.parse_args()
        data = json.load(open(args.file))
        #self.ground_body = RigidBody( data['bodies'] )


    def reset(self):
        self.contacts = []

    def update_display(self, show_contacts):
        if len(self.contacts) == 0 or not show_contacts:
            self.ps_contacts = ps.register_point_cloud("contacts", np.zeros((0,3)))
        else:
            # can only update the points if they have the same number :(
            pos = np.array([c.p for c in self.contacts])
            depth = np.array([c.d for c in self.contacts])
            normal = np.array([c.n for c in self.contacts])
            t1 = np.array([c.t1 for c in self.contacts])
            t2 = np.array([c.t2 for c in self.contacts])
            force_viz_scale = 2
            force = np.array([force_viz_scale*(c.lamb[0] * c.n + c.lamb[1] * c.t1 + c.lamb[2] * c.t2) for c in self.contacts])
            self.ps_contacts = ps.register_point_cloud("contacts", pos)
            self.ps_contacts.add_scalar_quantity("contact depth", depth, enabled=True)
            # self.ps_contacts.add_vector_quantity("contact normal", normal, enabled=True, radius=0.01, color=(0,0,1), vectortype='ambient')
            # self.ps_contacts.add_vector_quantity("contact t1", t1, enabled=True, radius=0.01, color=(1,0,0), vectortype='ambient')
            # self.ps_contacts.add_vector_quantity("contact t2", t2, enabled=True, radius=0.01, color=(0,1,0), vectortype='ambient')
            self.ps_contacts.add_vector_quantity("contact force", force, enabled=True, radius=0.01, color=(1,1,0), vectortype='ambient')

    # collision check with ground plane
    def check_ground( self, body ):
        # check if any vertex is below the ground plane
        # if so, compute penetration depth and normal and add to contact list
        vt = body.Vertices @ body.R.T + body.x
        for v in vt:
            if v[1] < 0: # y is up
                print("hi")
                #self.contacts.append( Contact( self.ground_body, body, v, np.array([0,1,0]), v[1] ) )

    def check_body_pair( self, body1, body2 ):
        # check if any vertex of one body is inside the other
        # NOTE: this is super gross because the signed distance function is expensive
        # thus we check the larger number of vertices with the smaller body
        # but WATCH OUT becaues this appears to be buggy in general.
        # For vericies inside the other body, compute penetration depth and normal
        v1t = body1.V @ body1.R.T + body1.x
        v2t = body2.V @ body2.R.T + body2.x
        if ( v1t.shape[0] > v2t.shape[0] ):
            S,I,C,N = igl.signed_distance( v1t, v2t, body2.F, return_normals=True )
            for i in range(len(S)):
                if S[i] < 0:
                    self.contacts.append( Contact( body1, body2, C[i], -N[i], -S[i] ) )
        else:
            S,I,C,N = igl.signed_distance( v2t, v1t, body1.F, return_normals=True )
            for i in range(len(S)):
                if S[i] < 0:
                    self.contacts.append( Contact( body2, body1, C[i], -N[i], -S[i] ) )

    def check( self, rigid_body_list ):
        self.contacts = []
        for i in range(len(rigid_body_list)):
            self.check_ground(rigid_body_list[i])
            for j in range(i+1, len(rigid_body_list)):
                self.check_body_pair(rigid_body_list[i], rigid_body_list[j])
        return len(self.contacts) > 0

    def PhysicalProcess(self,rigid_body_list,mu,gravity,h):
        for contact in self.contacts:
            elasticity = 1
            t1 = np.array([-contact.n[1],contact.n[0],0])
            t2 = np.array([contact.n[0]*contact.n[2],contact.n[1]*contact.n[2],contact.n[0]**2-contact.n[1]**2])

            t1 = t1/np.linalg.norm(t1)
            t2 = t2/np.linalg.norm(t2)

            r2 = contact.p - contact.body2.x

            if contact.body1 == self.ground_body:
                E1 = 0.5*(contact.body2.mass*np.dot(contact.body2.v,contact.body2.v)+np.dot(contact.body2.omega,np.matmul(contact.body2.J0,contact.body2.omega)))
                if (contact.body2.v[1]<0):
                    contact.body2.v[1] = 0

                #need to project via angles onto basis t1 n t2 not get scalars
                om = contact.body2.omega
                mag = np.linalg.norm(om)
                if mag!=0:
                    theta = np.arccos(np.dot(om,t2)/(mag))
                    phi = np.arccos(np.dot(om,contact.n)/mag)
                    newOm = -t1*mag*np.sin(phi)*np.sin(theta)-t2*mag*np.sin(phi)*np.cos(theta)+contact.n*mag*np.cos(phi)
                    contact.body2.omega = newOm + np.matmul(contact.body2.Jinv0,np.cross(-contact.body2.v,r2))
                else:
                    contact.body2.omega=np.matmul(contact.body2.Jinv0,np.cross(-contact.body2.v,r2))
                newVmag = 2*(E1-0.5*np.dot(contact.body2.omega,np.matmul(contact.body2.J0,contact.body2.omega)))/contact.body2.mass

                if (newVmag<0):
                    newVmag=0

                newVmag = np.sqrt(newVmag)
                oldnorm = np.linalg.norm(contact.body2.v)

                contact.body2.v = (contact.body2.v/oldnorm)*newVmag
                if contact.body2.x[1]<=0:
                    contact.body2.v[1]+=5

                continue
            '''

            v1 = contact.body1.v
            v2 = contact.body2.v
            w1 = contact.body1.omega
            w2 = contact.body2.omega


            #angular:
            r1 = contact.p-contact.body1.x

            dOmega1 = contact.body2.mass_inv*np.matmul(contact.body1.Jinv0,
                                                         np.cross(contact.body2.v+np.cross(contact.body2.omega,r2),r1))
            dOmega2 = contact.body1.mass_inv*np.matmul(contact.body2.Jinv0,
                                                         np.cross(contact.body1.v+np.cross(contact.body1.omega,r1),r2))
            contact.body1.omega = contact.body1.omega+dOmega1
            contact.body2.omega = contact.body2.omega+dOmega2

            w1Tilda = contact.body1.omega
            w2Tilda = contact.body2.omega
            m1 = contact.body1.mass
            m2 = contact.body2.mass
            v1 = np.linalg.norm(v1)
            v2 = np.linalg.norm(v2)

            E1 = 0.5*(m1*v1**2+m2*v2**2+np.dot(w1,np.matmul(contact.body1.J0,w1))+np.dot(w2,np.matmul(contact.body2.J0,w2)))
            RotTilda = 0.5*(np.dot(w1Tilda,np.matmul(contact.body1.J0,w1Tilda))+np.dot(w2Tilda,np.matmul(contact.body2.J0,w2Tilda)))

            #linear:
            x = E1-RotTilda

            v1Tilda = (m1*m1*v1+m1*m2*v2-np.sqrt(np.abs(-(m1*m2*(m1*m1*v1*v1+2*m1*m2*v1*v2-2*m1*x-2*m2*x+m2*m2*v2*v2)))))
            v1Tilda = v1Tilda/(m1*(m1+m2))

            v2Tilda = (m1*v1+m2*v2-m1*v1Tilda)/m2

            v1 = contact.body1.v
            mag = np.linalg.norm(v1)
            if mag != 0:
                theta = np.arccos(np.dot(v1,t2)/(mag))
                phi = np.arccos(np.dot(v1,contact.n)/mag)
                newV1 = t1*v1Tilda*np.sin(phi)*np.sin(theta)+t2*v1Tilda*np.sin(phi)*np.cos(theta)-contact.n*v1Tilda*np.cos(phi)

                contact.body1.v = newV1
            else:
                contact.body1.v = v1Tilda*contact.n

            v2 = contact.body2.v
            mag = np.linalg.norm(v2)
            if mag != 0:
                theta = np.arccos(np.dot(v2,t2)/mag)
                phi = np.arccos(np.dot(v2,contact.n)/mag)
                newV2 = t2*v2Tilda*np.sin(phi)*np.sin(theta)+t2*v2Tilda*np.sin(phi)*np.cos(theta)-contact.n*v2Tilda*np.cos(phi)
                contact.body2.v = newV2
            else:
                contact.body2.v = v2Tilda*contact.n

            contact.body1.v += h*gravity
            contact.body2.v += h*gravity
            relativeVel = contact.body2.v-contact.body1.v
            if (np.linalg.norm(relativeVel) <= 0.005):
                contact.body1.v += contact.n*5
                contact.body2.v -= contact.n*5
        '''
        return



    def process(self,rigid_body_list,mu,stepsize):
        k=0
        for contact in self.contacts:
            if contact.body1 != self.ground_body:
                velocities = np.array([contact.body1.v,contact.body1.omega,contact.body2.v,contact.body2.omega])
                velocities = velocities.flatten()
                Fext = np.array([contact.body1.force,contact.body1.torque,contact.body2.force,contact.body2.torque]).flatten()
            else:
                velocities = np.array([np.zeros(3),np.zeros(3),contact.body2.v,contact.body2.omega])
                velocities = velocities.flatten()
                Fext = np.array([np.zeros(3),np.zeros(3),contact.body2.force,contact.body2.torque]).flatten()
            if k==0:
                Cjacobian=contact.compute_jacobian()
                contact.Jacobian = Cjacobian
                MINV = contact.compute_inv_effective_mass()
                contact.massINV = MINV
            b = np.matmul(contact.Jacobian,velocities)
            A = np.matmul(contact.Jacobian,contact.massINV)
            FextTerm = np.matmul(A,Fext)
            if k==0:
                A = np.matmul(A,contact.Jacobian.T)
                D = np.diag(A)
                U = np.triu(A,1)
                L = np.tril(A)
                Linv = np.linalg.inv(L)
                contact.AUDL = [A,D,U,L,Linv]
            RHS = -b
            G = contact.AUDL[1].size
            LambdasOLD = contact.lamb
            LambdasNew = LambdasOLD
            for i in range(0,G):
                LambdasNew[i] = (1/contact.AUDL[1][i])*(RHS[i])
                for j in range(1,G):
                    LambdasNew[i] -= (1/D[i])*(contact.AUDL[0][i][j]*LambdasOLD[j])
            LambdasOLD=LambdasNew

            LambdaHi = np.zeros(LambdasOLD.size)
            for GSD in range(0,100):
                LambdasNew = np.matmul(contact.AUDL[4],-np.matmul(contact.AUDL[2],LambdasOLD)+RHS)

                #projecting
                for iterate in range(LambdasNew.size):
                    LambdaHi[iterate] = mu*LambdasNew[iterate-iterate%3]
                    LambdaLo = -LambdaHi
                final = np.zeros(LambdasNew.size)
                for newIterate in range(final.size):
                    if LambdaLo[newIterate]<LambdasNew[newIterate]:
                        temp = LambdasNew[newIterate]
                    else:
                        temp = LambdaLo[newIterate]
                    if temp<LambdaHi[newIterate]:
                        final[newIterate] = temp
                    else:
                        final[newIterate] = LambdaHi[newIterate]
                LambdasNew = final
                LambdasOLD=LambdasNew

            for l in range(LambdasNew.size):
                if LambdasNew[l]<0:
                    LambdasNew[l]=0
            dV = np.matmul(np.matmul(contact.massINV,contact.Jacobian.T),LambdasNew)

            magnitude = np.linalg.norm(dV[0:3])
            if magnitude != 0:
                arg =(np.dot(dV[0:3],contact.n)/magnitude)
                if arg>1:
                    phi = 0
                elif arg<-1:
                    phi = np.pi
                else:
                    phi = np.arccos(arg)
                dV[0:3]  = contact.n*magnitude*np.cos(phi)

            magnitude = np.linalg.norm(dV[6:9])
            if magnitude != 0:
                dV[6:9] = contact.n*np.dot(dV[6:9],contact.n)


            contact.lamb = LambdasNew
            Forces = np.matmul(Cjacobian.T,LambdasOLD)
            #Forces = np.array([Forces[2],Forces[0],Forces[1],Forces[5],Forces[3],Forces[4],Forces[8],Forces[6],Forces[7],Forces[11],Forces[9],Forces[10]])

            contact.body1.v += np.array([dV[0],dV[1],dV[2]])
            contact.body1.omega +=np.array([dV[3],dV[4],dV[5]])
            contact.body2.v += np.array([dV[6],dV[7],dV[8]])
            contact.body2.omega += np.array([dV[9],dV[10],dV[11]])
            #project the forces for each body into t1,t2,n coordinates
            t1 = np.array([-contact.n[1],contact.n[0],0])
            t2 = np.array([contact.n[0]*contact.n[2],contact.n[1]*contact.n[2],contact.n[0]**2-contact.n[1]**2])

            t1 = t1/np.linalg.norm(t1)
            t2 = t2/np.linalg.norm(t2)
            if contact.body1!=self.ground_body:
                ForcesBody1 = Forces[0:3]
                magnitude = np.linalg.norm(ForcesBody1)
                if magnitude != 0:
                    arg = np.dot(ForcesBody1,contact.n)/magnitude
                    if arg>1:
                        phi = 0
                    elif arg<-1:
                        phi = np.pi
                    else:
                        phi = np.arccos(arg)

                    arg = np.dot(ForcesBody1,t2)/(magnitude)
                    if arg>1:
                        theta = 0
                    elif arg<-1:
                        theta = np.pi
                    else:
                        theta = np.arccos(arg)

                    Fb1N  = -contact.n*magnitude*np.cos(phi)

                else:
                    Fb1T1 = 0
                    Fb1T2 = 0
                    Fb1N = 0
                Fb1FrictionMax = mu*np.linalg.norm(Fb1N)

            ForcesBody2 = Forces[6:9]
            magnitude = np.linalg.norm(ForcesBody2)
            if magnitude != 0:
                arg = np.dot(ForcesBody2,t2)/(magnitude)
                if arg>1:
                    theta = 0
                elif arg<-1:
                    theta = np.pi
                else:
                    theta = np.arccos(arg)
                arg = np.dot(ForcesBody2,contact.n)
                Fb2N  = -contact.n*arg
            else:

                Fb2N = 0
            Fb2FrictionMax = -mu*np.linalg.norm(Fb2N)

            if (self.ground_body != contact.body1):
                #project velocities onto t1,t2,n coordinates
                v1 = contact.body1.v
                magnitude = np.linalg.norm(v1)
                if magnitude != 0:
                    theta = np.arccos(np.dot(v1,t2)/(magnitude))
                    arg = np.dot(v1,contact.n)/magnitude
                    if arg>1:
                        phi = 0
                    elif arg<-1:
                        phi = np.pi
                    else:
                        phi = np.arccos(arg)
                    Vb1T1 = t1*magnitude*np.sin(phi)*np.sin(theta)
                    Vb1T2 = t2*magnitude*np.sin(phi)*np.cos(theta)
                    Vb1N  = -contact.n*magnitude*np.cos(phi)

                else:
                    Vb1T1 = 0*t1
                    Vb1T2 = 0*t2
                    Vb1N  = 0*t1

            v2 = contact.body2.v
            magnitude = np.linalg.norm(v2)
            if magnitude > 10**(-10) and magnitude<10^13:
                theta = np.arccos(np.dot(v2,t2)/magnitude)
                arg = np.dot(v2,contact.n)/magnitude
                if arg>=1:
                    phi = 0
                elif arg<=-1:
                    phi = np.pi
                elif arg==0:
                    phi = np.pi/2
                else:
                    phi = np.arccos(arg)
                Vb2T1 = t1*magnitude*np.sin(phi)*np.sin(theta)
                Vb2T2 = t2*magnitude*np.sin(phi)*np.cos(theta)
                Vb2N = -contact.n*magnitude*arg
            else:
                Vb2T1 = 0*t1
                Vb2T2 = 0*t2
                Vb2N = 0*t1


            if (np.linalg.norm(Vb2T1)>0):
                Vb2T1 -= contact.body2.mass_inv*(Fb2FrictionMax*np.sign(Vb2T1))*stepsize*t1
            if (np.linalg.norm(Vb2T2)>0):
                Vb2T2 -= contact.body2.mass_inv*(Fb2FrictionMax*np.sign(Vb2T2))*stepsize*t2
            if (self.ground_body!=contact.body1):
                if (np.linalg.norm(Vb1T1)>0):
                    Vb1T1 -= contact.body1.mass_inv*(Fb1FrictionMax*np.sign(Vb1T1))*stepsize*t1
                if (np.linalg.norm(Vb1T2)>0):
                    Vb1T2 -= contact.body1.mass_inv*(Fb1FrictionMax*np.sign(Vb1T2))*stepsize*t2
                contact.body1.v = Vb1T1+Vb1T2+Vb1N


            contact.body2.v = Vb2T1+Vb2T2+Vb2N
            r2 = contact.p - contact.body2.x

            if contact.body1 == self.ground_body:
                #need to project via angles onto basis t1 n t2 not get scalars
                om = contact.body2.omega
                mag = np.linalg.norm(om)
                if mag>10**(-10) and mag<10**10:
                    theta = np.arccos(np.dot(om,t2)/(mag))
                    phi = np.arccos(np.dot(om,contact.n)/mag)
                    newOm = -t1*mag*np.sin(phi)*np.sin(theta)-t2*mag*np.sin(phi)*np.cos(theta)+contact.n*mag*np.cos(phi)
                    contact.body2.omega = newOm + np.matmul(contact.body2.Jinv0,np.cross(-contact.body2.v,r2))
                else:
                    contact.body2.omega=np.matmul(contact.body2.Jinv0,np.cross(-contact.body2.v,r2))
            else:
                r1 = contact.p-contact.body1.x

                dOmega1 = contact.body2.mass_inv*np.matmul(contact.body1.Jinv0,
                                                           np.cross(contact.body2.v+np.cross(contact.body2.omega,r2),r1))
                dOmega2 = contact.body1.mass_inv*np.matmul(contact.body2.Jinv0,
                                                           np.cross(contact.body1.v+np.cross(contact.body1.omega,r1),r2))
                contact.body1.omega = contact.body1.omega+dOmega1
                contact.body2.omega = contact.body2.omega+dOmega2





        return