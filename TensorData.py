
import sympy
import numpy as np
import igl
import argparse
from sympy.utilities.lambdify import implemented_function

def DoEverything(fileName,breadDensity,spreadDensity,spreadThickness,spreadFaces):

    #make the command-line parser with arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, default="data/cube.obj")
    parser.add_argument("--density", type=float, default=1000, help="density in kg/m^3")
    parser.add_argument("--scale", nargs=3, type=float, default=(1,1,1),help="scale the mesh")
    parser.add_argument("--translate", nargs=3, type=float, default=(0,0,0),help="translate (after scale)")
    parser.add_argument("--test", type=int,help="run a numbered unit test")

    args, unknown = parser.parse_known_args()
    print(args.translate)

    if type(fileName) == str:
        V, _, _, F, _, _ = igl.read_obj(fileName)
    else:
        print(fileName)
        V, _, _, F, _, _ = igl.read_obj(fileName['file_root']+'.obj')
    V = V * args.scale
    V = V + args.translate

    #assumptions!!!!! THIS CODE ASSUMES THE FOLLOWING ARE INPUTS
    #spreadFaces = list of all faces with butter on them
    #spreadThickness = thickness of the spread on those faces
    #spreadDensity
    #breadDensity
    #

    #breadDensity = 1
    #spreadDensity = 1
    spreadThickness = 1
    spreadFaces = []

    #accounting for translation
    origin = np.array([args.translate[0],args.translate[1],args.translate[2]])
    for vertex in range(len(V)):
        V[vertex] = V[vertex]-origin

    #defining the symbols required for calculation
    x = sympy.Symbol("x")
    y = sympy.Symbol("y")
    z = sympy.Symbol("z")

    A0 = sympy.Symbol("A0")
    A1 = sympy.Symbol("A1")
    A2 = sympy.Symbol("A2")

    B0 = sympy.Symbol("B0")
    B1 = sympy.Symbol("B1")
    B2 = sympy.Symbol("B2")

    C0 = sympy.Symbol("C0")
    C1 = sympy.Symbol("C1")
    C2 = sympy.Symbol("C2")

    D0 = sympy.Symbol("D0")
    D1 = sympy.Symbol("D1")
    D2 = sympy.Symbol("D2")

    chi = sympy.Symbol("chi")
    nu = sympy.Symbol("nu")
    sigma = sympy.Symbol("sigma")

    mu = sympy.Symbol("mu")


    # In[1147]:


    #gets the volume of a tetrahedron using V = (1/6)A*(BxC)
    #A,B,C the vector coordinates of 3 vertices
    def volumeOfTetrahedron(npTxTMatrix):

        origin = np.array([args.translate[0],args.translate[1],args.translate[2]])
        Volume = (np.dot(npTxTMatrix[0]-origin,np.cross(npTxTMatrix[1]-origin,npTxTMatrix[2]-origin)))/6

        return Volume


    #Calls volumeOfTetrahedron on each tetrahedron of the mesh
    #adds up the volumes to get total volume

    def meshVolume(breadDensity):
        totalVolume = 0

        for face in F:
            tempMatrix = np.array([V[face[0]],V[face[1]],V[face[2]]])
            v = volumeOfTetrahedron(tempMatrix)
            totalVolume +=v

        totalMass = totalVolume*breadDensity

        return totalVolume,totalMass


    # In[1148]:


    #Constant density => average of position of each vertex
    def findWeightedsCOM(tetrahedron):

        origin = np.array([args.translate[0],args.translate[1],args.translate[2]])
        CM = tetrahedron[0]-4*origin+tetrahedron[1]+tetrahedron[2]+tetrahedron[3]
        CM = CM/4 + origin

        return CM

    def getObjectCOM(breadDensity):
        COM = np.zeros(3)
        origin = np.array([args.translate[0],args.translate[1],args.translate[2]])
        V1, _, _, F1, _, _ = igl.read_obj(args.file)
        vol,mass = meshVolume(breadDensity)
        for face in F1:

            tet1 = np.array([V[face[0]],V[face[1]],V[face[2]],[0,0,0]])
            tet2 = np.array([V[face[0]]+origin,V[face[1]]+origin,V[face[2]]+origin,origin])
            COM = COM + (findWeightedsCOM(tet1))*volumeOfTetrahedron(tet2)*args.density/mass

        return COM+origin


    # In[1149]:


    #getting the inertia tensor:

    #these equations were calculated via triple integration
    #using the equation for the entries of an inertia tensor
    #lambdified for faster run time (1000 times faster)

    EqJ00Lambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,mu],
                             A1**2*mu/60 + A1*B1*mu/60 + A1*C1*mu/60 + A1*D1*mu/60 + A2**2*mu/60 +
                             A2*B2*mu/60 + A2*C2*mu/60 + A2*D2*mu/60 + B1**2*mu/60 + B1*C1*mu/60 +
                             B1*D1*mu/60 + B2**2*mu/60 + B2*C2*mu/60 + B2*D2*mu/60 + C1**2*mu/60 +
                             C1*D1*mu/60 + C2**2*mu/60 + C2*D2*mu/60 + D1**2*mu/60 + D2**2*mu/60)


    EqJ11Lambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,mu],
                             A0**2*mu/60 + A0*B0*mu/60 + A0*C0*mu/60 + A0*D0*mu/60 + A2**2*mu/60 +
                             A2*B2*mu/60 + A2*C2*mu/60 + A2*D2*mu/60 + B0**2*mu/60 + B0*C0*mu/60 +
                             B0*D0*mu/60 + B2**2*mu/60 + B2*C2*mu/60 + B2*D2*mu/60 + C0**2*mu/60 +
                             C0*D0*mu/60 + C2**2*mu/60 + C2*D2*mu/60 + D0**2*mu/60 + D2**2*mu/60)

    EqJ22Lambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,mu],
                             A0**2*mu/60 + A0*B0*mu/60 + A0*C0*mu/60 + A0*D0*mu/60 + A1**2*mu/60 +
                             A1*B1*mu/60 + A1*C1*mu/60 + A1*D1*mu/60 + B0**2*mu/60 + B0*C0*mu/60 +
                             B0*D0*mu/60 + B1**2*mu/60 + B1*C1*mu/60 + B1*D1*mu/60 + C0**2*mu/60 +
                             C0*D0*mu/60 + C1**2*mu/60 + C1*D1*mu/60 + D0**2*mu/60 + D1**2*mu/60)

    EqJ01Lambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,mu],
                             -A0*A1*mu/60 - A0*B1*mu/120 - A0*C1*mu/120 - A0*D1*mu/120 - A1*B0*mu/120 -
                             A1*C0*mu/120 - A1*D0*mu/120 - B0*B1*mu/60 - B0*C1*mu/120 - B0*D1*mu/120 -
                             B1*C0*mu/120 - B1*D0*mu/120 - C0*C1*mu/60 -
                             C0*D1*mu/120 - C1*D0*mu/120 - D0*D1*mu/60)


    EqJ02Lambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,mu],
                             -A0*A2*mu/60 - A0*B2*mu/120 - A0*C2*mu/120 - A0*D2*mu/120 - A2*B0*mu/120 -
                             A2*C0*mu/120 - A2*D0*mu/120 - B0*B2*mu/60 - B0*C2*mu/120 - B0*D2*mu/120 -
                             B2*C0*mu/120 - B2*D0*mu/120 - C0*C2*mu/60
                             - C0*D2*mu/120 - C2*D0*mu/120 - D0*D2*mu/60)


    EqJ12Lambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,D0,D1,D2,mu],
                             -A1*A2*mu/60 - A1*B2*mu/120 - A1*C2*mu/120 - A1*D2*mu/120 - A2*B1*mu/120 -
                             A2*C1*mu/120 - A2*D1*mu/120 - B1*B2*mu/60 - B1*C2*mu/120 - B1*D2*mu/120 -
                             B2*C1*mu/120 - B2*D1*mu/120 - C1*C2*mu/60 - C1*D2*mu/120 -
                             C2*D1*mu/120 - D1*D2*mu/60)


    def Jtensor(A,B,C,D,breadDensity):

        J = np.zeros((3,3))

        #principle inertias:

        J[0][0] = EqJ00Lambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2],breadDensity)
        J[1][1] = EqJ11Lambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2],breadDensity)
        J[2][2] = EqJ22Lambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2],breadDensity)

        #off diagonal entries (makes use of symmetric property of the tensor)

        J[0][1] = EqJ01Lambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2],breadDensity)
        J[1][0] = J[0][1]
        J[0][2] = EqJ02Lambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2],breadDensity)
        J[2][0] = J[0][2]
        J[1][2] = EqJ12Lambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2],breadDensity)
        J[2][1] = J[1][2]

        #scale appropriately
        J= J*4

        #ignore floating point errors
        J = np.where(np.abs(J)>10**(-12),J,0)

        return J



    IzzSPREADLambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,mu],
                mu*(-(6*A1**2*(B0**2*(A0*B1 - A0*C1 - A1*B0 + A1*C0 + B0*C1 - B1*C0) +
                C0**2*(-A0*B1 + A0*C1 + A1*B0 - A1*C0 - B0*C1 + B1*C0))*(A0**2*B0**2 - 2*A0**2*B0*C0 + A0**2*C0**2 -
                2*A0*B0**2*C0 + 4*A0*B0*C0**2 - 2*A0*C0**3 + B0**2*C0**2 - 2*B0*C0**3 + C0**4)*(A0**3*B0**3 -
                3*A0**3*B0**2*C0 + 3*A0**3*B0*C0**2 - A0**3*C0**3 - 3*A0**2*B0**3*C0 + 9*A0**2*B0**2*C0**2 - 9*A0**2*B0*C0**3
                + 3*A0**2*C0**4 + 3*A0*B0**3*C0**2 - 9*A0*B0**2*C0**3 + 9*A0*B0*C0**4 - 3*A0*C0**5 - B0**3*C0**3 + 3*B0**2*C0**4 - 3*B0*C0**5 + C0**6) + 4*A1*(B0**3*(A0**2*B1**2 - 2*A0**2*B1*C1 + A0**2*C1**2 - 2*A0*B1**2*C0 + 4*A0*B1*C0*C1 - 2*A0*C0*C1**2 - A1**2*B0**2 + 2*A1**2*B0*C0 - A1**2*C0**2 + 2*A1*B0**2*C1 - 4*A1*B0*C0*C1 + 2*A1*C0**2*C1 - B0**2*C1**2 + 2*B0*C0*C1**2 + B1**2*C0**2 - 2*B1*C0**2*C1) + C0**3*(-A0**2*B1**2 + 2*A0**2*B1*C1 - A0**2*C1**2 + 2*A0*B1**2*C0 - 4*A0*B1*C0*C1 + 2*A0*C0*C1**2 + A1**2*B0**2 - 2*A1**2*B0*C0 + A1**2*C0**2 - 2*A1*B0**2*C1 + 4*A1*B0*C0*C1 - 2*A1*C0**2*C1 + B0**2*C1**2 - 2*B0*C0*C1**2 - B1**2*C0**2 + 2*B1*C0**2*C1))*(A0*B0 - A0*C0 - B0*C0 + C0**2)*(A0**3*B0**3 - 3*A0**3*B0**2*C0 + 3*A0**3*B0*C0**2 - A0**3*C0**3 - 3*A0**2*B0**3*C0 + 9*A0**2*B0**2*C0**2 - 9*A0**2*B0*C0**3 + 3*A0**2*C0**4 + 3*A0*B0**3*C0**2 - 9*A0*B0**2*C0**3 + 9*A0*B0*C0**4 - 3*A0*C0**5 - B0**3*C0**3 + 3*B0**2*C0**4 - 3*B0*C0**5 + C0**6) + (B0**4*(3*A0**3*B0**2*B1 - 3*A0**3*B0**2*C1 - 6*A0**3*B0*B1*C0 + 6*A0**3*B0*C0*C1 + A0**3*B1**3 - 3*A0**3*B1**2*C1 + 3*A0**3*B1*C0**2 + 3*A0**3*B1*C1**2 - 3*A0**3*C0**2*C1 - A0**3*C1**3 - 3*A0**2*A1*B0**3 + 9*A0**2*A1*B0**2*C0 - 9*A0**2*A1*B0*C0**2 + 3*A0**2*A1*C0**3 + 3*A0**2*B0**3*C1 - 9*A0**2*B0**2*B1*C0 + 18*A0**2*B0*B1*C0**2 - 9*A0**2*B0*C0**2*C1 - 3*A0**2*B1**3*C0 + 9*A0**2*B1**2*C0*C1 - 9*A0**2*B1*C0**3 - 9*A0**2*B1*C0*C1**2 + 6*A0**2*C0**3*C1 + 3*A0**2*C0*C1**3 + 6*A0*A1*B0**3*C0 - 18*A0*A1*B0**2*C0**2 + 18*A0*A1*B0*C0**3 - 6*A0*A1*C0**4 - 6*A0*B0**3*C0*C1 + 9*A0*B0**2*B1*C0**2 + 9*A0*B0**2*C0**2*C1 - 18*A0*B0*B1*C0**3 + 3*A0*B1**3*C0**2 - 9*A0*B1**2*C0**2*C1 + 9*A0*B1*C0**4 + 9*A0*B1*C0**2*C1**2 - 3*A0*C0**4*C1 - 3*A0*C0**2*C1**3 - A1**3*B0**3 + 3*A1**3*B0**2*C0 - 3*A1**3*B0*C0**2 + A1**3*C0**3 + 3*A1**2*B0**3*C1 - 9*A1**2*B0**2*C0*C1 + 9*A1**2*B0*C0**2*C1 - 3*A1**2*C0**3*C1 - 3*A1*B0**3*C0**2 - 3*A1*B0**3*C1**2 + 9*A1*B0**2*C0**3 + 9*A1*B0**2*C0*C1**2 - 9*A1*B0*C0**4 - 9*A1*B0*C0**2*C1**2 + 3*A1*C0**5 + 3*A1*C0**3*C1**2 + 3*B0**3*C0**2*C1 + B0**3*C1**3 - 3*B0**2*B1*C0**3 - 6*B0**2*C0**3*C1 - 3*B0**2*C0*C1**3 + 6*B0*B1*C0**4 + 3*B0*C0**4*C1 + 3*B0*C0**2*C1**3 - B1**3*C0**3 + 3*B1**2*C0**3*C1 - 3*B1*C0**5 - 3*B1*C0**3*C1**2) + C0**4*(-3*A0**3*B0**2*B1 + 3*A0**3*B0**2*C1 + 6*A0**3*B0*B1*C0 - 6*A0**3*B0*C0*C1 - A0**3*B1**3 + 3*A0**3*B1**2*C1 - 3*A0**3*B1*C0**2 - 3*A0**3*B1*C1**2 + 3*A0**3*C0**2*C1 + A0**3*C1**3 + 3*A0**2*A1*B0**3 - 9*A0**2*A1*B0**2*C0 + 9*A0**2*A1*B0*C0**2 - 3*A0**2*A1*C0**3 - 3*A0**2*B0**3*C1 + 9*A0**2*B0**2*B1*C0 - 18*A0**2*B0*B1*C0**2 + 9*A0**2*B0*C0**2*C1 + 3*A0**2*B1**3*C0 - 9*A0**2*B1**2*C0*C1 + 9*A0**2*B1*C0**3 + 9*A0**2*B1*C0*C1**2 - 6*A0**2*C0**3*C1 - 3*A0**2*C0*C1**3 - 6*A0*A1*B0**3*C0 + 18*A0*A1*B0**2*C0**2 - 18*A0*A1*B0*C0**3 + 6*A0*A1*C0**4 + 6*A0*B0**3*C0*C1 - 9*A0*B0**2*B1*C0**2 - 9*A0*B0**2*C0**2*C1 + 18*A0*B0*B1*C0**3 - 3*A0*B1**3*C0**2 + 9*A0*B1**2*C0**2*C1 - 9*A0*B1*C0**4 - 9*A0*B1*C0**2*C1**2 + 3*A0*C0**4*C1 + 3*A0*C0**2*C1**3 + A1**3*B0**3 - 3*A1**3*B0**2*C0 + 3*A1**3*B0*C0**2 - A1**3*C0**3 - 3*A1**2*B0**3*C1 + 9*A1**2*B0**2*C0*C1 - 9*A1**2*B0*C0**2*C1 + 3*A1**2*C0**3*C1 + 3*A1*B0**3*C0**2 + 3*A1*B0**3*C1**2 - 9*A1*B0**2*C0**3 - 9*A1*B0**2*C0*C1**2 + 9*A1*B0*C0**4 + 9*A1*B0*C0**2*C1**2 - 3*A1*C0**5 - 3*A1*C0**3*C1**2 - 3*B0**3*C0**2*C1 - B0**3*C1**3 + 3*B0**2*B1*C0**3 + 6*B0**2*C0**3*C1 + 3*B0**2*C0*C1**3 - 6*B0*B1*C0**4 - 3*B0*C0**4*C1 - 3*B0*C0**2*C1**3 + B1**3*C0**3 - 3*B1**2*C0**3*C1 + 3*B1*C0**5 + 3*B1*C0**3*C1**2))*(A0*B0 - A0*C0 - B0*C0 + C0**2)*(A0**2*B0**2 - 2*A0**2*B0*C0 + A0**2*C0**2 - 2*A0*B0**2*C0 + 4*A0*B0*C0**2 - 2*A0*C0**3 + B0**2*C0**2 - 2*B0*C0**3 + C0**4))*(A0**5 - 2*A0**4*B0 - 3*A0**4*C0 + A0**3*B0**2 + 6*A0**3*B0*C0 + 3*A0**3*C0**2 - 3*A0**2*B0**2*C0 - 6*A0**2*B0*C0**2 - A0**2*C0**3 + 3*A0*B0**2*C0**2 + 2*A0*B0*C0**3 - B0**2*C0**3) + (A0*B0 - A0*C0 - B0*C0 + C0**2)*(A0**2*B0**2 - 2*A0**2*B0*C0 + A0**2*C0**2 - 2*A0*B0**2*C0 + 4*A0*B0*C0**2 - 2*A0*C0**3 + B0**2*C0**2 - 2*B0*C0**3 + C0**4)*(A0**3*B0**3 - 3*A0**3*B0**2*C0 + 3*A0**3*B0*C0**2 - A0**3*C0**3 - 3*A0**2*B0**3*C0 + 9*A0**2*B0**2*C0**2 - 9*A0**2*B0*C0**3 + 3*A0**2*C0**4 + 3*A0*B0**3*C0**2 - 9*A0*B0**2*C0**3 + 9*A0*B0*C0**4 - 3*A0*C0**5 - B0**3*C0**3 + 3*B0**2*C0**4 - 3*B0*C0**5 + C0**6)*(3*A0**8*B1 - 3*A0**8*C1 - 3*A0**7*A1*B0 + 3*A0**7*A1*C0 - 3*A0**7*B0*B1 + 6*A0**7*B0*C1 - 9*A0**7*B1*C0 + 6*A0**7*C0*C1 + 17*A0**6*A1**2*B1 - 17*A0**6*A1**2*C1 + 3*A0**6*A1*B0**2 + 3*A0**6*A1*B0*C0 - 7*A0**6*A1*B1**2 - 6*A0**6*A1*C0**2 + 7*A0**6*A1*C1**2 - 3*A0**6*B0**2*C1 + 9*A0**6*B0*B1*C0 - 12*A0**6*B0*C0*C1 + A0**6*B1**3 + 9*A0**6*B1*C0**2 - 3*A0**6*C0**2*C1 - A0**6*C1**3 - 17*A0**5*A1**3*B0 + 17*A0**5*A1**3*C0 - 3*A0**5*A1**2*B0*B1 + 34*A0**5*A1**2*B0*C1 - 51*A0**5*A1**2*B1*C0 + 20*A0**5*A1**2*C0*C1 - 6*A0**5*A1*B0**2*C0 - 3*A0**5*A1*B0*B1**2 + 3*A0**5*A1*B0*C0**2 - 14*A0**5*A1*B0*C1**2 + 21*A0**5*A1*B1**2*C0 + 3*A0**5*A1*C0**3 - 4*A0**5*A1*C0*C1**2 + 6*A0**5*B0**2*C0*C1 + A0**5*B0*B1**3 - 9*A0**5*B0*B1*C0**2 + 6*A0**5*B0*C0**2*C1 + 2*A0**5*B0*C1**3 - 3*A0**5*B1**3*C0 - 3*A0**5*B1*C0**3 + 10*A0**4*A1**3*B0**2 + 17*A0**4*A1**3*B0*C0 - 27*A0**4*A1**3*C0**2 - 3*A0**4*A1**2*B0**2*B1 - 11*A0**4*A1**2*B0**2*C1 + 9*A0**4*A1**2*B0*B1*C0 - 40*A0**4*A1**2*B0*C0*C1 + 51*A0**4*A1**2*B1*C0**2 - 6*A0**4*A1**2*C0**2*C1 - 3*A0**4*A1*B0**2*B1**2 + 3*A0**4*A1*B0**2*C0**2 + 7*A0**4*A1*B0**2*C1**2 + 9*A0**4*A1*B0*B1**2*C0 - 3*A0**4*A1*B0*C0**3 + 8*A0**4*A1*B0*C0*C1**2 - 21*A0**4*A1*B1**2*C0**2 - 3*A0**4*B0**4*B1 + 3*A0**4*B0**4*C1 + A0**4*B0**2*B1**3 - 3*A0**4*B0**2*C0**2*C1 - A0**4*B0**2*C1**3 - 3*A0**4*B0*B1**3*C0 + 3*A0**4*B0*B1*C0**3 + 3*A0**4*B1**3*C0**2 + 5*A0**3*A1**3*B0**3 - 19*A0**3*A1**3*B0**2*C0 + 3*A0**3*A1**3*B0*C0**2 + 11*A0**3*A1**3*C0**3 + A0**3*A1**2*B0**3*B1 - 4*A0**3*A1**2*B0**3*C1 + 9*A0**3*A1**2*B0**2*B1*C0 + 8*A0**3*A1**2*B0**2*C0*C1 - 9*A0**3*A1**2*B0*B1*C0**2 + 12*A0**3*A1**2*B0*C0**2*C1 - 17*A0**3*A1**2*B1*C0**3 + 3*A0**3*A1*B0**5 - 3*A0**3*A1*B0**4*C0 + A0**3*A1*B0**3*B1**2 - 4*A0**3*A1*B0**3*C1**2 + 9*A0**3*A1*B0**2*B1**2*C0 - 4*A0**3*A1*B0**2*C0*C1**2 - 9*A0**3*A1*B0*B1**2*C0**2 + 7*A0**3*A1*B1**2*C0**3 + 3*A0**3*B0**5*B1 - 6*A0**3*B0**5*C1 + 9*A0**3*B0**4*B1*C0 - 6*A0**3*B0**4*C0*C1 + A0**3*B0**3*B1**3 - 3*A0**3*B0**2*B1**3*C0 + 3*A0**3*B0*B1**3*C0**2 - A0**3*B1**3*C0**3 + A0**2*A1**3*B0**4 - 11*A0**2*A1**3*B0**3*C0 + 15*A0**2*A1**3*B0**2*C0**2 - 5*A0**2*A1**3*B0*C0**3 - 7*A0**2*A1**2*B0**4*C1 - 3*A0**2*A1**2*B0**3*B1*C0 + 16*A0**2*A1**2*B0**3*C0*C1 - 9*A0**2*A1**2*B0**2*B1*C0**2 + 3*A0**2*A1**2*B0*B1*C0**3 - 3*A0**2*A1*B0**6 - 3*A0**2*A1*B0**5*C0 + 6*A0**2*A1*B0**4*C0**2 + 5*A0**2*A1*B0**4*C1**2 - 3*A0**2*A1*B0**3*B1**2*C0 + 4*A0**2*A1*B0**3*C0*C1**2 - 9*A0**2*A1*B0**2*B1**2*C0**2 + 3*A0**2*A1*B0*B1**2*C0**3 + 3*A0**2*B0**6*C1 - 9*A0**2*B0**5*B1*C0 + 12*A0**2*B0**5*C0*C1 - 9*A0**2*B0**4*B1*C0**2 + 3*A0**2*B0**4*C0**2*C1 + A0**2*B0**4*C1**3 - 3*A0**2*B0**3*B1**3*C0 + 3*A0**2*B0**2*B1**3*C0**2 - A0**2*B0*B1**3*C0**3 - 2*A0*A1**3*B0**5 + 4*A0*A1**3*B0**4*C0 + 3*A0*A1**3*B0**3*C0**2 - 5*A0*A1**3*B0**2*C0**3 + 2*A0*A1**2*B0**5*C1 + 4*A0*A1**2*B0**4*C0*C1 + 3*A0*A1**2*B0**3*B1*C0**2 - 12*A0*A1**2*B0**3*C0**2*C1 + 3*A0*A1**2*B0**2*B1*C0**3 + 6*A0*A1*B0**6*C0 - 3*A0*A1*B0**5*C0**2 + 2*A0*A1*B0**5*C1**2 - 3*A0*A1*B0**4*C0**3 - 8*A0*A1*B0**4*C0*C1**2 + 3*A0*A1*B0**3*B1**2*C0**2 + 3*A0*A1*B0**2*B1**2*C0**3 - 6*A0*B0**6*C0*C1 + 9*A0*B0**5*B1*C0**2 - 6*A0*B0**5*C0**2*C1 - 2*A0*B0**5*C1**3 + 3*A0*B0**4*B1*C0**3 + 3*A0*B0**3*B1**3*C0**2 - A0*B0**2*B1**3*C0**3 - A1**3*B0**6 + 4*A1**3*B0**5*C0 - 6*A1**3*B0**4*C0**2 + 3*A1**3*B0**3*C0**3 + 3*A1**2*B0**6*C1 - 8*A1**2*B0**5*C0*C1 + 6*A1**2*B0**4*C0**2*C1 - A1**2*B0**3*B1*C0**3 - 3*A1*B0**6*C0**2 - 3*A1*B0**6*C1**2 + 3*A1*B0**5*C0**3 + 4*A1*B0**5*C0*C1**2 - A1*B0**3*B1**2*C0**3 + 3*B0**6*C0**2*C1 + B0**6*C1**3 - 3*B0**5*B1*C0**3 - B0**3*B1**3*C0**3))/(12*(A0*B0 - A0*C0 - B0*C0 + C0**2)*(A0**2*B0**2 - 2*A0**2*B0*C0 + A0**2*C0**2 - 2*A0*B0**2*C0 + 4*A0*B0*C0**2 - 2*A0*C0**3 + B0**2*C0**2 - 2*B0*C0**3 + C0**4)*(A0**5 - 2*A0**4*B0 - 3*A0**4*C0 + A0**3*B0**2 + 6*A0**3*B0*C0 + 3*A0**3*C0**2 - 3*A0**2*B0**2*C0 - 6*A0**2*B0*C0**2 - A0**2*C0**3 + 3*A0*B0**2*C0**2 + 2*A0*B0*C0**3 - B0**2*C0**3)*(A0**3*B0**3 - 3*A0**3*B0**2*C0 + 3*A0**3*B0*C0**2 - A0**3*C0**3 - 3*A0**2*B0**3*C0 + 9*A0**2*B0**2*C0**2 - 9*A0**2*B0*C0**3 + 3*A0**2*C0**4 + 3*A0*B0**3*C0**2 - 9*A0*B0**2*C0**3 + 9*A0*B0*C0**4 - 3*A0*C0**5 - B0**3*C0**3 + 3*B0**2*C0**4 - 3*B0*C0**5 + C0**6)))

    IyySPREADLambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,mu],
                                                 -(mu*(3*A0**3*B1 - 3*A0**3*C1 - 3*A0**2*A1*B0 + 3*A0**2*A1*C0 + 3*A0**2*B0*B1 - 3*A0**2*C0*C1 - 3*A0*A1*B0**2 + 3*A0*A1*C0**2 + 3*A0*B0**2*B1 + 6*A0*B1*z**2 - 3*A0*C0**2*C1 - 6*A0*C1*z**2 + A1*B0**3 + 6*A1*B0*z**2 - A1*C0**3 - 6*A1*C0*z**2 - 4*B0**3*B1 + 3*B0**3*C1 - 3*B0**2*B1*C0 + 3*B0**2*C0*C1 - 3*B0*B1*C0**2 - 12*B0*B1*z**2 + 3*B0*C0**2*C1 + 6*B0*C1*z**2 + B1*C0**3 + 6*B1*C0*z**2)/12))
    IxxSPREADLambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,mu],
                                                 mu*(-17*A0**5*A1**2*B0**2*B1 + 17*A0**5*A1**2*B0**2*C1 + 34*A0**5*A1**2*B0*B1*C0 - 34*A0**5*A1**2*B0*C0*C1 - 17*A0**5*A1**2*B1*C0**2 + 17*A0**5*A1**2*C0**2*C1 + 7*A0**5*A1*B0**2*B1**2 - 7*A0**5*A1*B0**2*C1**2 - 14*A0**5*A1*B0*B1**2*C0 + 14*A0**5*A1*B0*C0*C1**2 + 7*A0**5*A1*B1**2*C0**2 - 7*A0**5*A1*C0**2*C1**2 - A0**5*B0**2*B1**3 - 6*A0**5*B0**2*B1*z**2 + A0**5*B0**2*C1**3 + 6*A0**5*B0**2*C1*z**2 + 2*A0**5*B0*B1**3*C0 + 12*A0**5*B0*B1*C0*z**2 - 2*A0**5*B0*C0*C1**3 - 12*A0**5*B0*C0*C1*z**2 - A0**5*B1**3*C0**2 - 6*A0**5*B1*C0**2*z**2 + A0**5*C0**2*C1**3 + 6*A0**5*C0**2*C1*z**2 + 13*A0**4*A1**3*B0**3 - 39*A0**4*A1**3*B0**2*C0 + 39*A0**4*A1**3*B0*C0**2 - 13*A0**4*A1**3*C0**3 + 3*A0**4*A1**2*B0**3*B1 - 34*A0**4*A1**2*B0**3*C1 + 28*A0**4*A1**2*B0**2*B1*C0 + 65*A0**4*A1**2*B0**2*C0*C1 - 65*A0**4*A1**2*B0*B1*C0**2 - 28*A0**4*A1**2*B0*C0**2*C1 + 34*A0**4*A1**2*B1*C0**3 - 3*A0**4*A1**2*C0**3*C1 + 3*A0**4*A1*B0**3*B1**2 + 14*A0**4*A1*B0**3*C1**2 - 6*A0**4*A1*B0**3*z**2 - 20*A0**4*A1*B0**2*B1**2*C0 - 31*A0**4*A1*B0**2*C0*C1**2 + 18*A0**4*A1*B0**2*C0*z**2 + 31*A0**4*A1*B0*B1**2*C0**2 + 20*A0**4*A1*B0*C0**2*C1**2 - 18*A0**4*A1*B0*C0**2*z**2 - 14*A0**4*A1*B1**2*C0**3 - 3*A0**4*A1*C0**3*C1**2 + 6*A0**4*A1*C0**3*z**2 + 14*A0**4*B0**3*B1**3 - 17*A0**4*B0**3*B1**2*C1 + 7*A0**4*B0**3*B1*C1**2 + 24*A0**4*B0**3*B1*z**2 - 3*A0**4*B0**3*C1**3 - 18*A0**4*B0**3*C1*z**2 - 13*A0**4*B0**2*B1**3*C0 + 3*A0**4*B0**2*B1**2*C0*C1 + 3*A0**4*B0**2*B1*C0*C1**2 - 42*A0**4*B0**2*B1*C0*z**2 + 4*A0**4*B0**2*C0*C1**3 + 24*A0**4*B0**2*C0*C1*z**2 + 2*A0**4*B0*B1**3*C0**2 + 3*A0**4*B0*B1**2*C0**2*C1 + 3*A0**4*B0*B1*C0**2*C1**2 + 12*A0**4*B0*B1*C0**2*z**2 - 5*A0**4*B0*C0**2*C1**3 + 6*A0**4*B0*C0**2*C1*z**2 + A0**4*B1**3*C0**3 - A0**4*B1**2*C0**3*C1 - A0**4*B1*C0**3*C1**2 + 6*A0**4*B1*C0**3*z**2 - 12*A0**4*C0**3*C1*z**2 - 8*A0**3*A1**3*B0**4 + 16*A0**3*A1**3*B0**3*C0 - 16*A0**3*A1**3*B0*C0**3 + 8*A0**3*A1**3*C0**4 + 3*A0**3*A1**2*B0**4*B1 + 17*A0**3*A1**2*B0**4*C1 - 12*A0**3*A1**2*B0**3*B1*C0 - 28*A0**3*A1**2*B0**3*C0*C1 - 2*A0**3*A1**2*B0**2*B1*C0**2 + 2*A0**3*A1**2*B0**2*C0**2*C1 + 28*A0**3*A1**2*B0*B1*C0**3 + 12*A0**3*A1**2*B0*C0**3*C1 - 17*A0**3*A1**2*B1*C0**4 - 3*A0**3*A1**2*C0**4*C1 + 3*A0**3*A1*B0**4*B1**2 - 7*A0**3*A1*B0**4*C1**2 + 12*A0**3*A1*B0**4*z**2 - 12*A0**3*A1*B0**3*B1**2*C0 + 20*A0**3*A1*B0**3*C0*C1**2 - 24*A0**3*A1*B0**3*C0*z**2 + 22*A0**3*A1*B0**2*B1**2*C0**2 - 22*A0**3*A1*B0**2*C0**2*C1**2 - 20*A0**3*A1*B0*B1**2*C0**3 + 12*A0**3*A1*B0*C0**3*C1**2 + 24*A0**3*A1*B0*C0**3*z**2 + 7*A0**3*A1*B1**2*C0**4 - 3*A0**3*A1*C0**4*C1**2 - 12*A0**3*A1*C0**4*z**2 - 31*A0**3*B0**4*B1**3 + 34*A0**3*B0**4*B1**2*C1 - 14*A0**3*B0**4*B1*C1**2 - 30*A0**3*B0**4*B1*z**2 + 3*A0**3*B0**4*C1**3 + 18*A0**3*B0**4*C1*z**2 + 8*A0**3*B0**3*B1**3*C0 + 28*A0**3*B0**3*B1**2*C0*C1 - 20*A0**3*B0**3*B1*C0*C1**2 + 24*A0**3*B0**3*B1*C0*z**2 + 14*A0**3*B0**2*B1**3*C0**2 - 12*A0**3*B0**2*B1**2*C0**2*C1 - 12*A0**3*B0**2*B1*C0**2*C1**2 + 48*A0**3*B0**2*B1*C0**2*z**2 + 10*A0**3*B0**2*C0**2*C1**3 - 48*A0**3*B0**2*C0**2*C1*z**2 - 8*A0**3*B0*B1**3*C0**3 - 4*A0**3*B0*B1**2*C0**3*C1 - 4*A0**3*B0*B1*C0**3*C1**2 - 48*A0**3*B0*B1*C0**3*z**2 + 24*A0**3*B0*C0**3*C1*z**2 + A0**3*B1**3*C0**4 + 2*A0**3*B1**2*C0**4*C1 + 2*A0**3*B1*C0**4*C1**2 + 6*A0**3*B1*C0**4*z**2 + 3*A0**3*C0**4*C1**3 + 6*A0**3*C0**4*C1*z**2 - A0**2*A1**3*B0**5 + 5*A0**2*A1**3*B0**4*C0 - 10*A0**2*A1**3*B0**3*C0**2 + 10*A0**2*A1**3*B0**2*C0**3 - 5*A0**2*A1**3*B0*C0**4 + A0**2*A1**3*C0**5 - A0**2*A1**2*B0**5*B1 - 4*A0**2*A1**2*B0**4*B1*C0 - 3*A0**2*A1**2*B0**4*C0*C1 + 14*A0**2*A1**2*B0**3*B1*C0**2 + 12*A0**2*A1**2*B0**3*C0**2*C1 - 12*A0**2*A1**2*B0**2*B1*C0**3 - 14*A0**2*A1**2*B0**2*C0**3*C1 + 3*A0**2*A1**2*B0*B1*C0**4 + 4*A0**2*A1**2*B0*C0**4*C1 + A0**2*A1**2*C0**5*C1 - A0**2*A1*B0**5*B1**2 - 6*A0**2*A1*B0**5*z**2 - 4*A0**2*A1*B0**4*B1**2*C0 - 3*A0**2*A1*B0**4*C0*C1**2 - 6*A0**2*A1*B0**4*C0*z**2 + 14*A0**2*A1*B0**3*B1**2*C0**2 + 12*A0**2*A1*B0**3*C0**2*C1**2 + 48*A0**2*A1*B0**3*C0**2*z**2 - 12*A0**2*A1*B0**2*B1**2*C0**3 - 14*A0**2*A1*B0**2*C0**3*C1**2 - 48*A0**2*A1*B0**2*C0**3*z**2 + 3*A0**2*A1*B0*B1**2*C0**4 + 4*A0**2*A1*B0*C0**4*C1**2 + 6*A0**2*A1*B0*C0**4*z**2 + A0**2*A1*C0**5*C1**2 + 6*A0**2*A1*C0**5*z**2 + 14*A0**2*B0**5*B1**3 - 17*A0**2*B0**5*B1**2*C1 + 7*A0**2*B0**5*B1*C1**2 + 12*A0**2*B0**5*B1*z**2 - A0**2*B0**5*C1**3 - 6*A0**2*B0**5*C1*z**2 + 47*A0**2*B0**4*B1**3*C0 - 65*A0**2*B0**4*B1**2*C0*C1 + 31*A0**2*B0**4*B1*C0*C1**2 + 30*A0**2*B0**4*B1*C0*z**2 - 4*A0**2*B0**4*C0*C1**3 - 24*A0**2*B0**4*C0*C1*z**2 - 52*A0**2*B0**3*B1**3*C0**2 - 2*A0**2*B0**3*B1**2*C0**2*C1 + 22*A0**2*B0**3*B1*C0**2*C1**2 - 96*A0**2*B0**3*B1*C0**2*z**2 - 10*A0**2*B0**3*C0**2*C1**3 + 48*A0**2*B0**3*C0**2*C1*z**2 + 14*A0**2*B0**2*B1**3*C0**3 + 14*A0**2*B0**2*B1**2*C0**3*C1 + 14*A0**2*B0**2*B1*C0**3*C1**2 + 48*A0**2*B0**2*B1*C0**3*z**2 + 2*A0**2*B0*B1**3*C0**4 - A0**2*B0*B1**2*C0**4*C1 - A0**2*B0*B1*C0**4*C1**2 + 12*A0**2*B0*B1*C0**4*z**2 - 9*A0**2*B0*C0**4*C1**3 - 18*A0**2*B0*C0**4*C1*z**2 - A0**2*B1**3*C0**5 - A0**2*B1**2*C0**5*C1 - A0**2*B1*C0**5*C1**2 - 6*A0**2*B1*C0**5*z**2 + 2*A0*A1**3*B0**5*C0 - 4*A0*A1**3*B0**4*C0**2 + 4*A0*A1**3*B0**2*C0**4 - 2*A0*A1**3*B0*C0**5 + 2*A0*A1**2*B0**5*B1*C0 - A0*A1**2*B0**4*B1*C0**2 - 3*A0*A1**2*B0**4*C0**2*C1 - 4*A0*A1**2*B0**3*B1*C0**3 + 4*A0*A1**2*B0**3*C0**3*C1 + 3*A0*A1**2*B0**2*B1*C0**4 + A0*A1**2*B0**2*C0**4*C1 - 2*A0*A1**2*B0*C0**5*C1 + 2*A0*A1*B0**5*B1**2*C0 + 12*A0*A1*B0**5*C0*z**2 - A0*A1*B0**4*B1**2*C0**2 - 3*A0*A1*B0**4*C0**2*C1**2 - 24*A0*A1*B0**4*C0**2*z**2 - 4*A0*A1*B0**3*B1**2*C0**3 + 4*A0*A1*B0**3*C0**3*C1**2 + 3*A0*A1*B0**2*B1**2*C0**4 + A0*A1*B0**2*C0**4*C1**2 + 24*A0*A1*B0**2*C0**4*z**2 - 2*A0*A1*B0*C0**5*C1**2 - 12*A0*A1*B0*C0**5*z**2 - 28*A0*B0**5*B1**3*C0 + 34*A0*B0**5*B1**2*C0*C1 - 14*A0*B0**5*B1*C0*C1**2 - 24*A0*B0**5*B1*C0*z**2 + 2*A0*B0**5*C0*C1**3 + 12*A0*B0**5*C0*C1*z**2 - A0*B0**4*B1**3*C0**2 + 28*A0*B0**4*B1**2*C0**2*C1 - 20*A0*B0**4*B1*C0**2*C1**2 + 30*A0*B0**4*B1*C0**2*z**2 + 5*A0*B0**4*C0**2*C1**3 - 6*A0*B0**4*C0**2*C1*z**2 + 24*A0*B0**3*B1**3*C0**3 - 12*A0*B0**3*B1**2*C0**3*C1 - 12*A0*B0**3*B1*C0**3*C1**2 + 24*A0*B0**3*B1*C0**3*z**2 - 24*A0*B0**3*C0**3*C1*z**2 - 13*A0*B0**2*B1**3*C0**4 - 4*A0*B0**2*B1**2*C0**4*C1 - 4*A0*B0**2*B1*C0**4*C1**2 - 42*A0*B0**2*B1*C0**4*z**2 + 9*A0*B0**2*C0**4*C1**3 + 18*A0*B0**2*C0**4*C1*z**2 + 2*A0*B0*B1**3*C0**5 + 2*A0*B0*B1**2*C0**5*C1 + 2*A0*B0*B1*C0**5*C1**2 + 12*A0*B0*B1*C0**5*z**2 - A1**3*B0**5*C0**2 + 3*A1**3*B0**4*C0**3 - 3*A1**3*B0**3*C0**4 + A1**3*B0**2*C0**5 - A1**2*B0**5*B1*C0**2 + 2*A1**2*B0**4*B1*C0**3 + A1**2*B0**4*C0**3*C1 - A1**2*B0**3*B1*C0**4 - 2*A1**2*B0**3*C0**4*C1 + A1**2*B0**2*C0**5*C1 - A1*B0**5*B1**2*C0**2 - 6*A1*B0**5*C0**2*z**2 + 2*A1*B0**4*B1**2*C0**3 + A1*B0**4*C0**3*C1**2 + 18*A1*B0**4*C0**3*z**2 - A1*B0**3*B1**2*C0**4 - 2*A1*B0**3*C0**4*C1**2 - 18*A1*B0**3*C0**4*z**2 + A1*B0**2*C0**5*C1**2 + 6*A1*B0**2*C0**5*z**2 + 14*B0**5*B1**3*C0**2 - 17*B0**5*B1**2*C0**2*C1 + 7*B0**5*B1*C0**2*C1**2 + 12*B0**5*B1*C0**2*z**2 - B0**5*C0**2*C1**3 - 6*B0**5*C0**2*C1*z**2 - 15*B0**4*B1**3*C0**3 + 3*B0**4*B1**2*C0**3*C1 + 3*B0**4*B1*C0**3*C1**2 - 30*B0**4*B1*C0**3*z**2 + 12*B0**4*C0**3*C1*z**2 + 6*B0**3*B1**3*C0**4 + 3*B0**3*B1**2*C0**4*C1 + 3*B0**3*B1*C0**4*C1**2 + 24*B0**3*B1*C0**4*z**2 - 3*B0**3*C0**4*C1**3 - 6*B0**3*C0**4*C1*z**2 - B0**2*B1**3*C0**5 - B0**2*B1**2*C0**5*C1 - B0**2*B1*C0**5*C1**2 - 6*B0**2*B1*C0**5*z**2))

    #denom = /(12*(A0**4*B0**2 - 2*A0**4*B0*C0 + A0**4*C0**2 - 2*A0**3*B0**3 + 2*A0**3*B0**2*C0 + 2*A0**3*B0*C0**2 - 2*A0**3*C0**3 + A0**2*B0**4 + 2*A0**2*B0**3*C0 - 6*A0**2*B0**2*C0**2 + 2*A0**2*B0*C0**3 + A0**2*C0**4 - 2*A0*B0**4*C0 + 2*A0*B0**3*C0**2 + 2*A0*B0**2*C0**3 - 2*A0*B0*C0**4 + B0**4*C0**2 - 2*B0**3*C0**3 + B0**2*C0**4)))
    #these are already multiplied by -1
    IxySPREADLambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,mu],
                                    mu*(14*A0**4*A1*B0*B1 - 14*A0**4*A1*B0*C1 - 14*A0**4*A1*B1*C0 + 14*A0**4*A1*C0*C1 - 3*A0**4*B0*B1**2 + 3*A0**4*B0*C1**2 + 3*A0**4*B1**2*C0 - 3*A0**4*C0*C1**2 - 14*A0**3*A1**2*B0**2 + 28*A0**3*A1**2*B0*C0 - 14*A0**3*A1**2*C0**2 + 6*A0**3*A1*B0**2*B1 + 14*A0**3*A1*B0**2*C1 - 20*A0**3*A1*B0*B1*C0 - 20*A0**3*A1*B0*C0*C1 + 14*A0**3*A1*B1*C0**2 + 6*A0**3*A1*C0**2*C1 - 3*A0**3*B0**2*B1**2 - 3*A0**3*B0**2*C1**2 + 6*A0**3*B0*B1**2*C0 + 6*A0**3*B0*C0*C1**2 - 3*A0**3*B1**2*C0**2 - 3*A0**3*C0**2*C1**2 + 3*A0**2*A1**2*B0**3 - 3*A0**2*A1**2*B0**2*C0 - 3*A0**2*A1**2*B0*C0**2 + 3*A0**2*A1**2*C0**3 + 6*A0**2*A1*B0**3*B1 - 12*A0**2*A1*B0**2*B1*C0 + 6*A0**2*A1*B0**2*C0*C1 + 6*A0**2*A1*B0*B1*C0**2 - 12*A0**2*A1*B0*C0**2*C1 + 6*A0**2*A1*C0**3*C1 - 20*A0**2*B0**3*B1**2 + 14*A0**2*B0**3*B1*C1 - 3*A0**2*B0**3*C1**2 + 9*A0**2*B0**2*B1**2*C0 + 6*A0**2*B0**2*B1*C0*C1 - 6*A0**2*B0**2*C0*C1**2 + 6*A0**2*B0*B1*C0**2*C1 + 3*A0**2*B0*C0**2*C1**2 - A0**2*B1**2*C0**3 - 2*A0**2*B1*C0**3*C1 - 6*A0**2*C0**3*C1**2 - A0*A1**2*B0**4 - 2*A0*A1**2*B0**3*C0 + 6*A0*A1**2*B0**2*C0**2 - 2*A0*A1**2*B0*C0**3 - A0*A1**2*C0**4 - 2*A0*A1*B0**4*B1 - 4*A0*A1*B0**3*B1*C0 + 6*A0*A1*B0**2*B1*C0**2 + 6*A0*A1*B0**2*C0**2*C1 - 4*A0*A1*B0*C0**3*C1 - 2*A0*A1*C0**4*C1 + 14*A0*B0**4*B1**2 - 14*A0*B0**4*B1*C1 + 3*A0*B0**4*C1**2 + 20*A0*B0**3*B1**2*C0 - 20*A0*B0**3*B1*C0*C1 + 6*A0*B0**3*C0*C1**2 - 9*A0*B0**2*B1**2*C0**2 - 12*A0*B0**2*B1*C0**2*C1 + 3*A0*B0**2*C0**2*C1**2 - 2*A0*B0*B1**2*C0**3 - 4*A0*B0*B1*C0**3*C1 + 12*A0*B0*C0**3*C1**2 + A0*B1**2*C0**4 + 2*A0*B1*C0**4*C1 + A1**2*B0**4*C0 - A1**2*B0**3*C0**2 - A1**2*B0**2*C0**3 + A1**2*B0*C0**4 + 2*A1*B0**4*B1*C0 - 2*A1*B0**3*B1*C0**2 - 2*A1*B0**2*C0**3*C1 + 2*A1*B0*C0**4*C1 - 14*B0**4*B1**2*C0 + 14*B0**4*B1*C0*C1 - 3*B0**4*C0*C1**2 + 6*B0**3*B1*C0**2*C1 - 3*B0**3*C0**2*C1**2 + 3*B0**2*B1**2*C0**3 + 6*B0**2*B1*C0**3*C1 - 6*B0**2*C0**3*C1**2 - B0*B1**2*C0**4 - 2*B0*B1*C0**4*C1))#/(24*(A0**2*B0 - A0**2*C0 - A0*B0**2 + A0*C0**2 + B0**2*C0 - B0*C0**2)))
    IxzSPREADLambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,mu],
                        -(mu*z*(-2*A0**2*B1 + 2*A0**2*C1 + 2*A0*A1*B0 - 2*A0*A1*C0 - 2*A0*B0*B1 + 2*A0*C0*C1 - A1*B0**2 + A1*C0**2 + 3*B0**2*B1 - 2*B0**2*C1 + 2*B0*B1*C0 - 2*B0*C0*C1 - B1*C0**2)/6))
    IzySPREADLambdify = sympy.utilities.lambdify([A0,A1,A2,B0,B1,B2,C0,C1,C2,mu],
                    mu*z*(-5*A0**3*A1*B0*B1 + 5*A0**3*A1*B0*C1 + 5*A0**3*A1*B1*C0 - 5*A0**3*A1*C0*C1 + A0**3*B0*B1**2 - A0**3*B0*C1**2 - A0**3*B1**2*C0 + A0**3*C0*C1**2 + 2*A0**2*A1**2*B0**2 - 4*A0**2*A1**2*B0*C0 + 2*A0**2*A1**2*C0**2 - 2*A0**2*A1*B0**2*B1 - 5*A0**2*A1*B0**2*C1 + 7*A0**2*A1*B0*B1*C0 + 7*A0**2*A1*B0*C0*C1 - 5*A0**2*A1*B1*C0**2 - 2*A0**2*A1*C0**2*C1 + 8*A0**2*B0**2*B1**2 - 5*A0**2*B0**2*B1*C1 + 2*A0**2*B0**2*C1**2 - 7*A0**2*B0*B1**2*C0 - 2*A0**2*B0*B1*C0*C1 - A0**2*B0*C0*C1**2 + 2*A0**2*B1**2*C0**2 + A0**2*B1*C0**2*C1 + 2*A0**2*C0**2*C1**2 + A0*A1**2*B0**3 - A0*A1**2*B0**2*C0 - A0*A1**2*B0*C0**2 + A0*A1**2*C0**3 + A0*A1*B0**3*B1 + A0*A1*B0**2*B1*C0 - 2*A0*A1*B0**2*C0*C1 - 2*A0*A1*B0*B1*C0**2 + A0*A1*B0*C0**2*C1 + A0*A1*C0**3*C1 - 6*A0*B0**3*B1**2 + 5*A0*B0**3*B1*C1 - A0*B0**3*C1**2 - 4*A0*B0**2*B1**2*C0 + 7*A0*B0**2*B1*C0*C1 - A0*B0**2*C0*C1**2 + 5*A0*B0*B1**2*C0**2 + A0*B0*B1*C0**2*C1 - 4*A0*B0*C0**2*C1**2 - A0*B1**2*C0**3 - A0*B1*C0**3*C1 - A1**2*B0**3*C0 + 2*A1**2*B0**2*C0**2 - A1**2*B0*C0**3 - A1*B0**3*B1*C0 + A1*B0**2*B1*C0**2 + A1*B0**2*C0**2*C1 - A1*B0*C0**3*C1 + 6*B0**3*B1**2*C0 - 5*B0**3*B1*C0*C1 + B0**3*C0*C1**2 - 4*B0**2*B1**2*C0**2 - 2*B0**2*B1*C0**2*C1 + 2*B0**2*C0**2*C1**2 + B0*B1**2*C0**3 + B0*B1*C0**3*C1)/(6*(A0**2*B0 - A0**2*C0 - A0*B0**2 + A0*C0**2 + B0**2*C0 - B0*C0**2)))


    def JtensorSPREAD(A,B,C,spreadDensity,spreadThickness):

        J = np.zeros((3,3))

        #principle inertias:
        J[2][2] = IzzSPREADLambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],spreadDensity)

        tempXX = IxxSPREADLambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],spreadDensity)
        J[0][0] = tempXX.subs(z,spreadThickness) - tempXX.subs(z,0)
        tempYY = ((IyySPREADLambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],spreadDensity)))
        J[1][1] = tempYY.subs(z,spreadThickness) - tempYY.subs(z,0)

        #off diagonal entries (makes use of symmetric property of the tensor)
        J[0][1] = IxySPREADLambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],spreadDensity)
        J[1][0] = J[0][1]
        tempXZ = IxzSPREADLambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],spreadDensity)
        J[0][2] = tempXZ.subs(z,spreadThickness) - tempXZ.subs(z,0)
        J[2][0] = J[0][2]
        tempZY = IzySPREADLambdify(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],spreadDensity)
        J[1][2] = tempZY.subs(z,spreadThickness) - tempZY.subs(z,0)
        J[2][1] = J[1][2]

        #scale appropriately
        J= J*4

        #ignore floating point errors
        J = np.where(np.abs(J)>10**(-12),J,0)

        return J*spreadThickness

    def getFullMeshInertia(spreadDensity,breadDensity,spreadThickness,spreadFaces):
        J=np.zeros((3,3))
        origin = np.array([args.translate[0],args.translate[1],args.translate[2]])

        #use the previous function to retrieve the inertia tensor for the whole object
        for face in F:
            J+=Jtensor(V[face[0]],V[face[1]],V[face[2]],origin,breadDensity)*(volumeOfTetrahedron(np.array([V[face[0]],V[face[1]],V[face[2]]]))*breadDensity)
        for spreaded in spreadFaces:
            A = np.array(V[face[0]])
            B = np.array(V[face[1]])
            C = np.array(V[face[2]])
            #print(np.cross(A-C,A-B))
            J += JtensorSPREAD(A,B,C,spreadDensity,spreadThickness)*spreadThickness*np.linalg.norm(np.cross(A-C,A-B))*spreadDensity/2
            ##WE WANT TO ADD THE INERTIA TENSORS FROM BUTTER TO INERTIA TENSOR FROM BREAD
        return J


    # In[1151]:


    def getMassMatrix(MeshCOM,mass,meshJ):
        MassMatrix = np.zeros((6,6))

        #fill out M*Id block
        MassMatrix[0][0]=1
        MassMatrix[1][1]=1
        MassMatrix[2][2]=1

        #fill out off diagonal blocks
        #making use of skew symmetric property

        #x-component
        MassMatrix[1][5]=MeshCOM[0]
        MassMatrix[2][4]=-MeshCOM[0]
        MassMatrix[5][1]=MeshCOM[0]
        MassMatrix[4][2]=-MeshCOM[0]

        #y-component
        MassMatrix[0][5]=-MeshCOM[1]
        MassMatrix[5][0]=-MeshCOM[1]
        MassMatrix[3][2]=MeshCOM[1]
        MassMatrix[2][3]=MeshCOM[1]

        #z-component
        MassMatrix[0][4]=MeshCOM[2]
        MassMatrix[4][0]=MeshCOM[2]
        MassMatrix[1][3]=-MeshCOM[2]
        MassMatrix[3][1]=-MeshCOM[2]

        #fill out lower diagonal block (inertia tensor spot)

        for row in range(0,3):
            for column in range(0,3):

                #divide by mass because we're about to multiply the whole matrix by the mass
                MassMatrix[row+3][column+3] = meshJ[row][column]/(mass)

        MassMatrix=MassMatrix*mass
        MassMatrix=np.round(MassMatrix,decimals = 3)

        return MassMatrix





    # In[1152]:


    def getInitialInfo(breadDensity,spreadDensity,spreadThickness,spreadFaces):
        CM = getObjectCOM(breadDensity)
        Vertices = V
        Faces = F
        Volume = meshVolume(breadDensity)
        I = getFullMeshInertia(spreadDensity,breadDensity,spreadThickness,spreadFaces)
        return [CM,Vertices,Faces,Volume,I]
    return (getInitialInfo(breadDensity,spreadDensity,spreadThickness,spreadFaces))
    #print(getInitialInfo(1,1,1,[10]))

    # In[ ]:





    # In[ ]:





