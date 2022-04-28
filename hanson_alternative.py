import math
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy import interpolate


class Flow:
    def __init__(self, density, velocity, a_0):
        self.density = density
        self.velocity = velocity
        self.a_0 = a_0
        self.M_x = velocity / a_0


class Propeller:
    def __init__(self, flow, RPM, bladeNumber, diameter, x, theta, y):
        self.RPM = RPM
        self.flow = flow
        self.bladeNumber = bladeNumber
        self.diameter = diameter
        self.radius = 0.5 * diameter
        self.x = x
        self.theta = theta #np.acos((cos(theta)*(1-self.M_x*self.M_x*sin(theta)*sin(theta))**0.5)+self.M_x*sin(theta)*sin(theta)
        self.M_t = (self.RPM * 2 * math.pi / 60 * diameter / 2) / flow.a_0
        self.bladePassingFrequency = (RPM / 60) * bladeNumber
        self.y = y
        self.thicknesstoCord = 0.12
        self.dopplerFactor = 1 #(1 - self.flow.M_x * math.cos(self.theta))

        with open(r"numerical results/Liftpoints.csv") as file:
            self.CL = np.genfromtxt(file, delimiter=',', dtype=float)
            x = np.arange(0.2267,1,(1-0.2267)/len(self.CL))
            self.Lift = sc.interpolate.interp1d(x, self.CL, kind='cubic',fill_value="extrapolate")
        with open(r"numerical results/Dragpoints.csv") as file:
            self.CD = np.genfromtxt(file, delimiter=',', dtype=float)
            x = np.arange(0.2267, 1, (1-0.2267) / len(self.CD))
            self.Drag = sc.interpolate.interp1d(x, self.CD, kind='cubic',fill_value="extrapolate")

        with open(r"numerical results/CP Interpolation/mChord.csv") as chordCoordsFile:
            chordCoords = np.genfromtxt(chordCoordsFile, delimiter=',',dtype=float)
            with open(r"numerical results/CP Interpolation/mSpan.csv") as spanCoordsFile:
                spanCoords = np.genfromtxt(spanCoordsFile, delimiter=',',dtype=float)
                with open(r"numerical results/CP Interpolation/mDown.csv") as downPressureFile:
                    downPressures = np.genfromtxt(downPressureFile, delimiter=',', dtype=float)
                    with open(r"numerical results/CP Interpolation/mUp.csv") as upPressureFile:
                        upPressures = np.genfromtxt(upPressureFile, delimiter=',', dtype=float)
                        #print(downPressures)
                        chordCoordArr2=chordCoords[0]
                        spanCoordArr2=spanCoords[:,0]

        self.Cpl = sc.interpolate.interp2d(chordCoordArr2, spanCoordArr2, downPressures, kind='cubic')
        self.Cpu = sc.interpolate.interp2d(chordCoordArr2, spanCoordArr2, upPressures, kind='cubic')

    # propeller properties
    def B_D(self, z):  # chord to diameter ratio
        r_o = z - 0.15
        if (r_o < 0.1):
            return 0
        B_D_o = 0.127 + 0.7404112848810135 * (r_o ** 1) + -4.673821287763687 * (r_o ** 2) + 19.02690171288593 * (
                    r_o ** 3) + -48.12645454430276 * (r_o ** 4) + 36.20381779945387 * (r_o ** 5) + 111.0096041870429 * (
                            r_o ** 6) + -311.5856431939494 * (r_o ** 7) + 309.1780229926285 * (
                            r_o ** 8) + -113.64512633245806 * (r_o ** 9)
        B_D = 0.5 * (B_D_o + 0.05)
        return B_D

    def twist(self,x):
        twist = -8971.329845161235 + 206093.7010106019 * x - 2.0708079035287888e6 * x ** 2 + 1.2100709270279987e7 * x ** 3 - 4.5793508020059824e7 * x ** 4 + 1.1806041107129405e8 * x ** 5 - 2.1196809140788797e8 * x ** 6 + 2.6550027408285195e8 * x ** 7 - 2.2773507697170436e8 * x ** 8 + 1.2759430429713346e8 * x ** 9 - 4.208253179063628e7 * x ** 10 + 6.197207001092323e6 * x ** 11
        return twist

    def angle_of_attack(self,z):#different from that group but not that different
        return (self.twist(z)-(np.arctan(self.flow.M_x**2/(self.M_r(z)**2))*180/np.pi))*np.pi/180

    def chord(self, z):
        c = self.B_D(z) * self.diameter
        return c

    def M_r(self, z):
        # print(self.M_t)
        return np.sqrt((self.flow.M_x ** 2) + (z ** 2 * self.M_t ** 2))

    # calculate the first constant in the large formula
    def calcFirstPart(self, m):
        omega = (self.RPM / 60) * 2 * math.pi
        # print((self.flow.density * (self.flow.a_0**2) * self.bladeNumber * math.sin(self.theta)))
        return -(self.flow.density * (self.flow.a_0 ** 2) * self.bladeNumber * math.sin(self.theta) * np.exp(
            1j * m * self.bladeNumber * (((omega * self.radius) / (self.flow.a_0)) - (0.5 * math.pi)))) / (
                           8 * math.pi * (self.y / self.diameter) * self.dopplerFactor)

    # k_x and k_y
    def k_x(self, z, m):
        # print(str(self.bladeNumber) + " " + str(self.B_D(z)) + " " + str(m) + " " + str(self.M_t) + " ")
        k_x = (2 * self.bladeNumber * m * self.B_D(z) * self.M_t) / (self.M_r(z) * self.dopplerFactor)
        # print(k_x)
        return k_x

    def k_y(self, z, m):
        return (2 * self.bladeNumber * m * self.B_D(z) * (self.flow.M_x - self.M_r(z) ** (2) * np.cos(self.theta))) / (
                    z * self.M_r(z) * self.dopplerFactor)

    # The psi's
    # method 1
    def psi_D(self, z, m):
        return math2.integration(self.psi_D_derivative, -0.49999, 0.5, 100, [z, m], 'Simpsons')

    def psi_L(self, z, m):
        #print(math2.integration(self.psi_L_derivative, -0.499999, 0.5, 100, [z, m], 'Simpsons'))
        return math2.integration(self.psi_L_derivative, -0.499999, 0.5, 100, [z, m], 'Simpsons')

    def psi_V(self, z, m):
        return math2.integration(self.psi_V_derivative, -0.499999, 0.5, 100, [z, m], 'Simpsons')

    # method 2
    def psi_D2(self, z, m):
        return math2.integrateRiemannSums(self.psi_D_derivative, -0.49999, 0.5, 100, [z, m])

    def psi_L2(self, z, m):
        return math2.integrateRiemannSums(self.psi_L_derivative, -0.49999, 0.5, 100, [z, m])

    def psi_V2(self, z, m):
        return math2.integrateRiemannSums(self.psi_V_derivative, -0.49999, 0.5, 100, [z, m])

    # Derivative of psi's 🐙🐙🐙🐙🐙🐙
    def psi_D_derivative(self, x, args):
        z = args[0]
        m = args[1]

        k_x = self.k_x(z, m)
        Drag_sec = self.normalisedDragDist(x,z)
        return Drag_sec * np.exp(1j * k_x * x)

    def psi_L_derivative(self, x, args):
        z = args[0]
        m = args[1]

        k_x = self.k_x(z, m)
        Lift_sec = self.normalisedLiftDist(x,z)
        return Lift_sec * np.exp(1j * k_x * x)

    def psi_V_derivative(self, x, args):
        z = args[0]
        m = args[1]
        k_x = self.k_x(z, m)
        return (self.thicknessDistNormalized(x+0.5)[0]-self.thicknessDistNormalized(x+0.5)[1]) * np.exp(1j * k_x * x)

    # p(t)🐙🐙🐙🐙🐙🐙🐙
    def pressure(self, m):
        k = self.calcFirstPart(m)
       # p_Vm = k * math2.integrateRiemannSums(self.p_Vm_derivative, 0, 1, 100, m)
       # p_Dm = k * math2.integrateRiemannSums(self.p_Dm_derivative, 0, 1, 100, m)
       # p_Lm = k * math2.integrateRiemannSums(self.p_Lm_derivative, 0, 1, 100, m)
        p_Vm2 = k * math2.integration(self.p_Vm_derivative, 0.2267, 1, 100, m, 'Simpsons')
        p_Dm2 = k * math2.integration(self.p_Dm_derivative, 0.2267, 1, 100, m, 'Simpsons')
        p_Lm2 = k * math2.integration(self.p_Lm_derivative, 0.2267, 1, 100, m, 'Simpsons')
        #print(p_Vm2,p_Lm2,p_Dm2)
        #p_mb = p_Vm + p_Dm + p_Lm
        p_mb2 = p_Vm2 + p_Dm2 + p_Lm2
        return 2* p_mb2

    # harmonic noise level
    def noise(self, m):
        # p = self.pressure(m)
        p = self.pressure(m)
        p_ref = 20 * 10 ** -6
        # print((p**2) / (p_ref**2))

        OASPL =20 * np.log10((p) / (p_ref))
        return abs(OASPL)

    def p_Vm_derivative(self, z, m):
        k_x = self.k_x(z, m)
        P_Vm_1 = (self.M_r(z)) ** 2 * math2.besselsFunc(m * self.bladeNumber,
                                                        m * self.bladeNumber * z * self.M_t * np.sin(self.theta) / (
                                                        self.dopplerFactor))
        P_Vm_2 = self.k_x(z, m) ** 2 * (self.thicknesstoCord) * self.psi_V(z, m)
        #P_Vm_2_alt = self.k_x(z, m) ** 2 * (self.thicknesstoCord) * self.psi_V2(z, m)
        return P_Vm_1 * P_Vm_2

    def p_Dm_derivative(self, z, m): #🐙🐙🐙🐙🐙🐙🐙
        self.dragArea = self.areaDrag(z)
        k_x = self.k_x(z, m)
        C_D = self.Drag(z)
        P_Dm_1 = (self.M_r(z)) ** 2 * math2.besselsFunc(m * self.bladeNumber,
                                                        m * self.bladeNumber * z * self.M_t * np.sin(self.theta) /
                                                        (self.dopplerFactor))
        P_Dm_2 = 1j * k_x * (C_D / 2) * self.psi_D(z, m)
        P_Dm_2_alt = 1j * k_x * (C_D / 2) * self.psi_D2(z, m)
        return P_Dm_1 * P_Dm_2

    def p_Lm_derivative(self, z, m):
        self.liftArea = self.areaLift(z)
        k_x = self.k_x(z, m)
        C_L = self.Lift(z)
        # print(C_L)
        P_Lm_1 = (self.M_r(z)) ** 2 * math2.besselsFunc(m * self.bladeNumber,
                                                        m * self.bladeNumber * z * self.M_t * np.sin(self.theta) /
                                                        (self.dopplerFactor))
        P_Lm_2 = 1j * k_x * (C_L / 2) * self.psi_L(z, m)
        P_Lm_2_alt = 1j * k_x * C_L / 2 * self.psi_L2(z, m)
        return P_Lm_1 * P_Lm_2



    def thicknessDist(self,x): #😖😖😖
        #NACA 4412 thicknes from x = [0,1]
        #approximate as x coordinates on the x-axis are used instead of the actual surface; may need changing 😖😖😖
        m = 0.04
        p = 0.4
        t = 0.12
        #y_c = 0
        #d_y_d_x = 0
        if(x<p):
            y_c = (2*p*x-x**2)*(m)/(p**2)
            d_y_d_x = (p-x) * (2*m)/(p**2)
        else:
            y_c = ((1-2*p) +2*p*x -x**2)*m/(1-p)**2
            d_y_d_x = (p - x) * (2 * m) / ((1- p) ** 2)

        theta = math.atan(d_y_d_x)

        y_t = 5*t*(0.2696*math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843*x**3 - 0.1015*x**4)

        y_u = y_c + y_t*math.cos(theta)
        y_l = y_c - y_t * math.cos(theta)
        return [y_u,y_l] # [upper,lower]

    def thicknessDistNormalized(self,x): #🤠🤠🤠
        return [self.thicknessDist(x)[0]/ 0.12,self.thicknessDist(x)[1] / 0.12] #0.12 is the height of the airfoil

    def thicknessDistDerivative(self, x):
        #theta_1 = 0
        #theta_2 = 0
        h = 0.001

        d_x_1 = (self.thicknessDist(x+h)[0]  - self.thicknessDist(x)[0])/(h)
        #theta_1 = math.atan(d_x_1)
        d_x_2 = (self.thicknessDist(x + h)[1]  -  self.thicknessDist(x)[1]) / (h)
        #theta_2 = math.atan(d_x_2)

        # if (x == 0.0):
        #     d_x_1 = math.inf
        #     d_x_2 = -math.inf
        #return  [theta_1,theta_2] #angle
        return [d_x_1,d_x_2] #slope

    def thicknessDistfull(self,z):
        return np.sin(self.angle_of_attack(z))

    def Cp_u(self, x, z):
        # to be determined

        return self.Cpu(x,z)

    def Cp_l(self, x, z): #😎😎😎😎
        # to be determined

        return self.Cpl(x,z)

    def N_chord(self,z, x):
        dN = self.Cp_l(x+0.5,z) - self.Cp_u(x+0.5,z)
        #dN = self.C_pl(x, z)*np.cos(self.thicknessDistDericative(x)[1]) - self.C_pu(x, z)*np.cos(self.thicknessDistDericative(x)[0])
        return dN
    def A_chord(self, z,x):
        #dA = Cp_u(x,z) *np.sin(self.thicknessDistDericative(x)[1]) - Cpl(x,z) * df_thickness(x)[1]
        #print(self.Cp_u(x+0.5,z))
        #print(self.thicknessDistDerivative(x+0.5)[0])
        dA = self.Cp_u(x+0.5,z) *self.thicknessDistDerivative(x+0.5)[0] - self.Cp_l(x+0.5,z) *self.thicknessDistDerivative(x+0.5)[1]
        #print(self.thicknessDistDerivative(x+0.5)[0])
        return dA

    def liftDist(self,x,z):
        dN = self.N_chord(z,x)
        dA =self.A_chord(z,x)
        dL = dN*np.cos(self.angle_of_attack(z)) - dA*np.sin(self.angle_of_attack(z))
        #Lift = self.CLfunction(x + 0.5) * 0.5 * self.flow.density * (self.M_r(x) * self.flow.a_0) ** 2 * self.chord(x)
        #print(dL)
        return dL

    def normalisedLiftDist(self, x, z):
        return self.liftDist(x, z) /self.liftArea

    def areaLift(self,z):
        area = math2.integration(self.liftDist, -0.5, 0.5, 100, z, 'Simpsons')
        return area

    def dragDist(self, x,z):
        dN = self.N_chord(z, x)
        dA = self.A_chord(z, x)
        dD = dN * np.sin(self.angle_of_attack(z)) + dA * np.cos(self.angle_of_attack(z))
        #print(dD)
        #drag = self.CDfunction(x + 0.5) * 0.5 * self.flow.density * (self.M_r(x) * self.flow.a_0) ** 2 * self.chord(x)
        return dD

    def areaDrag(self,z):
        area = math2.integration(self.dragDist, -0.5, 0.5, 100, z, 'Simpsons')
        return area

    def normalisedDragDist(self, x, z):
        #area = math2.integration(self.dragDist, -0.5, 0.5, 100, z, 'Simpsons')
        return self.dragDist(x, z) / self.dragArea

class Math2: #😵
    def integration(self, func, a, b, steps, args, type):
        h = (b - a) / steps
        if type == 'trapezium':
            area = 0
            for i in range(0, steps + 1):
                x = a + i * h
                x2 = a + (i + 1) * h
                area += (func(x, args) + func(x2, args)) * h / 2
            return area
        if type == 'Simpsons':
            area = 0
            if ((steps + 1) / 2) % 1 == 0:
                i = 0
                while i <= steps - 3:
                    x = a + i * h
                    x2 = a + (i + 1) * h
                    x3 = a + (i + 2) * h
                    area += h / 3 * (func(x, args) + 4 * func(x2, args) + func(x3, args))
                    i = i + 2
                area += (func(a + steps * h, args) + func(a + (steps - 1) * h, args)) * h / 2
                return area
            else:
                i = 0
                while i <= steps - 2:
                    x = a + i * h
                    x2 = a + (i + 1) * h
                    x3 = a + (i + 2) * h
                    area += h / 3 * (func(x, args) + 4 * func(x2, args) + func(x3, args))
                    i = i + 2
                return area

    def besselsFunc(self, N, x, M=20):  # http://hyperphysics.phy-astr.gsu.edu/hbase/Math/bessel.html
        # standard bessels function
        J_x = 0
        for m in range(0, M):
            J_x += ((-1) ** m * x ** (2 * m)) / (2 ** (2 * m + N) * math.factorial(m) * math.factorial(N + m))
        J_x *= x ** N
        return J_x

    def integrateRiemannSums(self, func, a, b, steps, args=[]):
        # integrate function using riemann sums
        h = (b - a) / steps
        # print(steps)
        y = 0
        for i in range(0, steps + 1):
            x = a + i * h
            if (args != []):
                y += func(x, args) * h
            else:
                y += func(x) * h
        return y


math2 = Math2()

mic = 0.075
#microphone spacing is 0.075m
miccoord=[[12*mic, math.atan(1.2/(12*mic)),1.2], #x,theta,y
           [9*mic, math.atan(1.2/(9*mic)),1.2],
           [6*mic, math.atan(1.2/(6*mic)),1.2],
           [3*mic, math.atan(1.2/(3*mic)),1.2],
           [0*mic, np.pi/2,1.2],
           [-3*mic, math.atan(1.2/(-3*mic)),1.2],
           [-6*mic, math.atan(1.2/(-6*mic)),1.2],
           [-9*mic, math.atan(1.2/(-9*mic)),1.2],
           [-12*mic, math.atan(1.2/(-12*mic)),1.2],]

# flow: density,velocity,a_0
flow = Flow(1.225, 8, 343)
# propeller: flow,RPM,bladeNumber,diameter,x,theta,y
#propeller = Propeller(flow, 8000, 2, 0.3, 0, math.pi / 2, 1.2)

results=[]
#yArr = []
for i in range(len(miccoord)):
    propeller = Propeller(flow, 8000, 2, 0.3, miccoord[i][0], miccoord[i][1], miccoord[i][2])
    print('Microphone '+str(i+1))
    noise = ['microphone '+str(i+1) +' x: ' + str(miccoord[i][0]) + " y: 1.2"]
    for j in range(1, 4):
        x = j
        print('Harmonic '+ str(x))
        y = propeller.noise(x)
        print(str(y) + " dB")
        noise.append(y)
    results.append(noise)
        #xArr.append(133.33 * x)
        #yArr.append(y)
print(results)
#plt.plot(xArr, yArr, "r+")
#plt.show()
#😀😁😂🤣😃😄😅😆😉😊😋😎😍😘🥰😗😙😚☺🙂🤗🤩🤔🤨😐😑😶🙄😏😣😥😮🤐😯😪😫🥱😴😌😛😜😝🤤😒😓😔😕🙃🤑😲☹🙁😖😞😟😤😢😭😦😧😨😩🤯😬😰😱🥵🥶😳🤪😵🥴😠😡🤬😷🤒🤕🤢🤮🤧😇🥳🥺🤠🤡🤥🤫🤭🧐🤓😈👿
#stuff to do 🤩:
#more accurate airfoil x-coordinate 🤢, probably unfeasable actually
