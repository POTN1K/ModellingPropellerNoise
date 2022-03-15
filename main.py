import math
import numpy as np

class Flow:
    def __init__(self,density,velocity,a_0):

        self.density = density
        self.velocity = velocity
        self.a_0 = a_0
        self.M_x = velocity/a_0


class Pressure:
    def __init__(self,P_Vm,P_Dm,P_Lm):
        self.P_Vm = P_Vm
        self.P_Dm = P_Dm
        self.P_Lm = P_Lm

class Propeller:
    def __init__(self,flow,RPM,bladeNumber,diameter,x,theta,y):
        self.RPM = RPM
        self.flow = flow
        self.bladeNumber = bladeNumber
        self.diameter = diameter
        self.x = x
        self.theta = theta
        self.M_t = propeller.RPM * 2 * math.pi / 60 * diameter / 2
        self.bladePassingFrequency = (RPM/60)*bladeNumber
        self.y = y



    #propeller properties
    def B_D(self,z): #chord to diameter ratio
        r_o = z-0.2
        if(r_o<0.1):
            return 0
        B_D_o = 0.127 + 0.7404112848810135 * (r_o ** 1) + -4.673821287763687 * (r_o ** 2) + 19.02690171288593 * (r_o ** 3) + -48.12645454430276 * (r_o ** 4) + 36.20381779945387 * (r_o ** 5) + 111.0096041870429 * (r_o ** 6) + -311.5856431939494 * (r_o ** 7) + 309.1780229926285 * (r_o ** 8) + -113.64512633245806 * (r_o ** 9)
        B_D = 0.5 * (B_D_o + 0.05)
        return B_D

    def chord(self, z):
        c = self.B_D(z) * self.diameter
        return c

    def thickness(self,z):
        return 0.12 * self.chord(z)

    def M_r(self,z):
        return np.sqrt((self.flow.M_x ** 2) + (z ** 2 * self.M_t ** 2))


    #calculate the first constant in the large formula
    def calcFirstPart(self,m):
        omega = (self.RPM/60)*2*math.pi
        return -(self.flow.)/()

    #k_x and k_y
    def k_x(self, z, m):
        k_x = 2*self.bladeNumber*m*self.B_D(z)*self.M_t/(self.M_r(z)*(1-self.flow.M_x*math.cos(self.theta)))
        return k_x

    def k_y(self,z, m):
        return (2 * self.bladeNumber * m * self.B_D(z) * (self.flow.M_x - self.M_r(z) ** (2) * np.cos(self.theta))) / (z * self.M_r(z) * (1 - self.flow.M_x * math.cos(self.theta)))



    #The psi's
     #method 1
    def psi_D(self, step,z,m):
        z_arr = np.arange(-0.5, 0.5, step)
        derivativepsi_D = self.psi_D_derivative(z_arr, [z,m])
        return Math2.integration(step, derivativepsi_V, 'Simpsons')

    def psi_L(self, step,z,m):
        z_arr = np.arange(-0.5, 0.5, step)
        derivativepsi_L = self.psi_L_derivative(z_arr,[z,m]z)
        return Math2.integration(step,derivativepsi_L,'Simpsons')

    def psi_V(self,step, z,m):
        z_arr = np.arange(-0.5, 0.5, step)
        derivativepsi_V = self.psi_V_derivative(z_arr, [z,m])
        return Math2.integration(step, derivativepsi_V, 'Simpsons')


     #method 2
    def psi_D2(self, z,m):
        return math2.integrateRiemannSums(self.psi_D_derivative,-0.5,0.5,100,[z,m])

    def psi_L2(self, z,m):
        return math2.integrateRiemannSums(self.psi_L_derivative,-0.5,0.5,100,[z,m])

    def psi_V2(self,z,m):
        return math2.integrateRiemannSums(self.psi_V_derivative,-0.5,0.5,100,[z,m])



    #Derivative of psi's
    def psi_D_derivative(self,x,args):
        z = args[0]
        m = args[1]
        k_x = self.k_x(z, m)
        Drag_sec = self.CDfunction(x) * 0.5 * self.flow.density * (self.M_r(x) * self.flow.a_0) ** 2 * self.chord(z)
        return Drag_sec * np.exp(1j*k_x*x)

    def psi_L_derivative(self,x,args):
        z = args[0]
        m = args[1]
        k_x = self.k_x(z, m)
        Lift_sec = self.CLfunction(x) * 0.5 * self.flow.density * (self.M_r(x) * self.flow.a_0) ** 2 * self.chord(z)
        return Lift_sec * np.exp(1j*k_x*x)

    def psi_V_derivative(self,x,args):
        z = args[0]
        m = args[1]
        k_x = self.k_x(z, m)
        return self.thickness(x) * np.exp(1j*k_x*x)



    def CLfunction(self,x):

        return 0

    def CDfunction(self,x):

        return 0

class Math2:
    def integration(self, step, array, type):
        if type=='trapezium':
            area=0
            for i in range(len(array)-1):
                area += (array[i]+ array[i+1])*step/2
            return area
        if type== 'Simpsons':
            area=0
            if (len(array)/2)%1==0:
                i=0
                while i<=len(array)-3:
                    area += step/3*(array[i]+4*array[i+1]+array[i+2])
                    i = i+2
                area += (array[-2]+array[-1])*step/2
                return area
            else:
                i=0
                while i<=len(array)-2:
                    area += step/3*(array[i]+4*array[i+1]+array[i+2])
                    i = i+2
                return area
    def besselsFunc(self,N,x,M=10): #http://hyperphysics.phy-astr.gsu.edu/hbase/Math/bessel.html
        #standard bessels function
        J_x = 0
        for m in range(0,M):
            J_x += ((-1)**m * x**(2*m))/( 2**(2*m + N) * math.factorial(m) * math.factorial(N+m))
        J_x *= x**N
        return J_x

    def integrateRiemannSums(func,a,b,steps,args=[]):
        #integrate function using riemann sums
        h = (b - a) / steps
        print(steps)
        y = 0
        for i in range(0, steps + 1):
            x = a + i * h
            y += func(x,args) * h
        return y

#Questions to ask Ragni
#1. propeller radius. Why does it start at 0.2?

flow = Flow()
math2 = Math2
propeller = Propeller()

