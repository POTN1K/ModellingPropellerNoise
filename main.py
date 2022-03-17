import math
import numpy as np

class Flow:
    def __init__(self,density,velocity,a_0):

        self.density = density
        self.velocity = velocity
        self.a_0 = a_0
        self.M_x = velocity/a_0


class Propeller:
    def __init__(self,flow,RPM,bladeNumber,diameter,x,theta,y):
        self.RPM = RPM
        self.flow = flow
        self.bladeNumber = bladeNumber
        self.diameter = diameter
        self.radius = 0.5 * diameter
        self.x = x
        self.theta = theta
        self.M_t = (self.RPM * 2 * math.pi / 60 * diameter / 2)/flow.a_0
        self.bladePassingFrequency = (RPM/60)*bladeNumber
        self.y = y



    #propeller properties
    def B_D(self,z): #chord to diameter ratio
        r_o = z-0.2

        B_D_o=[]
        if isinstance(r_o,float):
            i=r_o
            B_D_o=0.127 + 0.7404112848810135 * (i ** 1) + -4.673821287763687 * (i ** 2) + 19.02690171288593 * (i ** 3) + -48.12645454430276 * (i ** 4) + 36.20381779945387 * (i ** 5) + 111.0096041870429 * (i ** 6) + -311.5856431939494 * (i ** 7) + 309.1780229926285 * (i ** 8) + -113.64512633245806 * (i ** 9)
            return 0.5 * (B_D_o + 0.05)
        else:
            B_D_o = []
            for i in r_o:
                if(i<0.1):
                    B_D_o.append(0)
                else:
                    B_D_o.append(0.127 + 0.7404112848810135 * (i ** 1) + -4.673821287763687 * (i ** 2) + 19.02690171288593 * (i ** 3) + -48.12645454430276 * (i ** 4) + 36.20381779945387 * (i ** 5) + 111.0096041870429 * (i ** 6) + -311.5856431939494 * (i ** 7) + 309.1780229926285 * (i ** 8) + -113.64512633245806 * (i ** 9))
            B_D = 0.5 * (np.array(B_D_o) + 0.05)
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
        return -(self.flow.density * (self.flow.a_0**2) * self.bladeNumber * math.sin(self.theta) * np.exp(1j * m * self.bladeNumber * ( ( (omega*self.radius)/(self.flow.a_0) ) - (0.5*math.pi)) ) )/( 8 * math.pi * (self.y/self.diameter) * (1-self.flow.M_x*math.cos(self.theta)))

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
        return math2.integration(step, derivativepsi_D, 1)

    def psi_L(self, step,z,m):
        z_arr = np.arange(-0.5, 0.5, step)
        derivativepsi_L = self.psi_L_derivative(z_arr,[z,m])
        return math2.integration(step,derivativepsi_L,1)

    def psi_V(self,step, z,m):
        z_arr = np.arange(-0.5, 0.5, step)
        derivativepsi_V = self.psi_V_derivative(z_arr, [z,m])
        return math2.integration(step=step,array= derivativepsi_V,type=1)


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
        Drag_sec = self.CDfunction(x+0.5) * 0.5 * self.flow.density * (self.M_r(x) * self.flow.a_0) ** 2 * self.chord(z)
        return Drag_sec * np.exp(1j*k_x*x)

    def psi_L_derivative(self,x,args):
        z = args[0]
        m = args[1]
        k_x = self.k_x(z, m)
        Lift_sec = self.CLfunction(x+0.5) * 0.5 * self.flow.density * (self.M_r(x) * self.flow.a_0) ** 2 * self.chord(z)
        return Lift_sec * np.exp(1j*k_x*x)

    def psi_V_derivative(self,x,args):
        z = args[0]
        m = args[1]
        k_x = self.k_x(z, m)
        return self.thickness(x+0.5) * np.exp(1j*k_x*x)



    #p(t)
    def pressure(self,m,method):
        k = self.calcFirstPart(m)
        if method==0:
            step=0.01
            z=np.arange(0,1,step)
            form=0
            type=1
            deriV=self.p_Vm_derivative(z,m)[form]
            deriD=self.p_Dm_derivative(z,m)[form]
            deriL=self.p_Lm_derivative(z,m)[form]
            p_Vm = k * math2.integration(step, deriV, type)
            p_Dm = k * math2.integration(step, deriD, type)
            p_Lm = k * math2.integration(step, deriL, type)
            #print(p_Dm,p_Lm,p_Vm)
        else:
            p_Vm = k*math2.integrateRiemannSums(self.p_Vm_derivative,0,1,100,m,method)
            p_Dm = k*math2.integrateRiemannSums(self.p_Dm_derivative,0,1,100,m,method)
            p_Lm = k*math2.integrateRiemannSums(self.p_Lm_derivative,0,1,100,m,method)


        p_mb = p_Vm + p_Dm + p_Lm
        return p_mb


    #harmonic noise level
    def noise(self,m,method):
        p = self.pressure(m,method)*2
        print(p)
        p_ref = 20 * 10**-6
        OASPL = 10 * math.log10(p / (p_ref))
        return OASPL

    def p_Vm_derivative(self,z,m):
        k_x = self.k_x(z, m)
        P_Vm_1 = (self.M_r(z)) ** 2 * math2.besselsFunc(m*self.bladeNumber, m * self.bladeNumber * z * self.M_t*np.sin(self.theta)/(1-self.flow.M_x*np.cos(self.theta)))
        P_Vm_2 = self.k_x(z,m)**2*self.thickness(z)*self.psi_V(0.01,z,m)
        P_Vm_2_alt = self.k_x(z,m)**2*self.thickness(z)*self.psi_V2(z,m)

        return (P_Vm_1*P_Vm_2,P_Vm_1*P_Vm_2_alt)

    def p_Dm_derivative(self, z,m):
        k_x = self.k_x(z,m)
        C_D = self.CDfunction(z)
        P_Dm_1 = (self.M_r(z)) ** 2 * math2.besselsFunc(m*self.bladeNumber, m*self.bladeNumber * z * self.M_t*np.sin(self.theta)/(1-self.flow.M_x*np.cos(self.theta)))
        P_Dm_2 = 1j * k_x * (C_D/2) * self.psi_D(0.01,z,m)
        P_Dm_2_alt = 1j * k_x * (C_D/2) * self.psi_D2(z,m)
        return (P_Dm_1*P_Dm_2,P_Dm_1*P_Dm_2_alt)

    def p_Lm_derivative(self, z,m):
        k_x = self.k_x(z, m)
        C_L = self.CDfunction(z)
        P_Lm_1 = (self.M_r(z))**2*math2.besselsFunc(m*self.bladeNumber,m*self.bladeNumber*z*self.M_t*np.sin(self.theta)/(1-self.flow.M_x*np.cos(self.theta)))
        P_Lm_2 = 1j * k_x *C_L/2*self.psi_L(0.01,z,m)
        P_Lm_2_alt = 1j * k_x *C_L/2*self.psi_L2(z,m)
        return (P_Lm_1*P_Lm_2,P_Lm_1*P_Lm_2_alt)

    def CLfunction(self, x):
        Lift = 1
        return .2

    def CDfunction(self, x):
        return .0005

class Math2:
    def integration(self, step, array, type):
        if type==0:
            area=0
            for i in range(len(array)-1):
                area += (array[i]+ array[i+1])*step/2
            return area
        if type == 1:
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

    def integrateRiemannSums(self,func,a,b,steps,args=[],index=0):
        #integrate function using riemann sums
        h = (b - a) / steps
        #print(steps)
        y = 0
        for i in range(0, steps + 1):
            x = a + i * h
            y += func(x,args)[index] * h
        return y

#Questions to ask Ragni
#1. propeller radius. Why does it start at 0.2?

math2 = Math2()
#flow: density,velocity,a_0
flow = Flow(1.225,8,343)
#propeller: flow,RPM,bladeNumber,diameter,x,theta,y
propeller = Propeller(flow,8000,2,0.3,0,math.pi/4,1.2)
z=np.arange(0,np.pi,0.01)
print(propeller.noise(1,0))
