import math

class Flow:
    def __init__(self,density,velocity,a_0):

        self.density = density
        self.velocity = velocity
        self.a_0 = a_0
        self.M_x = velocity/a_0




class Propeller:
    def __init__(self,flow,RPM,bladeNumber,diameter,x):
        self.RPM = RPM
        self.flow = flow
        self.bladeNumber = bladeNumber
        self.diameter = diameter
        self.x = x
        self.theta = #whatever
        self.M_t = propeller.RPM * 2 * math.pi / 60 * diameter / 2
        self.bladePassingFrequency = (RPM/60)*bladeNumber

    def M_r(self,z):
        return math.sqrt((self.flow.M_x ** 2) + (z ** 2 * self.M_t ** 2))

    def k_x(self):
        return 0

    def k_y(self):
        return 0

    def psi_V(self):
        return 0

    def psi_D(self):
        return 0

    def psi_L(self):
        return 0

flow = Flow()

propeller = Propeller()