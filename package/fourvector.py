import math
import numpy as np

def eta_to_angle(eta):
    """ convert pseudorapidity to angle
    """
    try:
        return 2 * math.atan(math.exp(-eta))
    except Exception:
        print(eta)
        return 1.4

def invariant_mass(four_momentum):
    """ return the invariant mass of particles

    Parameters
    ----------
    four_momentum : list of len 4
    """
    energy, p_x, p_y, p_z = 0, 0, 0, 0
    for each in four_momentum:
        if each[0] == 0:
            each[0] = np.linalg.norm([each[1], each[2], each[3]])
        energy += each[0]
        p_x += each[1]
        p_y += each[2]
        p_z += each[3]
    mass = energy**2 - pow(np.linalg.norm([p_x, p_y, p_z]), 2)

    # detect negative mass
    if mass < 0:
        print('Warning: negative mass encontered in invariant_mass() ', mass)
        return 0
    return pow(mass, 0.5)

def phi(p_x, p_y):
    ''' return the polar angle of a transverse momentum
    '''
    if (p_x > 0 and p_y > 0) or (p_x < 0 and p_y > 0):
        return np.arctan(p_x / p_y)
    if p_x > 0 and p_y < 0:
        return math.pi - abs(np.arctan(p_x / p_y))
    if p_x < 0 and p_y < 0:
        return - (math.pi - abs(np.arctan(p_x / p_y)))
    if abs(p_x) < 1e-5:
        if p_y > 1e-5:
            return math.pi/2.
        if p_y < -1e-7:
            return -math.pi/2.
    if abs(p_y) < 1e-5:
        if p_x > 1e-5:
            return 0.
        if p_x < -1e-5:
            return math.pi
    raise ValueError('Error: Vector (0, 0) enconterd in phi(p_x, p_y).')

class FourVector():
    def __init__(self, datatype, p1, p2, p3, p4):
        if datatype == "ptphietam":
            x = p1 * math.sin(p2)
            y = p1 * math.cos(p2)
            z = p1 * math.tan(math.pi / 2.0 - eta_to_angle(p3))
            e = np.linalg.norm([x, y, z, p4])
            self.four_momentum = [e, x, y, z]
            self.m = p4

        elif datatype == "xyzE":
            pass
        else:
            raise ValueError("Error: undefined constructor.")

    def __getitem__(self, key):
        return self.four_momentum[key]

class Particles():
    def __init__(self, particles):
        self.allparticle = particles
    def invariant_mass(self):
        fourvecotor = []
        for each in self.allparticle:
            fourvecotor.append(each.four_momentum)
        return invariant_mass(fourvecotor)
    def transverse_mass(self):
        energy, p_x, p_y = 0, 0, 0
        for each in self.allparticle:
            energy += np.linalg.norm([each[1], each[2], each.m])
            p_x += each[1]
            p_y += each[2]
        mass = energy**2 - pow(np.linalg.norm([p_x, p_y]), 2)
        # detect negative mass
        if mass < 0:
            print('Warning: negative mass encontered in invariant_mass() ', mass)
            return 0
        return pow(mass, 0.5)

if __name__ == '__main__':
    a = FourVector("ptphietam", 1,0 , 0, 0)
    b = FourVector("ptphietam", 1, -math.pi/2, 0, 0)
    p = Particles([a,b])
    print(p.invariant_mass())
