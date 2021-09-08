import numpy as np
from .constants import *
import ufo
import time 

#Results are normalized by ring length?
def touschek(optics, ma, beam=ufo.Beam()): #line, twi0_h, twi0_v, emit_h, emit_v, ma, slices=5):
        #pos_ma = ma[ma[:, 0] >= 0.]
        #neg_ma = ma[ma[:, 0] <= 0.]

        #neg_ma[:, 0] = -neg_ma[:, 0]
        #neg_ma = np.flipud(neg_ma)

        # Crude integral sampling points
        count = 100
        u = np.linspace(0., 1., count) + 1./count

        e_count = beam.charge / elementary_charge
        rate = 0.

        t0 = time.time()
        #ALLL WRONG, MA MUST BE EVALUATED AT EVERY POINT. PERIOD!
        #cycle over optics and look at closest MA sample point
        for (idx, _), ax, ay, bx, by, dx, dpx in zip(optics.where, optics.ax, optics.ay, optics.bx, optics.by, optics.dx, optics.dpx):
            ma_idx = min(ma[:, 0], key=lambda x:abs(x - idx))

            gx = (1. + ax*ax) / bx
            sigma_x = np.sqrt(bx * beam.ex + (dx * beam.energy_spread)**2)
            sigma_y = np.sqrt(by * beam.ey) #+ (d_v * self.energy_spread)**2)

            H = gx * dx**2 + 2.*ax*dx*dpx + bx*dpx**2
            sigma_dx = beam.ex / sigma_x * np.sqrt(1. + H * beam.energy_spread**2 / beam.ex)

''' 

                        for ma in [neg_ma, pos_ma]:
                                tau_m = 0.
                                for dp, j in ma:
                                        action = 0.5 * gx * (dx*dp)**2
                                        if action > j:
                                                break
                                        tau_m = dp

                                tau = (tau_m / self.gamma / sigma_dx)**2
                                integral = np.sum((1./u - 0.5*np.log(1./u) - 1.0) * np.exp(-tau / u)) / count
                                rate += element.L / slices / (sigma_x * sigma_y * sigma_dx * tau_m**2) * integral
           
        print(time.time() - t0)
'''
'''
                for slice_count in range(slices):
                        alpha_h, beta_h, gamma_h, d_h, dp_h, _ = propagate_twiss(twi0_h, Mh)
                        alpha_v, beta_v, gamma_v, d_v, dp_v, _ = propagate_twiss(twi0_v, Mv)
                        sigma_x = np.sqrt(beta_h * emit_h + (d_h * self.energy_spread)**2)
                        sigma_y = np.sqrt(beta_v * emit_v + (d_v * self.energy_spread)**2)

                        H = gamma_h * d_h**2 + 2.*alpha_h*d_h*dp_h + beta_h*dp_h**2
                        sigma_dx = emit_h / sigma_x * np.sqrt(1. + H * self.energy_spread**2 / emit_h)
                        #sigma_dx = np.sqrt(gamma_h * emit_h / np.pi)

                        for ma in [neg_ma, pos_ma]:
                                tau_m = 0.
                                for dp, j in ma:
                                        action = 0.5 * gamma_h * (d_h*dp)**2
                                        if action > j:
                                                break
                                        tau_m = dp

                                tau = (tau_m / self.gamma / sigma_dx)**2
                                integral = np.sum((1./u - 0.5*np.log(1./u) - 1.0) * np.exp(-tau / u)) / count
                                rate += element.L / slices / (sigma_x * sigma_y * sigma_dx * tau_m**2) * integral
                        
                        Mh = np.matmul(ele_mat_h, Mh)
                        Mv = np.matmul(ele_mat_v, Mv)

        rate /= 2. #Adding positive and negative energy
        rate *= e_count / (8. * np.pi * self.gamma**3 * self.bunch_length * ring_length)
        rate *= speed_of_light * electron_radius**2

        return 1./rate
'''

