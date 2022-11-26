import numpy as np


class Mars:
    sigma = 5.67e-8
    
    def __init__(self, absorptivity_PV, absorptivity_dessert, emissivity_dessert, emissivity_PV, delta_time, f,
                 cp, T_Atmosphere, rho, x):
        self.absorptivity_PV = absorptivity_PV
        self.absorptivity_dessert = absorptivity_dessert
        self.emissivity_dessert = emissivity_dessert
        self.emissivity_PV = emissivity_PV
        self.delta_time = delta_time
        self.f = f
        self.cp = cp
        self.T_Atmosphere = T_Atmosphere
        self.rho = rho
        self.x = x
    
    def l_out(self, T, emissivity):
        return emissivity * self.sigma * np.power(T, 4)
    
    def H(self, T, r):
        return self.rho * self.cp * (T - self.T_Atmosphere) / r
    
    def G(self, T_prior, T_current):
        return self.x * (T_prior - T_current) / self.delta_time
    
    def system_PV(self, T_PV_current, T_SH_current, l_in, s_in, r_H_PV):
        return s_in * (1 - self.absorptivity_PV) + l_in + self.l_out(T_SH_current, self.emissivity_dessert) - \
               self.l_out(T_PV_current, self.emissivity_PV) - self.H(T_PV_current, r_H_PV)
    
    def system_under_PV(self, T_PV_current, T_SH_prior, T_SH_current, s_in, r_H_dessert):
        return self.f * s_in * (1 - self.absorptivity_dessert) + self.l_out(T_PV_current, self.emissivity_dessert) / 2 \
               - self.H(T_SH_current, r_H_dessert) - self.G(T_SH_prior, T_SH_current)
    
    def system(self, z, num_days, l_in, s_in, r_H_PV, r_H_dessert):
        T_PV = z[0:num_days]
        T_SH = z[num_days:2 * num_days]
        
        F = np.empty((2 * num_days))
        for i in range(0, num_days):
            F[i] = self.system_PV(T_PV[i], T_SH[i], l_in=l_in[i], s_in=s_in[i], r_H_PV=r_H_PV[i])
        
        F[num_days] = self.system_under_PV(T_PV[0], T_SH_prior=T_SH[num_days - 1], T_SH_current=T_SH[0],
                                           s_in=s_in[0], r_H_dessert=r_H_dessert[0])
        for i in range(1, num_days):
            F[num_days + i] = self.system_under_PV(T_PV[i], T_SH_prior=T_SH[i - 1], T_SH_current=T_SH[i],
                                                   s_in=s_in[i], r_H_dessert=r_H_dessert[i])
        
        return F
