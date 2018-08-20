# -*- coding: utf-8 -*-
"""
Created on Sun Aug  5 12:14:57 2018

@author: T.J.-LAB-PC
"""
import os
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import optimize
from cea_post import Read_datset
from tqdm import tqdm

class Single_tank():
    """ Simulate combustion of hybrid rocket
    
    Parameter
    ---------
    
    Return
    ---------
    """
    def __init__(self, cea_fldpath, dt, tb_init, Pc_init, eta, Pti, Ptf, Vox, Do, Cd, Lf, Dfi, rho_ox, rho_f, Dti, De, Pa, rn, theta):
        self.cea_fldpath = cea_fldpath
        self.dt = dt # time interval [s]
        self.tb_init = tb_init # initial firing duration time for iteration [s]
        self.Pc_init = Pc_init # initial chamber pressure for iteration [Pa]
        self.eta = eta
        self.Pti = Pti # initial tank pressure [Pa]
        self.Ptf = Ptf # final tank pressure [Pa]
        self.Vox = Vox # oxidizer volume [m^3]
        self.Do = Do
        self.Cd = Cd
        self.Lf = Lf # fuel length [m]
        self.Dfi = Dfi # initial fuel port diameter [m]
        self.rho_ox = rho_ox
        self.rho_f = rho_f
        self.Mox_fill = Vox*rho_ox # total oxidizer mass [kg]
        self.Dti = Dti # initial nozzle throat diameter [m]
        self.De = De # nozzle exit diameter [m]
        self.Pa = Pa # ambient pressure [Pa]
        self.rn = rn # nozzle throat regression rate [m]
        self.lmbd = (1 + np.cos( np.deg2rad(theta) ))/2 # nozzle coefficient. theta is nozzle opening half-angle [rad]
        self.func_cstr = Read_datset(cea_fldpath).gen_func("CSTAR")
        self.func_gamma = Read_datset(cea_fldpath).gen_func("GAMMAs_c")
        self.of_max = Read_datset(cea_fldpath).of.max()
        self.Pc = np.array([])
        self.Pt = np.array([])        
        self.mox = np.array([])
        self.Mox = np.array([])
        self.Mf = np.array([])
        self.Df = np.array([])
        self.Gox = np.array([])
        self.r = np.array([])
        self.Dt = np.array([])        
        self.mf = np.array([])
        self.Pe = np.array([])
        self.CF = np.array([])
        self.F = np.array([])
        self.I = np.array([])
        
    def exe_sim(self):
        """ Execute Simulation
        """
        self.iterat_tb(tb_init=self.tb_init, tb_min=1.0, tb_max=10.0, maxiter=100)        


    def iterat_tb(self, tb_init=5.0, tb_min=1.0, tb_max=10, maxiter=100):
        """ Iteration of tb: firing duration time
        
        Parameter
        ---------
        
        Return
        --------
        """
        self.num_iterat_tb = 0
        try:
            self.tb = optimize.newton(self.func_error_tb, tb_init, tol= 1.0e-2,maxiter=maxiter)
        except:
            self.tb = optimize.brentq(self.func_error_tb, tb_min, tb_max, xtol=1.0e-2, maxiter=maxiter, full_output=False)
        self.df = pd.DataFrame([], index=np.arange(0, self.tb+self.dt/2, self.dt))
        self.df["Pc"] = self.Pc
        self.df["Pt"] = self.Pt        
        self.df["mox"] = self.mox
        self.df["Mox"] = self.Mox
        self.df["mf"] = self.mf
        self.df["Mf"] = self.Mf
        self.df["Df"] = self.Df
        self.df["Gox"] = self.Gox
        self.df["r"] = self.r
        self.df["Dt"] = self.Dt
        self.df["Pe"] = self.Pe
        self.df["CF"] = self.CF
        self.df["F"] = self.F
        self.df["I"] = self.I
        return(self.df)
    
    def func_error_tb(self, tb):
        self.num_iterat_tb += 1
        if self.num_iterat_tb == 1:
            ordinal = "st"
        elif self.num_iterat_tb == 2:
            ordinal = "nd"
        elif self.num_iterat_tb == 3:
            ordinal = "rd"
        else:
            ordinal = "th"
        print("{}".format(self.num_iterat_tb) + ordinal +" iteration about firing duration time; tb.")
        tb_cal = self.func_tb(tb)
        diff = tb_cal - tb
        error = diff/tb_cal
        return(error)
        
    def func_tb(self, tb):
        self.Pc = np.array([])
        self.Pt = np.array([])
        self.mox = np.array([])
        self.Mox = np.array([])
        self.Mf = np.array([])
        self.Df = np.array([])
        self.Gox = np.array([])
        self.r = np.array([])
        self.Dt = np.array([])        
        self.mf = np.array([])
        self.Pe = np.array([])
        self.CF = np.array([])
        self.F = np.array([])
        self.I = np.array([])
        self._tmp_Mox_ = 0
        t = 0
        pbar = tqdm(total=tb/self.dt)
        while(self._tmp_Mox_ <= self.Mox_fill):
            Pc = self.iterat_Pc(t, tb)
            of = self._tmp_of_
            eps = np.power(self.De/self._tmp_Dt_, 2.0)
            Pe = self.iterat_Pe(of, Pc, eps)
            gamma = self.func_gamma(of, Pc)
            if Pe==0:
                Pe = self.Pa
            if Pc <= 0:
                Pc = Pe
            CF = np.sqrt((2*np.power(gamma, 2)/(gamma-1))*np.power(2/(gamma+1), (gamma+1)/(gamma-1))*(1-np.power(Pe/Pc,(gamma-1)/gamma))) + (Pe-self.Pa)*eps/Pc
            F = self.lmbd*CF*Pc*np.power(self._tmp_Dt_, 2)
            self.Pe = np.append(self.Pe, Pe)
            self.CF = np.append(self.CF, CF)        
            self.F = np.append(self.F, F)
            t_list = np.arange(0, t+self.dt/2, self.dt)
            I = integrate.simps(self.F, t_list)
            self.I = np.append(self.I, I)
            t = t + self.dt
#            print("t = {}s".format(round(t,3)))
            pbar.update(1) # show the progress bar
        tb = t - self.dt
        return(tb)
        
        
    
    def iterat_Pe(self, of, Pc, eps):
        """
        of: float
            O/F
        Pc: float
            Chamber pressure [Pa]
        eps: float
            Nozzle expansion ratio
        """
        if of<=0:
            of = 1.0e-2
        try:
            Pe = optimize.newton(self.func_error_eps, Pc/2, maxiter=100, tol=1.0e-3, args=(of, Pc, eps))
        except:
            Pe = optimize.brentq(self.func_error_eps, 1, Pc/2, maxiter=100, xtol=1.0e-3, full_output=False, args=(of, Pc,eps))
        return(Pe)
    
    def func_error_eps(self, Pe, of , Pc, eps):
        eps_cal = self.func_eps_cal(Pe, of, Pc)
        diff = eps_cal - eps
        error = diff/eps_cal
        return(error)
    
    def func_eps_cal(self, Pe, of, Pc):
        gam = self.func_gamma(of, Pc)
        if Pe == 0:
            Pe = self.Pa
        if Pc == 0:
            Pc = Pe
        eps_cal = np.power((2/(gam+1)), 1/(gam-1)) * np.power(Pc/Pe, 1/gam) / np.sqrt((gam+1)/(gam-1)*(1-np.power(Pe/Pc, (gam-1)/gam)))
        return(eps_cal)

    

    
    def iterat_Pc(self, t, tb, Pc_min=0.2e+6, maxiter=100):
        try:
            if t == 0:
                Pc_init = self.Pc_init
            else:
                Pc_init = self.Pc[-1]
            Pc = optimize.newton(self.func_error_Pc, Pc_init, maxiter=maxiter, tol=1.0e-3, args=(t,tb))
        except:
            Pc = optimize.brentq(self.func_error_Pc, Pc_min, self.Pti, xtol=1.0e-3, maxiter=maxiter, args=(t,tb))
        self.Pc = np.append(self.Pc, Pc)
        self.mox = np.append(self.mox, self._tmp_mox_)
        self.Mox = np.append(self.Mox, self._tmp_Mox_)
        self.Mf = np.append(self.Mf, self._tmp_Mf_)
        self.Df = np.append(self.Df, self._tmp_Df_)
        self.Gox = np.append(self.Gox, self._tmp_Gox_)
        self.r = np.append(self.r, self._tmp_r)
        self.mf = np.append(self.mf, self._tmp_mf_)
        self.Dt = np.append(self.Dt, self._tmp_Dt_)
        self.Pt = np.append(self.Pt, self._tmp_Pt_)
        return(Pc)

        
    def func_error_Pc(self, Pc, t, tb):
        Pc_cal = self.func_Pc(Pc, t, tb)
        diff = Pc - Pc_cal
        error = diff/Pc_cal
        return(error)
    
    def func_Pc(self, Pc, t, tb):
        mox = self.func_mox(t, Pc, tb)
        mf = self.func_mf(t, mox)
        if mf == 0 and mox>0:
            of = self.of_max
        elif mox==0 and mox==0:
            of = 0
        else:
            of = mox/mf
        if of > self.of_max:
            of = self.of_max
        self._tmp_of_ = of
        cstr = self.eta*self.func_cstr(of, Pc)
        At = self.func_At(t)
        Pc_cal = cstr*(mox + mf)/At
        return(Pc_cal)

    def func_At(self, t):
        Dt = self.Dti -2*self.rn*t
        At = np.pi*np.power(Dt, 2)/4
        self._tmp_Dt_ = Dt
        return(At)

    def func_mf(self, t, mox):
        Df = self.func_Df(t)
        r = self.func_regression(Df, mox)
        mf = self.Lf*Df*self.rho_f*r
        self._tmp_mf_ = mf
        return(mf)

    def func_Df(self, t):
        Mf = self.func_Mf(t)
        Df = np.sqrt(4*Mf/(np.pi*self.rho_f*self.Lf) + np.power(self.Dfi, 2))
        self._tmp_Df_ = Df
        return(Df)

    def func_regression(self, Df, mox):
        a = 1.310e-4 # regression coefficient [m^3/kg]
        n = 0.340 # oxidizer mass flux exponent
        Gox = self.func_Gox(Df, mox)
        r = a*np.power(Gox, n)
        self._tmp_r = r
        return(r)

    def func_Gox(self, Df, mox):
        Gox = 4*mox/(np.pi*Df)
        self._tmp_Gox_ = Gox
        return(Gox)



    def func_Mf(self, t):
        """ Calculate total fuel consumption Mf
        
        Paratmeter
        ----------
        t_list: 1d-arra like
            time liset by t
        mf: 1d-array like (need same size of t)
            fuel mass flow rate array
    
        Return
        ------
        Mox: float
            total oxidizer mass [kg] by t [s]
        """    
        t_list = np.arange(0, t-self.dt/2, self.dt)
        if t == 0:
            Mf = 0
        else:
#            tmp_mf = np.append(self.mf, mf)
#            Mf = integrate.simps(tmp_mf, t_list)
            Mf = integrate.simps(self.mf, t_list)
        self._tmp_Mf_ = Mf
        return(Mf)
    


    def func_mox(self, t, Pc, tb):
        """ Calculate oxidizer mass flow rate
        
        Parameter
        ---------
        Cd: float
            discharge coefficient [-]
        Do: float
            orifice diameter [mm]
        rho_ox: float
            oxidizer density [kg/m^3]
        t: float
            time [s]
        Pt: float
            tank pressure [Pa]
        Pc: float
            camber pressure [Pa]
        """
        Pt = self.func_Pt(t, tb)
        mox = self.Cd*(np.pi*np.power(self.Do, 2)/4)*np.sqrt(2*self.rho_ox*(Pt-Pc))
        Mox =self.func_Mox(t, mox) # calculate total oxidizer mass flow
        self._tmp_mox_ = mox
#        self.mox = np.append(self.mox, self._tmp_mox)
        return(mox)

    def func_Pt(self, t, tb):
        """ Calculate tank pressure
        
        Parameter
        ---------
        t: float
            time [s]
        tb: float
            firing duration time [s]
        
        Return
        -------
        Pt: float
            tank pressure at t [s], [Pa]
        """
        Pt = (self.Ptf - self.Pti)*t/tb + self.Pti # assuming linear changing on tank pressure
        self._tmp_Pt_ = Pt
        return(Pt)


    def func_Mox(self, t, mox):
        """ Calculate total oxidizer mass Mox
        
        Paratmeter
        ----------
        t_list: 1d-arra like
            time liset by t
        mox: 1d-array like (need same size of t)
            oxidizer mass flow rate array
    
        Return
        ------
        Mox: float
            total oxidizer mass [kg] by t [s]
        """
        if t==0:
            Mox = 0
        else:
            t_list = np.arange(0, t+self.dt/2, self.dt)
            tmp_mox = np.append(self.mox, mox)
            Mox = integrate.simps(tmp_mox, t_list)
            self._tmp_Mox_ = Mox
        return(Mox)


class Double_tank(Single_tank):
    def __init__(self, cea_fldpath, dt, tb_init, Pc_init, eta, Pti, Vox, Do, Cd, Lf, Dfi, rho_ox, rho_f, Dti, De, Pa, rn, theta):
        self.dt = dt # time interval [s]
        self.tb_init = tb_init # initial firing duration time for iteration [s]
        self.Pc_init = Pc_init # initial chamber pressure for iteration [Pa]
        self.eta = eta
        self.Pti = Pti # initial tank pressure [Pa]
#        self.Ptf = Ptf # final tank pressure [Pa]
        self.Vox = Vox # oxidizer volume [m^3]
        self.Do = Do
        self.Cd = Cd
        self.Lf = Lf # fuel length [m]
        self.Dfi = Dfi # initial fuel port diameter [m]
        self.rho_ox = rho_ox
        self.rho_f = rho_f
        self.Mox_fill = Vox*rho_ox # total oxidizer mass [kg]
        self.Dti = Dti # initial nozzle throat diameter [m]
        self.De = De # nozzle exit diameter [m]
        self.Pa = Pa # ambient pressure [Pa]
        self.rn = rn # nozzle throat regression rate [m]
        self.lmbd = (1 + np.cos( np.deg2rad(theta) ))/2 # nozzle coefficient. theta is nozzle opening half-angle [rad]
        self.func_cstr = Read_datset(cea_fldpath).gen_func("CSTAR")
        self.func_gamma = Read_datset(cea_fldpath).gen_func("GAMMAs_c")
    


if __name__ == "__main__":
#    cea_fldpath = os.path.join("cea_db", "N2O_PE", "csv_database")
    cea_fldpath = os.path.join("cea_db", "GOX_PE", "csv_database")
    dt = 0.01 # time interval [s]
    tb_init = 2.0 # initial firing duration time for iteration [s]
    Pc_init = 1.3e+6 # initial chamber pressure for iteration [Pa]
    eta = 0.8
    Pti = 4.8e+6 # initial tank pressure [Pa]
    Ptf = 2.0e+6 # final tank pressure [Pa]
    Vox = 440e-6 # oxidizer volume [m^3]
    Do = 3.0e-3
    Cd = 0.333
    Lf = 247e-3 # fuel length [m]
    Dfi = 35e-3 # initial fuel port diameter [m]
    rho_ox = 852.24 # oxidizer mass density [kg/m^3]
    rho_f = 1190 # fuel density [kg/m^3]
    Dti = 12.70e-3 # initial nozzle throat diameter [m]
    De = 22.0e-3 # nozzle exit diameter [m]
    Pa = 0.1013e-6 # ambient pressure [Pa]
    rn = 0.0 # nozzle throat regression rate [m]
    theta = 15 # nozzle opening half-angle [deg]
    
    inst = Single_tank(cea_fldpath, dt, tb_init, Pc_init, eta, Pti, Ptf, Vox, Do, Cd, Lf, Dfi, rho_ox, rho_f, Dti, De, Pa, rn, theta)
    inst.exe_sim()
    dat = inst.df
#    inst = double_tank(cea_fldpath, dt, tb_init, Pc_init, eta, Pti, Vox, Do, Cd, Lf, Dfi, rho_ox, rho_f, Dti, De, Pa, rn, theta)
