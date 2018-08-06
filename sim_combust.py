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

class Main():
    """ Simulate combustion of hybrid rocket
    
    Parameter
    ---------
    
    Return
    ---------
    """
    def __init__(self, cea_fldpath, dt, tb_init, Pc_init, eta, Pti, Ptf, Vox, Do, Cd, Lf, Dfi, rho_ox, rho_f, Dti, De, Pa, rn, theta):
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
        self.iterat_tb(tb_init=self.tb_init, tb_min=1.0, tb_max=10.0, maxiter=100)

    def iterat_tb(self, tb_init=5.0, tb_min=1.0, tb_max=10, maxiter=100):
        """ Iteration of tb: firing duration time
        
        Parameter
        ---------
        
        Return
        --------
        """
        try:
            self.tb = optimize.newton(self.func_error_tb, tb_init, maxiter=maxiter, tol=1.0e+2)
        except:
            self.tb = optimize.brentq(self.func_error_tb, tb_min, tb_max, maxiter=maxiter, xtol=1.0e+2, full_output=False)
        self.df = pd.DataFrame([], index=np.arange(0, self.tb+self.dt, self.dt))
        self.df["Pc"] = self.Pc
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
        return(self.df)
    
    def func_error_tb(self, tb):
        tb_cal = self.func_tb(tb)
        diff = tb_cal - tb
        error = diff/tb_cal
        return(error)
        
    def func_tb(self, tb):
        self.Pc = np.array([])
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
        self._tmp_Mox_ = 0
        t = 0
        while(self._tmp_Mox_ <= self.Mox_fill):
            Pc = self.iterat_Pc(t, tb)
            of = self._tmp_of_
            eps = np.power(self.De/self._tmp_Dt_, 2.0)
            Pe = self.iterat_Pe(of, Pc, eps)
            gamma = self.func_gamma(of, Pc)
            if Pe==0:
                Pe = self.Pa
            if Pc == 0:
                Pc = Pe
            CF = np.sqrt((2*np.power(gamma, 2)/(gamma-1))*np.power(2/(gamma+1), (gamma+1)/(gamma-1))*(1-np.power(Pe/Pc,(gamma-1)/gamma))) + (Pe-self.Pa)*eps/Pc
            F = self.lmbd*CF*Pc*np.power(self._tmp_Dt_, 2)
            t = t + self.dt
            print("t = {}s".format(t))
            self.Pe = np.append(self.Pe, Pe)
            self.CF = np.append(self.CF, CF)        
            self.F = np.append(self.F, F)
        tb = t - self.dt
        self.I = integrate.simps(self.F, np.arange(0,tb+self.dt/2, self.dt))
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
            Pe = optimize.newton(self.func_error_eps, Pc/2, maxiter=100, tol=1.0e+3, args=(of, Pc, eps))
        except:
            Pe = optimize.brentq(self.func_error_eps, 1, Pc/2, maxiter=100, xtol=1.0e+3, full_output=False, args=(of, Pc,eps))
        self.Pe=Pe
        return(self.Pe)
    
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

    

    
    def iterat_Pc(self, t, tb, Pc_min=0.2e+6, maxiter=50):
        Pc = optimize.brentq(self.func_error_Pc, Pc_min, self.Pti, xtol=1.0e-2, maxiter=maxiter, args=(t,tb))        
        self.Pc = np.append(self.Pc, Pc)
        self.mox = np.append(self.mox, self._tmp_mox_)
        self.Mox = np.append(self.Mox, self._tmp_Mox_)
        self.Mf = np.append(self.Mf, self._tmp_Mf_)
        self.Df = np.append(self.Df, self._tmp_Df_)
        self.Gox = np.append(self.Gox, self._tmp_Gox_)
        self.r = np.append(self.r, self._tmp_r)
        self.mf = np.append(self.mf, self._tmp_mf_)
        self.Dt = np.append(self.Dt, self._tmp_Dt_)
        return(Pc)

        
    def func_error_Pc(self, Pc, t, tb):
        Pc_cal = self.func_Pc(Pc, t, tb)
        diff = Pc_cal - Pc
        error = diff/Pc_cal
        self.Pc_pre = Pc_cal
        return(error)
    
    def func_Pc(self, Pc, t, tb):
        mox = self.func_mox(t, Pc, tb)
        mf = self.func_mf(t, mox)
        if mf == 0:
            of = 0
        else:
            of = mox/mf
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
        Mox =self.func_Mox(t)
        if Mox <= self.Mox_fill:
            mox = self.Cd*(np.pi*np.power(self.Do, 2)/4)*np.sqrt(2*self.rho_ox*(Pt-Pc))   
        else:
            mox = 0
        self._tmp_mox_ = mox
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
        return(Pt)


    def func_Mox(self, t):
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
            t_list = np.arange(0, t-self.dt/2, self.dt)
            Mox = integrate.simps(self.mox, t_list)
            self._tmp_Mox_ = Mox
        return(Mox)
    

    


if __name__ == "__main__":
    cea_fldpath = os.path.join("cea_db", "N2O_PE", "csv_database")
    dt = 0.01 # time interval [s]
    tb_init = 2.8 # initial firing duration time for iteration [s]
    Pc_init = 1.3e+6 # initial chamber pressure for iteration [Pa]
    eta = 0.8
    Pti = 4.8e+6 # initial tank pressure [Pa]
    Ptf = 2.6e+6 # final tank pressure [Pa]
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
    
    inst = Main(cea_fldpath, dt, tb_init, Pc_init, eta, Pti, Ptf, Vox, Do, Cd, Lf, Dfi, rho_ox, rho_f, Dti, De, Pa, rn, theta)
    