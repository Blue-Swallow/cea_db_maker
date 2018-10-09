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
    def __init__(self, cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pti, Ptf, Vt, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, Dti, De, Pa, rn, theta):
        self.cea_fldpath = cea_fldpath
        self.dt = dt # time interval [s]
        self.tb_init = tb_init # initial firing duration time for iteration [s]
        self.Pc_init = Pc_init # initial chamber pressure for iteration [Pa]
        self.a = a # pre-expotential coefficient of Gox [m^3/kg]
        self.n = n # expotential coefficient of Gox
        self.eta = eta
        self.Pti = Pti # initial tank pressure [Pa]
        self.Ptf = Ptf # final tank pressure [Pa]
        self.Vt = Vt # oxidizer volume [m^3]
        self.Do = Do
        self.Cd = Cd
        self.Np = Np # the number of fuel port
        self.Lf = Lf # fuel length [m]
        self.Dfi = Dfi # initial fuel port diameter [m]
        self.Dfo = Dfo # fuel outer diameter [m]
        self.rho_ox = rho_ox
        self.rho_f = rho_f
        self.Mox_fill = Vt*rho_ox # total oxidizer mass [kg]
        self.Dti = Dti # initial nozzle throat diameter [m]
        self.De = De # nozzle exit diameter [m]
        self.Pa = Pa # ambient pressure [Pa]
        self.rn = rn # nozzle throat regression rate [m]
        self.lmbd = (1 + np.cos( np.deg2rad(theta) ))/2 # nozzle coefficient. theta is nozzle opening half-angle [rad]
        self.func_cstr = Read_datset(cea_fldpath).gen_func("CSTAR")
        self.func_gamma = Read_datset(cea_fldpath).gen_func("GAMMAs_c")
        self.of_max = Read_datset(cea_fldpath).of.max()
#        self.Pc = np.array([])
#        self.Pt = np.array([])        
#        self.mox = np.array([])
#        self.Mox = np.array([])
#        self.Mf = np.array([])
#        self.Df = np.array([])
#        self.Gox = np.array([])
#        self.r = np.array([])
#        self.Dt = np.array([])        
#        self.mf = np.array([])
#        self.cstr = np.array([])
#        self.of = np.array([])
#        self.Pe = np.array([])
#        self.CF = np.array([])
#        self.F = np.array([])
#        self.I = np.array([])
        
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
#            pass
            self.tb = optimize.brentq(self.func_error_tb, tb_min, tb_max, xtol=1.0e-2, maxiter=maxiter, full_output=False)
        self.df = pd.DataFrame([], index=np.arange(0, self.tb+self.dt/2, self.dt))
        self.df["Pc"] = self.Pc
        self.df["Pt"] = self.Pt    
        self.df["Pg"] = self.Pg
        self.df["mox"] = self.mox
        self.df["Mox"] = self.Mox
        self.df["Vox"] = self.Vox
        self.df["vg"] = self.vg
        self.df["mf"] = self.mf
        self.df["Mf"] = self.Mf
        self.df["Df"] = self.Df
        self.df["Gox"] = self.Gox
        self.df["r"] = self.r
        self.df["Dt"] = self.Dt
        self.df["Pe"] = self.Pe
        self.df["CF"] = self.CF
        self.df["cstr"] = self.cstr
        self.df["of"] = self.of
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
        self.Pg = np.array([])
        self.Tt = np.array([])
        self.Tg = np.array([])
        self.mox = np.array([])
        self.Mox = np.array([])
        self.Vox = np.array([])
        self.vg = np.array([])
        self.Mf = np.array([])
        self.Df = np.array([])
        self.Gox = np.array([])
        self.r = np.array([])
        self.Dt = np.array([])        
        self.mf = np.array([])
        self.Pe = np.array([])
        self.CF = np.array([])
        self.cstr = np.array([])
        self.of = np.array([])
        self.F = np.array([])
        self.I = np.array([])
        self._tmp_Mox_ = 0
        t = 0
        pbar = tqdm(total=tb/self.dt)
        while(self._tmp_Mox_ <= self.Mox_fill):
            mox  = self.iterat_mox(t, tb)
            Pc = self._tmp_Pc_
            of = self._tmp_of_
            eps = np.power(self.De/self._tmp_Dt_, 2.0)
            Pe = self.iterat_Pe(of, Pc, eps)
            gamma = self.func_gamma(of, Pc)
            if Pe==0:
                Pe = self.Pa
            if Pc <= 0:
                Pc = Pe
            CF = np.sqrt((2*np.power(gamma, 2)/(gamma-1))*np.power(2/(gamma+1), (gamma+1)/(gamma-1))*(1-np.power(Pe/Pc,(gamma-1)/gamma))) + (Pe-self.Pa)*eps/Pc
            F = self.lmbd*CF*Pc*np.power(self._tmp_Dt_, 2)*np.pi/4
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



    def iterat_mox(self, t, tb, mox_min=0.0, maxiter=100):
        try:
            if t == 0:
                mox_init = 0
            else:
                mox_init = self.mox[-1]
            mox = optimize.newton(self.func_error_mox, mox_init, maxiter=maxiter, tol=1.0e-3, args=(t,tb))
        except:
            try:
                Pt_max = self.Pti
            except:
                Pt_max = self.Pgi
#            mox = optimize.brentq(self.func_error_mox, mox_min, np.sqrt(2*self.rho_ox*self.Pti), xtol=1.0e-3, maxiter=maxiter, args=(t,tb))
            mox = optimize.brentq(self.func_error_mox, mox_min, np.sqrt(2*self.rho_ox*Pt_max), xtol=1.0e-3, maxiter=maxiter, args=(t,tb))
        self.mox = np.append(self.mox, mox)        
        self._tmp_Mox_ = self.func_Mox(t, mox)
        self._tmp_Vox_ = self.func_Vox(t, mox)
        self.Pc = np.append(self.Pc, self._tmp_Pc_)
        self.of = np.append(self.of, self._tmp_of_)
        self.cstr = np.append(self.cstr, self._tmp_cstr_)
        self.Mox = np.append(self.Mox, self._tmp_Mox_)
        self.Vox = np.append(self.Vox, self._tmp_Vox_)
        self.vg = np.append(self.vg, self._tmp_vg_)
        self.Df = np.append(self.Df, self._tmp_Df_)
        self._tmp_Mf_ = self.func_Mf(t ,mox, self._tmp_Df_)
        self.Mf = np.append(self.Mf, self._tmp_Mf_)
        self.Gox = np.append(self.Gox, self._tmp_Gox_)
        self.r = np.append(self.r, self._tmp_r)
        self.mf = np.append(self.mf, self._tmp_mf_)
        self.Dt = np.append(self.Dt, self._tmp_Dt_)
        self.Pt = np.append(self.Pt, self._tmp_Pt_)
        self.Pg = np.append(self.Pg, self._tmp_Pg_)
        self.Tt = np.append(self.Tt, self._tmp_Tt_)
        self.Tg = np.append(self.Tg, self._tmp_Tg_)


    def func_error_mox(self, mox, t, tb):
        mox_cal = self.func_mox(mox, t, tb)
        diff = mox - mox_cal
#        error = diff/mox_cal
        error = diff/mox
        return(error)
    
    def func_mox(self, mox, t, tb):
        Pt = self.func_Pt(t, tb, mox)
        Pc = self.iterat_Pc(t, tb, mox)
        if Pt-Pc < 0:
#            mox_cal = -1*self.Cd*(np.pi*np.power(self.Do, 2)/4)*np.sqrt(2*self.rho_ox*np.abs(Pt-Pc))
            mox_cal = 0
        else:
            mox_cal = self.Cd*(np.pi*np.power(self.Do, 2)/4)*np.sqrt(2*self.rho_ox*(Pt-Pc))
        self._tmp_mox_ = mox
        return(mox_cal)
#        mf = self.func_mf(t, mox)
#        if mf == 0 and mox>0:
#            of = self.of_max
#        elif mox==0 and mox==0:
#            of = 0
#        else:
#            of = mox/mf
#        if of > self.of_max:
#            of = self.of_max
#        self._tmp_of_ = of
#        cstr = self.eta*self.func_cstr(of, Pc)
#        self._tmp_cstr_ = cstr
#        At = self.func_At(t)
#        Pc_cal = cstr*(mox + mf)/At
#        return(Pc_cal)
    

    
    def iterat_Pc(self, t, tb, mox, Pc_min=0.2e+6, maxiter=100):
        try:
            if t == 0:
                Pc_init = self.Pc_init
            else:
                Pc_init = self.Pc[-1]
            Pc = optimize.newton(self.func_error_Pc, Pc_init, maxiter=maxiter, tol=1.0e-3, args=(t,tb, mox))
        except:
            Pc = optimize.brentq(self.func_error_Pc, Pc_min, self.Pti, xtol=1.0e-3, maxiter=maxiter, args=(t,tb, mox))
        self._tmp_Pc_ = Pc
#        self.of = np.append(self.of, self._tmp_of_)
#        self.cstr = np.append(self.cstr, self._tmp_cstr_)
#        self.mox = np.append(self.mox, self._tmp_mox_)
#        self.Mox = np.append(self.Mox, self._tmp_Mox_)
#        self.Mf = np.append(self.Mf, self._tmp_Mf_)
#        self.Df = np.append(self.Df, self._tmp_Df_)
#        self.Gox = np.append(self.Gox, self._tmp_Gox_)
#        self.r = np.append(self.r, self._tmp_r)
#        self.mf = np.append(self.mf, self._tmp_mf_)
#        self.Dt = np.append(self.Dt, self._tmp_Dt_)
#        self.Pt = np.append(self.Pt, self._tmp_Pt_)
        return(Pc)

        
    def func_error_Pc(self, Pc, t, tb, mox):
        Pc_cal = self.func_Pc(Pc, t, tb, mox)
        diff = Pc - Pc_cal
        error = diff/Pc_cal
        return(error)
    
    def func_Pc(self, Pc, t, tb, mox):
#        mox = self.func_mox(t, Pc, tb)
        Df = self.iterat_Df(t, mox)
        mf = self.func_mf(t, mox, Df)
        if mf <= 0 and mox>0:
            of = self.of_max
        elif mox <= 0:
            of = 0
        else:
            of = mox/mf
        if of > self.of_max:
            of = self.of_max
        self._tmp_of_ = of
        cstr = self.eta*self.func_cstr(of, Pc)
        self._tmp_cstr_ = cstr
        At = self.func_At(t)
        Pc_cal = cstr*(mox + mf)/At
        return(Pc_cal)




    def iterat_Df(self, t, mox, maxiter=100):
        try:
            if t <= 0:
                Df_init = self.Dfi
            else:
                Df_init = self.Df[-1]
            Df = optimize.newton(self.func_error_Df, Df_init, maxiter=maxiter, tol=1.0e-3, args=(t, mox))
        except:
            if t <= 0:
                Df_init = self.Dfi
            else:
                Df = optimize.brentq(self.func_error_Df, self.Dfi, self.Dfo, xtol=1.0e-3, maxiter=maxiter, args=(t, mox))
        self._tmp_Df_ = Df
        return(Df)

    def func_error_Df(self, Df, t, mox):
        Df_cal = self.func_Df(Df, t, mox)
        diff = Df - Df_cal
        error = diff/Df_cal
        return(error)    
    
    def func_Df(self, Df, t, mox):
        dMf = self.func_Mf(t, mox, Df) - self.func_Mf(t-self.dt, mox, Df)
        if t <= 0:
            Df_cal = Df
        else:
            Df_cal = np.sqrt(4*dMf/(self.Np*np.pi*self.rho_f*self.Lf) + np.power(self.Df[-1], 2))
#        self._tmp_Df_ = Df
        return(Df_cal)


    def func_Mf(self, t, mox, Df):
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
        t_list = np.arange(0, t+self.dt/2, self.dt)
        if len(t_list) != len(self.mf):
            mf_list = np.append(self.mf, self.func_mf(t, mox , Df))
        else:
            mf_list = self.mf
        if t <= 0:
            Mf = 0
            t_list = np.array([0,])
        else:
#            tmp_mf = np.append(self.mf, mf)
#            Mf = integrate.simps(tmp_mf, t_list)
            Mf = integrate.simps(mf_list, t_list)
#        self._tmp_Mf_ = Mf
        return(Mf)


    def func_mf(self, t, mox, Df):
        if mox <=0:
            r = 0
        else:
            r = self.func_regression(Df, mox)
        mf = self.Np*self.Lf*np.pi*Df*self.rho_f*r
        self._tmp_mf_ = mf
        return(mf)

    def func_regression(self, Df, mox):
        Gox = self.func_Gox(Df, mox)
        r = self.a*np.power(Gox, self.n)
        self._tmp_r = r
        return(r)

    def func_Gox(self, Df, mox):
        Gox = 4*mox/(np.pi*np.power(Df,2))
        self._tmp_Gox_ = Gox
        return(Gox)


    def func_At(self, t):
        Dt = self.Dti -2*self.rn*t
        At = np.pi*np.power(Dt, 2)/4
        self._tmp_Dt_ = Dt
        return(At)


    


#==============================================================================
#     def func_mox(self, t, Pc, tb):
#         """ Calculate oxidizer mass flow rate
#         
#         Parameter
#         ---------
#         Cd: float
#             discharge coefficient [-]
#         Do: float
#             orifice diameter [mm]
#         rho_ox: float
#             oxidizer density [kg/m^3]
#         t: float
#             time [s]
#         Pt: float
#             tank pressure [Pa]
#         Pc: float
#             camber pressure [Pa]
#         """
#         Pt = self.func_Pt(t, tb)
#         mox = self.Cd*(np.pi*np.power(self.Do, 2)/4)*np.sqrt(2*self.rho_ox*(Pt-Pc))
#         Mox =self.func_Mox(t, mox) # calculate total oxidizer mass flow
#         self._tmp_mox_ = mox
# #        self.mox = np.append(self.mox, self._tmp_mox)
#         return(mox)
#==============================================================================

    def func_Pt(self, t, tb ,mox):
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


    def func_Mox(self, t , mox):
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
        t_list = np.arange(0, t+self.dt/2, self.dt)
        if len(t_list) != len(self.mox):
            mox_list = np.append(self.mox, mox)
        else:
            mox_list = self.mox
        if t==0:
            Mox = 0
        else:
            Mox = integrate.simps(mox_list, t_list)
        return(Mox)
    
    def func_Vox(self, t, mox):
        if t == 0:
            Vox = self.Vt
        else:
            Vox = self.Vox[-1] - (self.func_Mox(t, mox) - self.Mox[-1])/self.rho_ox
        self._tmp_Vox_ = Vox
        return(Vox)


class Double_tank(Single_tank):
    """ Class to simulate when the system has pressurize gas tank
    
    Parameter
    ----------    
    zeta: float
        combined pressure loss coefficient.  
        For example:
        in the case of inlet of pipe, the coefficient is 0.5 but it depends on inlet shape;
        in the case of outlet of pipe, the coefficient is nearly equal to 1.0;
        in the case of elbow, the coefficient is 1.129. 
    """
    
    def __init__(self, cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pgi, Vt, Vg, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, rho_g, mu_g, gamma_g, Dpg, Lpg, eps, zeta, Dti, De, Pa, rn, theta):
        self.dt = dt # time interval [s]
        self.tb_init = tb_init # initial firing duration time for iteration [s]
        self.Pc_init = Pc_init # initial chamber pressure for iteration [Pa]
        self.a = a # pre-expotential coefficient of Gox [m^3/kg]
        self.n = n # expotential coefficient of Gox
        self.eta = eta
#        self.Pti = Pti # initial tank pressure [Pa]
#        self.Ptf = Ptf # final tank pressure [Pa]
        self.Pgi = Pgi # initial gass tank pressure [Pa]
        self.Vt = Vt # oxidizer volume [m^3]
        self.Vg = Vg # pressurigze gass tank volume [m^3]
        self.Do = Do
        self.Cd = Cd
        self.Np = Np # the number of fuel port
        self.Lf = Lf # fuel length [m]
        self.Dfi = Dfi # initial fuel port diameter [m]
        self.Dfi = Dfi # initial fuel port diameter [m]
        self.rho_ox = rho_ox
        self.rho_f = rho_f
        self.rho_g = rho_g # density of pressurize gas [kg/m^3]
        self.mu_g = mu_g # viscosity of pressurize gas [Pa-s]
        self.gamma_g = gamma_g # specific heat ratio of pressurize gas [-]
        self.Dpg = Dpg # pipe diameter of pressurize gas [m]
        self.Lpg = Lpg # straight pipe length [m]
        self.eps = eps # equivalent sand roughness [m]
        self.zeta = zeta # combined pressure loss coefficient
        self.Mox_fill = Vt*rho_ox # total oxidizer mass [kg]
        self.Dti = Dti # initial nozzle throat diameter [m]
        self.De = De # nozzle exit diameter [m]
        self.Pa = Pa # ambient pressure [Pa]
        self.rn = rn # nozzle throat regression rate [m]
        self.lmbd = (1 + np.cos( np.deg2rad(theta) ))/2 # nozzle coefficient. theta is nozzle opening half-angle [rad]
        self.func_cstr = Read_datset(cea_fldpath).gen_func("CSTAR")
        self.func_gamma = Read_datset(cea_fldpath).gen_func("GAMMAs_c")
        self.of_max = Read_datset(cea_fldpath).of.max()

    def func_Pt(self, t, tb ,mox):
#        Vox_ = self._tmp_Vox_ # volume of oxidizer in tank when t=t_i
        Vox =  self.func_Vox(t, mox)# volume of oxidizer in tank when t=t_i+1
        if t==0:
            Pt = self.Pgi
        else:
            Pt = np.power((self.Vg+self.Vt-self.Vox[-1])/(self.Vg+self.Vt-Vox) , self.gamma_g) *self.Pt[-1]
        self._tmp_Pt_ = Pt
        return(Pt)
    
#==============================================================================
#     def func_Pc(self, Pc, t, tb):
#         rho_g = self.rho_g
#         mu_g = self.mu_g
#         Pg = self.iterat_Pg(Pc, rho_g, mu_g)
#         Vol_ox = self.iterat_Vol_ox(Pg, mu_g, rho_g)
#         Pt = Pg - self.func_dPg(Vol_ox, rho_g, mu_g)
#         self._tmp_Pt_ = Pt
#         mox = self.func_mox(Pg, Pc, Vol_ox, rho_g, mu_g)
#         Mox =self.func_Mox(t, mox) # calculate total oxidizer mass flow
#         mf = self.func_mf(t,mox)
#         if mf == 0 and mox>0:
#             of = self.of_max
#         elif mox==0 and mox==0:
#             of = 0
#         else:
#             of = mox/mf
#         if of > self.of_max:
#             of = self.of_max
#         self._tmp_of_ = of
#         cstr = self.eta*self.func_cstr(of, Pc)
#         self._tmp_cstr_ = cstr
#         At = self.func_At(t)
#         Pc_cal = cstr*(mox + mf)/At
#         return(Pc_cal)
# 
# 
#     def iterat_Pg(self, Pc, rho_g, mu_g):
#         """ Iteration of gas tank pressure [Pa]
#         """
#         Pg = optimize.newton(self.func_error_mox, 5.0e+6, args=(Pc, rho_g, mu_g))
#         return(Pg)
#         
#     def func_error_mox(self, Pg, Pc, rho_g, mu_g):
#         Vol_ox = self.iterat_Vol_ox(Pg, mu_g, rho_g) # oxidizer volume flow rate
#         mox = self.rho_ox*Vol_ox/self.dt # oxidizer mass flow rate assumed by volume flow rate
#         mox_cal = self.func_mox(Pg, Pc, Vol_ox, rho_g, mu_g) # calculated oxidizer mass flow rate
#         diff = mox_cal - mox
#         error = diff/mox_cal
#         return(error)
# 
#     def func_mox(self, Pg, Pc, Vol_g, rho_g, mu_g):
#         dPg = self.func_dPg(Vol_g, rho_g, mu_g)
#         mox = self.Cd*(np.pi*np.power(self.Do, 2)/4)*np.sqrt(2*self.rho_ox*((Pg-dPg) - Pc))
#         self._tmp_mox_ = mox
#         return(mox)
# 
#     def iterat_Vol_ox(self, Pg, mu_g, rho_g):
#         """ Iteration of Vt: oxidizer volume flow rate [m^3/s]
#         """
# #        Vol_ox = optimize.brentq(self.func_error_Vol_ox)
#         Vol_ox = optimize.newton(self.func_error_Vol_ox, 1.0e-3, args=(Pg, mu_g, rho_g))
#         return(Vol_ox)
#         
#     def func_error_Vol_ox(self, Vol_ox, Pg, mu_g, rho_g):
#         Vol_ox_cal = self.func_Vol_ox(Vt, Pg, mu_g, rho_g)
#         diff = Vol_ox_cal - Vol_ox
#         error = diff/Vol_ox_cal
#         return(error)
#         
#     def func_Vol_ox(self, Vol_ox, Pg, mu_g, rho_g):
#         """ Calculate the oxidize volume flow rate in dt [m^3]
#         
#         Parameter
#         ---------
#         Vol_ox: float
#             oxidizer volume flow rate in dt [m^3]
#         Pg: float
#             gas tank pressure [Pa]
#         mu_g: float
#             viscosity of pressurize gas [Pa-s]
#         rho_g: float
#             density of pressurize gas [kg/m^3]
#         """
#         Vol_g = Vol_ox # pressurize gas volume in dt; that is equal to oxidizer volume in dt
#         dPg = self.func_dPg(Vol_g, rho_g, mu_g) # pressure loss throug pipe
#         Vol_ox_cal = (self.Pgi - Pg)*self.Vg/(Pg - dPg) # assuming isothermal changing *This should considers isentropic condition?.
#         return(Vol_ox_cal)
#     
#     def func_dPg(self, Vol_g, rho_g, mu_g):
#         """ Calculate the pressure loss between gas tank and oxidizer tank
#         
#         Parameter
#         ----------
#         Vol_ox: float
#             oxidizer volume flow rate in dt [m^3]
#         mu_g: float
#             viscosity of pressurize gas [Pa-s]
#         rho_g: float
#             density of pressurize gas [kg/m^3]
#         """
#         Ap = np.pi*np.power(self.Dpg, 2)/4 # closs section area of pressurize gas pipe [m^2]
#         ug = Vol_g/(Ap*self.dt) # gas velocity in the pipe
#         dP_zeta = self.zeta*rho_g*np.power(ug,2) # pressure loss of inlet, outlet and bended tube etc...
#         Re_g = rho_g*ug*self.Dpg/mu_g
#         if Re_g > 4000: # turbulent flow
#             lmbd = optimize.newton(self._func_lmbd_, 0.03, args=(Re_g,))
#         else: # laminar flow
#             lmbd = 64/Re_g
#         dP_pipe = lmbd*self.Lpg/(2*self.Dpg) * rho_g*np.power(ug,2) # pressure loss through straight pipe
#         return(dP_zeta + dP_pipe)
#         
#     
#     def _func_lmbd_(self, lmbd, Re_g):
#         """ Iterate calculation of pipe friction factor: lmbd, using Cole-brook empirical formula
#         
#         Parameter
#         ------------
#         lmbd: float
#             pipe friction factor
#         Re_g: float
#             Reynolds number of pressurzie gas in the pipe
#         """
#         ans = 1/np.sqrt(lmbd) + 2*np.log10(self.eps/(3.71*self.Dpg) + 2.51/(Re_g*np.sqrt(lmbd)))
#         return(ans)
#==============================================================================



class Double_tank_regulator_tmp(Single_tank):
    """ Class to simulate when the system has pressurize gas tank and regulator
    
    Parameter
    ----------    
    zeta: float
        combined pressure loss coefficient.  
        For example:
        in the case of inlet of pipe, the coefficient is 0.5 but it depends on inlet shape;
        in the case of outlet of pipe, the coefficient is nearly equal to 1.0;
        in the case of elbow, the coefficient is 1.129. 
    """

    def __init__(self, cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pgi, Tgi, Pt_reg, cv_reg, Vt, Vg, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, Rg, rho_g, mu_g, gamma_g, cpg, Dpg, Lpg, eps, zeta, Dti, De, Pa, rn, theta):
        self.dt = dt # time interval [s]
        self.tb_init = tb_init # initial firing duration time for iteration [s]
        self.Pc_init = Pc_init # initial chamber pressure for iteration [Pa]
        self.a = a # pre-expotential coefficient of Gox [m^3/kg]
        self.n = n # expotential coefficient of Gox
        self.eta = eta
#        self.Pti = Pti # initial tank pressure [Pa]
#        self.Ptf = Ptf # final tank pressure [Pa]
        self.Pgi = Pgi # initial gass tank pressure [Pa]
        self.Tgi = Tgi # initial gas temperature [K]
        self.Pt_reg = Pt_reg # ragulate pressure [Pa]
        self.cv_reg = cv_reg # Cv value of regulator [m^2]
        self.Vt = Vt # oxidizer volume [m^3]
        self.Vg = Vg # pressurigze gass tank volume [m^3]
        self.Do = Do
        self.Cd = Cd
        self.Np = Np # the number of fuel port
        self.Lf = Lf # fuel length [m]
        self.Dfi = Dfi # initial fuel port diameter [m]
        self.Dfo = Dfo # fuel outer diameter [m]
        self.rho_ox = rho_ox
        self.rho_f = rho_f
        self.Rg = Rg
        self.rho_g = rho_g # density of pressurize gas [kg/m^3]
        self.mu_g = mu_g # viscosity of pressurize gas [Pa-s]
        self.gamma_g = gamma_g # specific heat ratio of pressurize gas [-]
        self.cpg = cpg # specific heat capacity [J/kg/K]
        self.Dpg = Dpg # pipe diameter of pressurize gas [m]
        self.Lpg = Lpg # straight pipe length [m]
        self.eps = eps # equivalent sand roughness [m]
        self.zeta = zeta # combined pressure loss coefficient
        self.Mox_fill = Vt*rho_ox # total oxidizer mass [kg]
        self.Dti = Dti # initial nozzle throat diameter [m]
        self.De = De # nozzle exit diameter [m]
        self.Pa = Pa # ambient pressure [Pa]
        self.rn = rn # nozzle throat regression rate [m]
        self.lmbd = (1 + np.cos( np.deg2rad(theta) ))/2 # nozzle coefficient. theta is nozzle opening half-angle [rad]
        self.func_cstr = Read_datset(cea_fldpath).gen_func("CSTAR")
        self.func_gamma = Read_datset(cea_fldpath).gen_func("GAMMAs_c")
        self.of_max = Read_datset(cea_fldpath).of.max()
        
    def func_Pt(self, t, tb ,mox):
        Vox =  self.func_Vox(t, mox) # volume of oxidizer in tank when t=t_i+1
        if t==0:
            self._tmp_Pg_ = self.Pgi
        else:
            pass
        if self._tmp_Pg_ > self.Pt_reg:
            Pt = self.Pt_reg
            Pg = self.func_Pg(t, Vox)
        else:
            Pt = np.power((self.Vg+self.Vt-self.Vox[-1])/(self.Vg+self.Vt-Vox) , self.gamma_g) *self.Pt[-1]
            Pg = Pt
        self._tmp_Pt_ = Pt
        self._tmp_Pg_ = Pg
        self.func_vg(t, mox, Pt)
        return(Pt)
        
    def func_Pg(self, t, Vox):
        ni = self.Pgi*self.Vg/(self.Rg*self.Tgi)
        Tt = np.power(self.Pgi/self.Pt_reg, (1-self.gamma_g)/self.gamma_g)*self.Tgi
        self._tmp_Tt_ = Tt
#        Vox = self.func_Vox(t, mox)
        Tg = self.Tgi - self.Pt_reg*(self.Vt - Vox)/self.cpg
        self._tmp_Tg_ = Tg
        Pg = (ni*self.Rg - self.Pt_reg*(self.Vt - Vox)/Tt)*(Tg/self.Vg)
        return(Pg)
        
    def func_vg(self, t, mox, Pt):
        if t == 0:
            dVox = self.func_Mox(t, mox)/self.rho_ox
        else:
            dVox = (self.func_Mox(t, mox) - self.Mox[-1])/self.rho_ox
        self._tmp_vg_ =  dVox # [m^3/s] gas mass flow volume
        return(self._tmp_vg_)
        


if __name__ == "__main__":
#    cea_fldpath = os.path.join("cea_db", "N2O_PE", "csv_database")
    cea_fldpath = os.path.join("cea_db", "LOX_PE", "csv_database")
    dt = 0.01 # time interval [s]
    tb_init = 10.0 # initial firing duration time for iteration [s]
    Pc_init = 3.0e+6 # initial chamber pressure for iteration [Pa]
#    a = 1.17063e-4 # regression coefficient [m^3/kg]
    a = 4.19347377e-5 # regression coefficient [m^3/kg]
#    n = 0.62 # oxidizer mass flux exponent
    n = 0.498 # oxidizer mass flux exponent
#    eta = 0.7
    eta = 0.8
    Pti = 5.0e+6 # initial tank pressure [Pa]
    Ptf = 1.5e+6 # final tank pressure [Pa]
    Vt = 15.0e-3 # oxidizer volume [m^3]
    Do = 1.0e-3
#    Cd = 0.82 *64
    Cd = 0.82 *56
#    Np = 1
#    Np = 5
    Np = 13
#    Lf = 820e-3 # fuel length [m]
    Lf = 1045e-3 # fuel length [m]
#    Dfi = 100e-3 # initial fuel port diameter [m]
#    Dfi = 20e-3 # initial fuel port diameter [m]
    Dfi = 2*13.867e-3 # initial fuel port diameter [m]
    Dfo = 200e-3 # fuel outer diameter [m]
    rho_ox = 1190.0 # oxidizer mass density [kg/m^3]
#    rho_f = 820 # fuel density [kg/m^3]
    rho_f = 950 # fuel density [kg/m^3]
    Dti = 44.0e-3 # initial nozzle throat diameter [m]
    De = 103.0e-3 # nozzle exit diameter [m]
    Pa = 0.1013e+6 # ambient pressure [Pa]
    rn = 0.0 # nozzle throat regression rate [m]
    theta = 10 # nozzle opening half-angle [deg]
#    inst = Single_tank(cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pti, Ptf, Vt, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, Dti, De, Pa, rn, theta)
#    inst.exe_sim()
#    dat = inst.df
    
    Pgi = 15.0e+6 # initial gas tank pressure [Pa]
    Vg = 3.4e-3 # pressurize gas tank volume [m^3]
    rho_g = 0.1786 # pressurize gas density [kg/m^3]
    mu_g = 0.0186e-3 # pressurize gas viscosity [Pa-s]
    gamma_g = 1.4 # specific heat ratioi of pressuirze gas [-]
    Dpg = (25.4/4)*1.0e-3 # pipe diameter of pressurize gas [m]
    Lpg = 200e-3 # straight pipe length [m]
    eps = 0.0015e-3 # equivalent sand roughness [m]
    zeta = 0.25 # combined pressure loss coefficient
#    inst2 = Double_tank(cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pgi, Vt, Vg, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, rho_g, mu_g, gamma_g, Dpg, Lpg, eps, zeta, Dti, De, Pa, rn, theta)
#    inst2.exe_sim()
#    dat2 = inst2.df
    
    Tgi = 300 # initial gas tank temperature [K]
    Pt_reg = 5.0e+6 # regulator set pressure [Pa]
    cv_reg = 0.5e-6 # Cv value of ragulator [m^2]
    Rg = 8.314/4.002602 # specific gas constant of pressurize gas [J/kg/K]
#    inst3 = Double_tank_regulator(cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pgi, Tgi, Pt_reg, cv_reg, Vt, Vg, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, Rg, rho_g, mu_g, gamma_g, Dpg, Lpg, eps, zeta, Dti, De, Pa, rn, theta)
#    inst3.exe_sim()
#    dat3 = inst3.df
    
    cpg = 5237 # specific heat capacity [J/kg/K]
    inst4 = Double_tank_regulator_tmp(cea_fldpath, dt, tb_init, Pc_init, a, n, eta, Pgi, Tgi, Pt_reg, cv_reg, Vt, Vg, Do, Cd, Np, Lf, Dfi, Dfo, rho_ox, rho_f, Rg, rho_g, mu_g, gamma_g, cpg, Dpg, Lpg, eps, zeta, Dti, De, Pa, rn, theta)
    inst4.exe_sim()
    dat4 = inst4.df
    
#    import matplotlib.pyplot as plt
#    Pc_range = np.arange(0.2e+6, 10.0e+6, 0.1e+6)
#    val = np.array([inst4.func_error_Pc(Pc, 0, 15, 0) for Pc in Pc_range])
#    plt.plot(Pc_range, val)
