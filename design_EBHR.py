# -*- coding: utf-8 -*-
"""
Created on Sun May 20 18:43:45 2018

@author: T.J.-LAB-PC
"""

import os
import warnings
import pandas as pd
import numpy as np
from scipy import optimize
from tqdm import tqdm

from execute import Read_datset



class Cui_input():
    """
    Class to attract information using CUI to calculate firing-test condition
    
    Class variable
    ----------
    self._langlist_: list ["jp","en",...]
        Contain initial two characters of language name.
    
    self._sntns_df_ : dict {outer_key: innerdict, ...}
        outer_key: string\n
            it is representative of the questionaire type. \n
        innerdict: string\n
            it is a questionaier sentense.

    self._lang_: string
        Selected language from the "_langlist_"
        
    self.ex_path: string
        Path containing experiment data
    
    self.ex_file: string
        The name of experimant data file (.csv)
        
    self.ex_param: dict {param_name: symbol}
        Parameter dictionary
        param_name: string; e.g., time [s], Oxidizer mass flow rate [kg/s], ...
        symbol: string; e.g., t, mox, ... 
        
    self.ex_df: DataFrame
        Data frame of measured parameter in an experiment
    
    self.input_param: dict{key, value}
        key: str, "mode"
        value: int, 1, 2, 3, 4, 5
        the value is reperesentative number,
        1: RT-1, Re-construction technique 1 
        2: RT-2, Re-construction technique 2
        3: RT-3, Re-construction technique 3 
        4: RT-4, Re-construction technique 4
        5: RT-5, Re-construction technique 5
        
        key: str, "Dt"
        value: float, nozzle throat diameter [m]
        
        key: str, "eps"
        value: float, nozzle expansion ratio [-]
        
        key: str, "Mf"
        value: float, fuel consumption [kg]
    
    self.cea_path: string
        Path containing the results of cea calculation
   
        
    """
    _langlist_ = ["jp","en"]
    _tmp_ = dict()
    _tmp_["oxid"] = {"jp": "\n\n計算オプション(0~2)を選択してください．\n例: 0: 全域平衡計算\n    1: 燃焼器内のみ平衡計算\n    2: スロートまで平衡計算",
                     "en": "\n\nPlease select option (0-2) of calculation.\ne.g.: 0: equilibrium during expansion\n      1: frozen after the end of chamber\n      2: frozen after nozzle throat"}

    def __init__(self):
        self._sntns_df_ = pd.DataFrame([], columns=self._langlist_)
        for i in self._tmp_:
            self._sntns_df_ = self._sntns_df_.append(pd.DataFrame(self._tmp_[i], index=[i]))
#        self._inp_lang_()
#        self._get_expath_()
        self._get_ceapath_()
#        self.cea_db = Read_datset(self.cea_path)
#        self.input_param = dict()
#        self._select_mode_()
#        self._input_nozzle_()
#        self._input_eps_()
#        self._input_consump_()
        

    def _inp_lang_(self):
        """
        Select user language
        """
        print("Please select language.\n{}".format(self._langlist_))
        lang = input()
        if lang in self._langlist_:
            self._lang_ = lang
        else:
            print("There is no such language set!")

    def _get_expath_(self):
        """
        Get the experiment folder path, file name and data and, create template file to contain an experimental data
            ex_path: string, folder path 
            ex_file: string, file name
            ex_df: data frame of experintal data
            ex_param: dictionary of experimental data
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        while(True):
            print("\nInput a Experiment Name (The Name of Folder Containing Experiment Data) >>")
            foldername = input()
            self.ex_path = os.path.join(cadir, foldername)
            if os.path.exists(self.ex_path):
                self.ex_file = "ex_dat.csv"
                file_path = os.path.join(self.ex_path, self.ex_file)
                p_name = ("time [s]", "Oxidizer mass flow rate [g/s]", "Thrust [N]", "Chamber pressure [MPaG]")
                symbol = ("t", "mox", "F", "Pc")
                self.ex_param = dict(zip(p_name, symbol))
                if os.path.exists(file_path):
                    self.ex_df = pd.read_csv(file_path,header=1, index_col=0)
                    self.ex_df.mox = self.ex_df.mox * 1.0e-3 #convert [g/s] to [kg/s]
                    self.ex_df.Pc = self.ex_df.Pc * 1.0e+6 + 0.1013 #convert [MPaG] to [Pa]
                    break
                else: # create template file
                    print("\nThere is no such a experiment data/n{}".format(file_path))
                    print("\nDo you want to make a template file ?\ny/n ?")
                    flag = input()
                    if flag == "y":
                        df = pd.DataFrame(self.ex_param, index=[0], columns=p_name)
                        df.to_csv(file_path, index= False)
                    elif flag == "n":
                        sys.exit()
            else:
                print("\nThere is no such a Folder\n{}".format(self.ex_path))           


    def _get_ceapath_(self):
        """
        Return the folder path cantaining the results of cea calculation.
            cea.path: string, folder path containing "out" folder
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        while(True):
            print("\nInput the Folder Name Containing Results of CEA >>")
            foldername = input()
            self.cea_path = os.path.join(cadir, foldername)
            self.cea_path = os.path.join(self.cea_path, "csv_database")
            if os.path.exists(self.cea_path):
                break
            else:
                print("There is no such a dataset folder/n{}".format(self.cea_path))
                
    def _select_mode_(self):
        """
        Select a calculation mode; RT-1,2,3,4,5,...
        """
        while(True):
            print("Please select calculation mode.\n")
            print("1: RT-1\n2: RT-2\n3: RT-3\n4: RT-4\n5: RT-5\n")
            mode = {1: "RT-1",
                    2: "RT-2",
                    3: "RT-3",
                    4: "RT-4",
                    5: "RT-5"}
            inp = int(input())
            if inp in mode.keys():
                self.input_param["mode"] = inp
                break
            else:
                print("There is no such a mode \"{}\"".format(inp))
                
    def _input_nozzle_(self):
        """
        Input the nozzle diameter [mm]
        """
        print("Please input nozzle throat diameter [mm]")
        self.input_param["Dt"] = float(input())*1.0e-3

    def _input_eps_(self):
        """
        Input the nozzle expansion ratio [-]
        """
        print("Please input nozzle expansion ratio")
        self.input_param["eps"] = float(input())
        
    def _input_consump_(self):
        """
        Input the fuel consumption [g]
        """
        print("Please input fuel consumption [g]")
        self.input_param["Mf"] = float(input())*1.0e-3








class Gen_excond_table:
    def __init__(self, d, N, Df, eta, rho_f, Rm, Tox, C1, C2, cea_path):
        self.a = 1 - N*np.power(d, 2)/np.power(Df, 2)
        self.Df = Df*1.0e-3
        self.eta = eta
        self.rho_f = rho_f
        self.Rm = Rm
        self.Tox = Tox
        self.C1 = C1
        self.C2 = C2
        self.cea_path = cea_path
        self.mox_range, self.Dt_range = self._input_range_() #generate the calculation range of mox[g/s] and Dt [mm]
#        self.table = self.gen_table()

    def _input_range_(self):
        print("\n\nInput the range of mox [g/s], oxidizer mass flow rate, where you want to generate the table." )
        print("\ne.g. If the range is 10.0 to 20.0 g/s and the interval is 1.0 g/s\n10.0 20.0 1.0")
        tmp = list(map(lambda x: float(x) ,input().split()))
        mox_range = np.arange(tmp[0], tmp[1], tmp[2])*1.0e-3 # generate range and convert [g/s] to [kg/s]

        print("\n\nInput the range of Dt [mm], nozzle throat diameter, where you want to generate the table." )
        print("\ne.g. If the range is 5.0 to 10.0 mm and the interval is 1.0 mm\n5.0 10.0 1.0")
        tmp = list(map(lambda x: float(x) ,input().split()))
        Dt_range = np.arange(tmp[0], tmp[1], tmp[2])*1.0e-3 # generate range and convert [mm] to [m]
        return(mox_range, Dt_range)

    def gen_table(self):
        self.df_Pc = pd.DataFrame({}, index=(self.mox_range*1.0e+3)) #convert the unit of mox to [g/s]
        self.df_Vox = pd.DataFrame({}, index=(self.mox_range*1.0e+3)) #convert the unit of mox to [g/s]
        self.df_Vf = pd.DataFrame({}, index=(self.mox_range*1.0e+3)) #convert the unit of mox to [g/s]
        self.df_mf = pd.DataFrame({}, index=(self.mox_range*1.0e+3)) #convert the unit of mox to [g/s]
        self.df_of = pd.DataFrame({}, index=(self.mox_range*1.0e+3)) #convert the unit of mox to [g/s]
        self.df_cstr = pd.DataFrame({}, index=(self.mox_range*1.0e+3)) #convert the unit of mox to [g/s]
        for Dt in tqdm(self.Dt_range, desc="Dt loop", leave=True):
            Pc = np.array([])
            Vox = np.array([])
            Vf = np.array([])
            mf = np.array([])
            of = np.array([])
            cstr = np.array([])
            for mox in tqdm(self.mox_range, desc="mox loop", leave=False):
                Pc_tmp = self.converge_Pc(mox, Dt)
                Vox_tmp = func_Vox(Pc_tmp, mox, self.Rm, self.Tox, self.Df, self.a)
                Vf_tmp = func_Vf(Vox_tmp, Pc_tmp, self.C1, self.C2)
                mf_tmp = func_mf(Vf_tmp, self.rho_f, self.Df, self.a)
                of_tmp = mox/mf_tmp
                func_cstr = gen_func_cstr(self.cea_path)
                cstr_tmp = func_cstr(of_tmp, Pc_tmp*1.0e-6)[0]
                Pc = np.append(Pc, Pc_tmp)
                Vox = np.append(Vox, Vox_tmp)
                Vf = np.append(Vf, Vf_tmp)
                mf = np.append(mf, mf_tmp)
                of = np.append(of, of_tmp)
                cstr = np.append(cstr, cstr_tmp)
            self.df_Pc[(Dt*1.0e+3)] = Pc*1.0e-6  #convert the unit of Pc to [MPa]
            self.df_Vox[(Dt*1.0e+3)] = Vox
            self.df_Vf[(Dt*1.0e+3)] = Vf*1.0e+3  #convert the unit of Vf to [mm/s]
            self.df_mf[(Dt*1.0e+3)] = mf*1.0e+3  #convert the unit of mox to [g/s]
            self.df_of[(Dt*1.0e+3)] = of
            self.df_cstr[(Dt*1.0e+3)] = cstr
        return(self.df_Pc)
        

    def converge_Pc(self, mox, Dt, Pc_init=0.1e+5, maxiter=100):
        print("Dt={:3.2f}mm, mox={:3.2f}g/s".format(Dt*1.0e+3, mox*1.0e+3))
        bracket_end = self._find_min_iterat_func_Pc_(mox, Dt)
        try:
            Pc = optimize.brentq(self._iterat_func_Pc_, Pc_init, bracket_end, maxiter=maxiter, xtol=1.0e-3, args=(mox, Dt), full_output=False)
        except RuntimeError:
            Pc = np.nan
            print("RuntimeError")
        except ValueError:
            Pc = np.nan
            print("ValueError")
            error_begin = self._iterat_func_Pc_(Pc_init, mox, Dt)
            error_end = self._iterat_func_Pc_(bracket_end, mox, Dt)
            if (error_begin*error_end>0):
                print("There is no solution")
                print("bracket begin: Pc={}, Pc_cal={}".format(Pc_init, error_begin))
                print("bracket end: Pc={}, Pc_cal={}".format(bracket_end, error_end))
        return(Pc)
    
    def _find_min_iterat_func_Pc_(self, mox, Dt):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            Pc_min = optimize.minimize_scalar(self._iterat_func_Pc_, args=(mox, Dt)).x
        return(Pc_min)
    
    def _iterat_func_Pc_(self, Pc, mox, Dt):
        Vox = func_Vox(Pc, mox, self.Rm, self.Tox, self.Df, self.a)
        Vf = func_Vf(Vox, Pc, self.C1, self.C2, n=1.0)
        of = mox/func_mf(Vf, self.rho_f, self.Df, self.a)
        func_cstr = gen_func_cstr(self.cea_path)
        Pc_cal = func_Pc_cal(of, Pc, mox, func_cstr, Dt, self.eta)
        diff = Pc_cal - Pc
        error = diff/Pc_cal
        return(error)


def func_Vox(Pc, mox, Rm, Tox, Df, a):
    Af = np.pi*np.power(Df, 2)/4
    Vox = mox*Rm*Tox/(Pc*Af*(1-a)) #[m/s]
    return(Vox)

def func_Vf(Vox, Pc, C1, C2, n=1.0):
    Vf = (C1/Vox + C2)*np.power(Pc, n) #[m/s]
    return(Vf)

def func_mf(Vf, rho_f, Df, a):
    Af = np.pi*np.power(Df, 2)/4
    mf = rho_f*Af*a*Vf
    return(mf)

def gen_func_cstr(cea_path):
    func = Read_datset(cea_path).gen_func("CSTAR")
    return(func)

def func_Pc_cal(of, Pc, mox, func_cstr, Dt, eta):
    cstr = func_cstr(of,Pc*1.0e-6)[0]
    At = np.pi*np.power(Dt, 2)/4
    Pc_cal = eta*cstr*mox*(1 + 1/of)/At
    return(Pc_cal)

if __name__ == "__main__":
    cea_path = Cui_input().cea_path
    d = 0.3 #[mm]
    N = 433 #[個]
    Df = 38 #[mm]
    eta = 0.86 #[-]
    rho_f = 1191 #[kg/m^3]
#    Rm = 259.8 #[J/kg/K] GOX
    Rm = 1.89e+2 #[J/kg/K] N2O
    Tox = 280 #[K]
#    C1 = 9.34e-8 # SI-unit
#    C2 = 2.46e-9 # SI-unit
#    C1 = 1.39e-7 # SI-unit
#    C2 = 1.61e-9 # SI-unit
    C1 = 6.76e-9 # SI-unit, N2O
    C2 = 5.36e-10 # SI-unit, N2O
    instance = Gen_excond_table(d, N, Df, eta, rho_f, Rm, Tox, C1, C2, cea_path)
    table = instance.gen_table()
    table.to_csv("Pc_vs_mox+Dt.csv")
    
# =============================================================================
#     import matplotlib.pylab as plt
#     mox = 10.0e-3 #[kg/s]
#     Dt = 9.0e-3 #[mm]
#     a = 1 - N*np.power(d, 2.0)/np.power(Df, 2.0)
#     Pc_range = np.arange(0.2, 5.0, 0.1)*1.0e+6
#     error = [instance._iterat_func_Pc_(x, mox, Dt) for x in Pc_range]
#     Vox = [func_Vox(x, mox, Rm, Tox, Df*1.0e-3, a) for x in Pc_range]
#     plt.plot(Pc_range*1.0e-6, Vox)
#     Vf = [func_Vf(Vox[i], Pc_range[i], C1, C2) for i in range(len(Pc_range))]
#     plt.plot(Pc_range*1.0e-6, Vf)
#     mf = [func_mf(Vf[i], rho_f, Df*1.0e-3, a) for i in range(len(Pc_range))]
#     plt.plot(Pc_range*1.0e-6, mf)
#     func_cstr = gen_func_cstr(cea_path)
#     cstr = [func_cstr(mf[i]/mox, Pc_range[i]*1.0e-6)[0] for i in range(len(Pc_range))]
#     plt.plot(Pc_range*1.0e-6, cstr)
#     of = mox/np.array(mf)
#     plt.plot(Pc_range*1.0e-6, of)
#     Pc_cal = np.array([func_Pc_cal(mox/mf[i], Pc_range[i], mox, func_cstr, Dt, eta) for i in range(len(Pc_range))])
# 
#     plt.plot(Pc_range*1.0e-6, Pc_cal*1.0e-6, label="Pc_cal")
#     plt.plot(Pc_range*1.0e-6, (Pc_cal-Pc_range)*1.0e-6, label="difference")
#     plt.plot(Pc_range*1.0e-6, error, label="error")
#     plt.plot(Pc_range*1.0e-6, (Pc_cal-Pc_range)/Pc_cal, label="difference divided by Pc_cal")
#     plt.plot(Pc_range*1.0e-6, (Pc_cal-Pc_range)/Pc_range, label="difference divided by Pc_range")
#     print(instance._find_min_iterat_func_Pc_(mox, Dt))
# =============================================================================

