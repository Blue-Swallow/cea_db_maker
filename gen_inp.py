# -*- coding: utf-8 -*-
"""
Generate NASA-CEA input file "*.inp"
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm

cond_name = "cond.txt"


def _getpath_():
    """
    Return the folder path cantaining "cond.txt"
    """
    cadir = os.path.dirname(os.path.abspath(__file__))
    print("Input a Case Name (Folder Name) >>")
    foldername = input()
    path = cadir + "/{}".format(foldername)
    return(path)


def read_cond(path):
    """
    Read condition file "cond.txt" containing information to generate "*.inp"
    
    Parameters
    ----------
    path : string, folder path to generate folder "inp"
    
    Returns
    -------
    oxid : string
        Containing oxidizer infromation: name, mass fraction, initial temperature
    
    fuel : string
        Containing fuel infromation: name, mass fraction, initial temperature
        
    dh : tuple (x, y)
        Standard enthalpy of formation, [kJ/mol]
        x : enthalpy of termanal molecule
        y : enthalpy of monomer molecule
        
    Pc : list [Pc_start, Pc_end, Pc_interval]
        Chamber pressure, [MPa]
    
    of : list [O/F_start, O/F_end, O/F_interval]
        Ratio of oxidizer to fuel, [-]
    
    n : list, [int, ...]
        Polymerization number
    
    elem : list containing 3elm-tuple [(X, a, b), ...]
        List of elements contained in fuel molecule
        X: string, symbol of element, e.g. "C", "O", "H"
        a: string, the number of element contained in a terminal molecule
        b: string, the number of element contained in a monomer molecule
    """
    fpath = path+"/{}".format(cond_name)
    file = open(fpath,"r")

    value_inp = {"oxid": 1,
                 "oxid_wt": 1,
                 "oxid_t,k": 1,
                 "fuel": 1,
                 "fuel": 1,
                 "fuel_wt": 1,
                 "fuel_t,k": 1}
    range_inp = ("h,kj/mol","O/F","n","Pc")

    line = str(0)
    elem = []
    while line:
        line = file.readline()
        dat = line.split()
        if(len(dat)==0 or dat[0]=="#"):
            pass
        elif(dat[0] in value_inp):
            value_inp[dat[0]]=dat[1]    #get out-put values & through into "value_out"
        elif(dat[0] in range_inp):
            if(dat[0]==range_inp[0]):
                dh = tuple(map(float,dat[1:]))
            elif(dat[0]==range_inp[1]):
                of = (float(dat[1]), float(dat[2]), float(dat[3]))
            elif(dat[0]==range_inp[2]):
                n = tuple(map(int,dat[1:]))
            else:
                Pc = (float(dat[1]), float(dat[2]), float(dat[3]))
        else:
            elem = elem + [(dat[0], dat[1], dat[2]),]
    file.close()
    
    oxid = "oxid={} wt={} t,k={}".format(value_inp["oxid"], value_inp["oxid_wt"], value_inp["oxid_t,k"])
    fuel = "fuel={} wt={} t,k={}".format(value_inp["fuel"], value_inp["fuel_wt"], value_inp["fuel_t,k"])
#    Pc = float(value_inp["Pc"])

    return(oxid,fuel,dh,Pc,of,n,elem)


def make_inp(path, option, of, Pc, oxid, fuel, h, elem, eps, n=""):
    """
    Write information in input file, "*.inp".
    
    Parameter
    ---------
    path : string
        Path where folders containing "*.inp" are created 
    option: string
        Calculation option, wtheher using equilibrium composition or frozen composition.
    of: float,
        O/F
    Pc: float,
        Camberpressure, [MPa]
    oxid: string,
        Information about oxidizer
    fuel: string
        Information about fuel
    h: float,
        Standard enthalpy of formation [kJ/mol]
    elem: string,
        Information about element and the number of each element
    eps: float,
        Area ratio of nozzle throat & exit, Ae/At
    n: int
        polimerization number
    """
    #check existence & create folder"inp"
#    fld_name = "inp_n={}".format(n) # folder name e.g. "n=100"

    if len(n) == 0:
        fld_name = "inp"
    else:
        fld_name = os.path.join("inp", "n={}".fromat(n))
        
    if os.path.exists(os.path.join(path, fld_name)):
        pass
    else:
        print("{},  {}".format(path, fld_name))
        os.makedirs(os.path.join(path, fld_name))

    inp_fname = "Pc_{:0>4}__of_{:0>4}.inp".format(round(Pc,1), round(of,1)) #.inp file name, e.g. "Pc=1.0_of=6.0.inp"
    file = open(path+"/"+fld_name+"/"+ inp_fname, "w")

    Pc = Pc * 10    #Pc:Chamber pressure [bar]
    prob = "case={} o/f={} rocket {} tcest,k=3800 p,bar={} sup,ae/at={}".format(inp_fname, round(of,3), option, round(Pc,3), round(eps,3))
    fuel = fuel + " h,kj/mol={} {} ".format(h, elem)
#    outp = "siunits short"
    outp = ""
    file.write("prob\n\t{}\nreact\n\t{}\n\t{}\noutput\n\t{}\nend\n".format(prob,oxid,fuel,outp))
    file.close()


def gen_all_cond(path):
    """
    Generate input file with respect to every condition
    
    Parameters
    ----------
    path: string
        Folder path where this function will make "inp" floder storing ".inp"

        
    """
    oxid,fuel,dh,Pc,of,n,elem  = read_cond(path)
    of = np.arange(of[0], of[1], of[2])
    Pc = np.arange(Pc[0], Pc[1], Pc[2])
    
    for k in range(len(n)):
        count = k
        h_input = dh[0] + dh[1]*n[k]
        elem_input = ""
        for i in range(len(elem)):
            elem_name = elem[i][0]
            elem_num = int(elem[i][1]) + int(elem[i][2])*n[k]
            elem_input = elem_input+"{} {} ".format(elem_name, elem_num)
        for i in range(np.size(Pc)):
            for j in range(np.size(of)):
                make_inp(path, of[j], Pc[i], oxid, fuel, h_input, elem_input, 1.0, n=n[k])
#                pass
    return(count)


class Cui_input():
    """
    Class to attract information through CUI to generate .inp file
    
    Parameters
    ----------
    langlist: list ["jp","en",...]
        Contain initial two characters of language name.
    
    sntns_df : dict {outer_key: innerdict, ...}
        outer_key: string\n
            it is representative of the questionaire type. \n
        innerdict: string\n
            it is a questionaier sentense.

    lang: string
        Selected language from the "langlist"
    
    oxid: string
        must use a symbol written in NASA RP-1311-P2 app.B
        
    fuel: string
        fuel name. It isn't necessary to select from NASA RP-1311-P2 app.B
    
    o_itemp: float
        Oxidizer initial temperature [K]

    f_itemp: float
        Fuel initial temperature [K]

    f_enthlpy: float
        Fuel satndard enthalpy of formation [kJ/mol]
    
    f_elem: dict, {"A": num, ...}
        A: string, symbol of contained element in fuel
        num: int, the number of element contained in 1 mol fuel [mol]
        
    eps: float
        Nozzle area ratio, Ae/At.
    """
    
    langlist = ["jp","en"]
    _tmp_ = dict()
#    _tmp_["oxid"] = {"jp": "酸化剤の種類と質量分率(%)を入力してください.\n*記号はNASA RP-1311-P2 app.B に準拠\n\n例)\nO2(L) 80, N20 20\n",
#                     "en": "Please input oxidizer name and its mass fraction (%).\n*Concerning spiecies name, please refer \"NASA RP-1311-P2 app.B.\"\n\ne.g.\nO2(L) 80, N20 20"}
    _tmp_["option"] = {"jp": "\n\n計算オプション(0~2)を選択してください．\n例: 0: 全域平衡計算\n    1: 燃焼器内のみ平衡計算\n    2: スロートまで平衡計算",
                     "en": "\n\nPlease select option (0-2) of calculation.\ne.g.: 0: equilibrium during expansion\n      1: frozen after the end of chamber\n      2: frozen after nozzle throat"}

    _tmp_["oxid"] = {"jp": "\n\n酸化剤の種類を入力してください.\n*記号はNASA RP-1311-P2 app.B に準拠\n例: O2(L)",
                     "en": "\n\nPlease input oxidizer name.\n*Concerning spiecies name, please refer \"NASA RP-1311-P2 app.B.\"\ne.g.: O2(L)"}

    _tmp_["fuel"] = {"jp": "\n\n燃料の名前を入力してください\n例: PMMA",
                     "en": "\n\nPlease input fuel name.\"\ne.g.: PMMA"}
    
    _tmp_["o_itemp"] = {"jp": "\n\n酸化剤の初期温度[K]を入力してください",
                     "en": "\n\nPlease input initial oxidizer temperature [K]"}
    
    _tmp_["f_itemp"] = {"jp": "\n\n燃料の初期温度[K]を入力してください",
                     "en": "\n\nPlease input initial fuel temperature [K]"}
    
    _tmp_["f_enthlpy"] = {"jp": "\n\n燃料の標準生成エンタルピ[kJ/mol]を入力してください",
                     "en": "\n\nPlease input standard enthalpy of formation [kJ/mol] respect to fuel"}
    
    _tmp_["f_elem"] = {"jp": "\n\n1molの燃料に含まれる元素とそのmol数(整数)を入力してください.\n例: C 5 H 2 O 6",
                     "en": "\n\nPlease input the element symbols and its mol number (integer) contained per 1 mol of fuel\ne.g.: C 5 H 2 O 6"}

    _tmp_["eps"] = {"jp": "\n\n開口比Ae/Atを入力してください.",
                    "en": "\n\nPlease input the area ratio, Ae/At."}
    
    _tmp_["of"] = {"jp": "\n\n計算するO/Fの範囲を入力してください.\n例) 0.5~10 を 0.1毎に計算する場合.\n0.5　10　0.1",
                "en": "\n\nPlease input the range of O/F where you want to calculate.\ne.g. If the range is 0.5 to 10 and the interval is 0.1\n0.5 10 0.1"}

    _tmp_["Pc"] = {"jp": "\n\n計算する燃焼室圧力[MPa]の範囲を入力してください.\n例) 0.5 MPa ~ 5.0 MPa を 0.1 MPa毎に計算する場合.\n0.5　5.0　0.1",
                "en": "\n\nPlease input the range of Chamber pressure [MPa] where you want to calculate.\ne.g. If the range is 0.5 to 5.0 MPa and the interval is 0.1 MPa\n0.5 5.0 0.1"}


    def __init__(self):
        self.sntns_df = pd.DataFrame([], columns=self.langlist)
        for i in self._tmp_:
            self.sntns_df = self.sntns_df.append(pd.DataFrame(self._tmp_[i], index=[i]))
        self._inp_lang_()
        self._inp_option_()
        self._inp_oxid_()
        self._inp_fuel_()
        self._inp_o_itemp_()
        self._inp_f_itemp_()
        self._inp_f_enthlpy_()
        self._inp_f_elem_()
        self._inp_eps_()
        self._inp_of_()
        self._inp_Pc_()

    def _inp_lang_(self):
        """
        Select user language
        """
        print("Please select language.\n{}".format(self.langlist))
        lang = input()
        if lang in self.langlist:
            self.lang = lang
        else:
            print("There is no such language set!")

    def _inp_option_(self):
        """
        Select a option of calculation: whether equilibrium or frozen composition.
        """
        print(self.sntns_df[self.lang]["option"])
        option = int(input())
        if option == 0:
            option = "equilibrium"
        elif option == 1:
            option = "frozen nfz=1" # frozen composition after the end of chamber
        elif option == 2:
            option = "frozen nfz=2" # frozen composition after nozzle throat
        else:
            print("Please re-confirm input integer!")
        self.option = option

    def _inp_oxid_(self):
        """
        Input oxidizer species
        """
        print(self.sntns_df[self.lang]["oxid"])
        oxid = input()
#        if oxid in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        self.oxid = oxid

    def _inp_fuel_(self):
        """
        Input fuel species
        """
        print(self.sntns_df[self.lang]["fuel"])
        fuel = input()
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        self.fuel = fuel

    def _inp_o_itemp_(self):
        """
        Input oxidizer initial temperature
        """
        print(self.sntns_df[self.lang]["o_itemp"])
        o_itemp = float(input())
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        self.o_itemp = o_itemp

    def _inp_f_itemp_(self):
        """
        Input fuel initial temperature
        """
        print(self.sntns_df[self.lang]["f_itemp"])
        f_itemp = float(input())
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        self.f_itemp = f_itemp
    
    def _inp_f_enthlpy_(self):
        """
        Input fuel standard enthalpy of formation
        """
        print(self.sntns_df[self.lang]["f_enthlpy"])
        f_enthlpy = float(input())
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        self.f_enthlpy = f_enthlpy
        
    def _inp_f_elem_(self):
        """
        Input element contained in fuel and its mol number contained in per 1 mol
        """
        print(self.sntns_df[self.lang]["f_elem"])
        self.f_elem = {}
        tmp = input().split()
        for i in range(len(tmp)):
            if i%2==1:
                self.f_elem[tmp[i-1]] = int(tmp[i])
                
    def _inp_eps_(self):
        """
        Input nozzle-area ratio
        """      
        print(self.sntns_df[self.lang]["eps"])
        self.eps = float(input())
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        
    def _inp_of_(self):
        """
        Input calculation range of O/F
        """      
        print(self.sntns_df[self.lang]["of"])
        self.of = list(map(lambda x: float(x) ,input().split()))
        

    def _inp_Pc_(self):
        """
        Input calculation range of chamber pressure, Pc.
        """      
        print(self.sntns_df[self.lang]["Pc"])
        self.Pc = list(map(lambda x: float(x) ,input().split()))

    def gen_all(self):
        """
        Generate input file with respect to every condition
        
        Parameters
        ----------
        path: string
            Folder path where this function will make "inp" floder storing ".inp"
    
        func: function    
        -----------------
        (oxid, fuel, dh, Pc, of, n, elem) = func(path)
         
        Parameters
            path: string\n
        Returns
            oxid: string, e.g."oxid=O2(L) wt=100 t,k=90.15"  \n
            fuel: string, e.g."fuel=PMMA wt=100 t,k=298.15"  \n
            dh: float  \n
            Pc: list, [start, end, interval], each element type is float  \n
            of: list, [start, end, interval], each element type is float  \n
            
        """
        path = _getpath_()
#        oxid,fuel,dh,Pc,of,n,elem  = read_cond(path)
        of = np.arange(self.of[0], self.of[1], self.of[2])
        Pc = np.arange(self.Pc[0], self.Pc[1], self.Pc[2])
        oxid = "oxid={} wt={} t,k={}".format(self.oxid, 100, self.o_itemp)
        fuel = "fuel={} wt={} t,k={}".format(self.fuel, 100, self.f_itemp)

        elem_input = ""
        for i in tqdm(self.f_elem):
            elem_input = elem_input+"{} {} ".format(i, self.f_elem[i])
        for i in range(np.size(Pc)):
            for j in range(np.size(of)):
                make_inp(path, self.option, of[j], Pc[i], oxid, fuel, self.f_enthlpy, elem_input, self.eps)

    
if __name__ == "__main__":
#    path = _getpath_()
#    test = gen_all(path)

    myclass = Cui_input()
    myclass.gen_all()
