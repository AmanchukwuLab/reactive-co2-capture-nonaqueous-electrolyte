from pyomo.environ import *
import pyomo.environ as pyo
import cyipopt
import itertools
import pandas as pd
from pyomo.util.infeasible import log_infeasible_constraints
import logging

#========this is a flash model that represents the vapor-liquid equilibrum for electrolyte systems==========
#=================specify the cases to be used==============================================================
lifetime_parameter = pd.read_csv('input_lifetime.csv')

#input_parameter = pd.read_csv('Analysis_two_10cases.csv')

life_case = case_number_to_be_changed

input_parameter = pd.read_csv('Analysis_two_updated.csv')
print(input_parameter)
case = 8   # this will be the only parameters to enter for choosing parameters to run
#definition of cases


anolyte_rou   = input_parameter.loc[case,'anolyte density']   #g/cm3
catholyte_rou = input_parameter.loc[case,'catholyte density'] #g/cm3
#====capturing reaction will be different==============================

"""
case = 0 if cathode as Sn, solvent as DMSO, amine as EOEA
case = 1 if cathode as Ag, solvent as Water,amine as EOEA
case = 2 if cathode as Sn, solvent as DMSO, amine as MEA
a potential case is cathode as Zn, solvent as DMSO, amine as EOEA
"""


# Create a concrete model
model = ConcreteModel()

# Define the sets
model.set_components    = pyo.Set(initialize=['electrolyte','Carb','amine','cation','anion','H2','CO'])
model.components_act    = pyo.Set(within=model.set_components, initialize=['electrolyte','amine'])
model.set_ions          = pyo.Set(within=model.set_components, initialize=['cation','anion'])
model.set_prodspecies   = pyo.Set(within=model.set_components, initialize=['H2','CO'])
model.set_nongas        = pyo.Set(within=model.set_components, initialize=['electrolyte','amine','cation','anion'])
model.set_nonions       = pyo.Set(within=model.set_components, initialize=['amine','H2','CO'])
model.set_electrolyte   = pyo.Set(within=model.set_components, initialize=['electrolyte'])

model.set_groups        = pyo.Set(initialize=['CH3','CH2NH2','COC','DMSO','cation','anion'])
model.set_CAT           = pyo.Set(initialize=['PIP','ELE','INS','CIV','BAL'])

# Define the sets
model.set_allcomponents = pyo.Set(initialize=['anolyte', 'catholyte','amine'])
model.set_rxcomponents  = pyo.Set(initialize=['Carb', 'amine', 'RNH3+', 'H+', 'OH-', 'H2O(l)', 'CO', 'H2','O2'])
model.set_reactions     = pyo.Set(initialize=['COER', 'HER'])
model.set_outcomponents = pyo.Set(initialize=['Carb', 'amine', 'CO', 'H2'])
model.set_allcomponents = pyo.Set(initialize=['Carb', 'amine', 'electrolyte','CO', 'H2'])

model.set_liqcomponents = pyo.Set(initialize=['electrolyte','Carb','amine'])

#EoEA: CCOCCN
#MEA:  C(CO)N
#H2O
#DMSO: 
#CS
#ClO4
#================the model is better and anything to start from=====================
#==========henry coefficent=========================================================
Henry = {}
Henry['electrolyte'] = 0
Henry['amine']       = 0
Henry['cation']      = 0
Henry['anion']       = 0


#==========this is for EoEA-DMSO-Water-Cation-Anion-H2-CO===========
#['CH3','CH2NH2','COC','DMSO','cation','anion','H2','CO'])
#['electrolyte','amine','cation','anion','H2','CO'])
yita = {}
yita['CH3',   'electrolyte'] = 1e-09
yita['CH2NH2','electrolyte'] = 1e-09
yita['COC',   'electrolyte'] = 1e-09
yita['DMSO',  'electrolyte'] = 1
yita['cation','electrolyte'] = 1e-09
yita['anion', 'electrolyte'] = 1e-09


yita['CH3',   'amine'] = 1
yita['CH2NH2','amine'] = 1
yita['COC',   'amine'] = 1
yita['DMSO',  'amine'] = 1e-09
yita['cation','amine'] = 1e-09
yita['anion', 'amine'] = 1e-09


yita['CH3',   'cation'] = 1e-09
yita['CH2NH2','cation'] = 1e-09
yita['COC',   'cation'] = 1e-09
yita['DMSO',  'cation'] = 1e-09
yita['cation','cation'] = 1
yita['anion', 'cation'] = 1e-09


yita['CH3',   'anion'] = 1e-09
yita['CH2NH2','anion'] = 1e-09
yita['COC',   'anion'] = 1e-09
yita['DMSO',  'anion'] = 1e-09
yita['cation','anion'] = 1e-09
yita['anion', 'anion'] = 1

#==========system parameters==================
PT = 100000        #par
T  = 298.15        #K
Psat = {}
Psat['electrolyte'] = 0.0008 * 100000
Psat['amine']       = 0.1 * 100000

VM = {}             #molar volume data is from chemspider
VM['electrolyte'] = 71.1/1000000    #m3/mol
VM['amine']       = 104.6/1000000   #m3/mol
VM['cation']      = 0
VM['anion']       = 0
#VM['H2']          = 0
#VM['CO']          = 0


#========================================================================

# Define decision variablesA
model.F_feed    = Var(within=NonNegativeReals, bounds=(0, 10e+09), initialize=0)
model.F_V       = Var(within=NonNegativeReals, bounds=(0, 10e+09), initialize=0)
model.F_L       = Var(within=NonNegativeReals, bounds=(0, 10e+09), initialize=0)

#=======================how to write this as components==================
model.cost_anolyte  = Var(within=NonNegativeReals, bounds=(0, 10e+07))


lb_rx = {'electrolyte': 1e-09, 'Carb': 1e-09, 'amine': 1e-09, 'cation': 1e-09, 'anion': 1e-09, 'CO': 1e-09, 'H2': 1e-09}
ub_rx = {'electrolyte': 1, 'Carb': 1, 'amine': 1, 'cation': 1, 'anion': 1, 'CO': 1, 'H2': 1}

lb_ry = {'electrolyte': 1e-09, 'Carb': 1e-09, 'amine': 1e-09, 'cation': 1e-09, 'anion': 1e-09, 'CO': 1e-09, 'H2': 1e-09}
ub_ry = {'electrolyte': 1, 'Carb': 1, 'amine': 1, 'cation': 1, 'anion': 1, 'CO': 1, 'H2': 1}

lb_ac = {'electrolyte': 1e-09,  'amine': 1e-09}
ub_ac = {'electrolyte': 10E+09, 'amine': 10E+09}


lb_vf = {'electrolyte': 0,   'amine': 0,   'cation': 0,   'anion': 0}
ub_vf = {'electrolyte': 100, 'amine': 100, 'cation': 100, 'anion': 100}


def fb_x(model, i):
    return (lb_rx[i], ub_rx[i])
    
def fb_y(model, i):
    return (lb_ry[i], ub_ry[i])

def fb_ac(model, i):
    return (lb_ac[i], ub_ac[i])
    
def fb_vfc(model, i):
    return (lb_vf[i], ub_vf[i])

#'electrolyte','amine','cation','anion','H2','CO'
feed_comp = {}
feed_comp['electrolyte'] =  input_parameter.loc[case,'electrolyte_comp']
feed_comp['amine']       =  input_parameter.loc[case,'amine_comp']
feed_comp['cation']      =  input_parameter.loc[case,'cation_comp']
feed_comp['anion']       =  input_parameter.loc[case,'anion_comp']
feed_comp['CO']          =  0
feed_comp['H2']          =  0

vfeed_comp = {}
vfeed_comp['electrolyte'] = 0
vfeed_comp['amine']       = 0
vfeed_comp['cation']      = 0
vfeed_comp['anion']       = 0
vfeed_comp['CO']          = input_parameter.loc[case,'CO_vcomp']
vfeed_comp['H2']          = input_parameter.loc[case,'H2_vcomp']

nnl={} #this is the lower bound
nnl['CH3',   'electrolyte'] = 1e-09
nnl['CH2NH2','electrolyte'] = 1e-09
nnl['COC',   'electrolyte'] = 1e-09
nnl['DMSO',  'electrolyte'] = 1e-09
nnl['cation','electrolyte'] = 1e-09
nnl['anion', 'electrolyte'] = 1e-09


nnl['CH3',   'amine'] = 1e-09
nnl['CH2NH2','amine'] = 1e-09
nnl['COC',   'amine'] = 1e-09
nnl['DMSO',  'amine'] = 1e-09
nnl['cation','amine'] = 1e-09
nnl['anion', 'amine'] = 1e-09


nnl['CH3',   'cation'] = 1e-09
nnl['CH2NH2','cation'] = 1e-09
nnl['COC',   'cation'] = 1e-09
nnl['DMSO',  'cation'] = 1e-09
nnl['cation','cation'] = 1e-09
nnl['anion', 'cation'] = 1e-09


nnl['CH3',   'anion'] = 1e-09
nnl['CH2NH2','anion'] = 1e-09
nnl['COC',   'anion'] = 1e-09
nnl['DMSO',  'anion'] = 1e-09
nnl['cation','anion'] = 1e-09
nnl['anion', 'anion'] = 1e-09



nnu={} #this is the upper bound
nnu['CH3',   'electrolyte'] = 10e+09
nnu['CH2NH2','electrolyte'] = 10e+09
nnu['COC',   'electrolyte'] = 10e+09
nnu['DMSO',  'electrolyte'] = 10e+09
nnu['cation','electrolyte'] = 10e+09
nnu['anion', 'electrolyte'] = 10e+09


nnu['CH3',   'amine'] = 10e+09
nnu['CH2NH2','amine'] = 10e+09
nnu['COC',   'amine'] = 10e+09
nnu['DMSO',  'amine'] = 10e+09
nnu['cation','amine'] = 10e+09
nnu['anion', 'amine'] = 10e+09


nnu['CH3',   'cation'] = 10e+09
nnu['CH2NH2','cation'] = 10e+09
nnu['COC',   'cation'] = 10e+09
nnu['DMSO',  'cation'] = 10e+09
nnu['cation','cation'] = 10e+09
nnu['anion', 'cation'] = 10e+09


nnu['CH3',   'anion'] = 10e+09
nnu['CH2NH2','anion'] = 10e+09
nnu['COC',   'anion'] = 10e+09
nnu['DMSO',  'anion'] = 10e+09
nnu['cation','anion'] = 10e+09
nnu['anion', 'anion'] = 10e+09

"""
nnu['CH3',   'H2'] = 10e+09
nnu['CH2NH2','H2'] = 10e+09
nnu['COC',   'H2'] = 10e+09
nnu['DMSO',  'H2'] = 10e+09
nnu['cation','H2'] = 10e+09
nnu['anion', 'H2'] = 10e+09


nnu['CH3',   'CO'] = 10e+09
nnu['CH2NH2','CO'] = 10e+09
nnu['COC',   'CO'] = 10e+09
nnu['DMSO',  'CO'] = 10e+09
nnu['cation','CO'] = 10e+09
nnu['anion', 'CO'] = 10e+09
"""

lb_gr ={}
ub_gr ={}
lb_gr['CH3']    = 1e-09
lb_gr['CH2NH2'] = 1e-09
lb_gr['COC']    = 1e-09
lb_gr['DMSO']   = 1e-09
lb_gr['cation'] = 1e-09
lb_gr['anion']  = 1e-09


ub_gr['CH3']    = 10e+09
ub_gr['CH2NH2'] = 10e+09
ub_gr['COC']    = 10e+09
ub_gr['DMSO']   = 10e+09
ub_gr['cation'] = 10e+09
ub_gr['anion']  = 10e+09



lbx_gr ={}
ubx_gr ={}
lbx_gr['CH3']    = 0
lbx_gr['CH2NH2'] = 0
lbx_gr['COC']    = 0
lbx_gr['DMSO']   = 0
lbx_gr['cation'] = 0
lbx_gr['anion']  = 0


ubx_gr['CH3']    = 1
ubx_gr['CH2NH2'] = 1
ubx_gr['COC']    = 1
ubx_gr['DMSO']   = 1
ubx_gr['cation'] = 1
ubx_gr['anion']  = 1



#===========interaction parametes come from http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
Inter = {}
Inter['CH3',   'CH3'] = 1e-09
Inter['CH2NH2','CH3'] = -30.480
Inter['COC',   'CH3'] = 21.49
Inter['DMSO',  'CH3'] = 50.490
Inter['cation','CH3'] = -701.34   #this interaction parameter
Inter['anion', 'CH3'] = 156.5 #1566.5


Inter['CH3',   'CH2NH2'] = 391.50
Inter['CH2NH2','CH2NH2'] = 1e-09
Inter['COC',   'CH2NH2'] = 21.49    #this interaction parameter need to be checked again, missing
Inter['DMSO',  'CH2NH2'] = 50.490
Inter['cation','CH2NH2'] = 300      #this interaction parameter need to be checked again, missing
Inter['anion', 'CH2NH2'] = -300     #this interaction parameter need to be checked again, missing


Inter['CH3',   'COC'] = 408.30
Inter['CH2NH2','COC'] = 1e-09        #this interaction parameter need to be checked again, missing
Inter['COC',   'COC'] = 1e-09        
Inter['DMSO',  'COC'] = -300     #this interaction parameter need to be checked again, missing
Inter['cation','COC'] = 300      #this interaction parameter need to be checked again, missing
Inter['anion', 'COC'] = -300     #this interaction parameter need to be checked again, missing


Inter['CH3',   'DMSO'] = 526.50
Inter['CH2NH2','DMSO'] = 874.19
Inter['COC',   'DMSO'] = 100      #this interaction parameter need to be checked again, missing
Inter['DMSO',  'DMSO'] = 1e-09
Inter['cation','DMSO'] = 300      #this interaction parameter need to be checked again, missing
Inter['anion', 'DMSO'] = -300     #this interaction parameter need to be checked again, missing


Inter['CH3',   'cation'] = 3398.2
Inter['CH2NH2','cation'] = 2000     #this interaction parameter need to be checked again, missing
Inter['COC',   'cation'] = -229.76  #-2297.6  #this parameter is based on CH3-CO        
Inter['DMSO',  'cation'] = 300      #this interaction parameter need to be checked again, missing
Inter['cation','cation'] = 1e-09      
Inter['anion', 'cation'] = 213.78   #this interaction parameter need to be checked again, missing


Inter['CH3',   'anion'] = -645.05   #-6450.5
Inter['CH2NH2','anion'] = 1e-09
Inter['COC',   'anion'] = 1e-09        #this interaction parameter need to be checked again, missing
Inter['DMSO',  'anion'] = 400      #this interaction parameter need to be checked again, missing
Inter['cation','anion'] = 781.78   
Inter['anion', 'anion'] = 1e-09    



Func_L = {}
Func_L['CH3',   'CH3'] = 1e-09
Func_L['CH2NH2','CH3'] = 1e-09
Func_L['COC',   'CH3'] = 1e-09
Func_L['DMSO',  'CH3'] = 1e-09
Func_L['cation','CH3'] = 1e-09   #this Func_Laction parameter
Func_L['anion', 'CH3'] = 1e-09


Func_L['CH3',   'CH2NH2'] = 1e-09
Func_L['CH2NH2','CH2NH2'] = 1e-09
Func_L['COC',   'CH2NH2'] = 1e-09    #this Func_Laction parameter need to be checked again, missing
Func_L['DMSO',  'CH2NH2'] = 1e-09
Func_L['cation','CH2NH2'] = 1e-09    #this Func_Laction parameter need to be checked again, missing
Func_L['anion', 'CH2NH2'] = 1e-09    #this Func_Laction parameter need to be checked again, missing


Func_L['CH3',   'COC'] = 1e-09
Func_L['CH2NH2','COC'] = 1e-09       #this Func_Laction parameter need to be checked again, missing
Func_L['COC',   'COC'] = 1e-09        
Func_L['DMSO',  'COC'] = 1e-09     #this Func_Laction parameter need to be checked again, missing
Func_L['cation','COC'] = 1e-09     #this Func_Laction parameter need to be checked again, missing
Func_L['anion', 'COC'] = 1e-09    #this Func_Laction parameter need to be checked again, missing



Func_L['CH3',   'DMSO'] = 1e-09
Func_L['CH2NH2','DMSO'] = 1e-09
Func_L['COC',   'DMSO'] = 1e-09      #this Func_Laction parameter need to be checked again, missing
Func_L['DMSO',  'DMSO'] = 1e-09
Func_L['cation','DMSO'] = 1e-09     #this Func_Laction parameter need to be checked again, missing
Func_L['anion', 'DMSO'] = 1e-09    #this Func_Laction parameter need to be checked again, missing



Func_L['CH3',   'cation'] = 1e-09
Func_L['CH2NH2','cation'] = 1e-09   #this Func_Laction parameter need to be checked again, missing
Func_L['COC',   'cation'] = 1e-09  #this parameter is based on CH3-CO        
Func_L['DMSO',  'cation'] = 1e-09      #this Func_Laction parameter need to be checked again, missing
Func_L['cation','cation'] = 1e-09      
Func_L['anion', 'cation'] = 1e-09   #this Func_Laction parameter need to be checked again, missing


Func_L['CH3',   'anion'] = 1e-09
Func_L['CH2NH2','anion'] = 1e-09
Func_L['COC',   'anion'] = 1e-09        #this Func_Laction parameter need to be checked again, missing
Func_L['DMSO',  'anion'] = 1e-09      #this Func_Laction parameter need to be checked again, missing
Func_L['cation','anion'] = 1e-09  
Func_L['anion', 'anion'] = 1e-09    



Func_U = {}
Func_U['CH3',   'CH3'] = 10e+08
Func_U['CH2NH2','CH3'] = 10e+08
Func_U['COC',   'CH3'] = 10e+08
Func_U['DMSO',  'CH3'] = 10e+08
Func_U['cation','CH3'] = 10e+08   #this Func_Uaction parameter
Func_U['anion', 'CH3'] = 10e+08


Func_U['CH3',   'CH2NH2'] = 10e+08
Func_U['CH2NH2','CH2NH2'] = 10e+08
Func_U['COC',   'CH2NH2'] = 10e+08    #this Func_Uaction parameter need to be checked again, missing
Func_U['DMSO',  'CH2NH2'] = 10e+08
Func_U['cation','CH2NH2'] = 10e+08    #this Func_Uaction parameter need to be checked again, missing
Func_U['anion', 'CH2NH2'] = 10e+08    #this Func_Uaction parameter need to be checked again, missing


Func_U['CH3',   'COC'] = 10e+08
Func_U['CH2NH2','COC'] = 10e+08        #this Func_Uaction parameter need to be checked again, missing
Func_U['COC',   'COC'] = 10e+08        
Func_U['DMSO',  'COC'] = 10e+08     #this Func_Uaction parameter need to be checked again, missing
Func_U['cation','COC'] = 10e+08     #this Func_Uaction parameter need to be checked again, missing
Func_U['anion', 'COC'] = 10e+08     #this Func_Uaction parameter need to be checked again, missing


Func_U['CH3',   'DMSO'] = 10e+08
Func_U['CH2NH2','DMSO'] = 10e+08
Func_U['COC',   'DMSO'] = 10e+08      #this Func_Uaction parameter need to be checked again, missing
Func_U['DMSO',  'DMSO'] = 10e+08
Func_U['cation','DMSO'] = 10e+08      #this Func_Uaction parameter need to be checked again, missing
Func_U['anion', 'DMSO'] = 10e+08     #this Func_Uaction parameter need to be checked again, missing


Func_U['CH3',   'cation'] = 10e+08
Func_U['CH2NH2','cation'] = 10e+08   #this Func_Uaction parameter need to be checked again, missing
Func_U['COC',   'cation'] = 10e+08  #this parameter is based on CH3-CO        
Func_U['DMSO',  'cation'] = 10e+08      #this Func_Uaction parameter need to be checked again, missing
Func_U['cation','cation'] = 10e+08      
Func_U['anion', 'cation'] = 10e+08   #this Func_Uaction parameter need to be checked again, missing


Func_U['CH3',   'anion'] = 10e+08
Func_U['CH2NH2','anion'] = 10e+08
Func_U['COC',   'anion'] = 10e+08        #this Func_Uaction parameter need to be checked again, missing
Func_U['DMSO',  'anion'] = 10e+08      #this Func_Uaction parameter need to be checked again, missing
Func_U['cation','anion'] = 10e+08   
Func_U['anion', 'anion'] = 10e+08     



def fb_vf(model, i, k):
    return (nnl[i,k], nnu[i,k])

def fb_gr(model, i):
    return (lb_gr[i], ub_gr[i])
    
def fb_gxy(model, i):
    return (lbx_gr[i], ubx_gr[i])

def func_lu(model, i, j):
    return (Func_L[i,j], Func_U[i,j])


#================the following parameters are used to determine group area fraction==============================
#=====http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
Q = {}
Q['CH3']   =  0.848
Q['CH2NH2']=  1.236
Q['COC']   =  0.240 
Q['DMSO']  =  2.472
Q['cation']=  0.9061
Q['anion'] =  2.0643

#================the following parameters are used to determine group volume fraction==============================
R = {}
R['CH3']   =  0.9011
R['CH2NH2']=  1.3692
R['COC']   =  0.6829
R['DMSO']  =  2.8266
R['cation']=  0.8611
R['anion'] =  2.2702

#================the molar mass====================================================================================
MM = {}
MM['electrolyte'] =  78.13
MM['amine']       =  89.084
MM['cation']      =  132.905
MM['anion']       =  100.46
#MM['H2']          =  2.016
#MM['CO']          =  28.01
#================the molar mass====================================================================================
Zc = {}
Zc['electrolyte'] =  0
Zc['amine']       =  0
Zc['cation']      =  1
Zc['anion']       =  1
#Zc['H2']          =  0
#Zc['CO']          =  0


#=========the following two equations determine parameters=====================================================
Xmp = {}
Qmp = {}
for i in model.set_groups:
      for j in model.set_nongas:
          Xmp[i,j] = yita[i,j] /sum(yita[k,j] for k in model.set_groups)


for i in model.set_groups:
      for j in model.set_nongas:
          Qmp[i,j] = Q[i] * Xmp[i,j] /sum(Q[k] * Xmp[k,j] for k in model.set_groups) + 1e-09
          


#==============parameter definition following model specification================================================
model.z_feed   =  pyo.Param(model.set_components, initialize=feed_comp, default=0)                   #known feed composition
PT             =  100000    #  with unit of par as pressure
model.Psat     =  pyo.Param(model.components_act, initialize=Psat, default=0)                        #known feed composition
#model.Xmp      =  pyo.Param(model.set_groups, model.set_components, initialize=yita, default=0)      #known feed composition
#model.Qmp      =  pyo.Param(model.set_groups, model.set_components, initialize=yita, default=0)

#model.Henry    =  pyo.Param(model.set_prodspecies, initialize=Henry, default=0)    
#'electrolyte','amine','cation','anion','H2','CO'
#================the following parameters are used to determine volume fraction==================================
#CO and H2 data is obtained from paper "Torli, M., Geer, L., Kontogeorgis, G.M. and Fosbøl, P.L., 2018. Industrial & Engineering Chemistry Research, 57(49), pp.16958-16977."
#ions are from paper "Macedo, E.A., Skovborg, P. and Rasmussen, P., 1990. Chemical Engineering Science, 45(4), pp.875-882."
#unconventional ions are from paper "Mohs, A. and Gmehling, J., 2013. Fluid Phase Equilibria, 337, pp.311-322."


model.X        =  Var(model.set_components,     within=NonNegativeReals, bounds=fb_x)      #liquid phase composition
model.Xm       =  Var(model.set_groups,     within=NonNegativeReals, bounds=fb_gxy)    #liquid phase composition
model.Y        =  Var(model.set_components,    within=NonNegativeReals, bounds=fb_y)      #vapor phase composition
model.AC       =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #activity coefficient
model.K        =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #activity coefficient
model.AC_SR    =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #short-range activity coefficient
model.AC_LR    =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #long-range activity coefficient
model.AC_C     =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #combinatorial activity coefficient
model.AC_R     =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #residual activity coefficient
model.phi      =  Var(model.set_groups, within=NonNegativeReals, bounds=fb_gr)          #related with functional groups
model.phip     =  Var(model.set_groups, model.set_nongas, within=NonNegativeReals, bounds=fb_vf)          #related with functional groups
 
model.van      =  Var(model.set_components, within=NonNegativeReals, bounds=(0, 10e+09))      #molecule volume parameters
model.qar      =  Var(model.set_components, within=NonNegativeReals, bounds=(0, 10e+09))      #molecule surface areas

model.A        =  Var(within=NonNegativeReals, bounds=(1e-09, 10e+09))      #Debye-Hückel parameter A
model.B        =  Var(within=NonNegativeReals, bounds=(1e-09, 10e+09))      #Debye-Hückel parameter B
model.volume   =  Var(model.set_components, within=NonNegativeReals, bounds=(0, 10e+09))      #Volume calculation
model.area     =  Var(model.set_components, within=NonNegativeReals, bounds=(0, 10e+09))      #Area calculation
model.surfacfrac = Var(model.set_groups, within=NonNegativeReals, bounds=(0, 10e+09))      #Area calculation

#Volume calculation
model.I        =  Var(within=NonNegativeReals, bounds=(1e-09, 10e+09))      #Ionic strength
model.Dmix     =  Var(within=NonNegativeReals, bounds=(0, 10e+09))      #Mixture density
model.Emix     =  Var(within=NonNegativeReals, bounds=(0, 10e+09))      #Mixture electric constant

model.vfrac    =  Var(model.components_act, within=NonNegativeReals, bounds=fb_ac)      #residual activity coefficient
model.feasible =  Var(within=NonNegativeReals, bounds=(0, 10e+09), initialize=0)
model.atermp   =  Var(model.set_groups, model.set_groups, within=NonNegativeReals, bounds=func_lu)    #related with functional groups
 



# Define the constraints
# why capital cost is much higher than operating cost?
# capital to operating cost is around 0.6/0.5 


#==================definition for the list of parameters==================
#==================define a csv file to keep all these parameters=========
#in order to do the sensitivity analysis, this section is to identify which parameters are case-dependent
#v is the charge density for reaction in EOEA, this should be a case-dependent input
#vv is the standard reversible reduction potential w.r.t SHE, this should be a case-dependent input with the symbol as ERev0, COER, HER, OER (OER is in athode)
#nnu is the stochimetric coefficient, reaction is always the same
#pp and MW are molecular weight, which are always constant


#AN_rou = 8.9        #g/cm3 as Ni, this is anode material density, which should be an input
#CN_rou = 7.3        #g/cm3 as Sn, this is cathod materail density, which should be an input
#anolyte_rou = 1.1   #g/cm3    anolyte density, DMSO, case dependent
#catholyte_rou = 1.1 #g/cm3    catholyte density, DMSO, case dependent

#MW_EoEA = 89        #molecular weight of EOEA
#cCap_RNH2   = 0.3   #M   amine concentration before co2 capture 0.3 M in experiment
#conver_RNH2 = 0.9   #conversion of amine, this equation requires clarification, for each cell, amine converted will be easier
#CO2_captured = 2637 #mol/L, co2 loading


#E_CatAg = -2        #this is a parameter to be supplied
#E_AnAg  = 2         #this is a parameter to be supplied
#jTot_Cat= 40.9      #this is a parameter to be supplied, total charge density
#E_Ag = 0.62         # V. SRP of Ag+/Ag w.r.t RHE, is this a constant?

#C_Catmat=226.47     #this is cost coefficient for cathode material
#C_Anmat=320.21488   #this is cost coefficient for anode material
#========================================================================

v = {}
v['COER']  = input_parameter.loc[case,'faradaic efficiency of COER']
v['HER']   = input_parameter.loc[case,'faradaic efficiency of HER']
vv = {}
vv['COER'] = input_parameter.loc[case,'standard reversible reduction potential of COER']
vv['HER']  = input_parameter.loc[case,'standard reversible reduction potential of HER']

z = {}
z['COER']  = 2
z['HER']   = 2 # number of electrons per reaction

nnu={} #why this coefficient is not the same as draft
nnu['COER','Carb']  = -1
nnu['COER','amine']  =  1
nnu['COER','RNH3+'] =  0
nnu['COER','H+']    =  -1
nnu['COER','OH-']   =  0
nnu['COER','H2O(l)']=  1
nnu['COER','CO']    =  1
nnu['COER','H2']    =  0
nnu['COER','O2']    =  0

nnu['HER','Carb']   =  0
nnu['HER','amine']   =  0
nnu['HER','RNH3+']  =  0
nnu['HER','H+']     =  -1
nnu['HER','OH-']    =  0
nnu['HER','H2O(l)'] =  0
nnu['HER','CO']     =  0
nnu['HER','H2']     =  1
nnu['HER','O2']     =  0

pp={}
pp['CO']=28.05
pp['H2']=2


#=====define these two reactions here======================
#Cathode reaction: 
#COER: RNHCOO-(Carb) + 2H+ + 2e- to RNH2 + CO + OH-
#HER:  2H+ + 2e- to H2 
#Anode reaction: 4OH- to O2 + 2H2O + 4e- 

F = 96485 # C /mol. Farraday's constant
R = 8.314 # J / mol K. Ideal gas constant
Troom= 298.15 # K
#=======Section for defining parameters====================
# define model parameters
Height  = 160        #cm
length  = 160        #cm

AN_area = Height * length #cm2
AN_rou  = input_parameter.loc[case,'anode density']          #g/cm3 as Ni

CN_area = Height * length #cm2
CN_rou  = input_parameter.loc[case,'cathode density']        #g/cm3 as Sn

Width_plastic = 0.6   #cm
Plastic_rou   = 2.2   #g/cm3 as plastic
Width_casing  = 0.5   #cm
d_ch          = 0.5   #channel depth, excluding plastic
rou_steel     = 7.96  #g/cm3



V_anolyte     = 15    #cm3
V_catholyte   = 15    #cm3
Ve_elyte      = 6.5   #cm/s
MW_EoEA       = input_parameter.loc[case,'amine molecular weight']
MW_ELYE       = input_parameter.loc[case,'electrolyte molecular weight']



V_Cat         = AN_area * d_ch
V_An          = V_Cat  # mL. Volume of compartments where electrolytes occupy, excluding plastic
VFlow         = Height * d_ch * Ve_elyte

A_V_Cat       = (AN_area * 1e-4) / (V_Cat * 1e-6) # 1/m
A_V_An        = (AN_area * 1e-4) / (V_An * 1e-6) # 1/m


tau           = V_Cat/VFlow # s. Residence time for electrochemical cell.

## Calculations
VFlow_Cat     = V_Cat / tau # mL / s
VFlow_An      = V_An / tau # mL / s 




#========process specification=============================
cCap_RNH2     = 0.3     #M   amine concentration before co2 capture
conver_RNH2   = 0.9     #conversion of amine, this equation requires clarification
CO2_captured  = 2637    #mol/L
AmineRequired = CO2_captured/conver_RNH2
solvent_price = 100

if (anolyte_rou < 1):
    AmineRequired = CO2_captured/conver_RNH2*2
    solvent_price = 0.53   #=======$/ton CO==================================
    
if (anolyte_rou > 1):   
    AmineRequired = CO2_captured/conver_RNH2   # number of cells required in parallel
    solvent_price = 2662   #=======$/ton CO==================================

"""
Ncell         = AmineRequired/(cCap_RNH2*VFlow/1000)
Ncell_stack   = 500   #number of cell per stack
Nstack        = Ncell/Ncell_stack
"""

model.Ncell         = Var(within=NonNegativeReals, bounds=(0, 10E+14), initialize=10E+14)
model.Nstack        = Var(within=NonNegativeReals, bounds=(0, 10E+14), initialize=10E+14)
#model.required_amine /(cCap_RNH2*VFlow/1000)

#model.Ncell.fix(100000)

Ncell_stack   = 500   #number of cell per stack
model.number_of_stack = Constraint(expr = model.Nstack == model.Ncell/Ncell_stack)


cincat = {}
cincat['Carb']   =  conver_RNH2*cCap_RNH2
cincat['amine']   =  (1 - conver_RNH2)*cCap_RNH2
cincat['RNH3+']  =  conver_RNH2*cCap_RNH2
cincat['H+']     =  conver_RNH2*cCap_RNH2
cincat['OH-']    =  0
cincat['H2O(l)'] =  0
cincat['CO']     =  0
cincat['H2']     =  0
cincat['O2']     =  0

#===================potential and charge specification========================
#===================faradaic efficiency of two reactions======================

model.FE_Raw  =  pyo.Param(model.set_reactions, initialize=v, default=0)
model.ERev0   =  pyo.Param(model.set_reactions, initialize=vv, default=0)
model.Z       =  pyo.Param(model.set_reactions, initialize=z, default=0)
model.nu      =  pyo.Param(model.set_reactions, model.set_rxcomponents, initialize=nnu, default=0)
model.cIn_Cat =  pyo.Param(model.set_rxcomponents, initialize=cincat, default=0)
model.MW      =  pyo.Param(model.set_prodspecies, initialize=pp, default=0)

Mass_AN       =  AN_area * AN_rou * 0.1
Mass_CN       =  CN_area * CN_rou * 0.1
Mass_plastic  =  CN_area * Plastic_rou * Width_plastic 

Mass_anolyte  =  anolyte_rou * V_An
Mass_catholyte=  catholyte_rou * V_Cat   #======g/mL * mL

print(catholyte_rou, V_Cat, Mass_catholyte)


Mass_Amine    =  MW_EoEA * cCap_RNH2 * 1e-3 * V_Cat
Flow_cell     =  Height * length * Ve_elyte


#============be careful about this value, cathode is always negative for the current experiment=================================
E_CatAg       =  input_parameter.loc[case,'relative potential between cathode and Ag (V)']    #this is a parameter to be supplied
E_AnAg        =  input_parameter.loc[case,'relative potential between anode and Ag (V)']     #this is a parameter to be supplied
jTot_Cat      =  input_parameter.loc[case,'total current density']  #this is a parameter to be supplied, total charge density
E_Ag          =  input_parameter.loc[case,'E-Ag']  # V. SRP of Ag+/Ag w.r.t RHE
vapor_pressure=  input_parameter.loc[case,'vapor pressure (bar)']

E_Cat         =  E_CatAg + E_Ag
E_An          =  E_AnAg + E_Ag

iTot          =  jTot_Cat * AN_area
jTot_An       =  iTot / AN_area
ERev0_OER     =  input_parameter.loc[case,'standard reversible reduction potential of OER']

DE            =  E_An - E_Cat # V
W_cell        =  DE * iTot * 1e-3     # electrical work with unit as W, 1 W = 1 J/s


model.total_W             = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.total_mass_AN       = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.total_mass_CN       = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.total_mass_plastic  = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.total_area_mem      = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.total_mass_casing   = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.total_mass_Amine    = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)

model.mass1  = Constraint(expr = model.total_W            == model.Ncell * W_cell)
model.mass2  = Constraint(expr = model.total_mass_AN      == model.Ncell * Mass_AN)
model.mass3  = Constraint(expr = model.total_mass_CN      == model.Ncell * Mass_CN)
model.mass4  = Constraint(expr = model.total_mass_plastic == Mass_plastic* model.Ncell) # g
model.mass5  = Constraint(expr = model.total_area_mem     == model.Ncell * AN_area)
model.mass6  = Constraint(expr = model.total_mass_casing  == 2 * AN_area * Width_casing * model.Nstack * rou_steel)
model.mass7  = Constraint(expr = model.total_mass_Amine   == model.Ncell * Mass_Amine)


"""
total_W       = Ncell*W_cell
total_mass_AN = Ncell*Mass_AN
total_mass_CN = Ncell*Mass_CN
total_mass_plastic = Mass_plastic * Ncell # g
total_area_mem     = Ncell*AN_area
total_mass_casing  = 2*AN_area*Width_casing*Nstack*rou_steel
total_mass_Amine   = Ncell*Mass_Amine
"""

#==========cost coefficent===================
C_Catmat   = input_parameter.loc[case,'cost of cathode material'] #====$/g======================
#C_Anmat    = 16500

C_Anmat    = input_parameter.loc[case,'cost of anode material']   #====$/g======================
C_Anolyte  = input_parameter.loc[case,'anolyte cost']     #====$/g======================
C_Catholyte= input_parameter.loc[case,'catholyte cost']   #====$/g======================
C_Amine    = input_parameter.loc[case,'amine_cost']       #====$/g======================
C_Pt       = 27
C_MEA      = 0.038
C_Ag       = 630.8e-3
C_Cu       = 7.4575e-3
C_Steel    = 1.5
C_Plastic  = 0.015  # $ / g. Cost of material
C_Elec     = 0.06 / (1000 * 3600) # $/J. Cost of electricity. $0.06/kWh
C_Mem      = 0.2777 # $/cm2. Cost of membrane for 183 um thickness


OL_CA      = 5
OL_AN      = 5
OL_plastic = 5
OL_casing  = 5
OL_membrane= 5

Plant_LT   = 20
Elect_LT   = lifetime_parameter.loc[life_case,'lifetime_value']


Finst      = 0.12
Fmain      = 0.06
AmineDegradationFraction = 0.05 # fraction that degrades over one year
#========================================================================


# Define decision variablesA
model.cap           = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.cop           = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.cost_AN       = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.cost_CA       = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.cost_Plastic  = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=0)
model.cost_Mem      = Var(within=NonNegativeReals, bounds=(0, 10e+20), initialize=461587692.307692)
model.cost_casing   = Var(within=NonNegativeReals, bounds=(0, 10e+20))
#=======================how to write this as components==================
model.cost_anolyte  = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_catholyte= Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_amine    = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_solvent  = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_material = Var(within=NonNegativeReals, bounds=(0, 10e+20))

model.cost_eleclyzer= Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_equip    = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_replace  = Var(within=NonNegativeReals, bounds=(0, 10e+20))

model.total_elec    = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_maintain = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_makeup   = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_recycle  = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_treat    = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_tvom     = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.cost_tfom     = Var(within=NonNegativeReals, bounds=(0, 10e+20))





model.eta_An        = Var(within=NonNegativeReals, bounds=(0, 10e+20))


lb_rx = {'COER': 0, 'HER': 0}
ub_rx = {'COER': 1000, 'HER': 1000}

lb_j = {'COER': 0, 'HER': 0}
ub_j = {'COER': 10E+07, 'HER': 10E+07}

lb_cat = {'COER': 0, 'HER': 0}
ub_cat = {'COER': 10E+07, 'HER': 10E+07}

lb_deeq = {'COER': 0, 'HER': 0}
ub_deeq = {'COER': 10E+07, 'HER': 10E+07}

lb_ve = {'COER': 0, 'HER': 0}
ub_ve = {'COER': 10E+07, 'HER': 10E+07}

lb_ee = {'COER': 0, 'HER': 0}
ub_ee = {'COER': 10E+07, 'HER': 10E+07}

lb_xi = {'COER': 0, 'HER': 0}
ub_xi = {'COER': 10E+07, 'HER': 10E+07}

lb_km = {}
ub_km = {}

lb_km['COER','Carb']  =  0
lb_km['COER','amine']  =  0
lb_km['COER','RNH3+'] =  0
lb_km['COER','H+']    =  0
lb_km['COER','OH-']   =  0
lb_km['COER','H2O(l)']=  0
lb_km['COER','CO']    =  0
lb_km['COER','H2']    =  0
lb_km['COER','O2']    =  0

lb_km['HER','Carb']   =  0
lb_km['HER','amine']   =  0
lb_km['HER','RNH3+']  =  0
lb_km['HER','H+']     =  0
lb_km['HER','OH-']    =  0
lb_km['HER','H2O(l)'] =  0
lb_km['HER','CO']     =  0
lb_km['HER','H2']     =  0
lb_km['HER','O2']     =  0

ub_km['COER','Carb']  =  10E+07
ub_km['COER','amine']  =  10E+07
ub_km['COER','RNH3+'] =  10E+07
ub_km['COER','H+']    =  10E+07
ub_km['COER','OH-']   =  10E+07
ub_km['COER','H2O(l)']=  10E+07
ub_km['COER','CO']    =  10E+07
ub_km['COER','H2']    =  10E+07
ub_km['COER','O2']    =  10E+07

ub_km['HER','Carb']    =  10E+07
ub_km['HER','amine']   =  10E+07
ub_km['HER','RNH3+']   =  10E+07
ub_km['HER','H+']      =  10E+07
ub_km['HER','OH-']     =  10E+07
ub_km['HER','H2O(l)']  =  10E+07
ub_km['HER','CO']      =  10E+07
ub_km['HER','H2']      =  10E+07
ub_km['HER','O2']      =  10E+07

lb_cout = {}
ub_cout = {}

lb_cout['Carb']   =  0
lb_cout['amine']  =  0
lb_cout['RNH3+']  =  0
lb_cout['H+']     =  0
lb_cout['OH-']    =  0
lb_cout['H2O(l)'] =  0
lb_cout['CO']     =  0
lb_cout['H2']     =  0
lb_cout['O2']     =  0

ub_cout['Carb']   =  10E+14
ub_cout['amine']  =  10E+14
ub_cout['RNH3+']  =  10E+14
ub_cout['H+']     =  10E+14
ub_cout['OH-']    =  10E+14
ub_cout['H2O(l)'] =  10E+14
ub_cout['CO']     =  10E+14
ub_cout['H2']     =  10E+14
ub_cout['O2']     =  10E+14


def fb(model, i): 
    return (lb_rx[i], ub_rx[i])

def fb_j(model, i):
    return (lb_j[i], ub_j[i])

def fb_cat(model, i):
    return (lb_cat[i], ub_cat[i])

def fb_deeq(model, i):
    return (lb_deeq[i], ub_deeq[i])

def fb_ve(model, i):
    return (lb_ve[i], ub_ve[i])

def fb_ee(model, i):
    return (lb_ee[i], ub_ee[i])

def fb_xi(model, i):
    return (lb_xi[i], ub_xi[i])

def fb_km(model, i, s):
    return (lb_km[i,s], ub_km[i,s])

def fb_cout(model, s):
    return (lb_cout[s], ub_cout[s])


ub_prod = {'electrolyte': 0,  'Carb': 0, 'amine': 0,   'cation': 0,   'anion': 0, 'CO': PT * (V_Cat - AN_area*0.1)/1000000/R/Troom*vfeed_comp['CO'], 'H2': PT * (V_Cat - AN_area*0.1)/1000000/R/Troom*vfeed_comp['H2']}

def prod_flow(model, i):
    return (ub_prod[i])

model.FE                = Var(model.set_reactions, within=NonNegativeReals, bounds=fb)
model.J                 = Var(model.set_reactions, within=NonNegativeReals, bounds=fb_j)
model.eta_CA            = Var(model.set_reactions, within=NonNegativeReals, bounds=fb_cat)
model.DE_Eq_OER         = Var(model.set_reactions, within=NonNegativeReals, bounds=fb_deeq)
model.VE                = Var(model.set_reactions, within=NonNegativeReals, bounds=fb_ve)
model.EE                = Var(model.set_reactions, within=NonNegativeReals, bounds=fb_ee)
model.Xi_Cat            = Var(model.set_reactions, within=NonNegativeReals, bounds=fb_xi)
model.km                = Var(model.set_reactions, model.set_rxcomponents, within=NonNegativeReals, bounds=fb_km)
model.cOut_Cat          = Var(model.set_rxcomponents, within=NonNegativeReals, bounds=fb_cout)
model.ProdFlow_Molar    = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.ProdFlow_Mass     = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.ProdFlow_Total    = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.product_MW        = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.XCO2              = Var(within=NonNegativeReals, bounds=(0, 1))
model.tol_use           = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.vapor_mass        = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.solvent_degrade   = Var(within=NonNegativeReals, bounds=(0, 10e+07))
model.vapor_mol         =  pyo.Param(model.set_components, initialize=prod_flow, default=0)


lb_cout = {}
ub_cout = {}

lb_cout['Carb']        =  0
lb_cout['amine']       =  0
lb_cout['electrolyte'] =  0
lb_cout['CO']          =  0
lb_cout['H2']          =  0

ub_cout['Carb']        =  10E+10
ub_cout['amine']       =  10E+10
ub_cout['electrolyte'] =  10E+10
ub_cout['CO']          =  10E+10
ub_cout['H2']          =  10E+10

def fb_cout(model, s):
    return (lb_cout[s], ub_cout[s])

model.Stream1           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)
model.Stream2           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)
model.Stream3           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)

model.Stream5           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)
model.Stream6           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)
model.Stream7           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)
model.Stream8           = Var(model.set_allcomponents,  within=NonNegativeReals, bounds=fb_cout)

#======================model this simulation problem as a feasibility problem==================================
model.con1 = Constraint(expr=model.feasible >= 0)


#=========this is the overall molar balance, inlet flowrate = sum of all outlet components and electrolyte = outlet vapor + outlet liquid====
#model.con1m = Constraint(expr= Mass_catholyte/MW_ELYE + Mass_Amine/MW_EoEA + sum(model.vapor_mol[i] for i in model.set_prodspecies) == model.F_V + model.F_L)
#=========the unit of odel.cIn_Cat and model.cOut_Cat are mol/L and VFlow_Cat is mL/s
model.con1m = Constraint(expr= sum(model.Stream2[i] for i in model.set_outcomponents) + model.Stream2['electrolyte'] == model.F_V + model.F_L)

#====the components of amine, CO and H2 can be vaporized and most of amine stay in liquid phase=========
model.component_bal = pyo.ConstraintList()
for i in model.set_nonions:                                                #model.set_nonions       = pyo.Set(within=model.set_components, initialize=['amine','H2','CO'])
         model.component_bal.add(model.Stream2[i] == model.F_V * model.Y[i] + model.F_L * model.X[i])

#====the components of electrolyte can be vaporized and most of it stay in liquid phase
model.electrolyte_bal = Constraint(expr = model.Stream2['electrolyte'] == model.F_V * model.Y['electrolyte'] + model.F_L * model.X['electrolyte'])


"""
#====the components of ions stay in liquid phase
#====ioinic amounts are based on total inlet flowrate and ions concentration===========
model.component_ion = pyo.ConstraintList()
for i in model.set_ions:                                                   #model.set_ions          = pyo.Set(within=model.set_components, initialize=['cation','anion'])
         model.component_ion.add((sum(model.cOut_Cat[i] * VFlow_Cat/1000 for i in model.set_outcomponents) + model.Stream2['electrolyte']) * feed_comp[i]  == model.F_L * model.X[i])
"""

#====equilibirum relationship for electrolyte and amine===============================         
model.equil_relate = pyo.ConstraintList()
for i in model.components_act:                                             #model.components_act   = pyo.Set(within=model.set_components, initialize=['electrolyte','amine'])  
         model.equil_relate.add(model.K[i] * model.X[i]  == model.Y[i])    #vapor-liquid phase equilibiurm is considered for nongas and nonions component

#====settting of compositions=========== 
model.component_sum_X = pyo.ConstraintList()
model.component_sum_X.add(sum(model.X[i] for i in model.set_nongas) == 1)  #model.set_nongas        = pyo.Set(within=model.set_components, initialize=['electrolyte','amine','cation','anion'])
model.X['CO'].fix(0)
model.X['H2'].fix(0)
model.Y['anion'].fix(0)
model.Y['cation'].fix(0)
model.Y['Carb'].fix(0)


model.component_sum_Y = pyo.ConstraintList()
model.component_Y_UP0 = pyo.ConstraintList()
model.component_Y_UP1 = pyo.ConstraintList()
model.component_Y_UP2 = pyo.ConstraintList()    
model.component_Y_UP3 = pyo.ConstraintList() 

CK_electrolyte   = input_parameter.loc[case,'K_electrolyte'] 
CK_amine         = input_parameter.loc[case,'K_amine']   

                            
model.component_sum_Y.add(sum(model.Y[i] for i in model.set_nonions) + model.Y['electrolyte'] == 1)     #the composition of nonions in vapor phase is summed to one


model.K['electrolyte'].fix(CK_electrolyte)
model.K['amine'].fix(CK_amine)



#=======this is to model a splitter based on the liquid outlet flowrate at the flash operation
#=======model.set_liqcomponents = pyo.Set(initialize=['electrolyte','Carb','amine'])===========
model.splitter = pyo.ConstraintList()
for i in model.set_liqcomponents:
      model.splitter  = Constraint(expr = model.F_L * model.X[i] == model.Stream5[i] + model.Stream6[i])
      
#====make-up electrolyte and amine is to compensate purge and vaporization=====================
model.makeup1  = Constraint(expr = model.Stream5['electrolyte'] + model.F_V * model.X['electrolyte'] == model.Stream7['electrolyte'])
model.makeup2  = Constraint(expr = model.Stream5['amine']       + model.F_V * model.Y['amine']       == model.Stream7['amine'])



#====amine degradation is modeled as purge and related with captured CO2=======================
#====electrolyte amount is related with amine concentration====================================
model.degrade1 = Constraint(expr = model.Stream5['electrolyte'] == model.Stream5['amine']/cCap_RNH2*anolyte_rou*1000/MW_ELYE)
model.degrade2 = Constraint(expr = model.Stream5['amine'] == 1.4 * (model.Stream1['Carb']-model.Stream2['Carb']) *(12+16*2)/1000000*1000/61.08)
model.degrade3 = Constraint(expr = model.Stream5['Carb'] == conver_RNH2/(1 - conver_RNH2) * model.Stream5['amine'])

model.balance1 = Constraint(expr = model.Stream2['Carb'] == model.Stream5['Carb'] + model.Stream6['Carb'])
model.balance3 = Constraint(expr = model.F_L * model.X['electrolyte'] == model.Stream5['electrolyte'] + model.Stream6['electrolyte'])

"""
=========DMSO as electrolyte==============
A = 5.23039	
B = 2239.161
C = -29.215

=========amine============================
#===338.6 to 444.1	4.29252	1408.873	-116.093	
A = 4.29252
B = 1408.873
C = -116.093
"""

"""
model.kamine       = Constraint(expr = model.K['amine'] == 10**(4.29252 - (1408.873 / (298.15 - 116.093))))
model.kelectrolyte = Constraint(expr = model.K['electrolyte'] == 10**(5.23039 - (2239.161 / (298.15 - 29.215))))


model.component_Y_UP0.add(model.K['amine'] >= 0.01)
model.component_Y_UP1.add(model.K['amine'] <= 0.99)     #the composition of nonions in vapor phase is summed to one
model.component_Y_UP2.add(model.K['electrolyte'] >= 1)     #the composition of nonions in vapor phase is summed to one
model.component_Y_UP3.add(model.K['electrolyte'] <= 100)     #the composition of nonions in vapor phase is summed to one
"""

"""
model.component_Y_UP1.add(model.Y['cation'] <= 1e-09)     #the composition of nonions in vapor phase is summed to one
model.component_Y_UP2.add(model.Y['anion'] <= 1e-09)     #the composition of nonions in vapor phase is summed to one
model.Y['cation'].fix(0)
model.Y['anion'].fix(0)
"""

#======================why capitla cost is so high==================================
model.con1  = Constraint(expr=model.cap - (1 + Finst)*(1/OL_AN*model.cost_AN + 1/OL_CA*model.cost_CA + 1/OL_plastic*model.cost_Plastic + 1/OL_membrane*model.cost_Mem + 1/OL_casing*model.cost_casing) == 0)


model.con1m = Constraint(expr=model.cop - (model.total_elec + model.cost_maintain + model.cost_makeup + model.cost_recycle + model.cost_treat) == 0)
model.con1a = Constraint(expr=model.cost_solvent   - (model.cost_anolyte + model.cost_catholyte + model.cost_amine) == 0)
model.con1b = Constraint(expr=model.cost_material  - (model.cost_CA + model.cost_AN + model.cost_Plastic + model.cost_Mem + model.cost_casing) == 0)

Learn_Rate = 0.13

model.elec_cost = Constraint(expr = model.cost_eleclyzer - Learn_Rate * model.cost_material == 0)
model.rpla_cost = Constraint(expr = model.cost_replace   - Learn_Rate * (model.cost_Mem + model.cost_CA) == 0)
model.equp_cost = Constraint(expr = model.cost_eleclyzer - model.cost_equip == 0)


model.con1c = Constraint(expr=model.total_elec      - C_Elec*model.total_W * 365 * 24 * 60 * 60 == 0)#this calculation convert the electricty cost unit to $/yr
model.con2  = Constraint(expr=model.cost_AN         - C_Anmat * model.total_area_mem/10000 == 0)
#/14.24*total_mass_AN/1000000*1000 == 0)
model.con3  = Constraint(expr=model.cost_CA         - C_Catmat * model.total_area_mem/10000 == 0)
#/11.52*total_mass_CN/1000000*1000 == 0)
model.con4  = Constraint(expr=model.cost_Plastic    - C_Plastic*model.total_mass_plastic == 0)
model.con5  = Constraint(expr=model.cost_Mem        - C_Mem*model.total_area_mem == 0)
model.con6  = Constraint(expr=model.cost_casing     - C_Steel*model.total_mass_casing == 0)


model.con7 = Constraint(expr=model.cost_anolyte    - C_Anolyte*Mass_anolyte*model.Ncell == 0)
model.con8 = Constraint(expr=model.cost_catholyte  - C_Catholyte*Mass_catholyte*model.Ncell == 0)


model.con9 = Constraint(expr = model.cost_amine      - C_Amine*model.total_mass_Amine == 0)

model.con10 = Constraint(expr = model.cost_maintain - (1/Plant_LT * model.cost_replace * (Plant_LT/Elect_LT - 1) + Fmain * model.cost_equip) == 0)

#model.con10 = Constraint(expr=model.cost_maintain  - Fmain*(1 + Finst)*model.cost_material == 0)   #the maintainence cost for each year.
#model.con11 = Constraint(expr=model.cost_makeup    - ((AmineDegradationFraction)*model.cost_amine + C_Amine*model.Ncell* MW_EoEA * model.F_V * model.Y['amine'] +  C_Anolyte * model.Ncell * MW_ELYE * model.F_V * model.Y['electrolyte']) == 0)

model.con11q = Constraint(expr=model.cost_makeup    - (978*model.Stream7['amine']/1000*61.08 + solvent_price * model.Stream7['electrolyte']/1000*MW_ELYE) * model.ProdFlow_Total  == 0)

#Alternative 1
#Life cycle assessment of post-combustion CO2 capture: A comparison between membrane separation and chemical absorption processes
#the above reference gives the pumping cost for co2 capture using MEA
#this unit cost is 33.8 kWh/t CO2 removed
#C_CO2_removed = 33.8/1000    #kWh/kg

#model.con12_recycle = Constraint(expr = model.cost_recycle   - (C_CO2_removed * C_Elec * (1000 * 3600) * (model.Stream1['Carb']  - model.Stream2['Carb']) * (61.08+12+16*2)/1000 * 365 * 24 * 60 * 60) == 0)  #this convert the cost to $/yr
#   kWh/kg * $/kWh * kg/s  = $/s

#Alternative 2
#====another source of determining the pumping cost is the following reference:
#====https://www.osti.gov/servlets/purl/1499969===============================
#====0.3 m3/s as velocity and 300 kPa ΔP as pressure difference===============
#====power (W) = 0.3 m3/s *  300,000 Pa / 0.8 = 112500 W = 112500 J/s
#====the recylce flowrate is 150,000 kmol/hr = 150,000/3600 kmol/s
#E_unit = 112500/(150,000/3600) J/kmol = 2700 J/kmol
E_unit = 2700

model.con12_recycle = Constraint(expr = model.cost_recycle   - (E_unit * C_Elec * (model.Stream6['Carb'] + model.Stream6['amine'] + model.Stream6['electrolyte'] + model.Stream7['amine'] + model.Stream7['electrolyte']) /1000 * 365 * 24 * 60 * 60) == 0)  #this convert the cost to $/yr


#=================the following cost parameters from NETL excel analysis sheet================
C_DMSO_treat = 0.145  #$0.145 per kg DMSO
C_Amine_treat= 0.145  #$0.145 per kg amine

model.con13_treat = Constraint(expr = model.cost_treat   - (C_DMSO_treat * model.Stream5['electrolyte']*MW_ELYE/1000 + C_Amine_treat * model.Stream5['amine']*61.8/1000) * 365 * 24 * 60 * 60 == 0)  #this convert the cost to $/yr

model.tvom_cal = Constraint(expr = model.cost_tvom - (model.total_elec + model.cost_makeup + model.cost_recycle + model.cost_treat) == 0)

f_LA = 0.06
f_SU = 0.15
f_SP = 0.15
f_TX = 0.02
f_IU = 0.01

model.C_LA     = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.C_SU     = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.C_SP     = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.C_TX     = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.C_IU     = Var(within=NonNegativeReals, bounds=(0, 10e+20))
model.BEC      = Var(within=NonNegativeReals, bounds=(0, 10e+20))


"""
model.c_la_cal = Constraint(expr = model.C_LA - f_LA * model.cost_equip == 0)
model.c_su_cal = Constraint(expr = model.C_SU - f_SU * model.C_LA == 0)
model.c_sp_cal = Constraint(expr = model.C_SP - f_SP * model.cost_maintain == 0)
model.c_tx_cal = Constraint(expr = model.C_TX - f_TX * model.epcc == 0)
model.c_iu_cal = Constraint(expr = model.C_IU - f_IU * model.epcc == 0)

model.tfom_cal = Constraint(expr = model.cost_tfom - (model.C_LA + model.C_SU + model.cost_maintain + model.C_SP + model.C_TX + model.C_IU) == 0)
"""

# Define the objective function
#model.obj = Objective(expr=model.feasible, sense=minimize)
model.obj = Objective(expr=(model.cap + model.cop)/model.ProdFlow_Total, sense=minimize)



#model.con11b = Constraint(expr = model.solvent_degrade == 2 * AmineDegradationFraction*Mass_anolyte*model.Ncell)
#model.con12= Constraint(expr=model.FE_COER == FE_Raw_COER / (FE_Raw_COER + FE_Raw_HER)) # normalizing in case FE inputs do not sum to 1
#model.con13= Constraint(expr=model.FE_HER == FE_Raw_HER / (FE_Raw_COER + FE_Raw_HER))



model.constraint_FE_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_FE_cal.add(model.FE[i] == model.FE_Raw) #[i]/sum(model.FE_Raw[ii] for ii in model.set_reactions))
model.constraint_J_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_J_cal.add(model.J[i]   == model.FE[i]*jTot_Cat)   #check this constraint in Tony's code
model.con12= Constraint(expr = model.eta_An      == E_An - ERev0_OER)


#model.con13= Constraint(expr = model.eta_An == E_An - ERev0_OER)
model.constraint_eta_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_eta_cal.add(model.eta_CA[i] == model.ERev0[i] - E_Cat)


model.constraint_deeq_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_deeq_cal.add(model.DE_Eq_OER[i] == ERev0_OER - model.ERev0[i])

model.constraint_ve_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_ve_cal.add(model.VE[i] == model.DE_Eq_OER[i]/DE)

model.constraint_ee_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_ee_cal.add(model.EE[i] == model.FE[i] * model.VE[i])


model.constraint_rate_cal = pyo.ConstraintList()
for i in model.set_reactions:
         model.constraint_rate_cal.add(model.Xi_Cat[i] == model.FE[i] * iTot / (model.Z[i] * F)) # mmol/s. extent of reaction



model.constraint_km_cal = pyo.ConstraintList()
n=2
mylist = model.set_rxcomponents
list_2 = itertools.islice(mylist, 2)
print(list_2)


"""
for i in model.set_reactions:
    for s in model.set_rxcomponents:
          if (model.nu[i,s] < 0):
              model.constraint_km_cal.add(model.km[i,s] == model.J[i] * (-model.nu[i,s]) / (model.Z[i] * F * model.cIn_Cat[s]) * 1e-2) # m/s
"""

"""
#====molar balance around the cell: Outlet = Inlet + Generation, mol/L as the unit
model.constraint_cout_cal = pyo.ConstraintList()
for s in model.set_rxcomponents:
    model.constraint_cout_cal.add(model.cOut_Cat[s] == model.cIn_Cat[s] + sum(model.nu[i,s] * model.Ncell * model.Xi_Cat[i] for i in model.set_reactions)/VFlow)
#=============this equation does not have relationship of DMSO==========================================================
"""
#===model.set_rxcomponents  = pyo.Set(initialize=['Carb', 'amine', 'RNH3+', 'H+', 'OH-', 'H2O(l)', 'CO', 'H2','O2'])


model.constraint_cout_cal = pyo.ConstraintList()
for s in model.set_outcomponents:   #===model.set_outcomponents = pyo.Set(initialize=['Carb', 'amine', 'CO', 'H2'])=====
    model.constraint_cout_cal.add(model.Stream2[s] == model.Stream1[s] + sum(model.nu[i,s] * model.Ncell * model.Xi_Cat[i]/1000 for i in model.set_reactions))



model.constraint_cell_cal1 = Constraint(expr = model.Stream2['Carb']   ==  model.Stream1['Carb']   -  model.Ncell * model.Xi_Cat['COER']/1000)
model.constraint_cell_cal2 = Constraint(expr = model.Stream2['amine']  ==  model.Stream1['amine']  +  model.Ncell/1000 * (1 * model.Xi_Cat['COER'] + 0 * model.Xi_Cat['HER'])) # - model.required_amine  * (1 + AmineDegradationFraction)) # + 2 * model.Xi_Cat['HER'])


model.constraint_cell_cal3 = Constraint(expr = model.Stream2['CO']  ==  0  +  model.Ncell * model.Xi_Cat['COER']/1000)
model.constraint_cell_cal4 = Constraint(expr = model.Stream2['H2']  ==  0  +  model.Ncell * model.Xi_Cat['HER']/1000)

"""
#158.5 ppm is the Carb concentration==================================================================
#production centration (mol/L) = 158.5 ppm / ((61.08+12+16*2) * 10^6) mol/L after carbon capture======
#model.stream1a  = Constraint(expr =  model.Stream1['electrolyte']/1000 == model.cIn_Cat['Carb']* VFlow_Cat/1000/(158.5 / ((61.08+12+16*2) * 10E+06))*anolyte_rou*1000/MW_ELYE)   #mol/s/(mol/L)*g/cm3/(g/mol)= L/s * g/cm3 / (g/mol) = 1000 g/s /(g/mol) = 1000 mol/s electrolyte molar flowrate
#model.stream1b  = Constraint(expr =  model.Stream1['electrolyte'] == model.cIn_Cat['amine'] * model.Ncell * VFlow/1000/cCap_RNH2*anolyte_rou*1000/MW_ELYE)
model.stream1c  = Constraint(expr =  model.Stream1['Carb']  == conver_RNH2/(1 - conver_RNH2) * model.Stream1['amine'])
#model.stream1d  = Constraint(expr =  model.Stream1['amine'] == model.Stream6['amine'] + model.Stream7['amine'])
"""

model.stream1b  = Constraint(expr =  model.Stream2['electrolyte'] == model.Stream2['amine']/cCap_RNH2*anolyte_rou*1000/MW_ELYE)
model.stream2a  = Constraint(expr =  model.Stream2['electrolyte'] == model.Stream1['electrolyte'])
#model.stream2c  = Constraint(expr =  model.Stream2['electrolyte'] == model.Stream8['electrolyte'])

#============relate stream 1 and stream 8============================
model.absorber1    = Constraint(expr = model.Stream1['Carb']   ==  model.Stream8['Carb']  + conver_RNH2 * model.Stream8['amine'])                           #====the unit for this equation is mol per second
model.absorber2    = Constraint(expr = model.Stream1['amine']  ==  model.Stream8['amine'] - conver_RNH2 * model.Stream8['amine'])     #====the unit for this equation is mol per second

model.stream8_define_1 = Constraint(expr = model.Stream8['Carb']       == model.Stream6['Carb'])
model.stream8_define_2 = Constraint(expr = model.Stream8['amine']      == model.Stream6['amine'] + model.Stream7['amine'])
model.stream8_define_2 = Constraint(expr = model.Stream8['electrolyte'] == model.Stream6['electrolyte'] + model.Stream7['electrolyte'])


#====determination of product molar flowrate: mmol/s = mol/L * mL/s
model.constraint_pflow_cal = pyo.ConstraintList()
model.constraint_pflow_cal.add(model.ProdFlow_Molar == model.Stream2['CO']) # mol / s

#====determination of product mass flowrate: mg/s = mol/L * mL/s * g/mol
model.constraint_pflow_mass = pyo.ConstraintList()
model.constraint_pflow_mass.add(model.ProdFlow_Mass == model.Stream2['CO'] * model.MW['CO']) # g / s
model.constraint_mass1      = Constraint(expr = model.ProdFlow_Mass == 1000)

#====determination of annual product mass flowrate
model.constraint_totalflow_mass = pyo.ConstraintList()
model.constraint_totalflow_mass.add(model.ProdFlow_Total == model.ProdFlow_Mass * 1e-6 * 365 * 24 * 60 * 60) # tonne / y


model.constraint_product_MW = pyo.ConstraintList()
model.constraint_product_MW.add(model.product_MW == model.ProdFlow_Mass/ model.ProdFlow_Molar) # g/mol


#=========determination of CO2 single-pass conversion=======
model.constraint_xco2 = pyo.ConstraintList()
model.constraint_xco2.add(model.XCO2 == 1 - model.Stream2['Carb'] / (model.Stream1['Carb']))  # []. conversion
#model.constraint_cout3_cal = Constraint(expr = model.Ncell * model.Xi_Cat['COER']/1000 == 0.749 * model.cIn_Cat['Carb'] *  model.Ncell * VFlow/1000)
model.conversion = Constraint(expr= model.XCO2 == 0.75)


#  object
solver = SolverFactory('scip')     # Change the solver name if you have a different solver installed
solver.options['max_iter']         = 10000
solver.options['limits/time']      = 50
solver.options['numerics/feastol'] = 1e-08
solver.options['limits/gap']       = 0.01

solver.solve(model, tee=True)
#log_infeasible_constraints(model, log_expression=True, log_variables=True)
#logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.INFO)


"""
df.loc[0, ['CAPEX']]          = model.cap.value/model.ProdFlow_Total.value
df.loc[0, ['OPEX']]           = model.cop.value/model.ProdFlow_Total.value
df.loc[0, ['cost_AN']]        = model.cost_AN.value/model.ProdFlow_Total.value * (1 + 0.12)/20
df.loc[0, ['cost_CA']]        = model.cost_CA.value/model.Ncell.value       #model.ProdFlow_Total.value * (1 + 0.12)/20
df.loc[0, ['cost_Plastic']]   = model.cost_Plastic.value/model.Ncell.value  #model.ProdFlow_Total.value * (1 + 0.12)/20
df.loc[0, ['cost_Mem']]       = model.cost_Mem.value/model.Ncell.value      #model.ProdFlow_Total.value * (1 + 0.12)/20
df.loc[0, ['cost_casing']]    = model.cost_casing.value/model.Ncell.value   #model.ProdFlow_Total.value * (1 + 0.12)/20

df.loc[0, ['total_elec']]     = model.total_elec.value/model.ProdFlow_Total.value
df.loc[0, ['cost_maintain']]  = model.cost_maintain.value/model.ProdFlow_Total.value
df.loc[0, ['cost_makeup']]    = model.cost_makeup.value/model.ProdFlow_Total.value
df.loc[0, ['cost_recycle']]   = model.cost_recycle.value/model.ProdFlow_Total.value
df.loc[0, ['cost_treat']]     = model.cost_treat.value/model.ProdFlow_Total.value

df.loc[0, ['Amine_makeup(kg/s)']]  = (model.Stream8['amine'].value)/1000*61.08 #  + (model.Ncell.value * model.F_V.value * model.Y['amine'].value)/1000*61.08

df.loc[0, ['prod_H2 (kg/s)']]      = model.Ncell.value*model.Xi_Cat['COER'].value/1000/1000 / model.FE['COER'].value*model.FE['HER'].value * 2
df.loc[0, ['prod_CO (kg/s)']]      = model.Ncell.value*model.Xi_Cat['COER'].value/1000/1000 * (12+16)

df.loc[0, ['Electricity(kJ/s)']]          = model.total_W.value /1000
df.loc[0, ['Reneration heat(kJ/hr)']]     = 0
df.loc[0, ['Compression work (kJ/hr)']]   = 0  #mol/s * KJ/mol * 3600 s/ hr
df.loc[0, ['prod_flow']]                  = model.ProdFlow_Total.value
df.loc[0, ['Number of cell']]             = model.Ncell.value
"""

df_netl = pd.DataFrame(columns=['Labor ($/yr)','Supervision ($/yr)','Maintenance ($/yr)', 'Supplies ($/yr)','Taxes ($/yr)','Insurance ($/yr)','BEC ($)','EPCC ($)','TPC ($)','TOC ($)','TFOM ($/yr)', 'TVOM ($/yr)'])

#=====provide values of A[ii] and beta[ii]=========================================================
Alpha = {}
beta  = {}

Alpha['PIP']=0.003596
Alpha['ELE']=0.0009276
Alpha['INS']=0.02446
Alpha['CIV']=0.03604
Alpha['BAL']=0.1111

beta['PIP'] =1.149
beta['ELE'] =1.173
beta['INS'] =0.9903
beta['CIV'] =1.067
beta['BAL'] =0.978

f_count = 0.09


BEC  = model.cost_equip.value + sum(Alpha[ii] * (model.cost_equip.value)**beta[ii] for ii in model.set_CAT)

EPCC = BEC*(1 + f_count) 

Cost_LA = f_LA * model.cost_equip.value

equip_value  = Fmain * model.cost_equip.value;
Cat_replace  = 1/Plant_LT * Learn_Rate * model.cost_CA.value  * (Plant_LT/Elect_LT - 1);
Mem_replace  = 1/Plant_LT * Learn_Rate * model.cost_Mem.value * (Plant_LT/Elect_LT - 1);

#=============postprocessing for determine NETL cost=======================================
df_netl.loc[0, ['Labor ($/yr)']]         = Cost_LA
df_netl.loc[0, ['Supervision ($/yr)']]   = f_SU * Cost_LA
df_netl.loc[0, ['Maintenance ($/yr)']]   = model.cost_maintain.value
df_netl.loc[0, ['Supplies ($/yr)']]      = f_SP * model.cost_maintain.value
df_netl.loc[0, ['Taxes ($/yr)']]         = f_TX * EPCC
df_netl.loc[0, ['Insurance ($/yr)']]     = f_IU * EPCC
df_netl.loc[0, ['BEC ($)']]              = BEC
df_netl.loc[0, ['EPCC ($)']]             = EPCC

f_proc      = 0.25
f_proj      = 0.2
f_preprod   = 0.02
f_inventory = 0.005
f_owner     = 0.15

f_TASC_TOC  = 1.07
f_FCF       = 0.06
f_CAP       = 0.907185

TPC = EPCC * (1 + f_proc + f_proj)
df_netl.loc[0, ['TPC ($)']]    = TPC

TOC = TPC * (1 + f_preprod + f_inventory + f_owner)

df_netl.loc[0, ['TOC ($)']]    = TPC * (1 + f_preprod + f_inventory + f_owner)


tfom = Cost_LA + f_SU * Cost_LA + model.cost_maintain.value + f_SP * model.cost_maintain.value + f_TX * EPCC + f_IU * EPCC
#model.cost_maintain - (1/Plant_LT * model.cost_replace * (Plant_LT/Elect_LT - 1) + Fmain * model.cost_equip

NETL_AN = (Learn_Rate * model.cost_AN.value * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_CA = (Learn_Rate * model.cost_CA.value * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_material = (Learn_Rate * (model.cost_Plastic.value + model.cost_casing.value) * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_mem      = (Learn_Rate * (model.cost_Mem.value) * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_BOP      = ((sum(Alpha[ii] * (model.cost_equip.value)**beta[ii] for ii in model.set_CAT)) * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

TOC_no_BOP = NETL_AN + NETL_CA  + Cat_replace/model.ProdFlow_Total.value + NETL_material + NETL_mem + Mem_replace/model.ProdFlow_Total.value
print(NETL_BOP, TOC_no_BOP, '=======================================bop====================================================')

df_netl.loc[0, ['Levelized TOC ($/Mt)']] = TOC_no_BOP #* f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

df_netl.loc[0, ['TFOM ($/yr)']] = tfom
df_netl.loc[0, ['TVOM ($/yr)']] = model.cost_tvom.value
df_netl.loc[0, ['Total equipment cost ($)']] = model.cost_equip.value

#df_netl.loc[0, ['H2 Profit ($/ton CO)']]     = R_H2/f_CAP
#df_netl.loc[0, ['Levelized TFOM ($/Mt)']]    = tfom/model.ProdFlow_Total.value
df_netl.loc[0, ['Levelized TVOM ($/Mt)']]    = model.cost_tvom.value/model.ProdFlow_Total.value

R_H2 = model.Ncell.value*model.Xi_Cat['COER'].value/1000/1000 / model.FE['COER'].value*model.FE['HER'].value * 2 *10 * 1000

print(R_H2, "Hydrogen cost")

df_netl.loc[0, ['H2 Profit ($/ton CO)']]     = R_H2/f_CAP
#df_netl.loc[0,['LCOP']] = (tfom + TOC * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value - R_H2/f_CAP
df_netl.loc[0,['LCOP with hydrogen revenue ($/Mt)']] = (model.cost_tvom.value)/model.ProdFlow_Total.value + TOC_no_BOP - R_H2/f_CAP
df_netl.loc[0,['LCOP without hydrogen revenue ($/Mt)']] = (model.cost_tvom.value)/model.ProdFlow_Total.value + TOC_no_BOP

df = pd.DataFrame(columns=['CAPEX','OPEX','cost_AN ($/ton CO)','cost_CA ($/ton CO)','cost_Material ($/ton CO)','cost_Mem ($/ton CO)','total_elec ($/ton CO)','cost_makeup + treat ($/ton CO)'])

NETL_AN = (Learn_Rate * model.cost_AN.value * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_CA = (Learn_Rate * model.cost_CA.value * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_material = (Learn_Rate * (model.cost_Plastic.value + model.cost_casing.value) * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_mem      = (Learn_Rate * (model.cost_Mem.value) * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

NETL_BOP      = ((sum(Alpha[ii] * (model.cost_equip.value)**beta[ii] for ii in model.set_CAT)) * (1 + f_count) * (1 + f_proc + f_proj) * (1 + f_preprod + f_inventory + f_owner) * f_TASC_TOC * f_FCF / f_CAP)/model.ProdFlow_Total.value

#df.loc[0, ['CAPEX']]          = model.cap.value/model.ProdFlow_Total.value
#df.loc[0, ['OPEX']]           = model.cop.value/model.ProdFlow_Total.value
df.loc[0, ['cost_AN ($/ton CO)']]        = NETL_AN
df.loc[0, ['cost_CA ($/ton CO)']]        = NETL_CA  + Cat_replace/model.ProdFlow_Total.value
df.loc[0, ['cost_Material ($/ton CO)']]  = NETL_material
df.loc[0, ['cost_Mem ($/ton CO)']]       = NETL_mem + Mem_replace/model.ProdFlow_Total.value
#df.loc[0, ['other BOP ($/ton CO)']]      = NETL_BOP
#model.ProdFlow_Total.value * (1 + 0.12)/20

#df = pd.DataFrame(columns=['CAPEX','OPEX','cost_AN','cost_CA','cost_Material','cost_Mem','other BOP','cost_maintain','other_maintain','total_elec','cost_makeup','cost_recylce','cost_treat','prod_flow'])

#df.loc[0, ['cost_maintain ($/ton CO)']]  = equip_value/model.ProdFlow_Total.value


other_maintain_cost = (tfom - model.cost_maintain.value)/model.ProdFlow_Total.value
#df.loc[0, ['other_maintain ($/ton CO)']] = other_maintain_cost

#df.loc[0, ['Other Cost ($/ton CO)']] = other_maintain_cost + NETL_BOP

df.loc[0, ['total_elec ($/ton CO)']]           = model.total_elec.value/model.ProdFlow_Total.value + model.cost_recycle.value/model.ProdFlow_Total.value #cell electricity and recycle electricity

df.loc[0, ['cost_makeup + treat ($/ton CO)']]  = model.cost_makeup.value/model.ProdFlow_Total.value + model.cost_treat.value/model.ProdFlow_Total.value
#df.loc[0, ['cost_recycle ($/ton CO)']]  = model.cost_recycle.value/model.ProdFlow_Total.value
#df.loc[0, ['cost_treat ($/ton CO)']]    = model.cost_treat.value/model.ProdFlow_Total.value

#df.loc[0, ['Other Cost ($/ton CO)']]  = other_maintain_cost + NETL_BOP  #other maintenance cost + balance of plant cost

"""
df.loc[0, ['BOP ($/ton CO)']]         = NETL_BOP
df.loc[0, ['Labor Cost ($/ton CO)']]  = Cost_LA/model.ProdFlow_Total.value
df.loc[0, ['Supervision ($/ton CO)']] = f_SU * Cost_LA/model.ProdFlow_Total.value
df.loc[0, ['Supplies ($/ton CO)']]    = f_SP * model.cost_maintain.value/model.ProdFlow_Total.value
df.loc[0, ['Taxes ($/ton CO)']]       = f_TX * EPCC/model.ProdFlow_Total.value
df.loc[0, ['Insurance ($/ton CO)']]   = f_IU * EPCC/model.ProdFlow_Total.value
"""



df.loc[0, ['Amine_makeup (kg/s)']]   = (model.Stream8['amine'].value)/1000*61.08 #  + (model.Ncell.value * model.F_V.value * model.Y['amine'].value)/1000*61.08

df.loc[0, ['prod_H2 (kg/s)']]        = model.Ncell.value*model.Xi_Cat['COER'].value/1000/1000 / model.FE['COER'].value*model.FE['HER'].value * 2
df.loc[0, ['prod_CO (kg/s)']]        = model.Ncell.value*model.Xi_Cat['COER'].value/1000/1000 * (12+16)

df.loc[0, ['Electricity (kJ/s)']]         = model.total_W.value /1000
df.loc[0, ['prod_flow (ton/yr)']]         = model.ProdFlow_Total.value
df.loc[0, ['H2 Profit ($/ton CO)']]       = R_H2/f_CAP
df.loc[0, ['Number of cell']]             = model.Ncell.value


df = df.transpose()
df_netl = df_netl.transpose()


#print(f"Optimal Solution: CAP = {model.cap.value}, OP = {model.cop.value}")
print(f"Optimal Solution: X_elect = {model.X['electrolyte'].value}, X_amine = {model.X['amine'].value}")
print(f"Optimal Solution: X_anion = {model.X['anion'].value}, X_canion = {model.X['cation'].value}")
print(f"Optimal Solution: Y_anion = {model.Y['anion'].value}, Y_canion = {model.Y['cation'].value}")
print(f"Optimal Solution: X_CO = {model.X['CO'].value}, X_H2 = {model.X['H2'].value}")
print(f"Optimal Solution: F_L = {model.F_L.value}, F_V = {model.F_V.value}, FY_CO = {model.F_V.value * model.Y['CO'].value}, FY_H2 = {model.F_V.value * model.Y['H2'].value}")
print(f"Optimal Solution: FY_elect = {model.F_V.value * model.Y['electrolyte'].value/tau}, FY_amine = {model.F_V.value * model.Y['amine'].value/tau}")
print(f"Optimal Solution: Y_elect = {model.Y['electrolyte'].value}, Y_amine = {model.Y['amine'].value}, Y_CO = {model.Y['CO'].value}, Y_H2 = {model.Y['H2'].value}")
print(f"Optimal Solution: vapor_CO = {model.vapor_mol['CO']}, vapor_H2 = {model.vapor_mol['H2']}")
print(f"Optimal Solution: k_elect = {model.K['electrolyte'].value}, k_amine = {model.K['amine'].value}")


# Print the optimal solution
print(f"Optimal Solution: CAP = {model.cap.value}, OP = {model.cop.value}")
print(f"Optimal Solution: X_COER = {model.Xi_Cat['COER'].value}, X_HER = {model.Xi_Cat['HER'].value}")
print(f"Optimal Solution: cost_AN = {model.cost_AN.value}, cost_CA = {model.cost_CA.value}, cost_Plastic= {model.cost_Plastic.value}, cost_Mem={model.cost_Mem.value}, cost_casing={model.cost_casing.value}")
print(f"Optimal Solution: total_elec = {model.total_elec.value}, cost_anolyte = {model.cost_anolyte.value}, cost_catholyte = {model.cost_catholyte.value}, cost_amine={model.cost_amine.value}, cost_maintain = {model.cost_maintain.value}, cost_makeup = {model.cost_makeup.value}")
#(OL_AN*model.cost_AN + OL_CA*model.cost_CA + OL_plastic*model.cost_Plastic + OL_membrane*model.cost_Mem + OL_casing*model.cost_casing)1

print(f"FE_COER = {model.FE['COER'].value}, FE_HER = {model.FE['HER'].value}, J_COER = {model.J['COER'].value}, J_HER = {model.J['HER'].value}")
print(f"etaAn = {model.eta_An.value}, eta_HER = {model.eta_CA['HER'].value}, eta_COER = {model.eta_CA['COER'].value}, DE_An={model.DE_Eq_OER['HER'].value}, DE_COER={model.DE_Eq_OER['COER'].value}")
print(f"VE_COER = {model.VE['COER'].value}, VE_HER = {model.VE['HER'].value}, EE_COER = {model.EE['COER'].value}, EE_HER = {model.EE['HER'].value}")
#print(f"Inlet = {model.cIn_Cat['Carb'].value}, Rate_HER = {model.Xi_Cat['HER'].value}")

print(f"X_CARB = {model.XCO2.value}")
print(f"Totoal_membrane_area = ", model.total_area_mem.value)
print(f"prod_mass = {model.ProdFlow_Total.value}")
print("number of cell", model.Ncell.value)
print("flowrate comparison", VFlow_Cat, VFlow)

print(model.F_L.value * model.X['electrolyte'].value, model.Stream5['electrolyte'].value, model.Stream6['electrolyte'].value)

print("CO flowrate",model.Stream2['CO'].value, model.Xi_Cat['COER'].value/1000)

df_2 = pd.DataFrame(columns=['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'])
#'H2O','Carb','amine','electrolyte','H2','CO',
#======the following table have stream flowrate unit of kg/s=========================
#======the unit of Cin_Cat[k] is mmol/L==============================================
#======the unit of VFlow is mL/s=====================================================

df_2.loc[0, 'S1'] = 0                                                                  #Water===========0 if the solvent is DMSO============
df_2.loc[1, 'S1'] = model.Stream1['Carb'].value  *(61.08+12+16*2)/1000     #carbomate
df_2.loc[2, 'S1'] = model.Stream1['amine'].value * 61.08/1000              #amine
df_2.loc[3, 'S1'] = model.Stream1['electrolyte'].value /1000 * MW_ELYE                                             #DMSO, mol/s/(mol/L)*g/L===========0 if the solvent is water============
df_2.loc[4, 'S1'] = 0
df_2.loc[5, 'S1'] = 0



df_2.loc[0, 'S2'] = 0                                                                 #Water===========0 if the solvent is DMSO============
df_2.loc[1, 'S2'] = model.Stream2['Carb'].value  /1000 * (61.08+12+16*2)              #carbomate
df_2.loc[2, 'S2'] = model.Stream2['amine'].value /1000 * 61.08                        #amine
df_2.loc[3, 'S2'] = model.Stream2['electrolyte'].value/1000*MW_ELYE                   #DMSO===========0 if the solvent is water============
df_2.loc[4, 'S2'] = model.Stream2['H2'].value /1000 * 2
df_2.loc[5, 'S2'] = model.Stream2['CO'].value /1000 * (12 + 16)


df_2.loc[0, 'S3'] = 0
df_2.loc[1, 'S3'] = model.Stream2['Carb'].value  /1000 * (61.08+12+16*2)
df_2.loc[2, 'S3'] = model.F_L.value * model.X['amine'].value * 61.08 /1000
df_2.loc[3, 'S3'] = model.F_L.value * model.X['electrolyte'].value * MW_ELYE /1000
df_2.loc[4, 'S3'] = 0 
df_2.loc[5, 'S3'] = 0



df_2.loc[0, 'S4'] = 0
df_2.loc[1, 'S4'] = 0
df_2.loc[2, 'S4'] = model.F_V.value * model.Y['amine'].value * 61.08/1000
df_2.loc[3, 'S4'] = model.F_V.value * model.Y['electrolyte'].value * MW_ELYE/1000
df_2.loc[4, 'S4'] = model.F_V.value * model.Y['H2'].value * 2/1000
df_2.loc[5, 'S4'] = model.F_V.value * model.Y['CO'].value * (12 + 16)/1000



df_2.loc[0, 'S5'] = 0
df_2.loc[1, 'S5'] = model.Stream5['Carb'].value  /1000 * (61.08+12+16*2)
df_2.loc[2, 'S5'] = model.Stream5['amine'].value * 61.08/1000
df_2.loc[3, 'S5'] = model.Stream5['electrolyte'].value * MW_ELYE/1000
df_2.loc[4, 'S5'] = 0 
df_2.loc[5, 'S5'] = 0


df_2.loc[0, 'S6'] = 0
df_2.loc[1, 'S6'] = model.Stream6['Carb'].value  /1000 * (61.08+12+16*2)
df_2.loc[2, 'S6'] = model.Stream6['amine'].value * 61.08/1000
df_2.loc[3, 'S6'] = model.Stream6['electrolyte'].value * MW_ELYE/1000
df_2.loc[4, 'S6'] = 0  
df_2.loc[5, 'S6'] = 0


df_2.loc[0, 'S7'] = 0
df_2.loc[1, 'S7'] = 0
df_2.loc[2, 'S7'] = model.Stream7['amine'].value * 61.08 /1000
df_2.loc[3, 'S7'] = model.Stream7['electrolyte'].value * MW_ELYE /1000
df_2.loc[4, 'S7'] = 0  
df_2.loc[5, 'S7'] = 0

df_2.loc[0, 'S8'] = 0
df_2.loc[1, 'S8'] = model.Stream6['Carb'].value  /1000 * (61.08+12+16*2)
df_2.loc[2, 'S8'] = (float(model.Stream6['amine'].value) * 61.08 + float(model.Stream7['amine'].value) * 61.08)/1000
df_2.loc[3, 'S8'] = (float(model.Stream6['electrolyte'].value * MW_ELYE) + float(model.Stream7['electrolyte'].value * MW_ELYE))/1000
df_2.loc[4, 'S8'] = 0  
df_2.loc[5, 'S8'] = 0

#cost times area or cost times mass
obj_val = model.obj()
print("Optimal Objective:", obj_val)
print("residential time", tau)
print(PT * (V_Cat - AN_area*0.1)/1000000/R/Troom, vfeed_comp['CO'])

print(model.Stream1['Carb'].value, model.Stream2['Carb'].value, model.Stream3['Carb'].value)# model.Stream4['Carb'].value)
print("===================================================================================")
print(model.Stream5['Carb'].value, model.Stream6['Carb'].value, model.Stream7['Carb'].value) 
print("CO flowrate", model.ProdFlow_Mass.value)
#print("recyle cost", model.cost_recycle.value, C_CO2_removed * C_Elec * 365 * 24 * 60 * 60 * (1000 * 3600) * (model.Stream1['Carb'].value  - model.Stream2['Carb'].value)*(61.08+12+16*2)/1000)
print("treat cost", model.cost_treat.value)
print("stream7 values", model.Stream7['amine'].value, model.Stream7['electrolyte'].value)
print(model.cost_makeup.value, model.ProdFlow_Total.value)



df.to_csv('./lifetime/cost_Analysis2_case_'+str(life_case)+'.csv')
df_2.to_csv('./lifetime/stream_Analysis2_table_'+str(life_case)+'.csv')
df_netl.to_csv('./lifetime/cost_netl_Analysis2_'+str(life_case)+'_main.csv')
