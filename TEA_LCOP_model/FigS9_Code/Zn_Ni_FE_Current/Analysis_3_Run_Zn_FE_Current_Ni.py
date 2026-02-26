from pyomo.environ import *
import pyomo.environ as pyo
import cyipopt
import itertools
import pandas as pd
from pyomo.util.infeasible import log_infeasible_constraints
import logging
import os
import multiprocessing
from joblib import Parallel, delayed

#========this is a flash model that represents the vapor-liquid equilibrum for electrolyte systems==========
#=================specify the cases to be used==============================================================
input_parameter = pd.read_csv('input_Current_FE_10K.csv') 
print(input_parameter)
#case = 6   # this will be the only parameters to enter for choosing parameters to run
#definition of cases


import pandas as pd
l = []
df_csv_append = pd.DataFrame()

os.system('mkdir Analysis3_Zn_FE_Current_Ni')

inputs = range(0, 10000)

def processInput(i):

    os.system('cp flash_analysis3_Zn_FE_Current_new_Ni.py flash_cal_analysis3_Zn_lifetime_case_'+str(i)+'.py')
    
    with open('flash_cal_analysis3_Zn_lifetime_case_'+str(i)+'.py', 'r') as file:
         filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('case_number_to_be_changed', ''+str(i)+'')

    # Write the file out again
    with open('flash_cal_analysis3_Zn_lifetime_case_'+str(i)+'.py', 'w') as file:
         file.write(filedata)
         
    os.system('python3 flash_cal_analysis3_Zn_lifetime_case_'+str(i)+'.py')
    os.system('rm flash_cal_analysis3_Zn_lifetime_case_'+str(i)+'.py')   
         
    input1 = [0]

    return input1
         
         
         
num_cores = 20 #multiprocessing.cpu_count()
print("numCores = " + str(num_cores))
results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)
print(results)
