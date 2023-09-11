#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 12:23:55 2023

@author: samuel
"""
import os
from bitstring import Bits
import matplotlib.pyplot  as plt
import glob
import re
import numpy as np

class ApproxCircuit(object):
    """
    An approximate adder or multiplier. Class created to carry the error, area, delay, power, etc metrics.
    """

    def __init__(self, Name, MAEp, MAE, WCEp, WCE, WCREp, EPp, MREp, MSE, PWR, AREA, DELAY):
        """
        """        
        self.Name = Name
        self.MAEp = MAEp
        self.MAE = MAE
        self.WCEp = WCEp
        self.WCE = WCE
        self.WCREp = WCREp
        self.EPp = EPp
        self.MREp = MREp
        self.MSE = MSE
        self.PWR = PWR
        self.AREA = AREA
        self.DELAY = DELAY


def CreateApproxCircuits():
   
    files =  glob.glob("all_circuits/*.v")
    #print(files)
    files.sort()
    
    Circuits = []
    
    for file in files:
        # Empty list to get the values
        metrics = []
        
        # Getting the contents of the files 
        with open(file, 'r') as ckt:

            for i in range(16):                                                         # Going through the first 16 lines of each file 
                line = next(ckt).strip()
                if (line[0] == '/' and line[1] == '/'):                                 # If they match the pattern
                  metrics.append(float(re.findall("(?<==)\s*(\d*\.*\d+)?", line)[0]))   # Get the metrics using the Regex
        
        Name = re.findall("([^\/]+$)",file)[0]
        Circuit = ApproxCircuit(Name,metrics[0],metrics[1],metrics[2],metrics[3],metrics[4],metrics[5],metrics[6],metrics[7],metrics[8],metrics[9],metrics[10])          
        
        Circuits.append(Circuit)
                       
    return Circuits           



def testbench_builder(individual):
    
    #individual.phenotype = "n 8'h00 8'h01"
    
    
    adders = individual.adders
    multipliers = individual.multipliers
    
    #pheno = individual.phenotype.split()
    
    coefficients = individual.coefs_real
    
    #coefficients = pheno[1:]
    taps = len(coefficients)
    
    # # Checking and processing for symmetry in coefficients
    # if pheno[0] == 's':
        
    #     original = coefficients
    #     reverse  = list(reversed(coefficients))
    
    #     coefficients = original + reverse
        
    # Initialising all the dynamic parts of the testbench
    dynamic1 = ""
    
    dynamic2 = "parameter taps = "+str(taps)+";\n"
    
    
    dynamic3= "FIR_Filter inst0(clk, en, "
    
    dynamic5 = "input real "
    
    dynamic6 =  "wire real "
    
    dynamic6b = "wire [N-1:0] "
    
    dynamic7 =  "DFF DFF0(clk,data_in,x1);\n"
    
    dynamic8 = ""
    
    dynamic9 = "assign real_add_final =  "
    
    i = 0
    for coef in coefficients:
        dynamic1 = dynamic1 + "parameter real c_"+str(i)+" = " + str(coef) + ";\n"
        dynamic3 = dynamic3 + "c_"+str(i)+", "
        dynamic5 = dynamic5 + " b"+str(i)+", "
        dynamic6 = dynamic6 + " Mul"+str(i)+", "
        dynamic6b = dynamic6b + "x"+str(i)+", "
        
        if i<(len(coefficients)-2):
            dynamic7 = dynamic7 + "DFF DFF"+str(i+1)+"(clk,x"+str(i+1)+",x"+str(i+2)+");\n"
        
        dynamic8 = dynamic8 + "assign Mul"+str(i)+" = $bitstoreal(x"+str(i)+") * b"+str(i)+";\n"
        dynamic9 = dynamic9 + "Mul"+str(i)+" + "
        i += 1
        
    # Wrapping the dynamic contents
    dynamic3 = dynamic3 + "data_in, data_out);"
    
    dynamic4 = dynamic3.replace("FIR_Filter inst0", "module FIR_Filter")
    dynamic4 = dynamic4.replace("c_", "b")
    
    dynamic5 = dynamic5 + ";\n"
    dynamic5 = dynamic5.replace(", ;", ";")
    
    dynamic6 = dynamic6 + ";\n"
    dynamic6 = dynamic6.replace(", ;", ";")
    dynamic6b = dynamic6b + ";\n"
    dynamic6b = dynamic6b.replace(", ;", ";")
    
    dynamic8 = dynamic8.replace("x0","data_in")    
    
    dynamic9 = dynamic9 + ";\n"
    dynamic9 = dynamic9.replace(" + ;", ";")
    
    
    # Getting the contents of the testbench template
    with open('generic_tb_template.v', 'r') as file:
        tb_template = file.read()
    
    tb_template = tb_template.replace("//DYNAMIC1","//DYNAMIC1\n" + dynamic1)
    tb_template = tb_template.replace("//DYNAMIC2","//DYNAMIC2\n" + dynamic2)
    tb_template = tb_template.replace("//DYNAMIC3","//DYNAMIC3\n" + dynamic3)
    tb_template = tb_template.replace("//DYNAMIC4","//DYNAMIC4\n" + dynamic4)
    tb_template = tb_template.replace("//DYNAMIC5","//DYNAMIC5\n" + dynamic5)
    tb_template = tb_template.replace("//DYNAMIC6","//DYNAMIC6\n" + dynamic6 + dynamic6b)
    tb_template = tb_template.replace("//DYNAMIC7","//DYNAMIC7\n" + dynamic7)
    tb_template = tb_template.replace("//DYNAMIC8","//DYNAMIC8\n" + dynamic8)
    tb_template = tb_template.replace("//DYNAMIC9","//DYNAMIC9\n" + dynamic9)
    
    return tb_template


def aggregate_fitness(individual):
    
    # Gets the list of used circuits on that specific individual 
    used_ckts = individual.used_ckts
    
    # Initiating the lists for all the metrics
    MAEp = []
    MAE = []
    WCEp = []
    WCE = []
    WCREp = []
    EPp = []
    MREp = []
    MSE = []
    PWR = []
    AREA = []
    DELAY = []
    
    # Looping through the used circuits and getting their metrics from the LUT
    for each_ckt in used_ckts:
       
       for circuit in lookup_table:
            
            if each_ckt == circuit.Name:
                
                MAEp.append(circuit.MAEp)
                MAE.append(circuit.MAE)
                WCEp.append(circuit.WCEp)
                WCE.append(circuit.WCE)
                WCREp.append(circuit.WCREp)
                EPp.append(circuit.EPp)
                MREp.append(circuit.MREp)
                MSE.append(circuit.MSE)
                PWR.append(circuit.PWR)
                AREA.append(circuit.AREA)
                DELAY.append(circuit.DELAY)
                #print(circuit.Name)

    # Error metrics
    individual.MAEp  = np.mean(MAEp) #  Mean Absolute Error (percentage)
    individual.MAE   = np.mean(MAE)  #  Mean Absolute Error
    individual.WCEp  = np.sum(WCEp)  #  Worst-Case Absolute Error (percentage)
    individual.WCE   = np.sum(WCE)   #  Worst-Case Absolute Error
    individual.WCREp = np.sum(WCREp) #  Worst-Case Relative Error (percentage)
    individual.EPp   = np.sum(EPp)   #  Error Probability (percentage)
    individual.MREp  = np.mean(MREp) #  Mean Relative Error (percentage)
    individual.MSE   = np.mean(MSE)  #  Mean Squared Error
    # Design metrics
    individual.PWR   = np.sum(PWR)   #  Power
    individual.AREA  = np.sum(AREA)  #  Area
    individual.DELAY = np.sum(DELAY) #  Delay
    

def build_grammar(individual):
    
    # Function that will build a grammar according to the size of the individual coming from GE
    coefficients = individual.coefs_real
    taps = len(coefficients)
    
    multipliers = "<multipliers> ::= <mul> <mul>"
    adders      = "<adders> ::= <adder>"
    
    for i in range(taps-2):
        multipliers = multipliers + " <mul>"
        adders      = adders + " <adder>"
     
    terminals = "<adder> ::= 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 \n<mul> ::=  15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 \n"
    
    bnf_file = "<circuit> ::=  <adders> <multipliers> \n" + adders + "\n" + multipliers + "\n" + terminals
    
    return bnf_file


def scale_coefficients(individual):
    
    # Function to scale the coefficients from real numbers [-1,1] to integers [-127,127]
    
    coefficients = individual.coefs_real
    
    scaled = np.rint(coefficients*127).astype(int)
    
    vhex = np.vectorize(hex)

    hex_param = vhex(scaled).astype(str)
    
    output = []
    for number in hex_param: 
        output.append("8'h" + number[-2:])
        
    return output 
    
    
################ MAIN


lookup_table = CreateApproxCircuits()  


# To run the golden individual:
os_command = 'iverilog -W all -g2012 ./aprox_in_FIR.sv -o ./FIR_tb.out && vvp FIR_tb.out'

output = os.popen(os_command).read().strip()

# Setting the "x" at the output of filter to zeros 
output = output.replace("x", "0")

# Split the string output into lines (as list)
out_list = output.splitlines()

int_array = []

# 4 is the number of taps in the filter, 16 is our wordlength from testbench
for line in out_list:
    # Get the LSB 16 bits  
    out = line[-16:]
    a = Bits(bin=out,length=16)
    int_array.append(a.int)

# Check result:
#np.binary_repr(-5760,width=16)

# Seeing output 
plt.plot(int_array, color='red', label="FIR")
plt.title('Output from FIR filter')
plt.show()