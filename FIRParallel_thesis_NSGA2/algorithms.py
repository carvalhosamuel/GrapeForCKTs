#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.

import random
import math
import numpy as np
import time
import warnings
import os
from bitstring import Bits
import string
import matplotlib.pyplot  as plt
from scipy.stats import pearsonr
from datetime import datetime
from functools import partial
from operator import attrgetter
from scipy.fftpack import fft, ifft

from deap import tools
from deap.benchmarks.tools import diversity, convergence, hypervolume

def varAnd(population, toolbox, cxpb, mutpb,
           bnf_grammar, codon_size, max_tree_depth, codon_consumption,
           genome_representation, max_genome_length):
    """Part of an evolutionary algorithm applying only the variation part
    (crossover **and** mutation). The modified individuals have their
    fitness invalidated. The individuals are cloned so returned population is
    independent of the input population.

    :param population: A list of individuals to vary.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param cxpb: The probability of mating two individuals.
    :param mutpb: The probability of mutating an individual.
    :returns: A list of varied individuals that are independent of their
              parents.

    """
    offspring = [toolbox.clone(ind) for ind in population]

    # Apply crossover and mutation on the offspring
    for i in range(1, len(offspring), 2):
        if random.random() < cxpb:
            offspring[i - 1], offspring[i] = toolbox.mate(offspring[i - 1],
                                                          offspring[i],
                                                          bnf_grammar, 
                                                          max_tree_depth, 
                                                          codon_consumption,
                                                          genome_representation,
                                                          max_genome_length)
            del offspring[i - 1].fitness.values, offspring[i].fitness.values

    for i in range(len(offspring)):
        offspring[i], = toolbox.mutate(offspring[i], mutpb,
                                       codon_size, bnf_grammar, 
                                       max_tree_depth, codon_consumption,
                                       max_genome_length)
        del offspring[i].fitness.values

    return offspring

class hofWarning(UserWarning):
    pass

def testbench_builder(individual):
    
    #individual.phenotype = "n 8'h00 8'h01"
    
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

def plotResponse(phenotype):
    
    # Target output obtained from "golden individual" (coefficients all set to 8'h20)
    # target = [14592, 15424, 16256, 15424, 14592, 12608, 10624, 8224, 5824, 2880, -64, -2880, -5696, -8224, -10752, -12544, -14336, -14624, -14912, -13792, -12672, -10464, -8256, -5760]
    
    # Target being a clean sine wave
    #target = [0, 3256, 6383, 9255, 11758, 13793, 15277, 16152, 16384, 15962, 14904, 13252, 11072, 8450,  5491,  2313, -957,  -4189, -7253, -10029, -12405, -14286, -15598, -16288, -16328, -15718, -14481, -12666, -10347, -7615, -4580,  -1362]
    
    # Target 64 samples
    target = [0,3223,6320,9171,11662,13698,15199,16105,16383,16019,15030,13453,11351,8805,5914,2793,-438,-3651,-6722,-9530,-11966,-13934,-15357,-16180,-16371,-15922,-14850,-13199,-11031,-8432,-5504,-2360,875,4077,7119,9883,12260,14159,15504,16243,16347,15813,14660,12934,10703,8054,5090,1926,-1312,-4499,-7511,-10228,-12547,-14374,-15640,-16295,-16312,-15693,-14459,-12661,-10368,-7670,-4672,-1491]

   
    # Getting the contents of the testbench template
    with open('template_header', 'r') as file:
        tb_template_header = file.read()
   
    # Getting the contents of the testbench template
    with open('template_footer', 'r') as file:
        tb_template_footer = file.read()
               
     
    # Creating the temporary testbench for this individual
    
    #Workaround using timestamp instead of random number, so results with/without plot match for the same seed
    dt = datetime.now()
    ts = datetime.timestamp(dt)
    ts = str(ts)
    ts = ts.replace('.', '')
    filename = 'temp_tb_' + ts + '.v'
    
    
    with open('./tmp/' + filename, "a") as testbench:
        
        # Writing the template header
        testbench.write(tb_template_header)
        
        # Phenotype
        testbench.write(phenotype)
       
        # Writing the template footer
        testbench.write(tb_template_footer)
        
    
    #TODO: CHECK EXECUTION TIME
    
    # Get fitness from call to simulator
    
    os_command = 'iverilog \'-g2005\' ./tmp/' + filename + ' -o ./tmp/a' + filename +'.out && vvp ./tmp/a' + filename +'.out' 
    
    # To run the golden individual:
    #os_command = 'iverilog -W all -g2012 ./FIR_tb.sv -o ./FIR_tb.out && vvp FIR_tb.out'
    
    output = os.popen(os_command).read().strip()
    
    # Setting the "x" at the output of filter to zeros 
    output = output.replace("x", "0")
    
    # Split the string output into lines (as list)
    out_list = output.splitlines()
    
    int_array = []
    
    # 16 is our wordlength from testbench
    for line in out_list:
        # Get the LSB 16 bits  
        out = line[-16:]
        a = Bits(bin=out,length=16)
        int_array.append(a.int)
    
    # Check result:
    #np.binary_repr(-5760,width=16)
    
    # Correlation as fitness score
    correlation, p_value = pearsonr(int_array, target)
    
    RMSE = np.sqrt(np.mean((np.array(int_array)-np.array(target))**2))
    
    # Seeing output 
    plt.plot(int_array, color='red', label="FIR")
    plt.plot(target, color='blue', label="Target")
    #plt.title('Best individual. Correlation:' +str(correlation))
    plt.title('Best: Correlation:' +str(correlation) + ' RMSE: ' + str(RMSE))
    plt.legend(loc="upper right")
    plt.show()
    
    # Delete temp testbench
    os.system('rm -f ./tmp/' + filename)
    
    # Delete temp .out files
    os.system('rm -f ./tmp/a' + filename +'.out')
    
    

def plotResponse_generic(individual,input_data,target):
  
    # Ran once to generate values below
    fft_target = fft(target)
    fft_target = np.abs(fft_target)[0:int(len(target)/2)]
    
         
     
    # Generating timestamp to save plots
    dt = datetime.now()
    ts = datetime.timestamp(dt)
    ts = str(ts)
    ts = ts.replace('.', '')
    
    # Time domain plot
    plt.figure()
    plt.plot(individual.time_response, color='red', label="FIR")
    plt.plot(target, color='blue', label="Target")
    plt.plot(input_data, color='green', label="Input")
    #plt.title('Best individual. Correlation:' +str(correlation))
    plt.title('Python Execution\n FFT Correlation:' +str(individual.fitness) + ' RMSE: ' + str(individual.RMSE))
    plt.legend(loc="upper right")
    plt.xlabel('Phenotype: ' +individual.phenotype)
    plt.savefig('./plots/python_timedomain'+ts+'.png')
    plt.show()
    
    # Frequency domain plot  
    plt.figure(figsize = (12, 6))
    plt.title('FFTs: Target vs Output')
    plt.subplot(121)    
    plt.stem(np.arange(len(fft_target)), fft_target, 'b', markerfmt=" ", basefmt="-b")
    plt.xlabel('Freq (Hz)')
    plt.ylabel('FFT Amplitude |X(freq)|')
    plt.suptitle('Python Execution FFTs: Target vs Output \n Phenotype: ' +individual.phenotype, fontsize=16)


    plt.subplot(122)
    plt.stem(np.arange(len(fft_target)), individual.freq_response, 'b', markerfmt=" ", basefmt="-b")
    plt.xlabel('Freq (Hz)')
    plt.ylabel('FFT Amplitude |X(freq)|')
    plt.savefig('./plots/python_freqdomain'+ts+'.png')
    plt.show()
    
def twos(val_str, bytes):
    import sys
    val = int(val_str, 2)
    b = val.to_bytes(bytes, byteorder=sys.byteorder, signed=False)                                                          
    return int.from_bytes(b, byteorder=sys.byteorder, signed=True)
    

def plotResponse_verilog(individual, input_data, target):
    
    # Getting the contents of the testbench template
    DUT = testbench_builder(individual)
                
    # Creating the temporary testbench for this individual
    
    #Workaround using timestamp instead of random number, so results with/without plot match for the same seed
    dt = datetime.now()
    ts = datetime.timestamp(dt)
    ts = str(ts)
    ts = ts.replace('.', '')
    filename = 'temp_tb_' + ts + '.v'
    
    
    with open('./tmp/' + filename, "a") as testbench:
        
     # Writing the template header
     testbench.write(DUT)
        
    
    #TODO: CHECK EXECUTION TIME
    
    # Get fitness from call to simulator
    
    os_command = 'iverilog \'-g2005\' ./tmp/' + filename + ' -o ./tmp/a' + filename +'.out && vvp ./tmp/a' + filename +'.out' 
    
    # To run the golden individual:
    #os_command = 'iverilog -W all -g2012 ./FIR_tb.sv -o ./FIR_tb.out && vvp FIR_tb.out'
    
    output = os.popen(os_command).read().strip()
    
    # Setting the "x" at the output of filter to zeros 
    output = output.replace("x", "0")
    
    # Split the string output into lines (as list)
    out_list = output.splitlines()
    
    #int_array = []
    double_array = []
  
    # 16 is our wordlength from testbench
    for line in out_list:
      # Get the LSB 16 bits  
      out = line
      #out = line[-16:]
      a = Bits(bin=out,length=64)
      #int_array.append(a.int)
      double_array.append(a.f)
    
    # output = []
    # for val in double_array:
    #     output.append(val/max(double_array))
    
    output = double_array
    
    #%varexp --plot int_array
    
    
    # Correlation as fitness score
    correlation, p_value = pearsonr(output, target)
    
    RMSE = np.sqrt(np.mean((np.array(output)-np.array(target))**2))
     
    # Seeing output 
    plt.figure()
    plt.plot(output, color='red', label="FIR")
    plt.plot(target, color='blue', label="Target")
    plt.plot(input_data, color='green', label="Input")
    #plt.title('Best individual. Correlation:' +str(correlation))
    plt.title('Verilog Execution\n Signal Correlation:' +str(correlation) + ' RMSE: ' + str(RMSE))
    plt.legend(loc="upper right")
    plt.xlabel('Phenotype: ' +individual.phenotype)
    plt.savefig('./plots/verilog_timedomain'+ts+'.png')
    plt.show()
          
    
    # FFT of target
    sr = 128
    X = fft(target)
    N = len(target)
    n = np.arange(N)
    T = N/sr
    freq = n/T
    
    # FFT of output    
    X2 = fft(output)    
    
    # sampling interval
    st = 1.0/sr
    t = np.arange(0,1,st)
      
    plt.figure(figsize = (12, 6))
    
    plt.subplot(121)    
    plt.stem(freq[0:int(len(freq)/2)], (2*np.abs(X)/N)[0:int(len(freq)/2)], 'b', markerfmt=" ", basefmt="-b")
    plt.xlabel('Freq (Hz)')
    plt.ylabel('FFT Amplitude |X(freq)|')
    plt.suptitle('Verilog Execution FFTs: Target vs Output \n Phenotype: ' +individual.phenotype, fontsize=16)
    #plt.xlim(0, 10)

    plt.subplot(122)
    plt.stem(freq[0:int(len(freq)/2)], (2*np.abs(X2)/N)[0:int(len(freq)/2)], 'b', markerfmt=" ", basefmt="-b")
    plt.xlabel('Freq (Hz)')
    plt.ylabel('FFT Amplitude |X(freq)|')
    plt.savefig('./plots/verilog_freqdomain'+ts+'.png')
    plt.show()
    
    # Delete temp testbench
    os.system('rm -f ./tmp/' + filename)
    
    # Delete temp .out files
    os.system('rm -f ./tmp/a' + filename +'.out')    



def ge_eaSimpleWithElitism(population, toolbox, cxpb, mutpb, ngen, elite_size, 
                bnf_grammar, codon_size, max_tree_depth, 
                max_genome_length=None,
                points_train=None, points_test=None, codon_consumption='eager', 
                report_items=None,
                genome_representation='list',
                stats=None, halloffame=None, 
                verbose=__debug__, showplots=False):
    """This algorithm reproduce the simplest evolutionary algorithm as
    presented in chapter 7 of [Back2000]_, with some adaptations to run GE
    on GRAPE.
    :param population: A list of individuals.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param cxpb: The probability of mating two individuals.
    :param mutpb: The probability of mutating an individual.
    :param ngen: The number of generation.
    :param elite_size: The number of best individuals to be copied to the 
                    next generation.
    :params bnf_grammar, codon_size, max_tree_depth: Parameters 
                    used to mapper the individuals after crossover and
                    mutation in order to check if they are valid.
    :param stats: A :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: Whether or not to log the statistics.
    :returns: The final population
    :returns: A class:`~deap.tools.Logbook` with the statistics of the
              evolution
    """
    
    logbook = tools.Logbook()
    
    if halloffame is None:
        if elite_size != 0:
            raise ValueError("You should add a hof object to use elitism.") 
        else:
            warnings.warn('You will not register results of the best individual while not using a hof object.', hofWarning)
            logbook.header = ['gen', 'invalid'] + (stats.fields if stats else []) + ['avg_length', 'avg_nodes', 'avg_depth', 'avg_used_codons', 'behavioural_diversity', 'structural_diversity', 'fitness_diversity', 'selection_time', 'generation_time']
    else:
        if halloffame.maxsize < 1:
            raise ValueError("HALLOFFAME_SIZE should be greater or equal to 1")
        if elite_size > halloffame.maxsize:
            raise ValueError("HALLOFFAME_SIZE should be greater or equal to ELITE_SIZE")         
        if points_test:
            logbook.header = ['gen', 'invalid'] + (stats.fields if stats else []) + ['fitness_test', 'best_ind_length', 'avg_length', 'best_ind_nodes', 'avg_nodes', 'best_ind_depth', 'avg_depth', 'avg_used_codons', 'best_ind_used_codons', 'behavioural_diversity', 'structural_diversity', 'fitness_diversity', 'selection_time', 'generation_time']
        else:
            logbook.header = ['gen', 'invalid'] + (stats.fields if stats else []) + ['best_ind_length', 'avg_length', 'best_ind_nodes', 'avg_nodes', 'best_ind_depth', 'avg_depth', 'avg_used_codons', 'best_ind_used_codons', 'behavioural_diversity', 'structural_diversity', 'fitness_diversity', 'selection_time', 'generation_time']

    start_gen = time.time()        
    # Evaluate the individuals with an invalid fitness
    for ind in population:
        if not ind.fitness.valid:
            ind.fitness.values = toolbox.evaluate(ind, points_train)
        
    valid0 = [ind for ind in population if not ind.invalid]
    valid = [ind for ind in valid0 if not math.isnan(ind.fitness.values[0])]
    if len(valid0) != len(valid):
        warnings.warn("Warning: There are valid individuals with fitness = NaN in the population. We will avoid them.")
    invalid = len(population) - len(valid0) #We use the original number of invalids in this case, because we just want to count the completely mapped individuals    
    
    list_structures = []
    if 'fitness_diversity' in report_items:
        list_fitnesses = []
    if 'behavioural_diversity' in report_items:
        behaviours = np.zeros([len(valid), len(valid[0].fitness_each_sample)], dtype=float)
    
    #for ind in offspring:
    for idx, ind in enumerate(valid):
        list_structures.append(str(ind.structure))
        if 'fitness_diversity' in report_items:
            list_fitnesses.append(str(ind.fitness.values[0]))
        if 'behavioural_diversity' in report_items:
            behaviours[idx, :] = ind.fitness_each_sample
            
    unique_structures = np.unique(list_structures, return_counts=False)  
    if 'fitness_diversity' in report_items:
        unique_fitnesses = np.unique(list_fitnesses, return_counts=False) 
    if 'behavioural_diversity' in report_items:
        unique_behaviours = np.unique(behaviours, axis=0)
    
    structural_diversity = len(unique_structures)/len(population)
    fitness_diversity = len(unique_fitnesses)/(len(points_train[1])+1) if 'fitness_diversity' in report_items else 0 #TODO generalise for other problems, because it only works if the fitness is proportional to the number of testcases correctly predicted
    behavioural_diversity = len(unique_behaviours)/len(population) if 'behavioural_diversity' in report_items else 0

    # Update the hall of fame with the generated individuals
    if halloffame is not None:
        halloffame.update(valid)
        best_ind_length = len(halloffame.items[0].genome) 
        best_ind_nodes = halloffame.items[0].nodes
        best_ind_depth = halloffame.items[0].depth
        best_ind_used_codons = halloffame.items[0].used_codons
        if not verbose:
            print("gen =", 0, ", Best fitness =", halloffame.items[0].fitness.values)
            
    
    length = [len(ind.genome) for ind in valid]
    avg_length = sum(length)/len(length)
    
    nodes = [ind.nodes for ind in valid]
    avg_nodes = sum(nodes)/len(nodes)
    
    depth = [ind.depth for ind in valid]
    avg_depth = sum(depth)/len(depth)
    
    used_codons = [ind.used_codons for ind in valid]
    avg_used_codons = sum(used_codons)/len(used_codons)
    
    end_gen = time.time()
    generation_time = end_gen-start_gen
        
    selection_time = 0
    
    if points_test:
        fitness_test = np.NaN
    
    record = stats.compile(population) if stats else {}
    if points_test: 
        logbook.record(gen=0, invalid=invalid, **record, 
                       fitness_test=fitness_test,
                       best_ind_length=best_ind_length, avg_length=avg_length, 
                       best_ind_nodes=best_ind_nodes,
                       avg_nodes=avg_nodes,
                       best_ind_depth=best_ind_depth,
                       avg_depth=avg_depth,
                       avg_used_codons=avg_used_codons,
                       best_ind_used_codons=best_ind_used_codons,
                       behavioural_diversity=behavioural_diversity,
                       structural_diversity=structural_diversity,
                       fitness_diversity=fitness_diversity,
                       selection_time=selection_time, 
                       generation_time=generation_time)
    else:
        logbook.record(gen=0, invalid=invalid, **record, 
                       best_ind_length=best_ind_length, avg_length=avg_length, 
                       best_ind_nodes=best_ind_nodes,
                       avg_nodes=avg_nodes,
                       best_ind_depth=best_ind_depth,
                       avg_depth=avg_depth,
                       avg_used_codons=avg_used_codons,
                       best_ind_used_codons=best_ind_used_codons,
                       behavioural_diversity=behavioural_diversity,
                       structural_diversity=structural_diversity,
                       fitness_diversity=fitness_diversity,
                       selection_time=selection_time, 
                       generation_time=generation_time)
    if verbose:
        print(logbook.stream)

    # Begin the generational process
    for gen in range(logbook.select("gen")[-1]+1, ngen + 1):
        start_gen = time.time()    
    
        # Select the next generation individuals
        start = time.time()    
        offspring = toolbox.select(population, len(population)-elite_size)
        end = time.time()
        selection_time = end-start
        # Vary the pool of individuals
        offspring = varAnd(offspring, toolbox, cxpb, mutpb,
                           bnf_grammar, codon_size, max_tree_depth, 
                           codon_consumption, genome_representation,
                           max_genome_length)

        # Evaluate the individuals with an invalid fitness
        for ind in offspring:
            if not ind.fitness.valid:
                ind.fitness.values = toolbox.evaluate(ind, points_train)
                
        #Update population for next generation
        population[:] = offspring
        #Include in the population the elitist individuals
        for i in range(elite_size):
            population.append(halloffame.items[i])
            
        valid0 = [ind for ind in population if not ind.invalid]
        valid = [ind for ind in valid0 if not math.isnan(ind.fitness.values[0])]
        if len(valid0) != len(valid):
            warnings.warn("Warning: There are valid individuals with fitness = NaN in the population. We will avoid in the statistics.")
        invalid = len(population) - len(valid0) #We use the original number of invalids in this case, because we just want to count the completely mapped individuals
        
        list_structures = []
        if 'fitness_diversity' in report_items:
            list_fitnesses = []
        if 'behavioural_diversity' in report_items:
            behaviours = np.zeros([len(valid), len(valid[0].fitness_each_sample)], dtype=float)
        
        for idx, ind in enumerate(valid):
            list_structures.append(str(ind.structure))
            if 'fitness_diversity' in report_items:
                list_fitnesses.append(str(ind.fitness.values[0]))
            if 'behavioural_diversity' in report_items:
                behaviours[idx, :] = ind.fitness_each_sample
                
        unique_structures = np.unique(list_structures, return_counts=False)  
        if 'fitness_diversity' in report_items:
            unique_fitnesses = np.unique(list_fitnesses, return_counts=False) 
        if 'behavioural_diversity' in report_items:
            unique_behaviours = np.unique(behaviours, axis=0)
        
        structural_diversity = len(unique_structures)/len(population)
        fitness_diversity = len(unique_fitnesses)/(len(points_train[1])+1) if 'fitness_diversity' in report_items else 0 #TODO generalise for other problems, because it only works if the fitness is proportional to the number of testcases correctly predicted
        behavioural_diversity = len(unique_behaviours)/len(population) if 'behavioural_diversity' in report_items else 0
        
        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(valid)
            best_ind_length = len(halloffame.items[0].genome)
            best_ind_nodes = halloffame.items[0].nodes
            best_ind_depth = halloffame.items[0].depth
            best_ind_used_codons = halloffame.items[0].used_codons
            if not verbose:
                #print("gen =", gen, ", Best fitness =", halloffame.items[0].fitness.values, ", Number of invalids =", invalid, " Phenotype = ",halloffame.items[0].phenotype)
                print("gen =", gen, ", Number of invalids =", invalid)
            if showplots:    
                plotResponse(halloffame.items[0].phenotype)
            if points_test:
                # if gen < ngen:
                #     fitness_test = np.NaN
                # else:
                #     #fitness_test = toolbox.evaluate(halloffame.items[0], points_test)[0]
                fitness_test = halloffame.items[0].fitness_test
            
        length = [len(ind.genome) for ind in valid]
        avg_length = sum(length)/len(length)
        
        nodes = [ind.nodes for ind in valid]
        avg_nodes = sum(nodes)/len(nodes)
        
        depth = [ind.depth for ind in valid]
        avg_depth = sum(depth)/len(depth)
        
        used_codons = [ind.used_codons for ind in valid]
        avg_used_codons = sum(used_codons)/len(used_codons)
        
        end_gen = time.time()
        generation_time = end_gen-start_gen
        
        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        if points_test: 
            logbook.record(gen=gen, invalid=invalid, **record, 
                       fitness_test=fitness_test,
                       best_ind_length=best_ind_length, avg_length=avg_length, 
                       best_ind_nodes=best_ind_nodes,
                       avg_nodes=avg_nodes,
                       best_ind_depth=best_ind_depth,
                       avg_depth=avg_depth,
                       avg_used_codons=avg_used_codons,
                       best_ind_used_codons=best_ind_used_codons,
                       behavioural_diversity=behavioural_diversity,
                       structural_diversity=structural_diversity,
                       fitness_diversity=fitness_diversity,
                       selection_time=selection_time, 
                       generation_time=generation_time)
        else:
            logbook.record(gen=gen, invalid=invalid, **record, 
                       best_ind_length=best_ind_length, avg_length=avg_length, 
                       best_ind_nodes=best_ind_nodes,
                       avg_nodes=avg_nodes,
                       best_ind_depth=best_ind_depth,
                       avg_depth=avg_depth,
                       avg_used_codons=avg_used_codons,
                       best_ind_used_codons=best_ind_used_codons,
                       behavioural_diversity=behavioural_diversity,
                       structural_diversity=structural_diversity,
                       fitness_diversity=fitness_diversity,
                       selection_time=selection_time, 
                       generation_time=generation_time)
                
        if verbose:
            print(logbook.stream)
            
    #print("Final population hypervolume FROM ALGORITHMS is %f" % hypervolume(population,[0,0]))
    
 
   
    return population, logbook

def selDoubleTournament(individuals, k, fitness_size, parsimony_size, fitness_first, fit_attr="fitness"):
    """Tournament selection which use the size of the individuals in order
    to discriminate good solutions. This kind of tournament is obviously
    useless with fixed-length representation, but has been shown to
    significantly reduce excessive growth of individuals, especially in GP,
    where it can be used as a bloat control technique (see
    [Luke2002fighting]_). This selection operator implements the double
    tournament technique presented in this paper.
    The core principle is to use a normal tournament selection, but using a
    special sample function to select aspirants, which is another tournament
    based on the size of the individuals. To ensure that the selection
    pressure is not too high, the size of the size tournament (the number
    of candidates evaluated) can be a real number between 1 and 2. In this
    case, the smaller individual among two will be selected with a probability
    *size_tourn_size*/2. For instance, if *size_tourn_size* is set to 1.4,
    then the smaller individual will have a 0.7 probability to be selected.
    .. note::
        In GP, it has been shown that this operator produces better results
        when it is combined with some kind of a depth limit.
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :param fitness_size: The number of individuals participating in each \
    fitness tournament
    :param parsimony_size: The number of individuals participating in each \
    size tournament. This value has to be a real number\
    in the range [1,2], see above for details.
    :param fitness_first: Set this to True if the first tournament done should \
    be the fitness one (i.e. the fitness tournament producing aspirants for \
    the size tournament). Setting it to False will behaves as the opposite \
    (size tournament feeding fitness tournaments with candidates). It has been \
    shown that this parameter does not have a significant effect in most cases\
    (see [Luke2002fighting]_).
    :param fit_attr: The attribute of individuals to use as selection criterion
    :returns: A list of selected individuals.
    .. [Luke2002fighting] Luke and Panait, 2002, Fighting bloat with
        nonparametric parsimony pressure
    """
    assert (1 <= parsimony_size <= 2), "Parsimony tournament size has to be in the range [1, 2]."

    def _sizeTournament(individuals, k, select):
        chosen = []
        for i in range(k):
            # Select two individuals from the population
            # The first individual has to be the shortest
            prob = parsimony_size / 2.
            ind1, ind2 = select(individuals, k=2)

            if ind1.RMSE > ind2.RMSE:
                ind1, ind2 = ind2, ind1
            elif ind1.RMSE == ind2.RMSE:
                # random selection in case of a tie
                prob = 0.5

            # Since size1 <= size2 then ind1 is selected
            # with a probability prob
            chosen.append(ind1 if random.random() < prob else ind2)

        return chosen

    def _fitTournament(individuals, k, select):
        chosen = []
        for i in range(k):
            aspirants = select(individuals, k=fitness_size)
            chosen.append(max(aspirants, key=attrgetter(fit_attr)))
        return chosen

    if fitness_first:
        tfit = partial(_fitTournament, select=selRandom)
        return _sizeTournament(individuals, k, tfit)
    else:
        tsize = partial(_sizeTournament, select=selRandom)
        return _fitTournament(individuals, k, tsize)
    
def selRandom(individuals, k):
    """Select *k* individuals at random from the input *individuals* with
    replacement. The list returned contains references to the input
    *individuals*.
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :returns: A list of selected individuals.
    This function uses the :func:`~random.choice` function from the
    python base :mod:`random` module.
    """
    return [random.choice(individuals) for i in range(k)]    