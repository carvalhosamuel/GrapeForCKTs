# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 15:21:08 2021

@author: allan
"""

import grape
import algorithms
from functions import not_, and_, or_, nand_, nor_
import os
import pandas as pd
import numpy as np
from deap import creator, base, tools
import random

import time
import matplotlib.pyplot  as plt
import multiprocessing as mp
import string
from bitstring import Bits
from scipy.stats import pearsonr


RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# Delete any temp testbench
os.system("rm -f ./tmp/temp_tb*")

# Starting timer to measure execution time 
t = time.time()


# Truth table (only used within the testbench, but loaded here for compatibility)
X_train = np.zeros([3,8], dtype=bool)
Y_train = np.zeros([2,8], dtype=bool)


# Loading grammar file
GRAMMAR_FILE = 'FIR_4taps.bnf'

f = open("grammars/" + GRAMMAR_FILE, "r")
print(f.read())
f.close() 

BNF_GRAMMAR = grape.Grammar(os.path.join("grammars", GRAMMAR_FILE))

#print(BNF_GRAMMAR.production_rules)

# Function for generating random strings for temp testbench names

def id_generator(size=8, chars=string.ascii_uppercase + string.digits + string.ascii_lowercase):
    return ''.join(random.choice(chars) for _ in range(size))

def Icarus_fitness_eval_FIR(individual,points):
    
    # Target output obtained from "golden individual" (coefficients all set to 8'h20)
    target = [6656, 8864, 11072, 12832, 14592, 15424, 16256, 15424, 14592, 12608, 10624, 8224, 5824, 2880, -64, -2880, -5696, -8224, -10752, -12544, -14336, -14624, -14912, -13792, -12672, -10464, -8256, -5760]

    
    if individual.invalid == True:
        return -99,
    
    if isinstance(individual.phenotype, str):
        
        # Getting the contents of the testbench template
        with open('individual_tb_FIR.sv', 'r') as file:
            tb_template = file.read()
               
        # Creating the temporary testbench for this individual
        filename = 'temp_tb_' + id_generator() + '.sv'
        
        with open('./tmp/' + filename, "a") as testbench:
            
            # Phenotype
            testbench.write(individual.phenotype)
           
            # Writing the template
            testbench.write(tb_template)
            
        
        #TODO: CHECK EXECUTION TIME
        
        # Get fitness from call to simulator
        
        os_command = 'iverilog \'-g2012\' ./tmp/' + filename + ' -o ./tmp/a' + filename +'.out && vvp ./tmp/a' + filename +'.out' 
        
        # To run the golden individual:
        #os_command = 'iverilog -W all -g2012 ./FIR_tb.sv -o ./FIR_tb.out && vvp FIR_tb.out'
        
        output = os.popen(os_command).read().strip()
        
        # Split the string output into lines (as list)
        out_list = output.splitlines()
        
        int_array = []
        
        # 4 is the number of taps in the filter, 16 is our wordlength from testbench
        for line in out_list[4:]:
            # Get the LSB 16 bits  
            out = line[-16:]
            a = Bits(bin=out,length=16)
            int_array.append(a.int)
        
        # Check result:
        #np.binary_repr(-5760,width=16)
        
        # Seeing output 
        # plt.plot(int_array, color='red', label="FIR")
        # plt.plot(target, color='blue', label="Target")
        # plt.title('Output from FIR filter vs Target output')
        # plt.legend(loc="upper right")
        # plt.show()
        
        # Pearson's r correlation as fitness score
        fitness, p_value = pearsonr(int_array, target)
        
        # Cross-correlation as fitness score
        #fitness = np.correlate(np.array(int_array),np.array(target))
        #fitness = float(fitness)
        
        # RMSE as fitness score (needs change to FitnessMax or FitnessMin)
        #fitness = np.sqrt(np.mean((np.array(int_array)-np.array(target))**2))
        
        # Delete temp testbench
        os.system('rm -f ./tmp/' + filename)
        
        # Delete temp .out files
        os.system('rm -f ./tmp/a' + filename +'.out')
        
    
        return fitness,
    else:
        return -99,    
    
    
# Now the regular DEAP/GRAPE settings and initialisations:

toolbox = base.Toolbox()

# define a single objective, maximising fitness strategy:
creator.create("FitnessMax", base.Fitness, weights=(1.0,))

creator.create('Individual', grape.Individual, fitness=creator.FitnessMax)

toolbox.register("populationCreator", grape.sensible_initialisation, creator.Individual) 
#toolbox.register("populationCreator", grape.random_initialisation, creator.Individual) 
#toolbox.register("populationCreator", grape.PI_Grow, creator.Individual) 

toolbox.register("evaluate", Icarus_fitness_eval_FIR)

# Tournament selection:
toolbox.register("select", tools.selTournament, tournsize=6)
#toolbox.register("select", tools.selLexicase)

# Single-point crossover:
toolbox.register("mate", grape.crossover_onepoint)

# Flip-int mutation:
toolbox.register("mutate", grape.mutation_int_flip_per_codon)

POPULATION_SIZE = 50
MAX_GENERATIONS = 50
P_CROSSOVER = 0.9
P_MUTATION = 0.01
ELITE_SIZE = round(0.01*POPULATION_SIZE) #it should be smaller or equal to HALLOFFAME_SIZE
HALLOFFAME_SIZE = max(round(0.01*POPULATION_SIZE),1) #it should be at least 1

CODON_CONSUMPTION = 'lazy'
GENOME_REPRESENTATION = 'list'
RANDOM_SEED = 0 #Pay attention that the seed is set up inside the loop of runs, so you are going to have similar runs

INIT_GENOME_LENGTH = 50 #used only for random initialisation
random_initilisation = False #put True if you use random initialisation

MAX_INIT_TREE_DEPTH = 100
MIN_INIT_TREE_DEPTH = 10

MIN_INIT_GENOME_LENGTH = 20
MAX_INIT_GENOME_LENGTH = 100

MAX_TREE_DEPTH = 100
MAX_WRAPS = 0
CODON_SIZE = 255

REPORT_ITEMS = ['gen', 'invalid', 'avg', 'std', 'min', 'max', 
                'fitness_test',
                'best_ind_length', 'avg_length', 
                'best_ind_nodes', 'avg_nodes', 
                'best_ind_depth', 'avg_depth', 
                'avg_used_codons', 'best_ind_used_codons', 
                'structural_diversity', 'fitness_diversity',
                'selection_time', 'generation_time']

N_RUNS = 1

# Mutex for the prints - needed to avoid clashes when printing on the console in parallel runs
lock = mp.Lock()

# Function containing the evolutionary loop (defined as function so we can use it with the multiprocessing library)
def evoloop(run):
    
    # Run number
    r = run
    
    # Setting seed here for repeatability
    random.seed(r)
    
    # create initial population (generation 0):
    if random_initilisation:
        population = toolbox.populationCreator(pop_size=POPULATION_SIZE, 
                                           bnf_grammar=BNF_GRAMMAR,                            
                                           min_init_genome_length=MIN_INIT_GENOME_LENGTH,
                                           max_init_genome_length=MAX_INIT_GENOME_LENGTH,
                                           max_init_depth=MAX_INIT_TREE_DEPTH,
                                           codon_size=CODON_SIZE,
                                           codon_consumption=CODON_CONSUMPTION,
                                           genome_representation=GENOME_REPRESENTATION
                                           )
    else:
        population = toolbox.populationCreator(pop_size=POPULATION_SIZE, 
                                           bnf_grammar=BNF_GRAMMAR, 
                                           min_init_depth=MIN_INIT_TREE_DEPTH,
                                           max_init_depth=MAX_INIT_TREE_DEPTH,
                                           codon_size=CODON_SIZE,
                                           codon_consumption=CODON_CONSUMPTION,
                                           genome_representation=GENOME_REPRESENTATION
                                            )
    
     # define the hall-of-fame object:
    hof = tools.HallOfFame(HALLOFFAME_SIZE)
    
    # prepare the statistics object:
    stats = tools.Statistics(key=lambda ind: ind.fitness.values)
    stats.register("avg", np.nanmean)
    stats.register("std", np.nanstd)
    stats.register("min", np.nanmin)
    stats.register("max", np.nanmax)
    
    # perform the Grammatical Evolution flow:
    population, logbook = algorithms.ge_eaSimpleWithElitism(population, toolbox, cxpb=P_CROSSOVER, mutpb=P_MUTATION,
                                              ngen=MAX_GENERATIONS, elite_size=ELITE_SIZE,
                                              bnf_grammar=BNF_GRAMMAR, 
                                              codon_size=CODON_SIZE, 
                                              max_tree_depth=MAX_TREE_DEPTH,
                                              points_train=[X_train, Y_train], 
                                              codon_consumption=CODON_CONSUMPTION,
                                              report_items=REPORT_ITEMS,
                                              genome_representation=GENOME_REPRESENTATION,
                                              stats=stats, halloffame=hof, verbose=False, showplots=False)
    
    import textwrap
    best = hof.items[0].phenotype
    
    # Printing outcome from this run
    lock.acquire()
    time.sleep(0.1)
    print("\n\nResults of run", r+1, "of",N_RUNS)
    print("Best individual: \n","\n".join(textwrap.wrap(best,80)))
    print("\nBest Fitness: ", hof.items[0].fitness.values[0])
    print("Depth: ", hof.items[0].depth)
    print("Length of the genome: ", len(hof.items[0].genome))
    print(f'Used portion of the genome: {hof.items[0].used_codons/len(hof.items[0].genome):.2f}')
    
    lock.release()
    
    #valListFitness.append(val_fitness)
    #bestPhenotipes.append(best)
    
    return logbook, r+1 , hof.items[0]



### Running the evoloop in parallel

# Defines our pool of processors using all available ones
pool = mp.Pool(mp.cpu_count())

# Runs the evoloop function in parallel, passing only r as the argument 
results = pool.starmap_async(evoloop, [([r]) for r in range(0, N_RUNS)]).get()

pool.close()

# Waits for all queued processes to reach here before resuming
pool.join() 


### Getting the objects from the runs

logbooks, runs, hallof = zip(*results)

 
elapsed = time.time() - t
print('Elapsed Time: ', elapsed)
print('========')

### Statistics

# # Genetic Programming is done (all runs) - plot statistics:
x_axis = np.arange(0, MAX_GENERATIONS+1)
# avgArray = np.array(avgListFitness)
# stdArray = np.array(stdListFitness)
# minArray = np.array(minListFitness)
# maxArray = np.array(maxListFitness)

# testArray = np.array(testListFitness)

max_Fitness_List = []
avg_Fitness_List = []

for Logbook in logbooks:
    
    max_Fitness_List.append(Logbook.select("max"))
    avg_Fitness_List.append(Logbook.select("avg"))


max_Fitness_List = np.array(max_Fitness_List)
avg_Fitness_List = np.array(avg_Fitness_List)


plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.title('Average Max Fitness FIR Filter - Tournament')
plt.errorbar(x_axis, max_Fitness_List.mean(0), yerr=max_Fitness_List.std(0),label="Avg", color="deepskyblue", )
plt.show()

# plt.xlabel('Generation')
# plt.ylabel('Fitness')
# plt.title('Average Fitness for FIR Filter - Tournament')
# #plt.ylim(7,16)
# plt.errorbar(x_axis, avg_Fitness_List.mean(0), yerr=avg_Fitness_List.std(0),label="Avg", color="deepskyblue")
# plt.show()


# Standard deviation of miminum fitness
max_Fitness_List[:,MAX_GENERATIONS].std(0)

# Average minimum fitness
max_Fitness_List[:,MAX_GENERATIONS].mean(0)



