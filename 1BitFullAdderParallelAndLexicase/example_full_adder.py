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

RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# Delete any temp testbench
os.system("rm -f temp_tb.sv")

# Starting timer to measure execution time 
t = time.time()


# Truth table (only used within the testbench, but loaded here for compatibility)
X_train = np.zeros([3,8], dtype=bool)
Y_train = np.zeros([2,8], dtype=bool)


# Loading grammar file
GRAMMAR_FILE = 'full_adder.bnf'

f = open("grammars/" + GRAMMAR_FILE, "r")
print(f.read())
f.close() 

BNF_GRAMMAR = grape.Grammar(os.path.join("grammars", GRAMMAR_FILE))

# Fitness evaluation using Icarus Verilog

def Icarus_fitness_eval(individual,points):
    
    if individual.invalid == True:
        return np.NaN,
    
    if isinstance(individual.phenotype, str):
        # Getting the contents of the testbench template
        with open('individual_tb.sv', 'r') as file:
            tb_template = file.read()
               
        # Creating the temporary testbench for this individual
        filename = "temp_tb.sv"
        with open(filename, "a") as testbench:
            
            # Writing the template
            testbench.write(tb_template)
            
            # Header
            testbench.write("module individual(input wire a, input wire b, input wire ci, output wire sum, output wire co);\n")
            
            # Phenotype
            testbench.write(individual.phenotype)
            
            # Footer
            testbench.write("\nendmodule ")
        
        #TODO: CHECK EXECUTION TIME
        
        # Get fitness from call to simulator
        fitness = float(os.popen("iverilog '-g2012' temp_tb.sv -o a.out && vvp a.out").read())
        
        # Delete temp testbench
        os.system("rm -f temp_tb.sv")
             
        return fitness,
    else:
        return 0,


toolbox = base.Toolbox()

# define a single objective, maximising fitness strategy:
creator.create("FitnessMax", base.Fitness, weights=(1.0,))

creator.create('Individual', grape.Individual, fitness=creator.FitnessMax)

#toolbox.register("populationCreator", grape.sensible_initialisation, creator.Individual) 
#toolbox.register("populationCreator", grape.random_initialisation, creator.Individual) 
toolbox.register("populationCreator", grape.PI_Grow, creator.Individual) 

toolbox.register("evaluate", Icarus_fitness_eval)

# Tournament selection:
toolbox.register("select", tools.selTournament, tournsize=6)

# Single-point crossover:
toolbox.register("mate", grape.crossover_onepoint)

# Flip-int mutation:
toolbox.register("mutate", grape.mutation_int_flip_per_codon)

POPULATION_SIZE = 250
MAX_GENERATIONS = 50
P_CROSSOVER = 0.9
P_MUTATION = 0.01
ELITE_SIZE = round(0.01*POPULATION_SIZE) #it should be smaller or equal to HALLOFFAME_SIZE
HALLOFFAME_SIZE = round(0.01*POPULATION_SIZE) #it should be at least 1

CODON_CONSUMPTION = 'eager'
RANDOM_SEED = 0 #Pay attention that the seed is set up inside the loop of runs, so you are going to have similar runs

INIT_GENOME_LENGTH = 50 #used only for random initialisation
random_initilisation = False #put True if you use random initialisation

MAX_INIT_TREE_DEPTH = 100
MIN_INIT_TREE_DEPTH = 10
MAX_TREE_DEPTH = 100
MAX_WRAPS = 0
CODON_SIZE = 255

N_RUNS = 3

for i in range(N_RUNS):
    print()
    print()
    print("Run:", i+1)
    print()
    
    random.seed(i) #Comment this line or set a different RANDOM_SEED each run if you want distinct results
    
    # create initial population (generation 0):
    if random_initilisation:
        population = toolbox.populationCreator(pop_size=POPULATION_SIZE, 
                                           bnf_grammar=BNF_GRAMMAR, 
                                           #init_genome_length=INIT_GENOME_LENGTH,
                                           min_init_depth=MIN_INIT_TREE_DEPTH,
                                           max_init_depth=MAX_INIT_TREE_DEPTH, 
                                           codon_size=CODON_SIZE,
                                           codon_consumption=CODON_CONSUMPTION
                                           )
    else:
        population = toolbox.populationCreator(pop_size=POPULATION_SIZE, 
                                           bnf_grammar=BNF_GRAMMAR, 
                                           min_init_depth=MIN_INIT_TREE_DEPTH,
                                           max_init_depth=MAX_INIT_TREE_DEPTH,
                                           codon_size=CODON_SIZE,
                                           codon_consumption=CODON_CONSUMPTION
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
                                              stats=stats, halloffame=hof, verbose=False)
    
    # import textwrap
    # best = hof.items[0].phenotype
    # print("Best individual: \n","\n".join(textwrap.wrap(best,80)))
    # print("\nTraining Fitness: ", hof.items[0].fitness.values[0])
    # print("Depth: ", hof.items[0].depth)
    # print("Length of the genome: ", len(hof.items[0].genome))
    # print(f'Used portion of the genome: {hof.items[0].used_codons/len(hof.items[0].genome):.2f}')
    
    # max_fitness_values, mean_fitness_values = logbook.select("max", "avg")
    # min_fitness_values, std_fitness_values = logbook.select("min", "std")
    # best_ind_length = logbook.select("best_ind_length")
    # avg_length = logbook.select("avg_length")

    # selection_time = logbook.select("selection_time")
    # generation_time = logbook.select("generation_time")
    # gen, invalid = logbook.select("gen", "invalid")
    # avg_used_codons = logbook.select("avg_used_codons")
    # best_ind_used_codons = logbook.select("best_ind_used_codons")
    
    # best_ind_nodes = logbook.select("best_ind_nodes")
    # avg_nodes = logbook.select("avg_nodes")

    # best_ind_depth = logbook.select("best_ind_depth")
    # avg_depth = logbook.select("avg_depth")

    # structural_diversity = logbook.select("structural_diversity") 
    
    # import csv
    # import random
    # r = random.randint(1,1e10)
    
    # header = ['gen', 'invalid', 'avg', 'std', 'min', 'max',
    #           'best_ind_length', 'avg_length', 
    #           'best_ind_nodes', 'avg_nodes', 
    #           'best_ind_depth', 'avg_depth', 
    #           'avg_used_codons', 'best_ind_used_codons', 
    #           'structural_diversity',
    #           'selection_time', 'generation_time']
    # with open("results/" + str(r) + ".csv", "w", encoding='UTF8', newline='') as csvfile:
    #     writer = csv.writer(csvfile, delimiter='\t')
    #     writer.writerow(header)
    #     for value in range(len(max_fitness_values)):
    #         writer.writerow([gen[value], invalid[value], mean_fitness_values[value],
    #                          std_fitness_values[value], min_fitness_values[value],
    #                          max_fitness_values[value], 
    #                          best_ind_length[value], 
    #                          avg_length[value], 
    #                          best_ind_nodes[value],
    #                          avg_nodes[value],
    #                          best_ind_depth[value],
    #                          avg_depth[value],
    #                          avg_used_codons[value],
    #                          best_ind_used_codons[value], 
    #                          structural_diversity[value],
    #                          selection_time[value], 
    #                          generation_time[value]])
    
    
### Statistics

# # Genetic Programming is done (all runs) - plot statistics:
x_axis = np.arange(0, MAX_GENERATIONS+1)

max_fitness_values, mean_fitness_values = logbook.select("max", "avg")
min_fitness_values, std_fitness_values = logbook.select("min", "std")
best_ind_length = logbook.select("best_ind_length")
avg_length = logbook.select("avg_length")
max_length = logbook.select("max_length")
selection_time = logbook.select("selection_time")
generation_time = logbook.select("generation_time")
gen, invalid = logbook.select("gen", "invalid")

header = ['gen', 'invalid', 'avg', 'std', 'min', 'max', 'best_ind_length', 'avg_length', 'max_length', 'selection_time', 'generation_time']

results = pd.DataFrame(list(zip(gen, invalid, mean_fitness_values, std_fitness_values, min_fitness_values, max_fitness_values, best_ind_length, avg_length, max_length, selection_time, generation_time)),
               columns = header)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 10000)
pd.set_option('display.colheader_justify', 'center')

#pd.display(results)

import textwrap
best = hof.items[0].phenotype

print("\n\nResults of run")
print("Best individual: \n","\n".join(textwrap.wrap(best,80)))
print("\nBest Fitness: ", hof.items[0].fitness.values[0])
print("Depth: ", hof.items[0].depth)
print("Length of the genome: ", len(hof.items[0].genome))
print(f'Used portion of the genome: {hof.items[0].used_codons/len(hof.items[0].genome):.2f}')


elapsed = time.time() - t
print('Elapsed Time: ', elapsed)
print('========')

plt.plot(max_fitness_values, color='red')
plt.xlabel('Generations')
plt.ylabel('Best Fitness')
#plt.title('Max and Average Fitness over Generations')
plt.title('Best Fitness over Generations')
plt.show()

gen = np.arange(0, MAX_GENERATIONS+1)

plt.xlabel('Generations')
plt.ylabel('Average Fitness')
plt.title('Average Fitness over Generations')
plt.errorbar(gen, mean_fitness_values, yerr=std_fitness_values,label="Best", color="Green")
plt.show()