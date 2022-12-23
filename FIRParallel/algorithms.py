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

from deap import tools

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

# Function for generating random strings for temp testbench names

def id_generator(size=8, chars=string.ascii_uppercase + string.digits + string.ascii_lowercase):
    return ''.join(random.choice(chars) for _ in range(size))


def plotResponse(phenotype):
    
    # Target output obtained from "golden individual" (coefficients all set to 8'h20)
    target = [6656, 8864, 11072, 12832, 14592, 15424, 16256, 15424, 14592, 12608, 10624, 8224, 5824, 2880, -64, -2880, -5696, -8224, -10752, -12544, -14336, -14624, -14912, -13792, -12672, -10464, -8256, -5760]

   
    # Getting the contents of the testbench template
    with open('individual_tb_FIR.sv', 'r') as file:
        tb_template = file.read()
           
    # Creating the temporary testbench for this individual
    filename = 'temp_tb_' + id_generator() + '.sv'
    
    with open('./tmp/' + filename, "a") as testbench:
        
        # Phenotype
        testbench.write(phenotype)
       
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
    
    # Correlation as fitness score
    correlation, p_value = pearsonr(int_array, target)
    
    # Seeing output 
    plt.plot(int_array, color='red', label="FIR")
    plt.plot(target, color='blue', label="Target")
    plt.title('Best individual. Correlation:' +str(correlation))
    plt.legend(loc="upper right")
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
                print("gen =", gen, ", Best fitness =", halloffame.items[0].fitness.values, ", Number of invalids =", invalid, " Phenotype = ",halloffame.items[0].phenotype)
            if showplots:    
                plotResponse(halloffame.items[0].phenotype)
            if points_test:
                if gen < ngen:
                    fitness_test = np.NaN
                else:
                    fitness_test = toolbox.evaluate(halloffame.items[0], points_test)[0]
            
            
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

    return population, logbook