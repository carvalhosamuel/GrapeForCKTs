#TODO LIST
# More noise, different waveforms

#### New Experiments Checklist:

# Experiment name: xSNR_(tourn | double)_(sym | none)
#exp_name = "10SNR_double_sym"

exp_name = "NSGA2_test"

# line 64: change input data file name    
# generic_tb_template: change data file name 

  
import grape
import algorithms
from functions import not_, and_, or_, nand_, nor_
import os
import pandas as pd
import numpy as np
from deap import creator, base, tools
from deap.benchmarks.tools import diversity, convergence, hypervolume
import random

import time
import matplotlib.pyplot  as plt
import multiprocessing as mp
import string
from bitstring import Bits
from scipy.stats import pearsonr
import chime
from scipy.fftpack import fft, ifft
import csv


RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# Delete any temp testbench
os.system("rm -f ./tmp/temp_tb*")

# Create the plots folder
os.system("mkdir ./plots")

# Starting timer to measure execution time 
t = time.time()


# Truth table (only used within the testbench, but loaded here for compatibility)
X_train = np.zeros([3,8], dtype=bool)
Y_train = np.zeros([2,8], dtype=bool)


# Loading grammar file| <coef> <coef> <extra> | <coef> <coef> <extra>
GRAMMAR_FILE = 'FIR_generic_sym.bnf'

f = open("grammars/" + GRAMMAR_FILE, "r")
print(f.read())
f.close() 

BNF_GRAMMAR = grape.Grammar(os.path.join("grammars", GRAMMAR_FILE))


# Loading the data

input_data = np.genfromtxt('Data/input_1.1SNR.data', delimiter=',')
fft_input = fft(input_data)
fft_input = np.abs(fft_input)[0:int(len(input_data)/2)]


target = np.genfromtxt('Data/target.data', delimiter=',')
fft_target = fft(target)
fft_target = np.abs(fft_target)[0:int(len(target)/2)]

# Original RMSE

Original_RMSE = np.sqrt(np.mean((np.array(input_data)-np.array(target))**2))
print("Original RMSE: " + str(Original_RMSE))

# Original correlation
original, p_value = pearsonr(fft_input, fft_target)

# Function for generating random strings for temp testbench names

def id_generator(size=8, chars=string.ascii_uppercase + string.digits + string.ascii_lowercase):
    return ''.join(random.choice(chars) for _ in range(size))




def fft_fitness_NSGA2(individual,points):
   
  ######
  # This function will use the np.convolve function to calculate fitness
  ######  

  if individual.invalid == True:
      individual.RMSE = 1000000
      return -99,-99
  
  if isinstance(individual.phenotype, str):
            
        
      # Getting the contents of the phenotype
      pheno = individual.phenotype.split()
      
      coefficients = pheno[1:]
      
      # Checking and processing for symmetry in coefficients
      if pheno[0] == 's':
          
          original = coefficients
          reverse  = list(reversed(coefficients))
      
          coefficients = original + reverse
                 
      coefficients = np.array(coefficients,dtype=float)
                
      # Getting the convolution between the coefficients and the input data
      output = np.convolve(input_data,coefficients, mode='valid')
      
      # Padding to match Verilog simulation
      output =np.delete(output,-1)
      output =np.delete(output,0)
      size_diff=len(input_data)-len(output)
      output = np.pad(output,pad_width=(size_diff,0), mode='constant')
     
      # Getting the spectrum of the output    
      fft_output = fft(output)
      fft_output = np.abs(fft_output)[0:int(len(output)/2)]
      
          
      # Pearson's r correlation as fitness score
      fitness, p_value = pearsonr(fft_target, fft_output)
      
      # Correlation in time domain
      #corr_time, p_v2 = pearsonr(target,output) 
      
      # RMSE as fitness score  
      RMSE = np.sqrt(np.mean((np.array(output)-np.array(target))**2))
      
      
      # Saving the processed data in the individual object 
      individual.coefs_real = coefficients
      individual.RMSE = RMSE
      individual.time_response = output
      individual.freq_response = fft_output          
      individual.fitness_test = RMSE
      
       
      return fitness,(Original_RMSE - RMSE)
  else:
      return -99,-99


def fft_fitness_icarus(individual,points):
    
  if individual.invalid == True:
      individual.RMSE = 1000000
      return -99,
  
  if isinstance(individual.phenotype, str):
            
      # Getting the contents of the testbench template
      DUT = algorithms.testbench_builder(individual)
                 
             
      # Creating the temporary testbench for this individual
      filename = 'temp_tb_' + id_generator() + '.v'
      
      with open('./tmp/' + filename, "a") as testbench:
          # Writing the template header
          testbench.write(DUT)
          
      
      #TODO: CHECK EXECUTION TIME
      
      # Get fitness from call to simulator
      os_command = 'iverilog \'-g2005\' ./tmp/' + filename + ' -o ./tmp/a' + filename +'.out && vvp ./tmp/a' + filename +'.out' 
      
      # To run the golden individual:
      #os_command = 'iverilog -W all -g2012 ./FIR_tb.v -o ./FIR_tb.out && vvp FIR_tb.out'
      
      output = os.popen(os_command).read().strip()
      
      # Setting the "x" at the output of filter to zeros 
      output = output.replace("x", "0")
      
      # Split the string output into lines (as list)
      out_list = output.splitlines()
      
      # Array to save the result in double numbers
      double_array = []
      
      # 64 is our wordlength from testbench
      for line in out_list:
          out = line
          a = Bits(bin=out,length=64)
          double_array.append(a.f)
          
      fft_output = fft(double_array)
      fft_output = np.abs(fft_output)[0:int(len(double_array)/2)]
      
      # Pearson's r correlation as fitness score      
      fitness, p_value = pearsonr(fft_target, fft_output)
      
      # RMSE as fitness score      
      RMSE = np.sqrt(np.mean((np.array(double_array)-np.array(target))**2))
     
      
      individual.icarus_RMSE = RMSE
      individual.icarus_fitness = fitness
      
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
creator.create("FitnessMax", base.Fitness, weights=(1.0,1.0))

creator.create('Individual', grape.Individual, fitness=creator.FitnessMax)

#toolbox.register("populationCreator", grape.sensible_initialisation, creator.Individual) 
toolbox.register("populationCreator", grape.random_initialisation, creator.Individual) 
#toolbox.register("populationCreator", grape.PI_Grow, creator.Individual) 

toolbox.register("evaluate", fft_fitness_NSGA2)

# Tournament selection:
#toolbox.register("select", tools.selTournament, tournsize=6)
#toolbox.register("select", tools.selLexicase)
toolbox.register("select", tools.selNSGA2)
#toolbox.register("select", algorithms.selDoubleTournament, fitness_size=2, parsimony_size=1.4, fitness_first=True)


# Single-point crossover:
toolbox.register("mate", grape.crossover_onepoint)

# Flip-int mutation:
toolbox.register("mutate", grape.mutation_int_flip_per_codon)

POPULATION_SIZE = 150
MAX_GENERATIONS = 150
P_CROSSOVER = 0.9
P_MUTATION = 0.1
ELITE_SIZE = round(0.01*POPULATION_SIZE) #it should be smaller or equal to HALLOFFAME_SIZE
HALLOFFAME_SIZE = max(round(0.01*POPULATION_SIZE),1) #it should be at least 1

CODON_CONSUMPTION = 'lazy'
GENOME_REPRESENTATION = 'list'
RANDOM_SEED = 0 #Pay attention that the seed is set up inside the loop of runs, so you are going to have similar runs

INIT_GENOME_LENGTH = 5 #used only for random initialisation
random_initilisation = True #put True if you use random initialisation

MAX_INIT_TREE_DEPTH = 100
MIN_INIT_TREE_DEPTH = 10

MIN_INIT_GENOME_LENGTH = 6
MAX_INIT_GENOME_LENGTH = 100

MAX_TREE_DEPTH = 100
MAX_WRAPS = 0
CODON_SIZE = 255

REPORT_ITEMS = ['gen',
 'invalid',
 'avg',
 'std',
 'min',
 'max',
 'fitness_test',
 'best_ind_length',
 'avg_length',
 'best_ind_nodes',
 'avg_nodes',
 'best_ind_depth',
 'avg_depth',
 'avg_used_codons',
 'best_ind_used_codons',
 'fitness_diversity',
 'structural_diversity',
 'evaluated_inds',
 'selection_time',
 'generation_time']

N_RUNS = 3

address = r'./Results/'+exp_name+'/'
  
command = 'mkdir ' + address +'' 
run_mkdir = os.popen(command).read().strip()

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
    #hof = tools.HallOfFame(HALLOFFAME_SIZE)
    hof = tools.ParetoFront()
    hof.maxsize = 20
    
    # prepare the statistics object:
    stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
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
                                              points_test=[X_train, Y_train],
                                              codon_consumption=CODON_CONSUMPTION,
                                              report_items=REPORT_ITEMS,
                                              genome_representation=GENOME_REPRESENTATION,
                                              stats=stats, halloffame=hof, verbose=False, showplots=False)
    
    # Printing outcome from this run
    lock.acquire()
    time.sleep(0.1)
    print("\n\nResults of run", r+1, "of",N_RUNS)
    print("\nMax FFT Correlation: ", np.nanmax(logbook.select("max")))
    print("Best RMSE: ", np.nanmin(logbook.select("fitness_test")))
    print("Final population hypervolume is %f" % hypervolume(population,[original,0]))
    
    lock.release()
    
    
    # # Calculate the Pareto fronts
    fronts = tools.emo.sortLogNondominated(population, len(population))
    
    
    # x = list()
    # y = list()
    
    # for ind in population:
    #     if ind.invalid != True:
    #         x.append(ind.fitness.values[0]) 
    #         y.append(ind.fitness.values[1]) 
    
    # plt.xlabel('FFT Correlation')
    # plt.ylabel('1/RMSE')
    # plt.title('Pareto Front of run ' + str(r+1))
    # plt.xlim(0.9, 1)
    # #plt.ylim(0, 20)
    # plt.scatter(x, y, marker="+", color='blue')
    # x = list(ind.fitness.values[0] for ind in fronts[0])
    # y = list(ind.fitness.values[1] for ind in fronts[0])
    # plt.scatter(x, y, marker="x", color='red')
    # plt.show()
    
    
    #valListFitness.append(val_fitness)
    #bestPhenotipes.append(best)
    
    
    ### CSV Reports

    max_fitness_values, mean_fitness_values = logbook.select("max", "avg")
    min_fitness_values, std_fitness_values = logbook.select("min", "std")
    fitness_test = logbook.select("fitness_test")
    
    best_ind_length = logbook.select("best_ind_length")
    avg_length = logbook.select("avg_length")

    selection_time = logbook.select("selection_time")
    generation_time = logbook.select("generation_time")
    gen, invalid = logbook.select("gen", "invalid")
    avg_used_codons = logbook.select("avg_used_codons")
    best_ind_used_codons = logbook.select("best_ind_used_codons")
    
    best_ind_nodes = logbook.select("best_ind_nodes")
    avg_nodes = logbook.select("avg_nodes")

    best_ind_depth = logbook.select("best_ind_depth")
    avg_depth = logbook.select("avg_depth")

    structural_diversity = logbook.select("structural_diversity") 
    fitness_diversity = logbook.select("fitness_diversity")     
    evaluated_inds = logbook.select("evaluated_inds") 
    
    
    header = REPORT_ITEMS
    
           
    with open(address + str(r) + ".csv", "w", encoding='UTF8', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(header)
        for value in range(len(max_fitness_values)):
            writer.writerow([gen[value], invalid[value], mean_fitness_values[value],
                             std_fitness_values[value], min_fitness_values[value],
                             max_fitness_values[value], 
                             fitness_test[value],
                             best_ind_length[value], 
                             avg_length[value], 
                             best_ind_nodes[value],
                             avg_nodes[value],    
                             best_ind_depth[value],
                             avg_depth[value],
                             avg_used_codons[value],
                             best_ind_used_codons[value], 
                             fitness_diversity[value],
                             structural_diversity[value],
                             evaluated_inds[value],
                             selection_time[value], 
                             generation_time[value]])
            
            

    return logbook, r+1 , hof.items[0], fronts



### Running the evoloop in parallel

# Defines our pool of processors using all available ones
pool = mp.Pool(mp.cpu_count())

# Runs the evoloop function in parallel, passing only r as the argument 
results = pool.starmap_async(evoloop, [([r]) for r in range(0, N_RUNS)]).get()

pool.close()

# Waits for all queued processes to reach here before resuming
pool.join() 


### Getting the objects from the runsint_array

logbooks, runs, hallof, pareto_fronts = zip(*results)

 


### Statistics

logbooks_bkp = logbooks
hallof_bkp = hallof

# GE is done (all runs) - plot statistics:
x_axis = np.arange(0, MAX_GENERATIONS+1)
# avgArray = np.array(avgListFitness)
# stdArray = np.array(stdListFitness)
# minArray = np.array(minListFitness)
# maxArray = np.array(maxListFitness)

# testArray = np.array(testListFitness)

max_Fitness_List = []
min_Fitness_List = []
avg_Fitness_List = []
RMSE_List = []


for Logbook in logbooks:
    
    max_Fitness_List.append(Logbook.select("max"))
    avg_Fitness_List.append(Logbook.select("avg"))
    min_Fitness_List.append(Logbook.select("min"))
    RMSE_List.append(Logbook.select("fitness_test"))



max_Fitness_List = np.array(max_Fitness_List)
avg_Fitness_List = np.array(avg_Fitness_List)
min_Fitness_List = np.array(min_Fitness_List)
RMSE_List = np.array(RMSE_List)


plt.xlabel('Generation')
plt.ylabel('Fitness (Correlation)')
plt.title('Average Max Fitness FIR Filter')
#plt.ylim([0.4, 1])
plt.errorbar(x_axis, max_Fitness_List.mean(0), yerr=max_Fitness_List.std(0),label="Avg", color="deepskyblue", )
plt.savefig(address + 'Fitness.png')
plt.show()

plt.xlabel('Generation')
plt.ylabel('RMSE')
plt.title('Average RMSE FIR Filter')
#plt.ylim([0.4, 1])
plt.errorbar(x_axis, RMSE_List.mean(0), yerr=RMSE_List.std(0),label="Avg", color="deepskyblue", )
plt.savefig(address + 'RMSE.png')
plt.show()

# Plotting the best individuals

best_Phenotypes_List = []
DUTS =[]    
flat_paretos =[]


count = 0
for run in pareto_fronts:
    for individual in run[0]: 
        best_Phenotypes_List.append(individual.phenotype)
        #algorithms.plotResponse_generic(individual,input_data,target)
        #fft_fitness_icarus(individual,[])
        flat_paretos.append(individual)
        ## Verilog execution
        #algorithms.plotResponse_verilog(individual,input_data,target)
        #DUTS.append(algorithms.testbench_builder(individual))
       
        count +=1
        #print(count)
        # algorithms.plotResponse_verilog(individual)

# DUTS =[]
# for individual in pareto_fronts[0][0]:
#     algorithms.plotResponse_verilog(individual,input_data,target)
#     DUTS.append(algorithms.testbench_builder(individual))




# Saving backup of the logbooks
import pickle
file_name = 'logbooks_'+exp_name+'.pkl'
with open(file_name, 'wb') as file:
    pickle.dump(logbooks_bkp, file)
    print(f'Object successfully saved to "{file_name}"')

# To restore    
# with open('./logbooks.pkl', 'rb') as f:
#     logbooks = pickle.load(f)

# Sound showing that code has finished running

command = 'mv ./plots ' + address +'plots' 
run_mkdir = os.popen(command).read().strip()

elapsed = time.time() - t
print('Elapsed Time: ', elapsed)
print('========')
chime.success()

#### Final Pareto Front (front from fronts)

# Calculate the Pareto fronts
final_fronts = tools.emo.sortLogNondominated(flat_paretos, len(flat_paretos))

x = list()
y = list()

for ind in flat_paretos:
    if ind.invalid != True:
        x.append(ind.fitness.values[0]) 
        y.append(ind.fitness.values[1]) 

plt.xlabel('FFT Correlation')
plt.ylabel('1/RMSE')
plt.title('Pareto Front of all the runs')
plt.xlim(original, 1)
plt.ylim(0, 1)
plt.scatter(x, y, marker="+", color='blue')
x = list(ind.fitness.values[0] for ind in final_fronts[0])
y = list(ind.fitness.values[1] for ind in final_fronts[0])
plt.scatter(x, y, marker="x", color='red')


plt.scatter(original, 0, marker="o", color='black')

plt.show()



 
    

