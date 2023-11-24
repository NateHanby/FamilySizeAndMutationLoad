import numpy as np
import model
import time
import csv
import os
import multiprocessing as mp
import datetime as dt
import secrets

previous_time = time.perf_counter()

# a simple timer function to test the speed of operations.
# this function should not be called in the final version.
def timer(text=None):

    global previous_time
    new_time = time.perf_counter()
    if text != None:
        print(text + ": " + str((new_time - previous_time)*10000//10) + "ms")
    
    previous_time = new_time

# an object to hold the model's settings that are being used for the current simulation.
class Settings:

    # returns the number of settings.
    def count(self):
        return len(self.keys)

    # gets a setting based on the setting name, using S["Setting Name"] syntax.
    def __getitem__(self, item):
        index = self.keys.index(item)
        datatype = self.datatypes[index]
        return datatype(self.values[index])

    # generates a list of settings should appear in the UI for the user to set them and then be used by the model.
    # the UI will automatically loop through all of these items and create boxes and labels for them.
    # this object is also used to store settings used by the command line version of the program. 
    def __init__(self):

        # used to retrieve the setting in the model's code. 
        # This is also the label of the setting in the user interface. 
        # if the key has a "(" character in it, the setting will automatically be displayed and grouped on the graphs
        # the symbol on the graph is the next character that follows "("
        self.keys = []

        # default values for every setting 
        self.values = []

        # the datatypes are the data type of each setting
        self.datatypes = []

        # whether the datatype can be changed if a second simulation is run with the same data
        self.mutable = []

        # the tooltips display when you mouse over the setting in the UI, explaining what it does.
        self.tooltips = []


        self.keys.append("Simulation Name")
        self.values.append("Enter Name")
        self.datatypes.append(str)
        self.mutable.append(True)
        self.tooltips.append("Simulation Name: \n" +
                            "After the simulation is run, data will be saved " +
                            "in the file system under the name entered here. " +
                            "If 'Keep Raw Data?' is set to 1, the data can be " +
                            "used as a starting point for future simulations, "
                            "and in all cases the output data can be graphed and "
                            "compared to other simulations. To find this data later, " 
                            "open the .np file with this simulation name in the " + 
                            "'Starter Data:' drop down box.")

        self.keys.append("Population Size (N)")
        self.values.append("10000")
        self.datatypes.append(int)
        self.mutable.append(True)
        self.tooltips.append("Population Size (N):\n" +
                        "The target number of individuals in each generational cohort. " +
                        "The population of the next generation will be this number or it will be " +
                        "the Legal Child Limit times the number of surviving members of one of the sexes, " +
                        "whatever number is lower. ")

        self.keys.append("Legal Child Limit (Z)")
        self.values.append("16")
        self.datatypes.append(int)
        self.mutable.append(True)
        self.tooltips.append("Legal Child Limit (Z): \n" +
                            "The maximum number of children a person can have. " +
                            "Can be lowered to simulate the effects of " +
                            "family-size restriction laws on the population.")
        
        self.keys.append("Marriage Flag (M)")
        self.values.append("1")        
        self.datatypes.append(bool)
        self.mutable.append(True)
        self.tooltips.append("Marriage Flag (M): \n" +
                            "Set if the population is monogamous.")
        
        # self.keys.append("Gene Linkage (X)")
        # self.values.append("1")        
        # self.datatypes.append(bool)
        # self.mutable.append(True)
        # self.tooltips.append("Gene Linkage (X): \n" +
        #                     "Set if genes can be linked--causing mutations at nearby loci to be likely to be inherited together. " +
        #                     "Crossover recombination is implemented. " +
        #                     "If this is disbled, the program will run slower and all mutations will assort completely independently. ")
        
        self.keys.append("Number of Generations (t)")
        self.values.append("25000")
        self.datatypes.append(int)
        self.mutable.append(True)
        self.tooltips.append("Number of Generations (t): \n" +
                        "The number of total generations that will be simulated.")

        self.keys.append("Deleterious Mutation Rate (U)")
        self.values.append("2.2")
        self.datatypes.append(float)
        self.mutable.append(True)
        self.tooltips.append("Deleterious Mutation Rate: \n" +
                        "The total number of deleterious mutations per diploid per generation.")
        
        self.keys.append("Genome Size (G)")
        self.values.append("20000")
        self.datatypes.append(int)
        self.mutable.append(False)
        self.tooltips.append("Genome Size: \n" +
                            "The number of distinct functionalities that can mutate. The total number of sets " +
                            "of loci in the genome, where each locus in the set will break the same functionality " +
                            "if it mutates as every other loci of the same set. " +
                            "Mutations within the same set of loci are not distinguished from each other. ")
        
        self.keys.append("Hard Selection Coefficient (s)")
        self.values.append(".014")
        self.datatypes.append(float)
        self.mutable.append(True)
        self.tooltips.append("Hard Selection Coefficient (s): \n" +
                "A proportional value that indicates how much the average heterozygous mutation reduces " +
                "an individuals chance of surviving to adulthood. This represents absolute fitness (" +
                " i.e. hard selection). This can be decreased to simulate the development of modern medicine.")


        self.keys.append("Soft Selection Coefficient (c)")
        self.values.append(".071")
        self.datatypes.append(float)
        self.mutable.append(True)
        self.tooltips.append("Competition Penalty (c): \n" +
                "For each child of the next generation, the average heterozygous mutation reduces each person's odds of "
                "being the parent of that child by this proportion." )
        
        
        self.keys.append("Non-Selective Death Chance (Q)")
        self.values.append("0")
        self.datatypes.append(float)
        self.mutable.append(True)
        self.tooltips.append("Non-selective Death Chance (Q):\n" +
                "The proportion of individuals that due purely from bad luck having nothing to do with " +
                "their genetics. If this is zero, all deaths are selective deaths." )

        self.keys.append("Mutation Variability (Y)")
        self.values.append("4.35")
        self.datatypes.append(float)
        self.mutable.append(False)
        self.tooltips.append("Mutation Variability (Y): \n" +
                            "A metric of how much variability there is in the deleteriousness of each mutation. " +
                            "The deleteriuosness of each mutation follows a gamma distribution with shape = 1 / K, " +
                            "and scale = K * s [for absolulte fitness metrics] and " + 
                            "scale = K * c [for competition selection]. If K is zero, all mutations have the " +
                            "same effect.")

        self.keys.append("Fully Recessive Chance (F)")
        self.values.append("0")
        self.datatypes.append(float)
        self.mutable.append(False)
        self.tooltips.append("Fully Recessive Chance : \n" +
                        "Proportion of deleterious mutations that are fully recessive. ")       
        
        self.keys.append("Dominance Coefficient (h)")
        self.values.append(".2")
        self.datatypes.append(float)
        self.mutable.append(False)
        self.tooltips.append("Dominance Coefficient : \n" +
                        "Proportion of the fitness penalty of the average heterozygous mutation to the same mutation in " +
                        "homozygous state. Note that this is just the average: the effect of a mutation in homozygoust state and the " +
                        "same mutation in heterozygous state are generated independently. Note that the s and c paramters are the average penalties " +
                        "for mutations in the heterozygous state. ")       
        
        self.keys.append("Fitness Overlap (r)")
        self.values.append("0.5")
        self.datatypes.append(float)
        self.mutable.append(False)
        self.tooltips.append("Fitness Overlap (r) : \n" +
                        "Proportion of mutations that have the same penalty weighing for both hard and soft selection." +
                        "For other mutations, weightings are generated independently")
        
        self.keys.append("Keep Raw Data?")
        self.values.append("1")
        self.datatypes.append(bool)
        self.mutable.append(True)
        self.tooltips.append("Save Data?: \n" +
                        "Flag that determines if the raw mutation data generated by the simulation is stored on the hard drive or not. " +
                        "This data will be needed to use the current simulation as a starting point for future simulations. " +
                        "If you do not plan to use the current simulation as a starting point for future simulations, " +
                        "then you can save hard drive space by setting this flag to anything other than 1")
        

# when the run simulation button is pressed in the User Interface
def start_model_from_UI(view):

    #create a new settings object specifically for the upcoming simulation
    s = Settings()

    # get the settings for this simulation from the current values of the user-entered fields.
    for i in range(s.count()):
        s.values[i] = view.SettingFields[i].get()

    # get the name of data file with the data from the previous simulation, or "New Simulation" if fresh.
    if view.dboxSelect.get() == "New Simulation":
        previous_simulation_name = "New Simulation"
    else:
        previous_simulation_name = view.dboxSelect.get_folder() + "/" + view.dboxSelect.get()

    # number of times this simulation should be replicated
    replication_count = int(view.tboxReplicationCount.get())

    # a list of multi processing objects that the simulations are running on. 
    threads = []

    # start the model with these settings. 
    # pass the UI as the view parameter so that the UI can be updated with the progress of the simulation. 
    for replication_number in range(replication_count):
        threads.append(start_model(s, previous_simulation_name, replication_number, replication_count, view))

    return threads

# when the start_script.py is run from the command line.
def start_model_from_csv(filename, super_replication=None):

    if super_replication is None:
        SRA = ""
    else:
        SRA = "-" + str(super_replication)

    # generates distinct random seeds for every simulation for the random number generators.
    # so that the pseudorandom simulation results aren't correlated. 
    seed_count = 0
    with open(filename) as csv_file:
        ignore_first_row_because_it_is_the_header_flag = True
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:

            if ignore_first_row_because_it_is_the_header_flag:
                ignore_first_row_because_it_is_the_header_flag = False
                continue

            #number of replications
            seed_count += int(row[1])
    base_seq = np.random.SeedSequence(entropy=secrets.randbits(128))
    random_seeds = base_seq.spawn(seed_count)
    current_seed = 0

    MAXIMUM_THREADS_THAT_SHOULD_BE_RUN_AT_ONCE = 5

    threads = []
    for i in range(MAXIMUM_THREADS_THAT_SHOULD_BE_RUN_AT_ONCE):
        threads.append(None)

    #0: not done. 1: done. -1: deferred until all 0's done. 
    progress_list = []

    first_pass = True

    while True:
        
        current_row = 0
        with open(filename) as csv_file:
            
            ignore_first_row_because_it_is_the_header_flag = True
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:

                if ignore_first_row_because_it_is_the_header_flag:
                    ignore_first_row_because_it_is_the_header_flag = False
                    continue

                if first_pass == False:
                    if progress_list[current_row] == 1:
                        current_row = current_row + 1
                        continue

                if row[0] == 'New Simulation':
                    previous_simulation_name = 'New Simulation'
                    if first_pass == True:
                        progress_list.append(1)
                    else:
                        progress_list[current_row] = 1
                else:
                    asdf = row[0] + SRA + '.npy'
                    
                    # skip if deferred
                    if first_pass == False and progress_list[current_row] == -1:
                        current_row = current_row + 1
                        continue

                    if os.path.exists('simulations/' + row[0] + SRA + '.npy'):
                        previous_simulation_name = 'simulations/' + row[0] + SRA + '.npy'
                        if first_pass == True:
                            progress_list.append(1)
                        else:
                            progress_list[current_row] = 1
                    else:
                        if first_pass == True:
                            if row[0] != row[2]:
                                progress_list.append(0)
                            else:
                                progress_list.append(-1)
                        current_row = current_row + 1
                        continue

                replication_count = int(row[1])
                s = Settings()
                for i in range(s.count()):
                

                    # need this because bool('FALSE') and bool('0') both return true!!!
                    if s.datatypes[i] == bool:
                        if row[i+2] in ('FALSE', 'False', 'F', '0', 'Off', 'OFF', 'No', 'NO'):
                            s.values[i] = False
                        else:
                            s.values[i] = True
                    else:
                        if s.keys[i] == "Simulation Name":
                            s.values[i] = s.datatypes[i](row[i+2] + SRA)
                        else:
                            s.values[i] = s.datatypes[i](row[i+2])

                thread_number = 0
                for replication_number in range(replication_count):
                    simulation_assigned_to_thread = False
                    while not simulation_assigned_to_thread:
                        if threads[thread_number] is not None:
                            threads[thread_number].join(timeout=.1)
                        if threads[thread_number] is None or not threads[thread_number].is_alive():
                            threads[thread_number] = start_model(s, previous_simulation_name, replication_number, 
                                                                 replication_count, seed=random_seeds[current_seed])
                            current_seed += 1
                            simulation_assigned_to_thread = True
                            print("(controller) Simulation " + row[2] + SRA + " Starting " + str(dt.datetime.now()))
                        thread_number += 1
                        if thread_number >= MAXIMUM_THREADS_THAT_SHOULD_BE_RUN_AT_ONCE:
                            thread_number = 0
                
                current_row = current_row + 1
            
            first_pass = False
            threads_still_alive = False
            for thread in threads:
                if thread is not None and thread.is_alive():
                    threads_still_alive = True
            
            if threads_still_alive == False:
                if -1 in progress_list:
                    for i in range(len(progress_list)):
                        if progress_list[i] == -1:
                            progress_list[i] = 0
                else:
                    break

            time.sleep(0.1)

    # wait for any remaining threads to complete. 
    for thread in threads:
        if thread is not None:
            thread.join()
def start_model(s, previous_simulation_name, replication_number, replication_count, view=None, seed=None):

    if previous_simulation_name != "New Simulation":

        retry = 0
        while retry < 2:
            try:
                # load data from the file system about the previous simulation
                with open(previous_simulation_name, 'rb') as f:
                    np.load(f)  # don't need the log data right now
                    np.load(f, allow_pickle=True)  # don't need the list of log items right now.
                    np.load(f),
                    girls_mutations = np.load(f)
                    boys_mutations = np.load(f)
                    gene_surv_pen_ratios = np.load(f)
                    gene_hoz_surv_pen_ratios = np.load(f)
                    gene_comp_pen_ratios = np.load(f)
                    gene_hoz_comp_pen_ratios = np.load(f)
                    fully_recessive = np.load(f)
                    genome_map = np.load(f)
                
                break
            except:
                    if view is not None:
                        view.show_msg("Error: Could not load raw data from " + previous_simulation_name)
                        view.completed_simulations.value += 1
                        view.message_code.value = 0
                        return
                    else:
                        if retry == 0:
                            retry = 1
                            time.sleep(60)
                        else:
                            retry = 2
                            print("Error: Could not load raw data from " + previous_simulation_name)
                            return

    else:
        girls_mutations = None
        boys_mutations = None 
        gene_surv_pen_ratios = None
        gene_hoz_surv_pen_ratios = None
        gene_comp_pen_ratios = None
        gene_hoz_comp_pen_ratios = None
        fully_recessive = None 
        genome_map = None

    if replication_count == 1:
        thread_num = None
    else:
        thread_num = replication_number

    if view is None:
        thread = mp.Process(target=model.run_simulation, args=([s, girls_mutations, boys_mutations, gene_surv_pen_ratios, gene_hoz_surv_pen_ratios, 
                            gene_comp_pen_ratios, gene_hoz_comp_pen_ratios, fully_recessive, genome_map, thread_num, seed]))
    else:
        thread = mp.Process(target=model.run_simulation, args=([s, girls_mutations, boys_mutations, gene_surv_pen_ratios, gene_hoz_surv_pen_ratios, 
                            gene_comp_pen_ratios, gene_hoz_comp_pen_ratios, fully_recessive, genome_map, thread_num, seed,
                            view.simulations_stopped, view.log_updates, view.message_code, view.simulation_name_id, view.completed_simulations]))
    thread.start()
    
    return thread



def save_data(S, G, girls, boys, log, thread=None):

    if thread is None:
        savefilename = S["Simulation Name"] + '.npy'
    else:
        savefilename = S["Simulation Name"] + '_' + str(thread) + '.npy'

    with open('simulations/' + savefilename, 'wb') as f:
        np.save(f, log)
        np.save(f, SimulationLog.get_log_items())
        np.save(f, np.array(S.values))
        if S["Keep Raw Data?"] == 1:
            np.save(f, girls.mutations)
            np.save(f, boys.mutations)
            np.save(f, G.GENE_SURV_PEN_RATIOS)
            np.save(f, G.GENE_HOZ_SURV_PEN_RATIOS)
            np.save(f, G.GENE_COMP_PEN_RATIOS)
            np.save(f, G.GENE_HOZ_COMP_PEN_RATIOS)
            np.save(f, G.FULLY_RECESSIVE)
            np.save(f, G.GENOME_MAP)

class SimulationLog:

    @staticmethod
    def get_log_items():

        # list of items to be displayed in the log section of the user interface. 
        # this list is saved after every simulation and is then used by the UI to determine what statistics to display
        # note that old save files can perfectly well use different versions of the below list.
        # changing this list will affect the log items for future simulations only.
        # if the list is modified, the corresponding calculations in the below record_data method should also be modified. 
        # simulations can be compared and graphed together if they have a log item with the same string value. 
        LogItems = [("Expected Child Mortality %", float),
                    ("Average Survival Fitness", float),
                    ("Average Reproduction Chance", float),
                    ("Average Competition Weighting", float),
                    ("Survivor Competition Weighting", float),
                    ("Effective # of Survival Mutations", float),
                    ("Effective # of Competition Mutations", float),
                    ("Total People", int),
                    ("Deleterious Mutations Per Person", float),
                    ("Roughly Fixed Mutations", int),
                    ("Neutral Locus Variance", float),
                    ("Estimated Effective Population", float)]

        return LogItems
    
    @staticmethod
    def get_log_tooltip(index):
        
        LogItems = SimulationLog.get_log_items()
        tooltip = ""

        
        if LogItems[index][0] == "Expected Child Mortality %":
            tooltip = "Expected Child Mortality %: \n" +\
                      "Proportion of people who were expected to die before reproduction based on their genetics and the non-selective death chance. Can differ slightly from the actual percent that died."

        if LogItems[index][0] == "Average Survival Fitness":
            tooltip = "Average Survival Fitness: \n" +\
                      "Proportion of people who would be expected to survive to reproductive age from this generation based upon " +\
                      "their genetics, the Hard Selection Coefficient (s) setting."
            
        if LogItems[index][0] == "Average Reproduction Chance":
            tooltip = "Average Reproduction Chance: \n" +\
                      "Proportion of people this generation (out of all people, including those who died before reproductive age) who had " +\
                      "at least one child (whether or not that child survived)."
        
        if LogItems[index][0] == "Average Competition Weighting":
            tooltip = "Average Competition Weighting: \n" +\
                      "The average of a weighted proportion which represents each person's chances, relative to living people of the same sex and the same generation " +\
                      "of being the one to parent any particular child of the next generation. This average includes people who did not survive to " +\
                      "reproductive age. Since this is a metric of soft selection, it can be arbitrarily low without the population dying off. " +\
                      "it is useful to compare this valule to different populations or different generations to determine how soft fitness " +\
                      "changes over time. "
        
        if LogItems[index][0] == "Survivor Competition Weighting":
            tooltip = "Survivor Competition Weighting: \n" +\
                      "The average of a weighted proportion which represents each person's chances, relative to living people of the same sex and the same generation " +\
                      "of being the one to parent any particular child of the next generation. This average does not include people who did not survive to " +\
                      "reproductive age. Since this is a metric of soft selection, it can be arbitrarily low without the population dying off. " +\
                      "it is useful to compare this valule to different populations or different generations to determine how soft fitness " +\
                      "changes over time. "
            
        if LogItems[index][0] == "Effective # of Survival Mutations":
            tooltip = "Effective # of Survival Mutations: \n" +\
                      "A metric of genetic hard fitness which indicates the number of times the Hard Selection Coefficient (s) proportion should be multiplied together " +\
                      "to determine the average person's fitness according to hard selection. Equivalent to k in the common w(k) = (1-s)^k model of fitness. " +\
                      "this average is inclusive of people who did not survive. This represents a measure of genetic fitness which " +\
                      "is separate from the value of s, since s may change depending on the environment."
            
        if LogItems[index][0] == "Effective # of Competition Mutations":
            tooltip = "Effective # of Competition Mutations: \n" +\
                      "A metric of genetic soft fitness which indicates the number of times the Competition Penalty (c) proportion should be multiplied together " +\
                      "to determine the average person's fitness according to soft selection. Note that c is usually labeled as s in the literature, but in this model " +\
                      "s is used to represent hard selection only. Equivalent to k in the W = (1-c)^k model of fitness. " +\
                      "this average is inclusive of people who did not survive. This represents a measure of genetic fitness which " +\
                      "is separate from the value of c, since c may change depending on the environment. "
            
        if LogItems[index][0] == "Total People":
            tooltip = "Total People: \n" +\
                      "the total number of people that survived to adulthood this generation. "
            
        if LogItems[index][0] == "Deleterious Mutations Per Person":
            tooltip = "Deleterious Mutations Per Person: \n" +\
                      "Average number of deleterious mutations in the genome of each person. A heterozygous allele counts as one mutation, " +\
                      "and a homozygous allele counts as two."
            
        if LogItems[index][0] == "Roughly Fixed Mutations":
            tooltip = "Roughly Fixed Mutations: \n" +\
                      "number of deleterious mutations that are completely fixed in a 100 person sample of people from this population. " +\
                      "if deleterious mutations are fixing at a high rate, that is a good indication that the effective population size " +\
                      "may be unhealthily low in this population."
            
        if LogItems[index][0] == "Neutral Locus Variance":
            tooltip = "Neutral Locus Variance: \n" +\
                      "the variance in a neutral locus that changes by an amount with a variance of 1 and a mean of zero " +\
                      "every generation. The Neutral Locus Variance is extremely noisy but averaging it over thousands of generations gives a rough " +\
                      "estimate of the effective population size."
        
        if LogItems[index][0] == "Estimated Effective Population":
            tooltip = "Estimated Effective Population: \n" +\
                      "An very rough estimate of effective population size, based on the variance in a neutral locus."

        return tooltip

    @staticmethod
    def record_data(S, log, men, women, generation_number):

        # calculations for the data that is reported to the user interface and the log file. This data can be graphed later. 
        if men.person_count == 0 and women.person_count == 0:
            log[generation_number]["Expected Child Mortality %"] = 0
            log[generation_number]["Average Survival Fitness"] = 0
            log[generation_number]["Average Reproduction Chance"] = 0
            log[generation_number]["Average Competition Weighting"] = 0
            log[generation_number]["Survivor Competition Weighting"] = 0
            log[generation_number]["Effective # of Survival Mutations"] = 0
            log[generation_number]["Effective # of Competition Mutations"] = 0
            log[generation_number]["Total People"] = 0
            log[generation_number]["Deleterious Mutations Per Person"] = 0
            log[generation_number]["Roughly Fixed Mutations"] = 0
            log[generation_number]["Neutral Locus Variance"] = 0
            log[generation_number]["Estimated Effective Population"] = 0
        else:
            total_people = sum(men.alive_flags) + sum(women.alive_flags)
            survival_chances = np.average(np.concatenate((men.survival_chances, women.survival_chances)))
            survival_fitness = np.average(np.concatenate((men.survival_fitness, women.survival_fitness)))
            log[generation_number]["Expected Child Mortality %"] = (1 - survival_chances)*100
            log[generation_number]["Average Survival Fitness"] = survival_fitness
            log[generation_number]["Average Reproduction Chance"] = np.average(np.where(np.concatenate((men.child_count,women.child_count)) >0, 1, 0))
            log[generation_number]["Average Competition Weighting"] = np.average(np.concatenate((men.competition_chances, women.competition_chances)))
            log[generation_number]["Survivor Competition Weighting"] = np.average(np.concatenate((men.competition_chances[np.where(men.alive_flags,True,False)],  
                                                                                                  women.competition_chances[np.where(women.alive_flags,True,False)])))
            log[generation_number]["Effective # of Survival Mutations"] = np.average(np.concatenate((men.surv_pens, women.surv_pens)))
            log[generation_number]["Effective # of Competition Mutations"] = np.average(np.concatenate((men.attr_pens, women.attr_pens)))
            log[generation_number]["Total People"] = total_people
            log[generation_number]["Deleterious Mutations Per Person"] = (np.sum(women.mutations) + np.sum(men.mutations))/(len(men.alive_flags) + len(women.alive_flags))
            log[generation_number]["Roughly Fixed Mutations"] = np.sum(np.all((women.mutations[:100,:,:]), axis=(0, 2)))
            log[generation_number]["Neutral Locus Variance"] = np.var(np.concatenate((men.neutral_locus, women.neutral_locus)))
            if generation_number >= 19999:
                log[generation_number]["Estimated Effective Population"] = np.average(log[10000: generation_number]["Neutral Locus Variance"])/2
            else:
                log[generation_number]["Estimated Effective Population"] = float("NaN")

        




    @staticmethod
    def create_log_array(S):
        dtype = SimulationLog.get_log_items()

        # Create arrays to hold Future Log Data for this simulation.
        log = np.stack(np.array(np.zeros((S["Number of Generations (t)"]+1), dtype)))
        return log