import numpy as np
import controller as ctr
import time
import datetime as dt
import secrets


# a simple timer function to test the speed of operations.
# this function should not be called in the final version.
previous_time = time.perf_counter()
def timer(text=None):

    global previous_time
    new_time = time.perf_counter()
    if text != None:
        print(text + ": " + str((new_time - previous_time)*10000//10) + "ms")
    
    previous_time = new_time


# the entry point of the simulation.
# this function is called by the "Run Simulation" button on the User Interface.
def run_simulation(S, girls_mutations, boys_mutations, gene_surv_pen_ratios, gene_hoz_surv_pen_ratios, gene_comp_pen_ratios, gene_hoz_comp_pen_ratios,
                   fully_recessive, genome_map, thread=None, seed=None, simulations_stopped=None, log_updates=None, message_code=None,
                   simulation_name_id=None, completed_simulations=None):
    try:

        if seed == None:
            rng = np.random.default_rng(seed=np.random.SeedSequence().entropy)
        else:
            rng = np.random.default_rng(seed=seed)

        if log_updates == None:
            if thread is None:
                thread_text = ""
            else:
                thread_text = "_" + str(thread)
            print("(model) Simulation " + S["Simulation Name"] + thread_text + " Starting " + str(dt.datetime.now()))
            
        # First, genome meta data is generated or imported from a previous simulation.
        # The args contain the genomic meta data (such as the deleteriousness of each possible mutation)
        # If the args are None, then this is a new simulation and fresh genomic meta data is generated according to the settings.
        G = Genome(S, rng, gene_surv_pen_ratios, gene_hoz_surv_pen_ratios, gene_comp_pen_ratios, gene_hoz_comp_pen_ratios, fully_recessive, genome_map)

        # if this is a new simulation, generate an empty set of mutations for each person.
        # otherwise, the mutation data contained in the args that was imported from a previous simulation will be used as a starting point.  
        if girls_mutations is None:
            girls_mutations = np.zeros((S["Population Size (N)"]//2, S["Genome Size (G)"], 2))
            boys_mutations = np.zeros((S["Population Size (N)"]//2, S["Genome Size (G)"], 2))

        # create the first generation of boys and girls based on the mutation data, the simulation settings, S, and the genomic metadata, G.
        girls = Population(S, rng, G, mutations=girls_mutations)
        boys = Population(S, rng, G, mutations=boys_mutations)

        #create some arrays to store log information about this simulation to be used later to record statistics and make graphs. 
        log = ctr.SimulationLog.create_log_array(S)

        # main loop.
        # loops once per generation
        # +1 is added because the very first iteration of this loop is from the previous simulation or from scratch, 
        # so it shouldn't count as a generation from the user's perspective.
        for generation_number in range(int(S["Number of Generations (t)"])+1):

            # some members of the new generation will die due to bad genetics or luck
            girls.chance_of_dying(rng)
            boys.chance_of_dying(rng)

            # the adults of the new generation are the children of the previous generation.
            women = girls
            men = boys

            # determines how many children each person will have. Genetic competitive fitness increases odds of reproducing.
            # returns lists of the parents of all future boys and girls for the next generation
            girls_mothers, girls_fathers, boys_mothers, boys_fathers = assignment_of_children(S, rng, women, men)

            # creates new children and their mutation data, including new mutations and mutations inherited from parents.
            girls = Population(S, rng, G, girls_mothers, girls_fathers, women.mutations, men.mutations, women.neutral_locus, men.neutral_locus)
            boys = Population(S, rng, G, boys_mothers, boys_fathers, women.mutations, men.mutations, women.neutral_locus, men.neutral_locus)

            # record data about this generation in the log arrays
            ctr.SimulationLog.record_data(S, log, men, women, generation_number)

            # program interacts with the user interface after every loop.
            # this step is skipped if the simulation was run from a CSV file. 
            if (simulations_stopped is not None):

                # if the cancel simulation button has been pressed or the program has been closed, abort simulation without saving data. 
                if simulations_stopped.value == 1:
                    if message_code.value < 2: message_code.value = 2
                    completed_simulations.value +=1
                    return
                if simulations_stopped.value == 2:
                    return
                
                # these fields are visible to multiple threads, changing them will update the user interface
                if message_code.value < 1: message_code.value = 1
                log_updates[0] = generation_number
                log_updates[1] = S["Number of Generations (t)"]
                log_updates[2] = simulation_name_id
                log_updates[3:] = log[generation_number]

        # when the simulation is finished, record all log data in a .np data file. 
        # this data file can be used to generate graphs or to run additional simulations using the completed simulation as a starting point. 
        ctr.save_data(S, G, girls, boys, log, thread)

        if log_updates == None:
            print("(model) Simulation " + S["Simulation Name"] + thread_text + " Complete " + str(dt.datetime.now()))
    except MemoryError:
        if simulations_stopped is None:
            print("Error: Insufficient memory. At least one simulation aborted.")
        else:
            message_code.value = 3
    

    if completed_simulations is not None:
        completed_simulations.value += 1


# This class holds constant genomic meta data, for example: the fitness effect of each possible mutation that could occur.
# Note: records for individual mutations for each person are not stored here, they are in the mutation array of the Population objects.
# this class stores data ABOUT possible mutations, not the mutations themselves.
class Genome:
 
    def __init__(self, S, rng, GENE_SURV_PEN_RATIOS=None, GENE_HOZ_SURV_PEN_RATIOS=None, GENE_COMP_PEN_RATIOS=None, GENE_HOZ_COMP_PEN_RATIOS=None,
                 FULLY_RECESSIVE=None, GENOME_MAP=None):

        # number of unique chromosomes in the human genome (haploid set)
        self.CHROMOSOME_COUNT = 23

        # the average number of crossover events in the human genome
        self.CROSSOVERS_PER_GENOME = 45.5

        # Estimate of how many genes each human chromosome has. Source: NIH. This info is outdated, though. 
        # I am not bothering to track Y-chromosomes -- the 23rd chromosome is treated as an autosome. The estimate here is for the X chromosome. 
        CHROMOSOME_SIZES = [3000, 2500, 1900, 1600, 1700, 1900, 1800, 1400, 1400, 1400, 2000, 1600, 800, 1200, 1200, 1300, 1600, 600, 1700, 900, 400, 800, 1400]

        #if this genome is not being imported from a previous simulation:
        if GENE_SURV_PEN_RATIOS is None:

            # Randomly generate the fitness penalties and dominance of each mutation according to disttribution specified in the settings
            # SURV: Hard Selection Coefficient Ratios. COMP: Soft Selection. HOZ: Penalty in Homozygous form. 

            # if Mutation Variability (Y) is zero, all mutations have the same proportional fitness penalties.
            if S["Mutation Variability (Y)"] == 0:
                self.GENE_SURV_PEN_RATIOS = np.ones(S["Genome Size (G)"])
                self.GENE_COMP_PEN_RATIOS = np.ones(S["Genome Size (G)"])
                self.GENE_HOZ_SURV_PEN_RATIOS = np.ones(S["Genome Size (G)"])/(S["Dominance Coefficient (h)"])
                self.GENE_HOZ_COMP_PEN_RATIOS = np.ones(S["Genome Size (G)"])/(S["Dominance Coefficient (h)"])
            else:
                # variability of fitness effects of absolute-fitness-reducing mutations follows a gamma distribution with shape = 1/Y and scale = Y
                self.GENE_SURV_PEN_RATIOS = rng.gamma(1/S["Mutation Variability (Y)"], S["Mutation Variability (Y)"], size=S["Genome Size (G)"])

                # for mutations in homozygous state, the fitness effect is uncorrelated to the effect in heterozygous state, 
                # but on average it is greater, according to the Dominance Coefficient (h) setting. 
                self.GENE_HOZ_SURV_PEN_RATIOS = self.GENE_SURV_PEN_RATIOS + rng.gamma(1/S["Mutation Variability (Y)"], 
                                                                S["Mutation Variability (Y)"], size=S["Genome Size (G)"])*(1/S["Dominance Coefficient (h)"]-1)

                # the Fitness Overlap (r) is the percentage of mutations where the fitness effects between hard and soft selection are perfectly correlated. 
                # the rest of the mutation's competition fitness penalties are generated independently along a gamma distribution. 
                correlation_rolls = rng.random(size=S["Genome Size (G)"])
                correlation_rolls = np.where(correlation_rolls < S["Fitness Overlap (r)"], 1, 0)

                self.GENE_COMP_PEN_RATIOS = np.where(correlation_rolls, 
                                                     self.GENE_SURV_PEN_RATIOS,
                                                     rng.gamma(1/S["Mutation Variability (Y)"], S["Mutation Variability (Y)"], size=S["Genome Size (G)"]) )

                self.GENE_HOZ_COMP_PEN_RATIOS = np.where(correlation_rolls, 
                                                     self.GENE_HOZ_SURV_PEN_RATIOS,
                                                     self.GENE_COMP_PEN_RATIOS + rng.gamma(1/S["Mutation Variability (Y)"], S["Mutation Variability (Y)"], size=S["Genome Size (G)"]) 
                                                     *(1/S["Dominance Coefficient (h)"]-1))

            # Genes can be fully recessive according to a random factor R.
            # A random number is generated between 0 and 1. If this number is less than R, the gene is fully recessive. 
            # if R is 1, all genes are fully recessive. 
            dominance_rolls = rng.random(size=S["Genome Size (G)"])
            self.FULLY_RECESSIVE = np.where(dominance_rolls < S["Fully Recessive Chance (F)"], 1, 0)

            # For the purposes of recombination and segregation, a GENOME_MAP is needed -- this will track the locations in the genome that 
            # segregate because they are at the boundaries of chromosomes. Using the above proportional list, a list of mutation 
            # identifiers will be created which record the location in the array of genes the first gene on every chromosome. 
            genome_map = []

            # the first gene of the first chromosome is always the 0th gene in the list of genes. 
            genome_current_gene_number = 0

            # the first gene of the first chromosome is added to the genome map, marking the starting point of that chromosome.
            genome_map.append(0)
            
            # take the list of the lengths (in number of genes) of each chromosome.
            for chromosome_gene_count in CHROMOSOME_SIZES:

                # the starting locations of the other chromosomes are added to the genome map, based on the list of lengths of each chromosome.
                genome_current_gene_number += chromosome_gene_count
                genome_map.append(genome_current_gene_number)

            #convert the list to a numpy array
            genome_map = np.array(genome_map)

            #rescales the genome_map based on the number of genes in the current simulation.
            self.GENOME_MAP = genome_map * float(S["Genome Size (G)"]) // genome_current_gene_number

        # alternatively, import genome data from a previous simulation 
        # so that multiple simulations can be run using the same randomly generated genome.
        else:
            self.GENE_SURV_PEN_RATIOS = GENE_SURV_PEN_RATIOS
            self.GENE_HOZ_SURV_PEN_RATIOS = GENE_HOZ_SURV_PEN_RATIOS
            self.GENE_COMP_PEN_RATIOS = GENE_COMP_PEN_RATIOS
            self.GENE_HOZ_COMP_PEN_RATIOS = GENE_HOZ_COMP_PEN_RATIOS
            self.FULLY_RECESSIVE = FULLY_RECESSIVE
            self.GENOME_MAP = GENOME_MAP


# Population Object
# this stores data for the people in the population, including data about what mutations they carry. 
class Population:

    def __init__(self, S, rng, G, mothers=None, fathers=None,
                 maternal_mutations=None, paternal_mutations=None,
                 maternal_neutral_locus=None, paternal_neutral_locus=None, mutations=None):

        # Loading the first generation of a simulation
        if mothers is None:

                # load mutation data for the new population
                self.mutations = mutations

                # the number of people who we have mutation data for
                self.person_count = len(mutations[:, 0, 0])

                # It is not worthwhile to pull this data from the previous simulation. 
                # due to the noisy nature of genetic drift, estimates of effective population are only attempted on simulations with >1000 generations.
                self.neutral_locus = np.zeros((self.person_count, 2), np.float)

        # If this is not the first generation, the new generation of people inherits genes from the previous generation of the same simulation
        else:

            # initialize data
            self.person_count = len(mothers)
            self.mothers = np.array(mothers)
            self.fathers = np.array(fathers)
            self.attr_pens_index = None
            self.mig_pens_index = None

            # create blank arrays for each new person-- inherited genes/mutations will later be inserted.
            self.mutations = np.ones((self.person_count, S["Genome Size (G)"], 2), np.uint8)
            self.neutral_locus = np.zeros((self.person_count, 2), np.float)

            # total number of segregating segments of the genome
            # each segment of each chromosome has a 50-50 chance of entering the gamete. 
            # used for genetic recombination. 
            segs = int(G.CHROMOSOME_COUNT + G.CROSSOVERS_PER_GENOME * 2)

            # total number of possible crossover events in the whole genome. 
            # note that the number of possible crossovers is double the number of actual crossovers, since 
            # each segment has a 50-50 chance of separating from the previous segment. 
            crosses = int(G.CROSSOVERS_PER_GENOME * 2)

            # this records which segments of the two parental chromosomes for each parent of each person will enter into the gamete.
            segregation_array = rng.integers(2, size=(self.person_count, segs, 2), dtype=np.uint8)

            # this records which parental chromosome of each parent will donate the neutral locus to the child
            # the neutral locus is used to calculate effective population size.
            segregation_array_neutral = rng.integers(2, size=(self.person_count, 2), dtype=np.uint8)


            # each person inherits genes from their parents.
            for person_id in range(self.person_count):


                # retrieves the index of the parents in the parental genomic data from the lists of parents. 
                # This is necessary because the parental mutation/genomic data is indexed differently than the lists of parents. 
                mother_id = self.mothers[person_id]
                father_id = self.fathers[person_id]

                # using the above index, looks up the value neutral locus in the genomic data of each parent.
                # each person has two values for the neutral locus: one from each parent. 
                # based on the random value in segregation_array_neutral, this picks a one of the two neutral loci from each of the parents. 
                self.neutral_locus[person_id, 0] = maternal_neutral_locus[mother_id, segregation_array_neutral[person_id][0]]
                self.neutral_locus[person_id, 1] = paternal_neutral_locus[father_id, segregation_array_neutral[person_id][1]]


                # # if gene linkage is turned off, every locus segregates independently
                # if S["Gene Linkage (X)"] != 1:
                    
                #     maternal_rand_array = rng.integers(0,2,S["Genome Size (G)"])
                #     paternal_rand_array = rng.integers(0,2,S["Genome Size (G)"])

                #     for gene_id in range(S["Genome Size (G)"]):
                #         self.mutations[person_id, gene_id, 0] = maternal_mutations[mother_id, gene_id, maternal_rand_array[gene_id]]
                #         self.mutations[person_id, gene_id, 1] = paternal_mutations[father_id, gene_id, paternal_rand_array[gene_id]]

                # # otherwise, the genome is divided into segments which can cross over. 
                # else:

                # for each new person, randomly determine locations of breaks in the parental genomes where crossovers occur.
                maternal_crossover_map = rng.integers(S["Genome Size (G)"], size=crosses)
                paternal_crossover_map = rng.integers(S["Genome Size (G)"], size=crosses)

                # locations of crossover breaks are appended to the list of locations where chromosomes begin and end
                maternal_crossover_map = np.sort(np.concatenate((G.GENOME_MAP, maternal_crossover_map)).astype(int))
                paternal_crossover_map = np.sort(np.concatenate((G.GENOME_MAP, paternal_crossover_map)).astype(int))


                # loop once for each segregating genome segment in the person's genome
                for genome_segment in range(segs):

                    # determines which of the two homologous chromosomes produces the segment that enters the egg
                    mother_chromosome = segregation_array[person_id][genome_segment][0]

                    # retrieving the gene loci where the current genome_segment starts and ends
                    # as determined by random recombination of the mother's genes. 
                    m_start = maternal_crossover_map[genome_segment]
                    m_end = maternal_crossover_map[genome_segment + 1]

                    # updating list of inherited genes with genes from the new person's mother
                    self.mutations[person_id, m_start:m_end, 0] = maternal_mutations[mother_id, m_start:m_end, mother_chromosome]

                    # determines which of the two homologous chromosomes produces the segment that enters the sperm
                    father_chromosome = segregation_array[person_id][genome_segment][1]

                    # retrieving the gene loci where the current genome_segment starts and ends
                    # as determined by random recombination of the father's genes. 
                    p_start = paternal_crossover_map[genome_segment]
                    p_end = paternal_crossover_map[genome_segment + 1]

                    # updating list of genes with genes from the new person's father
                    self.mutations[person_id, p_start:p_end, 1] = paternal_mutations[father_id, p_start:p_end, father_chromosome]

            # calculate odds that each gene site in the paternal genome will have a deleterious mutation.
            mu = S["Deleterious Mutation Rate (U)"]/S["Genome Size (G)"]

            # generate a random number for each gene to compare to the odds of mutation calculated above
            rarray = rng.random(size=(self.person_count, S["Genome Size (G)"]))

            # if the random number generated is below a threshold, a new random deleterious allele is generated.
            # the deleterious allele is represented as a 1. 0 represents non-deleteriuos alleles.
            # in the current model, all deleterious alleles are identical with the same effect, although this would be easy to change:
            # simply replace the "1" below with some function to produce a range of values. 
            # to represent additional alleles that could be behave differently or be tracked for whatever reason.
            new_mutations = rarray < mu

            # de-novo mutations and mutations inherited from parents are combined
            # if a site that already has an inherited mutation mutates further, the higher numbered allele is used.
            # in this model, the de-novo mutations always appear from the paternal lineage. Not that it makes much difference. 
            self.mutations[:,:,1] = np.maximum(self.mutations[:,:,1], new_mutations)

       
        # randomly alter the value of the neutral loci that were inherited 
        # The variance in this number is used to calculate effective population size
        self.neutral_locus = self.neutral_locus + rng.normal(0, 1, (self.person_count, 2))

        # make lists of which alleles are lower or higher in dominance priority, based on allele number.
        low_variants = np.minimum(self.mutations[:, :, 0], self.mutations[:, :, 1])
        high_variants = np.maximum(self.mutations[:, :, 0], self.mutations[:, :, 1])

        # calculate the fitness penalties for homozygous recessive mutations
        hoz_surv_pens = np.where(low_variants > 0, G.GENE_HOZ_SURV_PEN_RATIOS, 0)
        hoz_attr_pens = np.where(low_variants > 0, G.GENE_HOZ_COMP_PEN_RATIOS, 0)


        
        # calculate the fitness penalties for heterozygous mutations
        heteroz_surv_pens = np.where(high_variants > 0, np.where(low_variants == 0, (1-G.FULLY_RECESSIVE) * G.GENE_SURV_PEN_RATIOS, 0),0)
        heteroz_attr_pens = np.where(high_variants > 0, np.where(low_variants == 0, (1-G.FULLY_RECESSIVE) * G.GENE_COMP_PEN_RATIOS, 0),0)

        # adds together the fitness penalties of the dominant and recessive mutations.
        self.surv_pens = np.sum(heteroz_surv_pens, axis=1) + np.sum(hoz_surv_pens, axis=1)

        # calculate odds of survival based on:
        # the selective effect of the average gene (s) in the current simulation.
        # the non-selective death chance (Q) 
        # and the total fitness penalties of the person: the surv_pens, calculated above based on the person's genetics. 
        self.survival_fitness = np.power(1-S["Hard Selection Coefficient (s)"], self.surv_pens)
        self.survival_chances = (1 - S["Non-Selective Death Chance (Q)"]) * self.survival_fitness


        # do the same for the competitive fitness metric.
        self.attr_pens = np.sum(heteroz_attr_pens, axis=1) + np.sum(hoz_attr_pens, axis=1)
        self.competition_chances = np.power(1-S["Soft Selection Coefficient (c)"], self.attr_pens)


        # ordered list used to weight the chances of each person having offspring, based on their competition_chances and whether they are eligible
        # these lists are created and maintained by the generate_odds_lists function. they are used to quickly pick a person at random based on their weighting.
        # there are two lists, one for all people, used to pick who will have the next child,
        # and another for only unmarried people, used to pick the new spouse of the person picked from the first list.
        self.odds_list = []
        self.unmarried_odds_list = []

        # the above lists are indexed differently from the main population object because not everyone is eligible.
        # these indexes link the selection in the odds_list to the location of the person in the population object list. 
        self.odds_list_index = []
        self.unmarried_odds_list_index = []

        # stores data for each new person to represent whether they are still alive, who they have married, and how many children they have.
        # for now these are set to default values, because no one from this new generation has died, married, or had children of their own yet.
        self.alive_flags = np.ones(self.person_count)
        self.spouse_ids = (np.zeros(self.person_count) - 1).tolist()
        self.child_count = np.zeros(self.person_count)

        # the total number of married couples that could be produced by this generation, based on the number of surviving women and men. 
        # Tracked so the matchmaking algorithm doesn't try to find a spouse for someone when there is no one left to marry.
        self.marriage_limit = None

        # the total number of people in this population that are married. Used to check against the marriage limit.
        self.marriage_count = 0

    def chance_of_dying(self, rng):

        # flag people as dead if a random number was generated that was lower than their chance of dying.
        self.alive_flags = np.where(rng.random(self.person_count) < self.survival_chances, 1, 0)

    def pick_first_person(self, S, rng):

        return self.pick_person(S, rng, self.odds_list, self.odds_list_index, "first")

    def pick_second_person(self, S, rng, spouse_id, opposite_sexed_people):

        person_id = -1

        # find the person who is the spouse of the picked person, if anyone is. 
        if spouse_id in self.spouse_ids:
            person_id = self.spouse_ids.index(spouse_id)

        # otherwise find a spouse for them from the list of eligible people. 
        else:

            person_id = self.pick_person(S, rng, self.unmarried_odds_list, self.unmarried_odds_list_index, "second")

            # document the new marriage if marriage is enabled
            if person_id != -1 and S["Marriage Flag (M)"] == 1:
                opposite_sexed_people.spouse_ids[spouse_id] = person_id
                self.spouse_ids[person_id] = spouse_id

                self.marriage_count += 1
                opposite_sexed_people.marriage_count += 1

                # check if we have reached the limit for the number of possible marriages--if so, regenerate the odds lists. 
                if self.marriage_count >= self.marriage_limit:
                    self.generate_odds_lists(S)
                    opposite_sexed_people.generate_odds_lists(S)

        return person_id
    
    def pick_person(self, S, rng, odds_list, odds_list_index, order):

        person_id = -1

        failed_attempts = 0
        if len(odds_list) > 0:
            while person_id == -1:

                # generate a random number between 0 and the maximum value of the list. 
                r = rng.random()*odds_list[-1]

                # pick the person from the odds_list whose range in the weighted list is hit by the above random number
                index = np.searchsorted(odds_list, r)

                # lookup the person in the odds_list_index to find them in the population object. 
                person_id = odds_list_index[index]

                # check if they are over the legal child limit or if this is a second pick that is already married. 
                if self.child_count[person_id] >= S["Legal Child Limit (Z)"] or (self.spouse_ids[person_id] != -1 and order == "second"):
                    person_id = -1
                    failed_attempts += 1

                    # regenerate the lists from scratch if it takes too many tries to pick a valid person
                    if failed_attempts >= 5:

                        failed_attempts = 0

                        # regenerate the odds list if too many failed attempts. 
                        odds_list, odds_list_index = self.generate_odds_list(S, order)
                        
                        # if the regenerated odds list has no people in it, stop looping.
                        if len(odds_list) == 0:
                            break

        return person_id
    
    def generate_odds_lists(self, S):
        self.odds_list, self.odds_list_index = self.generate_odds_list(S, "first")
        self.unmarried_odds_list, self.unmarried_odds_list_index = self.generate_odds_list(S, "second")

    def generate_odds_list(self, S, order):

        total_attr = 0

        odds_list = []
        odds_list_index = []

        # for each person
        for i in range(self.person_count):

            #if the person is alive
            if self.alive_flags[i] == 1:

                # if the person is not already a parent of too many children 
                if (self.child_count[i] < S["Legal Child Limit (Z)"]):

                    # for the first list, picks everyone of the sex, or only married people if there are no new possible marriages. 
                    # for the second list, only pick unmarried people. 
                    if (self.spouse_ids[i] != -1 or self.marriage_count < self.marriage_limit) and order == "first" or (self.spouse_ids[i] == -1 and order == "second"):

                        # add them to the list
                        total_attr += self.competition_chances[i]
                        odds_list.append(total_attr)
                        odds_list_index.append(i)
        
        return odds_list, odds_list_index

def assignment_of_children(S, rng, women, men):

    girls_mothers = []
    girls_fathers = []
    boys_mothers = []
    boys_fathers = []

    # determines the maximum number of people who could possibly be married, based on the number of living men and women. 
    # this is needed for the generate_odds_lists functions to make sure it doesn't list unmarried people as options if there is no one left for them to marry.
    women.marriage_limit = np.minimum(np.sum(men.alive_flags),np.sum(women.alive_flags))
    men.marriage_limit = women.marriage_limit

    #creates lists to be used to pick random people weighted by their competition fitness. 
    women.generate_odds_lists(S)
    men.generate_odds_lists(S)

    # create children until the target population is reached.
    while len(girls_mothers) + len(boys_mothers) < S["Population Size (N)"]:
	
        # flip a coin to decide whether the woman or man has the initiative in producing this offspring.
        coin_flip = rng.random() > 0.5

        if coin_flip == True:

            # pick a random woman from the list of women weighted by their competition fitness.
            mother_id = women.pick_first_person(S, rng)
			
            # if the woman is already married, and marriage is enabled
            # this will pick her husband, otherwise, it will choose an unmarried man from the weighted list.
            father_id = men.pick_second_person(S, rng, mother_id, women)

        else:

            # pick a random man from the list of men weighted by their competition fitness. 
            father_id = men.pick_first_person(S, rng)

            # if the man is already married, and marriage is enabled
            # this will pick his wife, otherwise it will choose an unmarried woman from the weighted list. 
            mother_id = women.pick_second_person(S, rng, father_id, men)

        # if there are no more eligible women or men--they have all reached the child count limit--then the loop will break early 
        # in this case, the number of people this generation will be less than the target Population Size (N). 
        if mother_id == -1 or father_id == -1:
            break

        #add this child to the count of children for each parent
        men.child_count[father_id] += 1
        women.child_count[mother_id] += 1


        #testing litters why does this improve selection???
        #for i in range(5):

        #flip a coin to decide the gender of the child
        coin_flip_gender = rng.random() > 0.5

        if coin_flip_gender == True:

            # it's a girl!
            girls_mothers.append(mother_id)
            girls_fathers.append(father_id)

        else:

            # it's a boy!
            boys_mothers.append(mother_id)
            boys_fathers.append(father_id)

    return girls_mothers, girls_fathers, boys_mothers, boys_fathers