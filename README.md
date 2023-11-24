# FamilySizeAndMutationLoad
Code for a Computer Model of Population Genetics used to Estimate the Genetic Effects of Family Size Restriction Laws

model.py: the meat of the population genetics simulator code is here. 
view.py: graphical user interface for windows. 
script_start.py: a simple command line interface for running scripts outside of the UI. Takes as an input parameter a file of the same format as "simulations_to_run.csv".
controller.py: called by view.py or script_start.py, starts the model. Can handle multiple simulations running in parallel by using multiprocessing.
general_toolbox.py: a couple of minor functions used by view.py.
