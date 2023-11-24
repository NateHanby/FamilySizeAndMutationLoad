import controller as controller

import sys
 

if __name__ == '__main__':

    if len(sys.argv) == 1:
        filename = "simulations_to_run.csv"
    else:
        filename = sys.argv[1]

    if len(sys.argv) == 3:
        super_replication = int(sys.argv[2])
    else:
        super_replication = None
    # Creates And Starts The User Interface
    controller.start_model_from_csv(filename, super_replication=super_replication)
