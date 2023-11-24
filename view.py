
import multiprocessing as mp
import threading as td
import time as time
import os

import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg)
matplotlib.use('TkAgg')
import numpy as np
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from scipy import stats

import controller as controller
import general_toolbox as tlbx


class UserInterface:

    def __init__(self):

        # *** MULTI-PROCESSING CODE ***
        # this section of code relates to variables that can safely be passed between the different processes and threads running simultaneously.
        # there are several threads and processes:
        # #1: The main thread, which holds the user interface -- the current thread.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        # #2: the check_for_updates() thread, which runs continuously and waits for updates from the simulations to update the user interface.
        # #3: Each simulation runs on its own process. 

        # simulations_stopped lets the model know that the program has closed or that the cancel simulations button has been pressed.
        # possible values: 0: not canceled. 1: canceled or completed. 2: canceled and program closed. 
        self.simulations_stopped = mp.Value('i', 1)

        # the message_code corresponds with messages that can be displayed by the UI if something happens during the simulation.
        # higher messages are given priority.
        # 0: nothing to report
        # 1: Ordinary Log Data Reported
        # 2: Simulation Cancellation confirmed
        # 3: Memory Error
        self.message_code = mp.Value('i', 0)

        # this holds a list of simulation names of simulations that are running.
        # each simulation is given a unique simulation_name_id which indexes to the name of the simulation here. 
        # better to do it this way because making thread-safe string parameters is overly complicated. Each simulation knows its name_id, but not its name. 
        # the name_id of each simulation is sent back as log_updates[2] 
        self.simulation_names = []
        self.simulation_name_id = 0

        # used to determine whether any simulations are still running. 
        self.simulations_count = mp.Value('i',0)
        self.completed_simulations = mp.Value('i',0)

        # the log_updates multiprocessing array to hold all the log item data, plus generation_number, total_generations, and simulation_name_id
        # so that the log items can be updated from the model as it is running
        log_updates_length = 3 + len(controller.SimulationLog.get_log_items())
        self.log_updates = mp.Array('d', [0]*log_updates_length)

        # Create the tkinter user interface and name the window
        self.root = tk.Tk()
        self.root.title("Nate's Population Genetics Model")

        # *** MOUSEOVER FRAME ***
        # This area on the top left will display text describing any UI feature that is moused over.
        self.MouseOverFrame = tk.Frame(self.root)
        self.MouseOverFrame.grid(row=1, column=0)
        self.lblMouseOverDescription = tk.Label(self.root, text="Mouse over an element for a description.",
                                                justify="left", width=50, wraplength=250)
        self.lblMouseOverDescription.grid(row=0, column=0, sticky="nw")



        # *** SETTINGS FRAME *** 
        # settings for the model are displayed and can be changed here by the user.
        # IMPORTANT: if you add additional settings to the model, do not add them here -- add them in the list of settings in the controller.py Settings class.
        # Any new settings created in the Settings class will be added to the UI automatically by the code below. 
        self.SettingsFrame = tk.Frame(self.root)
        self.SettingsFrame.grid(row=0, column=1)
        self.SettingsFrame.columnconfigure(0, minsize=100)
        self.SettingsFrame.columnconfigure(1, minsize=200)

        # retrieves a list of all settings and their default values and their tooltips.
        # creates an entry field in the UI for every setting in the list of settings.
        # the list is in the same order as the list of settings.
        self.s = controller.Settings()
        self.SettingFields = []
        for i in range(self.s.count()):
            SettingField = self.UIField(self.SettingsFrame, self.s.keys[i], self.s.values[i], self.s.tooltips[i], self.s.datatypes[i], 
                                        self.lblMouseOverDescription)
            self.SettingFields.append(SettingField)


        # *** CONTROL PANEL Frame ***
        # frame for buttons to operate the application.
        # messages from the running simulations are also displayed here. 
        self.controlPanelFrame = tk.Frame(self.root)
        self.controlPanelFrame.grid(row=0, column=2)
        self.controlPanelFrame.columnconfigure(1, minsize=200)

        # dboxSelect: the user may select data from a previously saved simulation to use as starter data or to display.
        self.dboxSelect = self.OpenFileBox(self, self.controlPanelFrame, self.lblMouseOverDescription, 
                                           "Use Data: \n" 
                                           +"Select data from a previous simulation to graph or to view their logs or to use as a starting point for new simulations ", 
                                           row=0, column=0, )

        # When the Run Simulations button is pressed, the model begins doing its thing, based on the input settings.
        def start_model(): 
            self.simulation_names.append(self.SettingFields[self.s.keys.index("Simulation Name")].get())
            self.simulations_stopped.value = 0
            self.show_msg("Simulation Starting...")
            self.cancel_simulations_button.grid()
            self.simulations_count.value += int(self.tboxReplicationCount.get())
            controller.start_model_from_UI(self)
            self.simulation_name_id +=1
        self.run_simulation_button = tk.Button(self.controlPanelFrame, text="Run Simulations:", command=start_model)
        self.run_simulation_button.grid(row=1, column=0, sticky='e')

        def display_mouseover_run_simulation(event): 
            self.lblMouseOverDescription.config(text="Run Simulations Button and Replication Number Entry: \n"
                                                     "Click this button to start new simulations; their progress will be periodically displayed below. " +
                                                     "Once each simuation completes, a new .npy file will be created in the simulations folder. " +
                                                     "Multiple simulations with the same settings can be run in parallel by increasing the " +
                                                     "number in the box here. Clicking this button again when a simulation is already running " +
                                                     "will start another simulation in parallel.")
        self.run_simulation_button.bind("<Enter>", display_mouseover_run_simulation)

        # creates a text validator to make sure that non-int characters can't be input into int fields. 
        int_validator = self.root.register(tlbx.isint)

        # This textbox allows the user to run multiple replications of simulations simulataneously with the same settings.
        self.tboxReplicationCount = tk.Entry(self.controlPanelFrame, width=10, validate='all', validatecommand=(int_validator,'%P'))
        self.tboxReplicationCount.insert(0,"1")
        self.tboxReplicationCount.grid(row=1, column=1, sticky='w')
        self.tboxReplicationCount.bind("<Enter>", display_mouseover_run_simulation)

        # This label will display the progress of a currently running simulation or the simulation searched in the log.
        self.lblProgress = tk.Label(self.controlPanelFrame, text="")
        self.lblProgress.grid(row=2, column=0, columnspan=2)

        # Button to display the log data for the selected previously-run simulation and the indicated generation.
        # if multiple files are selected, their values will be averaged. 
        self.view_logs_button = tk.Button(self.controlPanelFrame, text="View Logs for Generation #:",
                                          command=self.open_log)
        self.view_logs_button.grid(row=3, column=0, sticky='e')
        def display_mouseover_view_logs(event): 
            self.lblMouseOverDescription.config(text="View Logs Button and Generation Number Entry: \n"
                                                     "Click this button to view data logs for previously run simulations. First, " +
                                                     "select the simulation(s) you wish to view the logs for from the 'Use Data:' box " +
                                                     "next, enter the generation number that you wish to see data for into the box here. " +
                                                     "if multiple simulations were selected at once from the file system using the 'Use Data:' box " +
                                                     "the values reported below will be the average values for all selected simulations.")
        self.view_logs_button.bind("<Enter>", display_mouseover_view_logs)

        # a frame specifically for holding the two log gen number boxes. 
        self.logGenNumberFrame = tk.Frame(self.controlPanelFrame)
        self.logGenNumberFrame.grid(row=3, column=1)

        # When the View Logs button is pressed, it will show data for generation indicated in this textbox. 
        self.tboxLogGenNumber_from = tk.Entry(self.logGenNumberFrame, width=5, validate='all', validatecommand=(int_validator,'%P'))
        self.tboxLogGenNumber_from.insert(0,"1")
        self.tboxLogGenNumber_from.grid(row=0, column=0, sticky='w')
        self.tboxLogGenNumber_from.bind("<Enter>", display_mouseover_view_logs)

        self.lblLogGenNumberTo = tk.Label(self.logGenNumberFrame, text="to")
        self.lblLogGenNumberTo .grid(row=0, column=1)

        self.tboxLogGenNumber_to = tk.Entry(self.logGenNumberFrame, width=5, validate='all', validatecommand=(int_validator,'%P'))
        self.tboxLogGenNumber_to.insert(0,"1")
        self.tboxLogGenNumber_to.grid(row=0, column=2, sticky='w')
        self.tboxLogGenNumber_to.bind("<Enter>", display_mouseover_view_logs)

        # Button clears all data on the graph so that you can create fresh new graphs. 
        def clear_graph():

            self.graphAdditions = []
            self.graphCanvas.get_tk_widget().pack_forget()
            self.clear_graph_button.grid_remove()
        self.clear_graph_button = tk.Button(self.controlPanelFrame, text="Clear Graph Data",
                                            command=clear_graph)
        self.clear_graph_button.grid(row=4, column=0)
        def display_mouseover_clear_graph(event): 
            self.lblMouseOverDescription.config(text="Clear Graph Button: \n"
                                                     "Click this button to clear all data on the graph. Graph data is not automatically " +
                                                     "cleared when you switch the selected files using the 'Use Data:' box. Clicking any of the " +
                                                     "graph buttons next to the log data fields when different files are selected will add more " +
                                                     "and more lines to the graphs until this button is pressed to clear the data.")
        self.clear_graph_button.bind("<Enter>", display_mouseover_clear_graph)
        
        # this button is invisible for now
        self.clear_graph_button.grid_remove()

        # Button cancels any currently running simulations
        def cancel_simulations():
            self.simulations_stopped.value = 1
            self.cancel_simulations_button.grid_remove()
            self.show_msg("Simulation cancelling... please wait...")
        self.cancel_simulations_button = tk.Button(self.controlPanelFrame, text="Cancel Simulations",
                                                   command=cancel_simulations)
        self.cancel_simulations_button.grid(row=4, column=1)
        def display_mouseover_cancel_simulations(event): 
            self.lblMouseOverDescription.config(text="Cancel Simulations Button: \n"
                                                     "Click this to cancel any currently running simulations. Note that this may take a few seconds " +
                                                     "to process as the simulations are running on separate threads. If you click the run simulations button " +
                                                     "too quickly after pressing this, the simulations will not be cancelled.")
        self.cancel_simulations_button.bind("<Enter>", display_mouseover_cancel_simulations)

        # this button is invisible for now
        self.cancel_simulations_button.grid_remove()



        # *** LOG LABELS FRAME ***
        # a frame for reporting data from the model to the user interface.
        # this frame contains information about what is going on in the model, displayed to the user
        # the data could come from either a currently running simulation or the logs from a previously run simulation.
        self.logLabelsFrame = tk.Frame(self.controlPanelFrame)
        self.logLabelsFrame.grid(row=5, column=0, columnspan=2)

        # retrieves a list of statistics that were logged by the model
        self.LogItems = None

        # creates a label for each statistic that was logged by the model and puts it in the UI.
        # also creates a "graph" button for each statistics that will display that statistic to the graph.
        self.LogLabels = []
        for i in range(len(controller.SimulationLog.get_log_items())):
            lbl = self.UILogLabel(self, self.logLabelsFrame, index=i)
            self.LogLabels.append(lbl)

        # generates the all the data tables that will be presented for this simulation in the thesis.
        def clipboard_report():



            data = None
            s = controller.Settings()
            if self.dboxSelect.get() in ("New Simulation", "Multiple Files"):
                self.show_msg("Select an individual file from 'Use Data:' to view logs.")
            else:

                # load currently selected simulation:
                with open(self.dboxSelect.SelectedValues[0], 'rb') as f:
                    np.load(f)
                    np.load(f, allow_pickle=True)
                    s.values = np.load(f)

                file_list = os.listdir("simulations\\")

                base_name = s["Simulation Name"]

                base_name_index = base_name.find("preindustrial")

                if base_name_index != -1:

                    results_array = np.zeros((5,7))
                    replication_count_array = np.zeros(5)
                    stats_list_Z3 = []
                    stats_list_Z16 = []

                    base_name = base_name[0:base_name_index]

                    for file in file_list:
                        
                        file_grid_column = -1

                        if file.startswith(base_name + "preindustrial"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 0

                        elif file.startswith(base_name + "modern_Z16"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 1

                            
                        
                        elif file.startswith(base_name + "modern_Z5"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 2
                            

                        elif file.startswith(base_name + "modern_Z4"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 3
                        
                        elif file.startswith(base_name + "modern_Z3"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 4
                        
                        if file_grid_column >= 0:

                            index_esc = log_list[:,0].tolist().index("Expected Survival Chance")
                            results_array[0,file_grid_column] +=(1 - data[-1][index_esc]) * 100

                            index_asp = log_list[:,0].tolist().index("Average Survival Penalties")
                            results_array[1,file_grid_column] += data[-1][index_asp]

                            index_acw = log_list[:,0].tolist().index("Average Competition Weighting")
                            results_array[2,file_grid_column] += data[-1][index_acw]

                            index_acp = log_list[:,0].tolist().index("Average Competition Penalties")
                            results_array[3,file_grid_column] += data[-1][index_acp]

                            index_dmpp = log_list[:,0].tolist().index("Deleterious Mutations Per Person")
                            results_array[4,file_grid_column] += data[-1][index_dmpp]

                            if file_grid_column == 1:
                                stats_list_Z16.append([(1 - data[-1][index_esc]) * 100, data[-1][index_asp], data[-1][index_acw],
                                                    data[-1][index_acp], data[-1][index_dmpp]])

                            if file_grid_column == 4:
                                stats_list_Z3.append([(1 - data[-1][index_esc]) * 100, data[-1][index_asp], data[-1][index_acw],
                                                    data[-1][index_acp], data[-1][index_dmpp]])

                            replication_count_array[file_grid_column] += 1

                            file_grid_column = -1
                    
                    for i in range(5):
                        for j in range(5):
                            results_array[i,j] = results_array[i,j] / replication_count_array[j]

                    results_array[0,5] = results_array[0,4] / results_array[0,1]
                    
                    for i in range(4):
                        results_array[i+1,5] = (results_array[i+1,4] - results_array[i+1,1]) / results_array[i+1,1] * 100

                    for i in range(5):
                        stats_list_Z3 = np.array(stats_list_Z3)
                        stats_list_Z16 = np.array(stats_list_Z16)
                        results_array[i,6] = stats.ttest_ind(stats_list_Z3[:,i], stats_list_Z16[:,i]).pvalue

                    self.root.clipboard_clear()
                    for i in range(5):
                        for j in range(7):
                            if j == 5 and i > 0:
                                if results_array[i,j] >= 0:
                                    self.root.clipboard_append("+ ")
                                else:
                                    self.root.clipboard_append("- ")
                                    results_array[i,j] = -results_array[i,j]
                            if i == 0 and j == 5:
                                self.root.clipboard_append("* ")
                            if i == 4 and j < 5:
                                self.root.clipboard_append('{:.4g}'.format(results_array[i,j]))
                            elif j == 6:
                                self.root.clipboard_append('{:.2g}'.format(results_array[i,j]))
                            else:
                                self.root.clipboard_append('{:.3g}'.format(results_array[i,j]))
                            if i == 0 and j != 5 and j !=6 or j == 5 and i != 0:
                                self.root.clipboard_append("%")
                            self.root.clipboard_append("\t")
                        self.root.clipboard_append('\n')

            
                base_name_index = base_name.find("Equilibrium_Test")
                if base_name_index != -1:

                    results_array = np.zeros((6,6))
                    replication_count_array = np.zeros(6)
                    k_stats_list = [[],[],[],[],[],[]]
                    w_stats_list = [[],[],[],[],[],[]]

                    for file in file_list:
                    
                        file_grid_column = -1
                        
                        if file.startswith("s25_Equilibrium_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 0

                        
                        if file.startswith("s25_U2_Equilibrium_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 1

                        
                        if file.startswith("s10_Equilibrium_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 2

                        if file.startswith("c25_Equilibrium_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 3

                        
                        if file.startswith("c25_U2_Equilibrium_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 4

                        
                        if file.startswith("c10_Equilibrium_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 5

                        if file_grid_column >= 0:
                            
                            replication_count_array[file_grid_column] += 1
                            
                            if file_grid_column < 3:
                                index_fitness = log_list[:,0].tolist().index("Expected Survival Chance")
                            else:
                                index_fitness = log_list[:,0].tolist().index("Average Competition Weighting")
                            index_mutations = log_list[:,0].tolist().index("Deleterious Mutations Per Person")

                            subtotal_fitness = 0
                            subtotal_mutations = 0

                            for i in range(10):
                                subtotal_fitness += data[-i-1][index_fitness]
                                subtotal_mutations += data[-i-1][index_mutations]
                            
                            results_array[0,file_grid_column] += subtotal_fitness/10
                            results_array[1,file_grid_column] += subtotal_mutations/10

                            w_stats_list[file_grid_column].append(subtotal_fitness/10)
                            k_stats_list[file_grid_column].append(subtotal_mutations/10)

                        file_grid_column = -1

                    for i in range(2):
                        for j in range(6):
                            results_array[i,j] = results_array[i,j] / replication_count_array[j]

                    w_theoretical_values_list = [0.36787944, 0.13533528, 0.36787944, 0.36787944, 0.13533528, 0.36787944 ]
                    k_theoretical_values_list = [4, 8, 10, 4, 8, 10]
                    for i in range(6):
                        w_stats_list[i]= np.array(w_stats_list[i])
                        k_stats_list[i] = np.array(k_stats_list[i])
                        results_array[2,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).pvalue
                        results_array[3,i] = stats.ttest_1samp(k_stats_list[i], k_theoretical_values_list[i]).pvalue
                        results_array[4,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).confidence_interval()[0]
                        results_array[5,i] = stats.ttest_1samp(k_stats_list[i], k_theoretical_values_list[i]).confidence_interval()[0]

                    self.root.clipboard_clear()
                    for i in range(4):
                        for j in range(6):
                            self.root.clipboard_append('{:.4g}'.format(results_array[i,j]))
                            if i in (0,1):
                                self.root.clipboard_append(" +/- " + '{:.4f}'.format(results_array[i,j] - results_array[i+4,j]))
                            self.root.clipboard_append("\t")
                        self.root.clipboard_append('\n')

                    
                base_name_index = base_name.find("Another_Test")
                if base_name_index != -1:

                    results_array = np.zeros((5,5))
                    replication_count_array = np.zeros(5)
                    w_stats_list = [[],[],[],[],[]]

                    for file in file_list:
                    
                        file_grid_column = -1
                        
                        if file.startswith("c25_Y1_Another_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 0

                        
                        if file.startswith("c25_M1_Another_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 1
                        
                        if file.startswith("c25_F1_Another_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 2

                        if file.startswith("c13_s13_Another_Test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 3

                        if file_grid_column >= 0:
                            
                            replication_count_array[file_grid_column] += 1

                            subtotal_fitness = 0
                            subtotal_mutations = 0

                            for i in range(10):
                                if file_grid_column == 3:
                                    index_surv_comp = log_list[:,0].tolist().index("Survivor Competition Weighting")
                                    index_surv_chance = log_list[:,0].tolist().index("Expected Survival Chance")
                                    subtotal_fitness += data[-1-i][index_surv_comp] * data[-1][index_surv_chance]
                                else:
                                    index_fitness = log_list[:,0].tolist().index("Average Competition Weighting")
                                    subtotal_fitness += data[-1-i][index_fitness]
                                
                                index_mutations = log_list[:,0].tolist().index("Deleterious Mutations Per Person")
                                subtotal_mutations += data[-1-i][index_mutations]

                            results_array[0,file_grid_column] += subtotal_fitness/10
                            results_array[1,file_grid_column] += subtotal_mutations/10

                            w_stats_list[file_grid_column].append(subtotal_fitness/10)


                    for i in range(3):
                        for j in range(4):
                            results_array[i,j] = results_array[i,j] / replication_count_array[j]

                    w_theoretical_values_list = [0.36787944, 0.36787944, 0.60653066, 0.36787944]
                    for i in range(4):
                        w_stats_list[i]= np.array(w_stats_list[i])
                        results_array[2,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).pvalue
                        results_array[3,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).confidence_interval()[0]

                    self.root.clipboard_clear()
                    for i in range(3):
                        for j in range(4):
                            self.root.clipboard_append('{:.4g}'.format(results_array[i,j]))
                            if i == 0:
                                self.root.clipboard_append(" +/- " + '{:.4f}'.format(results_array[i,j] - results_array[3,j]))
                            self.root.clipboard_append("\t")
                        self.root.clipboard_append('\n')

                # reporting the data for the simulations attempting to repicate the data from Lesecque et al. (2012)
                base_name_index = base_name.find("Lesecque_repl")
                if base_name_index != -1:

                    results_array = np.zeros([3,3])
                    replication_count_array = np.zeros(3)
                    w_stats_list = [[],[],[]]

                    for file in file_list:

                        file_grid_column = -1

                        if file.startswith("U1_Lesecque_repl"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 0

                        if file.startswith("U10_Lesecque_repl"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 1
                        
                        if file.startswith("U10_Y1_Lesecque_repl"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 2
                            
                        if file_grid_column >= 0:
                            
                            replication_count_array[file_grid_column] += 1

                            index_repr = log_list[:,0].tolist().index("Average Reproduction Chance")
                            subtotal_repr = 0

                            for i in range(10):
                                subtotal_repr += data[-1-i][index_repr]

                            results_array[file_grid_column, 0] += subtotal_repr/10

                            w_stats_list[file_grid_column].append(subtotal_repr/10)
                        
                    for i in range(3):
                        results_array[0,i] = results_array[i,0] / replication_count_array[i]

                    w_theoretical_values_list = [0.837, 0.658, 0.658]
                    for i in range(3):
                        w_stats_list[i]= np.array(w_stats_list[i])
                        results_array[1,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).pvalue
                        results_array[2,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).confidence_interval()[0]
                    
                    self.root.clipboard_clear()



                    for i in range(2):                    
                        for j in range(3):
                            self.root.clipboard_append('{:.4g}'.format(results_array[i,j]))
                            if i == 0:
                                self.root.clipboard_append(" +/- " + '{:.4f}'.format(results_array[i,j] - results_array[2,j]))
                            self.root.clipboard_append("\t")
                        self.root.clipboard_append('\n')

                base_name_index = base_name.find("effective_pop_test")
                if base_name_index != -1:

                    results_array = np.zeros([2,3])
                    replication_count_array = np.zeros(3)
                    w_stats_list = [[],[],[]]

                    for file in file_list:

                        file_grid_column = -1

                        if file.startswith("U0_effective_pop_test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 0
                        
                        if file.startswith("U2.2_effective_pop_test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 1                      
                            
                        if file.startswith("C25_effective_pop_test"):
                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)
                            file_grid_column = 2
                            
                        if file_grid_column >= 0:
                            
                            replication_count_array[file_grid_column] += 1

                            index_Ne = log_list[:,0].tolist().index("Estimated Effective Population")
                            results_array[0,file_grid_column] += data[-1][index_Ne]
                            w_stats_list[file_grid_column].append(data[-1][index_Ne])
                    
                    
                    w_theoretical_values_list = [500, 500, 500]
                    for i in range(3):
                        results_array[0,i] = results_array[0,i] / replication_count_array[i]
                        results_array[1,i] = stats.ttest_1samp(w_stats_list[i], w_theoretical_values_list[i]).confidence_interval()[0]
                    
                    self.root.clipboard_clear()
                    for i in range(3):
                        self.root.clipboard_append('{:.4g}'.format(results_array[0,i]))
                        self.root.clipboard_append(" +/- " + '{:.4f}'.format(results_array[0,i] - results_array[1,i]))
                        self.root.clipboard_append("\t")

                base_name_index = base_name.find("mort_test")
                if base_name_index != -1:

                    results_array = np.zeros(shape=(5,21))
                    replication_count_array = np.zeros(shape=(5,21))

                    for file in file_list:

                        file_grid_column = -1

                        if file.find("mort_test") != -1:


                            with open("simulations\\" + file, 'rb') as f:
                                data = np.load(f)
                                log_list = np.load(f, allow_pickle=True)
                                s.values = np.load(f)

                            for i in range(21):

                                if file.startswith("s15_c"+str(i)+"_mort_test"):
                                    file_grid_column = i
                                    file_row = 0

                                elif file.startswith("s10_c"+str(i)+"_mort_test"):
                                    file_grid_column = i
                                    file_row = 1
                                
                                elif file.startswith("s5_c"+str(i)+"_mort_test"):
                                    file_grid_column = i
                                    file_row = 2
                                
                                elif file.startswith("s1_c"+str(i)+"_mort_test"):
                                    file_grid_column = i
                                    file_row = 3


                                elif file.startswith("s1_c"+str(i)+"_Z3_mort_test"):
                                    file_grid_column = i
                                    file_row = 4

                            
                        if file_grid_column >= 0:
                            
                            replication_count_array[file_row, file_grid_column] += 1

                            index_total_people = log_list[:,0].tolist().index("Total People")
                            if data[-1][index_total_people] == 0:
                                exp_child_mort = float('nan')
                            else:
                                index_mort = log_list[:,0].tolist().index("Expected Child Mortality %")
                                exp_child_mort = data[-1][index_mort]

                            results_array[file_row, file_grid_column] += exp_child_mort
                        
                    for i in range(5):
                        for j in range(21):
                            results_array[i,j] = results_array[i,j] / replication_count_array[i,j] /100
                    
                    self.root.clipboard_clear()
                    for i in range(5):
                        for j in range(21):
                            self.root.clipboard_append('{:.7g}'.format(results_array[i,j]))
                            self.root.clipboard_append("\t")
                        self.root.clipboard_append("\n")

                        


        self.clipboard_report_button = tk.Button(self.controlPanelFrame, text = "Clipboard Report", command=clipboard_report)
        self.clipboard_report_button.grid(row=6, column=0)
        def display_mouseover_make_csv(event):
            self.lblMouseOverDescription.config(text="Clipboard Report: \n"
                                            "select a 'preindustrial' named simulation of the same format as what is in simulatios_to_run.csv file, " +\
                                            "this adds a report of the summarized data of this and all associated 'modern' simulations to the clipboard, " +\
                                            "to be pasted into word to generate the datatables used in the thesis.")
        self.clipboard_report_button.bind("<Enter>", display_mouseover_make_csv)


        # *** GRAPH FRAME ***
        # this frame holds a chart that can display data as line graphs.
        self.graphFrame = tk.Frame()
        self.graphFrame.grid(row=1, columnspan=3)

        # holds the graph object
        self.graphFigure = Figure(figsize=(8, 6), dpi=100)
        self.graphCanvas = FigureCanvasTkAgg(self.graphFigure, self.graphFrame)

        # list of statistics that are currently being displayed on the graph.
        # whenever a "graph" button is clicked, a note to graph its statistic will be added to this list.
        # since the program has just started, this list is currently empty.
        self.graphAdditions = []

        # when the user tries to close the program, open a dialog box to confirm
        def on_closing():
            if messagebox.askokcancel("Close Dialog", "End Program?"):
                self.simulations_stopped.value = 2
                self.root.destroy()
        self.root.protocol("WM_DELETE_WINDOW", on_closing)

        # a continuous loop thread checks the simulations running on different multiprocesses for updates to the UI.
        def check_for_updates():
            # this loop runs indefinitely until the program closes.
            while self.simulations_stopped.value != 2:
                time.sleep(.2)
                if self.simulations_stopped.value !=2 and self.message_code.value == 1:
                    self.message_code.value = 0
                    self.display_log_from_model(self.log_updates)
                    time.sleep(3)
                if self.simulations_stopped.value != 2 and self.message_code.value == 2 and self.simulations_count.value == self.completed_simulations.value:
                    self.show_msg("Simulation cancelled.")
                    self.message_code.value = 0
                    time.sleep(3)
                if self.simulations_stopped.value != 2 and self.message_code.value == 3:
                    self.show_msg("Error: Insufficient memory. At least one simulation aborted.")
                    self.message_code.value = 0
                    time.sleep(3)
                # disable the cancel simulations button if all simulations have finished. 
                # also need to check simulations_stopped (to see if the program is closed) before attempting to grid_remove()
                if self.simulations_stopped.value != 2 and self.simulations_count.value == self.completed_simulations.value:
                    # self.simulations_stopped.value = 1
                    self.cancel_simulations_button.grid_remove()

        # start the thread that communicates with the simulations running on different processes. 
        update_checker = td.Thread(target=check_for_updates)
        update_checker.start()

        # now that the details of the UI are worked out, the mainloop is ready to start. 
        self.root.mainloop()

    # Creates a field in the UI that the user can enter a setting into.
    # Both the tkinter label for the field and its tkinter textbox are created together through this class.
    class UIField:

        # a list of all UIFields that can be cycled through with "tab" or "return" key in the UI.
        FieldList = []

        def __init__(self, loc, lblText, defaultText, mouseover, datatype, mouseoverlabel):

            # gets the ordinal for the current field to be used as the current row and also for cycling between fields
            n = len(self.FieldList)

            # creates the label for the field with tkinter
            lbl = tk.Label(loc, text=lblText, anchor="e")
            lbl.grid(row=n, column=0, sticky="e")

            # creates the textbox for the field with tkinter. Or a checkbox if the field is a bool
            if datatype == bool:
                self.checked_flag = tk.IntVar()
                self.box = tk.Checkbutton(loc, variable=self.checked_flag, onvalue=1, offvalue=0)
            else:
                self.box = tk.Entry(loc, width=20)

            # do not allow the user to input non-numbers into numeric fields. 
            if datatype == float:
                float_validator = loc.register(tlbx.isfloat)
                self.box.config(validate='all', validatecommand=(float_validator,'%P'))
            elif datatype == int:
                int_validator = loc.register(tlbx.isint)
                self.box.config(validate='all', validatecommand=(int_validator, '%P'))

            self.set(defaultText)
            self.box.grid(row=n, column=1, sticky='w')


            # adds the current field to the static list of fields--used for cycling between fields with tab or return.
            self.FieldList.append(self.box)

            # this event is called when the user presses tab or enter--it changes focus to the next field in the list.
            def next_field(event):
                # checks if this is the last field on the UI. if it is, set focus to the first field in the UI.
                if len(self.FieldList) == n + 1:
                    self.FieldList[0].focus_set()
                # otherwise, sets focus to the next field in the list.
                else:
                    self.FieldList[n + 1].focus_set()
                # this return statement prevents the user from putting literal tabs or enters into any of the setting fields.
                return "break"
            self.box.bind("<Return>", next_field)
            self.box.bind("<Tab>", next_field)

            # when the user mouseovers a field, a description appears in lblMouseoverDescription
            def display_mouseover(event): mouseoverlabel.config(text=mouseover)
            if mouseover is not None:
                self.box.bind("<Enter>", display_mouseover)
                lbl.bind("<Enter>", display_mouseover)

        # functions to get or set the data entered into of the textbox or checkbox of this field
        def get(self):
            if type(self.box) == tk.Entry:
                return self.box.get().strip()
            if type(self.box) == tk.Checkbutton:
                return self.checked_flag.get()
        
        def set(self, input):
            if type(self.box) == tk.Entry:
                self.box.delete(0,tk.END)
                self.box.insert(0,input)
            if type(self.box) == tk.Checkbutton:
                if input == '1':
                    self.checked_flag.set(1)
                else:
                    self.checked_flag.set(0)

    
    # a class to create labels that display statistics from the model
    class UILogLabel:

        def __init__(self, UI, loc, index):

            # creates a label to hold the text describing the data
            # when initiated, this label is blank--it does not display until a model is run or "view logs" is pressed
            self.lblText = tk.Label(loc, text="", anchor="e")
            self.lblText.grid(row=index, column=0, sticky="e")

            # creates a label to hold the data itself
            self.lblValue = tk.Label(loc, text="", anchor="e")
            self.lblValue.grid(row=index, column=1, sticky="e")

            # creates the button to add this statistic to the graph
            def add_to_graph(): UI.add_to_graph(index)
            self.graphBtn = tk.Button(loc, text="graph", command=add_to_graph)
            font = tk.font.Font(size=7)
            self.graphBtn['font'] = font
            self.graphBtn.grid(row=index, column=2)

            # the button is only displayed if a simulation is opened from the drop down. 
            self.graphBtn.grid_remove()

            # when the user mouseovers a field, a description appears in lblMouseoverDescription
            def display_mouseover(event): UI.lblMouseOverDescription.config(text=controller.SimulationLog.get_log_tooltip(index))
            self.lblText.bind("<Enter>", display_mouseover)
            self.lblValue.bind("<Enter>", display_mouseover)
            self.graphBtn.bind("<Enter>", display_mouseover)

        # a function to set the text and values of this label
        def set_text(self, lbl, value):
            self.lblText.config(text=lbl)
            self.lblValue.config(text=value)
            self.graphBtn.grid_remove()

    # Creates a widget to select simulation data from the file system or to use new data.
    class OpenFileBox:
        def __init__(self, UI, loc, mouseoverlabel, mouseover, row, column):

            # a label describes the widget to the user.
            lblSelect = tk.Label(loc, text="Use Data: ")
            lblSelect.grid(row=row, column=column, sticky='e')

            # these strings will be displayed in the dropdown box.
            # The user can create a new simulation or retrieve data from an existing one from the file system.
            defaultSelection = "New Simulation"
            openFileString = "Open File..."
            defaultFolder = '\\simulations'

            # the filetype for numpy data files to be retrieved by the widget.
            filetypes = (("numpy data files", "*.npy"),)

            # a string variable holds what is currently displayed in the dropdown box.
            self.dboxS = tk.StringVar()
            self.dboxS.set(defaultSelection)
            self.currentFolder = defaultFolder
            self.SelectedValues = None

            # if the user selects the "Open File..." option, a popup will appear that opens the filesystem.
            def open_file(*args):

                #re-enable all disabled fields in the settings for now
                for settingfield in UI.SettingFields:
                    settingfield.box.config(state='normal')

                # only do anything if "Open File..." is selected from the dropdown box.
                if self.dboxS.get() == openFileString:
                    fullpaths = filedialog.askopenfilenames(initialdir=os.getcwd()+defaultFolder,
                                                  filetypes=filetypes)

                    # sets the string variable that is displayed in the dropdown box to show the selected file's name.
                    if fullpaths == "":
                        self.dboxS.set(defaultSelection)
                    else:
                        self.SelectedValues = fullpaths
                        # if one file was selected, display its filename.
                        if len(fullpaths) == 1:
                            self.dboxS.set(fullpaths[0].split("/")[-1])
                            self.currentFolder = fullpaths[0].split("/")[-2]

                            #open the old file to retrieve the settings of the previous simulation
                            with open(self.currentFolder + "/" + self.dboxS.get(), 'rb') as f:
                                np.load(f)  # don't need the log data right now
                                np.load(f, allow_pickle=True)  # don't need the list of log items right now.
                                previous_sim_setting_values = np.load(f)

                            #for each setting that isn't mutable, disable its entry field and set its value to the same as the previous simulation
                            S = controller.Settings()
                            for i in range(len(previous_sim_setting_values)):
                                if not S.mutable[i]:
                                    UI.SettingFields[i].set(previous_sim_setting_values[i])
                                    UI.SettingFields[i].box.config(state='disabled')

                                

                        # if multiple files were selected, display "Multiple Files"
                        else:
                            self.dboxS.set("Multiple Files")
                        
                        UI.display_graph_buttons()
                        UI.root.update()
                    self.dboxSelect.update()

            # creates a dropdown box to select between two options: A new object or an object from the file system.
            self.dboxSelect = tk.OptionMenu(loc, self.dboxS, *[defaultSelection, openFileString], command=open_file)
            self.dboxSelect.grid(row=row, column=column+1, sticky='w')

            # when the user mouseovers a field, a description appears in lblMouseoverDescription
            def display_mouseover(event): mouseoverlabel.config(text=mouseover)
            if mouseover is not None:
                self.dboxSelect.bind("<Enter>", display_mouseover)
                lblSelect.bind("<Enter>", display_mouseover)

        # returns the filename of what is currently selected or "New Simulation" or "Multiple Files"
        def get(self):
            return self.dboxS.get()
        
        # returns the foldername of what is currently selected. 
        def get_folder(self):
            return self.currentFolder

    # displays the logged data from a previously run simulation for a specified generation
    def open_log(self):

        # opens the log data for the selected simulation file(s).
        log = []
        s = controller.Settings()
        if self.dboxSelect.get() == "New Simulation":
            self.show_msg("Select a file from 'Use Data:' to view logs.")
        else:
            for i in range(len(self.dboxSelect.SelectedValues)):
                with open(self.dboxSelect.SelectedValues[i], 'rb') as f:
                    log.append(np.load(f))
                    log_list = np.load(f, allow_pickle=True)
                    s.values = np.load(f)

            # total number of generations recorded for the selected simulation in the log file.
            number_of_generations = len(log[0])-1

            # retrieves the specified generation from the user-entered Generation Number textbox.
            generation_number_from = int(self.tboxLogGenNumber_from.get())

            if self.tboxLogGenNumber_to.get() == "":
                generation_number_to = None
            else:
                generation_number_to = int(self.tboxLogGenNumber_to.get())

            # for each label in the list of log labels, set its text to display the data from the simulation file
            self.display_log(log, log_list, generation_number_from, number_of_generations, s["Simulation Name"], generation_number_to=generation_number_to )


        # update the user interface with the newly displayed log data.
        self.root.update()

    # instead of displaying logged data in the log statistics section, display an error message.
    def show_msg(self, message):

        # the first log label displays the error instead of whatever it normally displays.
        self.LogLabels[0].set_text(message, "")

        # the rest of the log labels display nothing.
        for i in range(len(self.LogLabels)-1):
            self.LogLabels[i+1].set_text("", "")

        self.lblProgress.config(text="")
        self.root.update()

    # sends log data to the line graph to display.
    def add_to_graph(self, index):

        # first, clear the graph
        self.graphFigure.clear()

        # add the name of the data file that is currently selected to the list of data to be graphed
        for filename in self.dboxSelect.SelectedValues:
            self.graphAdditions.append(filename)

        # this statement enforces the uniqueness of items in the list of items to be graphed.
        # if any item is listed twice, it will only be listed once after this statement.
        self.graphAdditions = list(set(self.graphAdditions))

        # creates the graph plot where the data will be displayed, retrieving its axes
        axes = self.graphFigure.add_subplot()

        # the text associated with the statistic to be graphed.
        # also used as a key to retreive the statistic from the data file.
        logItemLabel = self.LogItems[index][0]

        s = controller.Settings()

        data = []
        settings = []

        for i in range(len(self.graphAdditions)):

            # retrieves data from whatever is selected in the "Use Data: " dropdown box.
            try:
                # opens the log data for the selected simulation file.
                with open(self.graphAdditions[i], 'rb') as f:
                    data.append(np.load(f))
                    np.load(f, allow_pickle=True)
                    settings.append(np.load(f))

            # if it doesn't work, display error message and end function early.
            except FileNotFoundError:
                self.show_msg("Log file '" + self.graphAdditions[i] + "' does not exist.")
                return
        
        setsOfData = []
        for i in range(len(data)):
            newGrouping = True
            
            # TODO: think of a cleaner, less hacky way to do this
            settings[i][0] = "the same"  # Make all simulation names the same before checking if arrays are equal
            settings[i][-1] = "0" #also set the "keep data?" flag to the same. 

            for currentset in setsOfData:
                # to do: figure out how to set the labels here.
                if np.array_equal(currentset[0], settings[i]):
                    currentset.append(data[i])
                    newGrouping = False
                    break
            if newGrouping:
                setsOfData.append([settings[i], data[i]])

        #creates an array, one item for each setting, indicating which settings are distinct. 
        DistinctnessArray = []
        for i in range(len(setsOfData[0][0])):
            distinct = False
            for j in range(len(setsOfData)):
                for k in range(len(setsOfData)-j):
                    if setsOfData[j][0][i] != setsOfData[j+k][0][i] and s.keys[i].find("(") != -1:
                        distinct = True
            DistinctnessArray.append(distinct)

        titleSubheader = ""
        newLineTimer = 1
        for i in range(len(DistinctnessArray)):
            if not DistinctnessArray[i]:
                if titleSubheader != "":
                    titleSubheader = titleSubheader + "     "
                newLineTimer += 1
                if s.keys[i].find("(") != -1:
                    item_symbol = s.keys[i][s.keys[i].find("(")+1]
                    titleSubheader = titleSubheader + str(item_symbol) + ": " + setsOfData[0][0][i]

        sortData = []
        for i in range(len(setsOfData)):
            # total number of generations recorded for the selected simulation in the log file.
            number_of_generations = len(setsOfData[i][1])

            totals = np.zeros(len(setsOfData[i][1]))
            for j in range(len(setsOfData[i])-1):
                asdf = logItemLabel
                jkl = index
                totals = totals + setsOfData[i][j+1][logItemLabel]

            averageddata = totals/(len(setsOfData[i])-1)

            sortData.append(averageddata[number_of_generations-1])

            labelText = "Reps: " + str(len(setsOfData[i])-1)
            for j in range(len(DistinctnessArray)):
                if DistinctnessArray[j]:
                    item_symbol = s.keys[j][s.keys[j].find("(")+1]

                    white_space_length = 9 - len(str(setsOfData[i][0][j]))
                    
                    if white_space_length <= 0:
                        white_space = ""
                    else:
                        white_space = " " * white_space_length

                    labelText = labelText + white_space + str(item_symbol) + ": " + setsOfData[i][0][j]

            # plots the data -- the y-axis is the statistic being graphed and the x-axis is the generation.
            axes.plot(np.arange(number_of_generations), averageddata, label=labelText)

        # get handles and labels
        handles, labels = axes.get_legend_handles_labels()

        # specify order of items in legend
        order = []
        for i in range(len(sortData)):
            highestcurrent = max(sortData)
            index = sortData.index(highestcurrent)
            sortData[index] = float("-inf")
            order.append(index)

        # add legend to plot
        # axes.legend()

        # sets the title of the graph and labels the axes.
        self.graphFigure.suptitle(logItemLabel + ' by Generation', fontsize=16, y=0.99)
        axes.set_title('\n' + titleSubheader, fontdict={'fontsize': 7})
        axes.set_ylabel(logItemLabel)
        axes.set_xlabel("Generation", loc="left")

        labelsInLegend = len(setsOfData)
        box = axes.get_position()
        axes.set_position([box.x0, box.y0 + box.height * 0.03*labelsInLegend,
                         box.width, box.height * (1-0.03*labelsInLegend)])

        # Put a legend below current axis
        axes.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper right', bbox_to_anchor=(1.0, -0.05),
                  fancybox=True, shadow=True)

        # redraws the graph.
        self.graphCanvas.draw_idle()

        # now that there is data in the graph, displays it to the user.
        self.graphCanvas.get_tk_widget().pack()

        # creates toolbar widget for the graph
        # self.graphNavigationToolbar.pack(side=tk.LEFT)

        # enables the clear graph button
        self.clear_graph_button.grid()

        # tells the user interface to display with the new graph that was created.
        self.root.update()

    def display_log_from_model(self, log_updates):
        generation_number = log_updates[0]
        number_of_generations = log_updates[1]
        simulation_name = self.simulation_names[int(self.log_updates[2])]
        log = log_updates[3:]
        self.display_log(log, controller.SimulationLog.get_log_items(), generation_number, number_of_generations, simulation_name, fromModel=True)


    # display data in the logged fields
    def display_log(self, log, log_list, generation_number, number_of_generations, simulation_name, fromModel=False, generation_number_to=None):

        self.LogItems = log_list

        for i in range(len(self.LogItems)):

            # if not from modelfind the average value for the data in the selected log files
            if fromModel==False:
                data = 0
                for j in range(len(log)):
                        if generation_number_to is None: 
                            data += log[j][generation_number][i]
                        else:
                            for k in range(generation_number_to + 1 - generation_number):
                                data += log[j][generation_number+k][i]
                if generation_number_to is None:
                    data = data/len(log)
                else:
                    data = data/len(log)/(generation_number_to +1 - generation_number)
            else:
                data = log[i]

            if self.LogItems[i][1] == float:
                self.LogLabels[i].set_text(self.LogItems[i][0] + ": ", '{:.6g}'.format(data))
            else:
                text = self.LogItems[i][0] + ": ", self.LogItems[i][1](data)
                self.LogLabels[i].set_text(self.LogItems[i][0] + ": ", self.LogItems[i][1](data))

            if fromModel:
                self.LogLabels[i].graphBtn.grid_remove()
            else:
                self.LogLabels[i].graphBtn.grid()

        # update the progress label
        self.lblProgress.config(text=simulation_name + ": Generation " + str(int(generation_number)) + " of " + str(int(number_of_generations)))
        self.root.update()

    def display_graph_buttons(self):
        self.LogItems = controller.SimulationLog.get_log_items()
        for i in range(len(self.LogItems)):
            self.LogLabels[i].set_text(self.LogItems[i][0], "")
            self.LogLabels[i].graphBtn.grid()

if __name__ == '__main__':

    # Creates And Starts The User Interface
    UI = UserInterface()