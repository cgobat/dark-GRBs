"""
 Title: "Automating the Graphing of Beta_OX vs. Beta_X"

 Copyright (C) 2020 David Fitzpatrick

 From: "Analyzing Optically-Dark Short Gamma Ray Bursts"
 (1) David Fitzpatrick, (2) Professor Alexander van der Horst, Ph.D.

 1. Georgetown University, Department of Physics, 37 and O Streets NW, Washington D.C. 20057
 2. The George Washington University, Department of Physics, 725 21 Street NW, Washington
    D.C. 20052

  I hereby grant to Georgetown University and its agents the non-exclusive, worldwide right
  to reproduce, distribute, display and transmit my thesis in such tangible and electronic
  formats as may be in existence now or developed in the future. I retain all ownership
  rights to the copyright of the thesis including the right to use it in whole or in part
  in future works. I agree to allow the Georgetown University Department of Physics to
  serve as the institutional repository of my thesis and to make it available to the
  Georgetown University community through its website. I certify that the version that I
  have submitted is the same version that was approved by my senior research advisor.

 Description: Following the calculation of Beta_OX from the C++ program titled
 "Automating the Calculation of Beta_OX", this program reads in a file containing
 GRB IDs with paired X-ray and optical data pairs and ensuing parameters (including
 Beta_OX) and creates graphs of Beta_OX vs. Beta_X in order to identify optically
 dark GRBs by either the Jakobsson method (Jakobsson et al. 2004) or the Van der
 Horst method (Van der Horst et al. 2009).  The program does this by providing the
 user with a menu of options from which to choose.

"""

import numpy as np, pandas as pd, os, matplotlib, easygui
from matplotlib import rc, pyplot as plt
from scipy import integrate
from pylab import *

# for error bar caps
matplotlib.rcParams.update({'errorbar.capsize': 2})

# create a class of GRBs
class GRB:
    # define constructor for GRB objects
    def __init__(self, ID, dtX, dtO, del_t, BetaX, upper_sigmaX, lower_sigmaX, BetaOX,
                 upper_sigmaOX, lower_sigmaOX, D_Jakobsson, D_vanderHorst):
        self.ID = ID # define attribute for GRB ID
        self.dtX = dtX # define attribute for X-Ray temporal extent
        self.dtO = dtO # define attribute for optical temporal extent
        self.del_t = del_t # define attribute for X-Ray and optical temporal separation
        self.BetaX = BetaX # define attribute for X-Ray spectral index
        self.upper_sigmaX = upper_sigmaX # define attribute for upper bound of standard deviation of optical-to-X-Ray spectral index
        self.lower_sigmaX = lower_sigmaX # define attribute for lower bound of standard deviation of optical-to-X-Ray spectral index
        self.BetaOX = BetaOX # define attribute for optical-to-X-Ray spectral index
        self.upper_sigmaOX = upper_sigmaOX # define attribute for upper bound of standard deviation of optical-to-X-Ray spectral index
        self.lower_sigmaOX = lower_sigmaOX # define attribute for lower bound of standard deviation of optical-to-X-Ray spectral index
        self.D_Jakobsson = D_Jakobsson # define attribute for darkness distance according to Jakobsson criteria
        self.D_vanderHorst = D_vanderHorst # define attribute for darkness distance according to van der Horst criteria

# A function which loads in a file containing Terse Beta_OX data for GRBs and
# assigns them to corresponding attributes of a GRB object, all of which are then
# loaded into a list from which a graphing function can pull desired data
def load_file(filename):

    pd_csv = pd.read_csv(filename)

    # initialize list for GRB objects
    GRB_list = []

    # run through list containing lines read from file
    for line in pd_csv.values:
        # split csv and assign to corresponding variables
        ID, dtX, dtO, dt, BetaX, SigmaX_u, SigmaX_l, BetaOX, SigmaOX_u, SigmaOX_l = line
        # create and append GRB object from file contents
        GRB_list.append(GRB(ID, dtX, dtO, dt, BetaX, SigmaX_u, SigmaX_l, BetaOX, SigmaOX_u, SigmaOX_l, "", ""))

    GRB_counter = len(GRB_list)

    easygui.codebox(msg="Total number of lines (or number of GRBs) read from file: "+str(len(pd_csv))+"\nTotal number of GRBs in list: "+str(GRB_counter), # notify user of total number of lines read from file
                    title="File contents",
                    text=pd_csv.to_string()) # display formatted data table

    # return list of GRB objects
    return GRB_list

# a function to create graphs of Beta_OX vs. Beta_X
# @param: list of GRB objects, name of read-in file
def graph(GRB_list, parsed_filename, graph_title, y_or_n_delB, image_name):

    # remove file extension from filename
    parsed_filename = parsed_filename.split("\\")[-1]

    # initialize lists used to graph
    Beta_X_list = []
    Beta_OX_list = []
    Beta_X_upper_list = []
    Beta_X_lower_list = []
    Beta_OX_upper_list = []
    Beta_OX_lower_list = []

    # run through list containing GRB objects
    for z in GRB_list:
        if float(z.BetaOX) != 0:
            # append spectral index parameters to appropriate lists
            # convert from read-in string to float for data manipulation
            Beta_X_list.append(-1 * float(z.BetaX))
            Beta_X_upper_list.append(float(z.upper_sigmaX))
            Beta_X_lower_list.append(float(z.lower_sigmaX))
            Beta_OX_list.append(-1 * float(z.BetaOX))
            Beta_OX_upper_list.append(float(z.upper_sigmaOX) + delta_beta_ox_t)
            Beta_OX_lower_list.append(float(z.lower_sigmaOX) + delta_beta_ox_t)
        
    # create figure
    fig = plt.figure(figsize=(10,9))
    # formulate axes
    ax = fig.add_axes((0.1, 0.4, 0.8, 0.5))

    # set title
    if y_or_n_delB == 'Y':
        title = ax.set_title(graph_title + parsed_filename + " + \u0394\u03B2\u2092\u2093")
    else:
        title = ax.set_title(graph_title + parsed_filename)

    title.set_position([0.5, 1.05])
    # set x-axis
    ax.set_xlabel('\u03B2\u2093')
    ax.set_ylabel('\u03B2\u2092\u2093')
    plt.xticks()
    plt.yticks()

    # create scatter plot of data
    plt.errorbar(Beta_X_list, Beta_OX_list,
                xerr=np.array([Beta_X_lower_list, Beta_X_upper_list]),
                yerr=np.array([Beta_OX_lower_list, Beta_OX_upper_list]),
                fmt='go', ecolor='k', capthick=2)

    # graph Beta_OX = Beta_X
    x = np.linspace(0,10,1000)
    plt.plot(x, x, linestyle=':', color='red', label="\u03B2\u2092\u2093 = \u03B2\u2093")
    # graph Beta_OX = Beta_X - 0.5
    plt.plot(x,x-0.5,linestyle='-.', color='brown', label="\u03B2\u2092\u2093 = \u03B2\u2093 - 0.5")
    # graph Beta_OX = 0.5
    plt.axhline(y=0.5, linestyle='--', color='orange', label="\u03B2\u2092\u2093 = 0.5")

    # set bounds for x-axis
    plt.xlim(0.2, 3)
    # set bounds for y-axis
    plt.ylim(-0.4, 1.4)
    # make legend
    plt.legend()

    # show plot
    plt.show()

    # save file
    if y_or_n_delB == 'Y':
        # create string for graph's filename
        graph_filename = image_name+parsed_filename+"_w_delBeta.png"
    else:
        # create string for graph's filename
        graph_filename = image_name + parsed_filename + ".png"

    # save figure as PNG
    fig.savefig(f"./Required Files/Generated Files (Python)/{graph_filename}", dpi=600)

    return parsed_filename

# a function to determine if a GRB is optically dark
# @param: list of GRB objects, filename without file extension, Y or N to graph data, Y or N to del beta, image name, user defined ID
def determine_dark_Jakobsson(GRB_list, parsed_filename, yes_or_no_graph, yes_or_no_delB, image_name, user_defined_ID):
    
    # initialize list of dark GRBs
    dark_GRBs_list_Jakobsson = []

    # run through list of GRB objects
    for y in GRB_list:
        # calculate Jakobsson distance
        D_Jakobsson = 0.5 + float(y.BetaOX) - float(y.upper_sigmaOX) - delta_beta_ox_t
        #check if Jakobsson distance is positive (yielding dark)
        if D_Jakobsson > 0 and float(y.BetaOX) != 0:
            # set D_Jakobsson for dark GRB
            y.D_Jakobsson = D_Jakobsson
            # append GRB object to list of Jakobsson dark GRBs
            dark_GRBs_list_Jakobsson.append(y)

    # check that there exists a dark GRB according to the Jakobbson method
    if len(dark_GRBs_list_Jakobsson) != 0:
        if yes_or_no_graph == "Y":
            # graph dark bursts only
            graph(dark_GRBs_list_Jakobsson, parsed_filename, "\u03B2\u2092\u2093 vs. \u03B2\u2093 (" + user_defined_ID + "Jak Dark): ", yes_or_no_delB, image_name)

            # print out optically dark bursts
            easygui.codebox(msg="List of Optically-Dark (Jakobsson) GRBs:",
                            title="Optically dark GRBs",
                            text=pd.DataFrame([[l.ID, l.dtX, l.dtO, l.del_t, l.BetaX, l.upper_sigmaX, l.lower_sigmaX, l.BetaOX, l.upper_sigmaOX, l.lower_sigmaOX] for l in dark_GRBs_list_Jakobsson],
                                columns=["ID Number", "\u0394t\u2093 [hr]", "\u0394t\u2092 [hr]", "\u0394t [hr]", "\u03B2\u2093", "\u03C3\u2093_Up", "\u03C3\u2093_Low", "\u03B2\u2092\u2093",
                                "\u03C3\u2092\u2093_Up", "\u03C3\u2092\u2093_Low"]).to_string())

            # write data of optically dark bursts to file
            writer = []
            # write column headers to file
            writer.append(["ID Number", "dt_x [hr]", "dt_o [hr]", "dt [hr]",
                            "Beta_x", "sigma_x_Up", "sigma_x_Low", "Beta_ox",
                            "sigma_ox_Up", "sigma_ox_Low"])
            # run through list of dark GRBs according to the Van der Horst method
            for q in dark_GRBs_list_Jakobsson:
                # write GRB object attributes for each dark GRB object to file
                writer.append([q.ID, q.dtX, q.dtO, q.del_t, q.BetaX, q.upper_sigmaX,
                                q.lower_sigmaX, q.BetaOX, q.upper_sigmaOX, q.lower_sigmaOX])
                                    
            if yes_or_no_delB == 'Y':
                pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + "_w_delBeta.csv",header=False,index=False)
            else:
                pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + ".csv",header=False,index=False)

            print("\nNumber of Jakobsson Dark GRBs from " + parsed_filename + ": ", len(dark_GRBs_list_Jakobsson), "\n")

    # return list of optically dark bursts
    return dark_GRBs_list_Jakobsson

# a function to determine if a GRB is optically dark
# @param: list of GRB objects, filename without file extension, Y or N to graph data, Y or N to del beta, image name, user defined ID
def determine_dark_vanderHorst(GRB_list, parsed_filename, yes_or_no_graph, yes_or_no_delB, image_name, user_defined_ID):
    
    # initialize list of dark GRBs
    dark_GRBs_list_vanderHorst = []

    # run through list of GRB objects
    for y in GRB_list:
        # calculate van der Horst distance
        D_vanderHorst= ( -1*float(y.BetaX) - ( float(y.lower_sigmaX) )  + float(y.BetaOX)
                         - ( float(y.upper_sigmaOX) + delta_beta_ox_t ) - 0.5 )/np.sqrt(2)

        # check if van der Horst distance is positive (yielding dark)
        if -1*float(y.BetaX) - 0.5 > -1*float(y.BetaOX) + ( float(y.upper_sigmaOX) + delta_beta_ox_t ) and -1*float(y.BetaOX) + 0.5 < -1*float(y.BetaX) - float(y.lower_sigmaX) and float(y.BetaOX) != 0:
            # set D_vanderHorst for dark GRB
            y.D_vanderHorst = D_vanderHorst
            # append GRB object to list of van der Horst dark GRBs
            dark_GRBs_list_vanderHorst.append(y)

    # check that there exists a dark GRB according to the Van der Horst method
    if len(dark_GRBs_list_vanderHorst) != 0:

        if yes_or_no_graph == "Y":
            # graph dark bursts only
            graph(dark_GRBs_list_vanderHorst, parsed_filename, "\u03B2\u2092\u2093 vs. \u03B2\u2093 (" + user_defined_ID + "VdH Dark): ", yes_or_no_delB, image_name)

            # print out optically dark bursts
            easygui.codebox(msg="List of Optically-Dark GRBs:",
                            title="Optically dark GRBs",
                            text=pd.DataFrame([[l.ID, l.dtX, l.dtO, l.del_t, l.BetaX, l.upper_sigmaX, l.lower_sigmaX, l.BetaOX, l.upper_sigmaOX, l.lower_sigmaOX] for l in dark_GRBs_list_vanderHorst],
                                columns=["ID Number", "\u0394t\u2093 [hr]", "\u0394t\u2092 [hr]", "\u0394t [hr]", "\u03B2\u2093", "\u03C3\u2093_Up", "\u03C3\u2093_Low", "\u03B2\u2092\u2093",
                                "\u03C3\u2092\u2093_Up", "\u03C3\u2092\u2093_Low"]).to_string())

            # write data of optically dark bursts to file
            writer = []
            # write column headers to file
            writer.append(["ID Number", "dt_x [hr]", "dt_o [hr]", "dt [hr]",
                            "Beta_x", "sigma_x_Up", "sigma_x_Low", "Beta_ox",
                            "sigma_ox_Up", "sigma_ox_Low"])
            # run through list of dark GRBs according to the van der Horst method
            for q in dark_GRBs_list_vanderHorst:
                # write GRB object attributes for each dark GRB object to file
                writer.append([q.ID, q.dtX, q.dtO, q.del_t, q.BetaX, q.upper_sigmaX,
                                q.lower_sigmaX, q.BetaOX,q.upper_sigmaOX, q.lower_sigmaOX])

            if yes_or_no_delB == 'Y':
                # write data of optically dark bursts to file
                pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + "_w_delBeta.csv",header=False,index=False)
            else:
                # write data of optically dark bursts to file
                pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + ".csv",header=False,index=False)


            print("\nNumber of Van der Horst Dark GRBs from", parsed_filename, ":", len(dark_GRBs_list_vanderHorst), "\n")

    # return list of optically dark bursts
    return dark_GRBs_list_vanderHorst

# a function to determine optically-darkest GRB per unique GRB ID using the Jakobsson method
# @param: list of dark GRB objects by Jakobsson method, filename without file extension, Y or N to graph data, Y or N to del beta, image name, user defined ID
def determine_darkest_Jakobsson(dark_GRBs_list_Jakobsson, parsed_filename, yes_or_no_graph, yes_or_no_delB, image_name, user_defined_ID):
    # initialize counter
    counter = 0
    # initialize list of darkest GRBs
    darkest_GRBs_Jakobsson = []

    # check that there exists a dark GRB according to the Jakobbson method
    if len(dark_GRBs_list_Jakobsson) != 0:

        # run through dark_GRBs_list
        for z in dark_GRBs_list_Jakobsson:

            # initialize new_ID to be the GRB ID read from list
            new_ID = z.ID
            # initialize new_darkest to be the D read from list
            new_darkest_D = z.D_Jakobsson
            # initialize GRB object corresponding to darkest pairing for its unique GRB ID
            new_darkest_GRB = z

            # check if this is the first data point in dark_GRBs_list
            if counter == 0:
                # assign old_ID to first ID in list
                old_ID = z.ID
                # assign old_darkest_D to first D_Jakobsson in list
                old_darkest_D = z.D_Jakobsson
                # assign old_darkest_GRB to first GRB object in list
                old_darkest_GRB = z
            # check if new ID has been reached
            if new_ID == old_ID:
                # check if this pairing for the same GRB ID has larger D
                # (which would make it "darker")
                if new_darkest_D > old_darkest_D:
                    # re-assign darkest D_Jakobsson
                    old_darkest_D = new_darkest_D
                    # re-assign darkest GRB
                    old_darkest_GRB = new_darkest_GRB
            else:
                darkest_GRBs_Jakobsson.append(old_darkest_GRB)
                # re-assign old_darkest to D of different GRB ID just read in from list
                old_darkest_D = z.D_Jakobsson
                # re-assign darkest GRB to that corresponding to new
                # GRB ID just read in from list
                old_darkest_GRB = z

            # re-assign old_ID to GRB ID of GRB just read in from list
            old_ID = z.ID

            # increment dark counter
            counter+=1

        # add last darkest GRB to list of darkest GRBs
        darkest_GRBs_Jakobsson.append(old_darkest_GRB)

        # check that there exists a dark GRB according to the Jakobbson method
        if len(darkest_GRBs_Jakobsson) != 0:

            if yes_or_no_graph == "Y":
                # graph dark bursts only
                graph(darkest_GRBs_Jakobsson, parsed_filename, "\u03B2\u2092\u2093 vs. \u03B2\u2093 (" + user_defined_ID + "Jak Darkest): ", yes_or_no_delB, image_name)

                easygui.codebox(msg="List of Optically-Darkest GRBs by Jakobsson Criteria:",
                                title="Optically darkest GRBs",
                                text=pd.DataFrame([[l.ID, l.dtX, l.dtO, l.del_t, l.BetaX, l.upper_sigmaX, l.lower_sigmaX, l.BetaOX, l.upper_sigmaOX, l.lower_sigmaOX] for l in darkest_GRBs_Jakobsson],
                                    columns=["ID Number", "\u0394t\u2093 [hr]", "\u0394t\u2092 [hr]", "\u0394t [hr]", "\u03B2\u2093", "\u03C3\u2093_Up", "\u03C3\u2093_Low", "\u03B2\u2092\u2093",
                                    "\u03C3\u2092\u2093_Up", "\u03C3\u2092\u2093_Low"]).to_string())

                easygui.msgbox(msg="Number of Darkest Jakobsson Dark GRBs from " + parsed_filename + ": " + str(len(darkest_GRBs_Jakobsson)), title="Jakobsson results")

                # write data of optically dark bursts to file
                writer = []
                # write column headers to file
                writer.append(["ID Number", "dt_x [hr]", "dt_o [hr]", "dt [hr]",
                                "Beta_x", "sigma_x_Up", "sigma_x_Low", "Beta_ox",
                                "sigma_ox_Up", "sigma_ox_Low"])
                # run through list of dark GRBs according to the Jakobbson method
                for q in darkest_GRBs_Jakobsson:
                    # write GRB object attributes for each dark GRB object to file
                    writer.append([q.ID, q.dtX, q.dtO, q.del_t, q.BetaX, q.upper_sigmaX,
                                    q.lower_sigmaX, q.BetaOX,q.upper_sigmaOX, q.lower_sigmaOX])

                if yes_or_no_delB == 'Y':
                    # write data of optically dark bursts to file
                    pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + "_w_delBeta.csv",header=False,index=False)
                else:
                    # write data of optically dark bursts to file
                    pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + ".csv",header=False,index=False)

    # return list of optically dark bursts
    return darkest_GRBs_Jakobsson

# a function to determine optically-darkest GRB per unique GRB ID using the Van der Horst method
# @param: list of VdH dark GRB objects, filename without file extension, Y or N to graph data, Y or N to del beta, image name, user defined ID
def determine_darkest_vanderHorst(dark_GRBs_list_vanderHorst, parsed_filename, yes_or_no_graph, yes_or_no_delB, image_name, user_defined_ID):
    # initialize counter
    counter = 0
    # initialize list of darkest GRBs
    darkest_GRBs_vanderHorst = []

    # check that there exists a dark GRB according to the Van der Horst method
    if len(dark_GRBs_list_vanderHorst) != 0:

        # run through dark_GRBs_list
        for z in dark_GRBs_list_vanderHorst:

            # initialize new_ID to be the GRB ID read from list
            new_ID = z.ID
            # initialize new_darkest to be the D read from list
            new_darkest_D = z.D_vanderHorst
            # initialize GRB object corresponding to darkest pairing for its unique GRB ID
            new_darkest_GRB = z

            # check if this is the first data point in dark_GRBs_list
            if counter == 0:
                # assign old_ID to first ID in list
                old_ID = z.ID
                # assign old_darkest_D to first D_vanderHorst in list
                old_darkest_D = z.D_vanderHorst
                # assign old_darkest_GRB to first GRB object in list
                old_darkest_GRB = z
            # check if new ID has been reached
            if new_ID == old_ID:
                # check if this pairing for the same GRB ID has larger D
                # (which would make it "darker")
                if new_darkest_D > old_darkest_D:
                    # re-assign darkest D_vanderHorst
                    old_darkest_D = new_darkest_D
                    # re-assign darkest GRB
                    old_darkest_GRB = new_darkest_GRB
            else:
                darkest_GRBs_vanderHorst.append(old_darkest_GRB)
                # re-assign old_darkest to D of different GRB ID just read in from list
                old_darkest_D = z.D_vanderHorst
                # re-assign darkest GRB to that corresponding to new GRB ID
                # just read in from list
                old_darkest_GRB = z

            # re-assign old_ID to GRB ID of GRB just read in from list
            old_ID = z.ID

            # increment dark counter
            counter += 1

        # add last darkest GRB to list of darkest GRBs
        darkest_GRBs_vanderHorst.append(old_darkest_GRB)

        # check that there exists a dark GRB according to the Van der Horst method
        if len(darkest_GRBs_vanderHorst) != 0:

            if yes_or_no_graph == "Y":
                # graph dark bursts only
                graph(darkest_GRBs_vanderHorst, parsed_filename, "\u03B2\u2092\u2093 vs. \u03B2\u2093 (" + user_defined_ID + "VdH Darkest): ", yes_or_no_delB, image_name)

                easygui.codebox(msg=f"Number of Darkest Van der Horst Dark GRBs from {parsed_filename}: {len(darkest_GRBs_vanderHorst)}\nList of Optically-Darkest GRBs by Van der Horst Criteria:",
                                title="Optically darkest GRBs",
                                text=pd.DataFrame([[l.ID, l.dtX, l.dtO, l.del_t, l.BetaX, l.upper_sigmaX, l.lower_sigmaX, l.BetaOX, l.upper_sigmaOX, l.lower_sigmaOX] for l in darkest_GRBs_vanderHorst],
                                    columns=["ID Number", "\u0394t\u2093 [hr]", "\u0394t\u2092 [hr]", "\u0394t [hr]", "\u03B2\u2093", "\u03C3\u2093_Up", "\u03C3\u2093_Low", "\u03B2\u2092\u2093",
                                    "\u03C3\u2092\u2093_Up", "\u03C3\u2092\u2093_Low"]).to_string())
                
                # write data of optically dark bursts to file
                writer = []
                # write column headers to file
                writer.append(["ID Number", "dt_x [hr]", "dt_o [hr]", "dt [hr]",
                                "Beta_x", "sigma_x_Up", "sigma_x_Low", "Beta_ox",
                                "sigma_ox_Up", "sigma_ox_Low"])
                # run through list of dark GRBs according to the Van der Horst method
                for q in darkest_GRBs_vanderHorst:
                    # write GRB object attributes for each dark GRB object to file
                    writer.append([q.ID, q.dtX, q.dtO, q.del_t, q.BetaX, q.upper_sigmaX,
                                    q.lower_sigmaX, q.BetaOX, q.upper_sigmaOX, q.lower_sigmaOX])
                                    
                if yes_or_no_delB == 'Y':
                    pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + "_w_delBeta.csv",header=False,index=False)
                else:
                    pd.DataFrame(writer).to_csv(image_name + "dt_" + str(delta_t_beta) + ".csv",header=False,index=False)

        # return list of optically dark bursts
    return darkest_GRBs_vanderHorst

# function to isolate all data points of a particular, user-defined GRB ID
# @params: list of all loaded GRB IDs, parsed filename, user-defined GRB ID
# @returns: list of GRB objects corresponding to user-defined GRB ID
def isolate_ID(GRB_list, parsed_filename, user_defined_ID):

    # create list for GRB objects with user's ID of interest
    user_defined_ID_list = []
    # initialize counter for number of GRBs with user's ID of interest
    interest_counter = 0

    # run through list of GRB objects
    for a in GRB_list:
        # check if ID from full list of GRB objects matches the user's ID of interest
        if user_defined_ID == a.ID:

            # append GRB to list of GRB objects with user's ID of interest
            user_defined_ID_list.append(a)
    
    # check if no data point with user's ID of interest exists
    if len(user_defined_ID_list) > 0:
        # print out GRB object attributes for each GRB object in list of GRB objects with user's ID of interest
        easygui.codebox(msg=f"Found {len(user_defined_ID_list)} data points matching {user_defined_ID}.",
                        title="Selected GRB data",
                        text=pd.DataFrame([[l.ID, l.dtX, l.dtO, l.del_t, l.BetaX, l.upper_sigmaX, l.lower_sigmaX, l.BetaOX, l.upper_sigmaOX, l.lower_sigmaOX] for l in user_defined_ID_list],
                                            columns=["ID Number","\u0394t\u2093 [hr]","\u0394t\u2092 [hr]","\u0394t [hr]","\u03B2\u2093","\u03C3\u2093_Up","\u03C3\u2093_Low","\u03B2\u2092\u2093","\u03C3\u2092\u2093_Up","\u03C3\u2092\u2093_Low"]).to_string())
    else:
        # notify user that no data point exists for ID of interest
        easygui.msgbox("\nThere is no data point with GRB ID " + user_defined_ID + " in file " + parsed_filename + ".csv.  " + "Please try again.\n")
    
    # return list of optically dark bursts
    return user_defined_ID_list

# function to print out user's choices
# @params: Variable designating which menu list to print
# @returns: user's argument
def user_choice(fork_in_the_road,selected_ID):

    # check if user is choosing from main menu
    if fork_in_the_road == "main":
        # print out choices for user
        argument = easygui.choicebox(msg="Please choose from one of the below options.",
                                    title="Mode selection",
                                    choices=["1: Graph all data from loaded file.",
                                            "2: Graph those GRB pairings that are optically-dark according to the Jakobbson method.",
                                            "3: Graph those GRB pairings that are optically-dark according to the Van der Horst method.",
                                            "4: Graph only the GRB pairings that are the darkest for their unique GRB ID according to the Jakobbson method.",
                                            "5: Graph only the GRB pairings that are the darkest for their unique GRB ID according to the Van der Horst method.",
                                            "6: Graph all data points for a particular GRB ID.",
                                            "Q: Quit."])[0]
    # check if user is choosing from submenu
    elif fork_in_the_road == "unique":
        # print out choices for user
        argument = easygui.choicebox(msg="Please choose from one of the subsection choices below for your selected GRB ID.",
                                    title="Graph selection",
                                    choices=[f"A: Graph all data points for GRB {selected_ID}.",
                                            f"B: Graph those data points for GRB {selected_ID} that are optically-dark according to the Jakobbson method.",
                                            f"C: Graph those data points for GRB {selected_ID} that are optically-dark according to the Van der Horst method.",
                                            f"D: Graph only the data point that is the darkest for GRB {selected_ID} according to the Jakobbson method.",
                                            f"E: Graph only the data point that is the darkest for GRB {selected_ID} according to the Van der Horst method.",
                                            "R: Return to the main menu."])[0]

    # return's user's choice
    return argument

# determine execution mode and run main function
if __name__ == '__main__':
    # ask user for name of desired file to open
    # assign user input to variable for name of file
    file_name = easygui.fileopenbox(msg="Select the Paired Data file you wish to load.")

    # remove file extension from filename
    parsed_filename = os.path.splitext(file_name)[0]
    print(parsed_filename,"\n",file_name)
    # call function for loading in file
    GRB_list = load_file(file_name)
    delta_t_beta = int(file_name.split("_")[-2])

    # ask user if they wish to include delta beta due to temporal separation
    del_Beta_Y_N = easygui.ynbox("Include \u0394\u03B2 due to temporal separation?", "\u0394\u03B2 inclusion", ("Yes", "No"))
    if del_Beta_Y_N:
        del_Beta_Y_N = "Y"
    else:
        del_Beta_Y_N = "N"

    # define global variable for delta beta_ox due to temporal separation
    global delta_beta_ox_t
    if del_Beta_Y_N == "Y":
        delta_beta_ox_t = np.log10(1 + (delta_t_beta / 100))
    else:
        delta_beta_ox_t = 0

    # loop menu until user wishes to quit
    while True:
        # assign variable for user's choice
        # call function to display user's options
        argument = user_choice("main",None)

        # check if user wishes to graph all data from the loaded file
        if argument == '1':
            # call function for graphing Beta_OX parameters
            parsed_filename = graph(GRB_list, parsed_filename, "\u03B2\u2092\u2093 vs. \u03B2\u2093: ",del_Beta_Y_N, "Beta_OX_Graph_ALL-")
        # check if user wishes to graph those GRB pairings that are optically-dark according to the Jakobbson method
        elif argument == '2':
            # call function to determine if burst is optically dark using Jakobsson method
            # assign Y or N to graph to Yes
            check_J_dark = len(determine_dark_Jakobsson(GRB_list, parsed_filename, "Y", del_Beta_Y_N, "Jak_Dark-",""))
            # check if there are no dark GRBs according to the Jakobbson method for the user defined ID
            if check_J_dark == 0:
                # notify user
                print("\nThere are no dark GRBs according to the Jakobbson method for GRB " + user_defined_ID + ".")
        # check if user wishes to graph those GRB pairings that are optically-dark according to the Van der Horst method
        elif argument == '3':
            # call function to determine if burst is optically dark using Van der Horst method assign Y or N to graph to Yes
            check_vdH_dark = len(determine_dark_vanderHorst(GRB_list, parsed_filename, "Y", del_Beta_Y_N, "vdH_Dark-", ""))
            # check if there are no dark GRBs according to the Van der Horst method for the user defined ID
            if check_vdH_dark == 0:
                # notify user
                print("\nThere are no dark GRBs according to the Van der Horst method for GRB " + user_defined_ID + ".")
        # Check if user wishes to graph only the GRB pairings that are the darkest for their unique GRB ID according to the Jakobbson method
        elif argument == '4':
            # call function to determine if burst is optically dark using Jakobsson method
            # assign Y or N to graph to No
            dark_GRBs_list_Jakobsson = determine_dark_Jakobsson(GRB_list, parsed_filename, "N", del_Beta_Y_N, "Jak_Dark-", "")
            # call function to determine darkest burst per unique GRB ID by Jakobsson method
            check_J_dark = len(determine_darkest_Jakobsson(dark_GRBs_list_Jakobsson, parsed_filename, "Y", del_Beta_Y_N, "Jak_Darkest-", ""))
            # check if there are no dark GRBs according to the Jakobbson method for the user defined ID
            if check_J_dark == 0:
                # notify user
                print("\nThere are no dark GRBs according to the Jakobbson method for GRB " + user_defined_ID + ".")
        # Check if user wishes to graph only the GRB pairings that are the darkest for their unique GRB ID according to the Van der Horst method
        elif argument == '5':
            # call function to determine if burst is optically dark using Van der Horst method
            # assign Y or N to graph to No
            dark_GRBs_list_vanderHorst = determine_dark_vanderHorst(GRB_list, parsed_filename, "N", del_Beta_Y_N, "vdH_Dark-", "")
            # call function to determine darkest burst per unique GRB ID by the
            # Van der Horst method
            check_vdH_dark = len(determine_darkest_vanderHorst(dark_GRBs_list_vanderHorst, parsed_filename, "Y", del_Beta_Y_N, "vdH_Darkest-", ""))
            # check if there are no dark GRBs according to the Van der Horst method for the user defined ID
            if check_vdH_dark == 0:
                # notify user
                print("\nThere are no dark GRBs according to the Van der Horst method for GRB "+ user_defined_ID + ".")
        # check if user wishes to graph all data points for a particular GRB ID
        elif argument == '6':

            # ask user for their GRB ID of interest
            user_defined_ID = easygui.choicebox(msg="Please select the GRB ID of interest.", title="GRB selection", choices=list(set([grb.ID for grb in GRB_list])))

            # pass variable containing user's GRB of interest to function to create a
            # list of all data points with that unique ID and assign to list
            user_defined_ID_list = isolate_ID(GRB_list, parsed_filename, user_defined_ID)

            # check if no GRB object matches user's ID of interest
            if len(user_defined_ID_list) != 0:

                # loop menu until user wishes to quit
                while True:
                    # assign variable for user's choice
                    # call function to display user's options
                    sub_choice = user_choice("unique", user_defined_ID)
                    # check if user wishes to graph all data points with the GRB ID of interest
                    if sub_choice == 'A':
                        # graph all data points with user's ID of interest
                        graph(user_defined_ID_list, parsed_filename,
                              "\u03B2\u2092\u2093 vs. \u03B2\u2093 (" + user_defined_ID + " All): ", del_Beta_Y_N, "ALL_" + user_defined_ID + "-")
                    # check if user wishes to graph those data points with the GRB ID of interest that are optically-dark according to the Jakobbson method
                    elif sub_choice == 'B':
                        # call function to determine if burst is optically dark using the Jakobsson method
                        # assign Y or N to graph to Yes
                        check_J_dark = len(determine_dark_Jakobsson(user_defined_ID_list, parsed_filename, "Y", del_Beta_Y_N, "Jak_Dark_" + user_defined_ID + "-",
                                            user_defined_ID + " "))
                        # check if there are no dark GRBs according to the Jakobbson method for the user defined ID
                        if check_J_dark == 0:
                            # notify user
                            print("\nThere are no dark GRBs according to the Jakobbson method for GRB " + user_defined_ID + ".")
                    # check if user wishes to graph those data points with the GRB ID of interest that are optically-dark according to the Van der Horst method
                    elif sub_choice == 'C':
                        # call function to determine if burst is optically dark using the Van der Horst method
                        # assign Y or N to graph to Yes
                        check_vdH_dark = len(determine_dark_vanderHorst(user_defined_ID_list, parsed_filename, "Y",del_Beta_Y_N, "vdH_Dark_" + user_defined_ID + "-",
                                                user_defined_ID + " "))
                        # check if there are no dark GRBs according to the Van der Horst method for the user defined ID
                        if check_vdH_dark == 0:
                            # notify user
                            print("\nThere are no dark GRBs according to the Van der Horst method for GRB " + user_defined_ID + ".")
                    # check if user wishes to graph only the data point that is the darkest for the GRB ID of interest according to the Jakobbson method
                    elif sub_choice == 'D':
                        # call function to determine if burst is optically dark using Jakobsson method
                        # assign Y or N to graph to No
                        dark_GRBs_list_J_user_defined = determine_dark_Jakobsson(user_defined_ID_list, parsed_filename, "N", del_Beta_Y_N, "null", "")

                        # call function to determine darkest burst per unique GRB ID by Jakobsson criteria
                        check_J_dark = len(determine_darkest_Jakobsson(dark_GRBs_list_J_user_defined, parsed_filename, "Y", del_Beta_Y_N,
                                            "Jak_Darkest_" + user_defined_ID + "-", user_defined_ID + " ") )
                        # check if there are no dark GRBs according to the Jakobbson method for the user defined ID
                        if check_J_dark == 0:
                            # notify user
                            print("\nThere are no dark GRBs according to the Jakobbson method for GRB " + user_defined_ID + ".")
                    # check if user wishes to graph only the data point that is the darkest for the GRB ID of interest according to the Jakobbson method
                    elif sub_choice == 'E':
                        # call function to determine if burst is optically dark using Van der Horst method
                        # assign Y or N to graph to No
                        dark_GRBs_list_v_user_defined = determine_dark_vanderHorst(user_defined_ID_list, parsed_filename, "N", del_Beta_Y_N, "null", "")
                        # call function to determine darkest burst per unique GRB ID by Van der Horst method
                        check_vdH_dark = len(determine_darkest_vanderHorst(dark_GRBs_list_v_user_defined, parsed_filename, "Y", del_Beta_Y_N,
                                                "vdH_Darkest_" + user_defined_ID + "-", user_defined_ID + " "))
                        # check if there are no dark GRBs according to the Van der Horst method for the user defined ID
                        if check_vdH_dark == 0:
                            # notify user
                            print("\nThere are no dark GRBs according to the Van der Horst method for GRB " + user_defined_ID + ".")
                    # check if user wishes to return to the main menu
                    elif sub_choice == 'R':
                        break
                    # check if user does not type one of the available options
                    else:
                        print("\nYour response did not match one of the available options. Please try again.\n")
        elif argument == 'Q':
            break
        # check if user does not type one of the available options
        else:
            print("\nYour response did not match one of the available options. Please try again.\n")