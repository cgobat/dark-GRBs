"""
Title: "Automating the Calculation of Beta_OX"
 
Copyright (C) 2020 David Fitzpatrick
 
From: "Analyzing Optically-Dark Short Gamma Ray Bursts"
(1) David Fitzpatrick, (2) Professor Alexander van der Horst, Ph.D.

1. Georgetown University, of Physics, 37 and O Streets NW, D.C. 20057
2. The George Washington University, of Physics, 725 21 Street NW, D.C. 20052
 
I hereby grant to Georgetown University and its agents the non-exclusive, worldwide
right to reproduce, distribute, and transmit my thesis in such tangible and
electronic formats as may be in existence now or developed in the future. I retain all
ownership rights to the copyright of the thesis including the right to use it in whole
or in part in future works. I agree to allow the Georgetown University Department of
Physics to serve as the institutional repository of my thesis and to make it available
to the Georgetown University community through its website. I certify that the version
that I have submitted is the same version that was approved by my senior research
advisor.

Description: In an effort to automate the calculation of the optical to X-Ray spectral
index (Beta_OX) of well-documented Gamma Ray Bursts (GRBs), program loads in
multiple files containing disparate GRB characteristics (in order: X-ray flux data,
Beta_X data, flux data, optical telescope filter data; fields correspond to
Tables in Fong et al. 2015) and pairs burst measurements based on ID number and
user-defined temporal separation between optical and X-Ray measurements.  The fully-
populated GRBs, parameters which include a calculated value for Beta_OX, then
written to .csv files for further analysis.
"""

#import necessary modules
import numpy as np, pandas as pd
from easygui import fileopenbox, diropenbox

#initialize and declare necessary global constants

#set allowed temporal separation [hr]
DT_PERCENT_DIF = 0.
#set wavelength for X-Rays in nm
#self is based off of an energy of 1 keV and the fact that lambda = hf
FREQUENCY_XRAY = 2.415e+17

#initialize global variable for number of successful pairs after optical
#data is loaded
optical_pairs = 0

#initialize variable for total number of possible pairings; that is, if
#there was 100% accuracy for pairing, would be how many pairs
total_possible_pairings = 0

'''*********************************************************************
                      BEGIN CLASS DEFINITIONS
**********************************************************************'''

class GRB:
    # A constructor that creates a GRB beginning with
    # the GRB ID, X-Ray temporal separation, X-Ray
    # exposure time, X-Ray uncertainty
    def __init__(self, id, dt_X, ExpT_X, Fx, sig_x):
        #create constructor by the following assignment statements
        self.GRB_ID = id
        self.dt_XRay = dt_X
        self.ExpT_XRay = Fx
        self.F_x = Fx
        self.sigma_x = sig_x

        #initialize factors loaded after construction to 0 or None
        self.Beta_X = 31415926535
        self.Beta_X_lower_sigma = 31415926535
        self.Beta_X_upper_sigma = 31415926535
        self.dt_Opt = 0
        self.telescope = None
        self.instrument = None
        self.filter = None
        self.ExpT_Opt = 0
        self.F_o = 0
        self.sigma_o = 0
        self.References_Opt = None
        self.References_XRay = None
        self.frequency_XRay = FREQUENCY_XRAY
        self.Beta_OX = 0
        self.sigma_OX_lower = 0
        self.sigma_OX_upper = 0
        #initialize optical frequency to -1
        #used to check for full population for the Trial print and write function
        self.frequency_Opt = -1

    # Prints attributes of particular GRB neatly.  Items should
    # print in the order listed for private data members below
    def report(self):
        print(self.GRB_ID, self.dt_XRay, self.ExpT_XRay, self.F_x, self.sigma_x, self.Beta_X, self.Beta_X_upper_sigma, self.Beta_X_lower_sigma, self.dt_Opt, self.telescope, self.instrument, self.filter, self.ExpT_Opt, self.F_o, self.sigma_o, self.frequency_XRay, self.frequency_Opt, self.Beta_OX, self.sigma_OX_upper, self.sigma_OX_lower)

    # Accessor and Observers for each field of the X-Ray data.
    def set_Beta_X(self, b):
        #set X-Ray spectral flux density
        self.Beta_X = b
    def set_Beta_X_upper_sigma(self, u):
        #set upper bound of X-Ray spectral flux density
        self.Beta_X_upper_sigma = u
    def set_Beta_X_lower_sigma(self, l):
        #set lower bound of X-Ray spectral flux density
        self.Beta_X_lower_sigma = l

    # Accessor and Observers for each field of the Optical data
    def set_dt_Opt(self, dt):
        #set optical dt of GRB
        self.dt_Opt= dt

    def set_telescope(self, tel):
        #set optical telescope name
        self.telescope = tel

    def set_instrument(self, i):
        #set optical instrument name
        self.instrument = i

    def set_filter(self, f):
        #set filter name
        self.filter = f

    def set_Exp_Opt(self, e):
        #set optical exposure time
        self.ExpT_Opt = e

    def set_F_o(self, f):
        #set optical flux density
        self.F_o = f

    def set_sigma_o(self, s):
        #set optical flux density standard deviation
        self.sigma_o = s

    def set_References_Opt(self, r):
        #set optical references
        self.References_Opt = r

    def set_frequency_XRay(self, f):
        #set X-Ray frequency
        self.frequency_XRay = f

    def set_frequency_Opt(self, wa):
        #set optical frequency
        self.frequency_Opt = wa

    def set_Beta_OX(self, b):
        #set Beta_OX
        self.Beta_OX = b

    def set_sigma_OX_upper(self, s):
        #set upper bound on Beta_OX uncertainty
        self.sigma_OX_upper = s

    def set_sigma_OX_lower(self, s):
        #set lower bound on Beta_OX uncertainty
        self.sigma_OX_lower = s

# A class used for determining total number of possible outcomes
# between X-Ray and optical data
class Possibility:
    # A constructor consisting of the GRB ID number and its
    # corresponding multiplicity
    def __init__(self, id, mult):
        self.ID = id
        self.multiplicity = mult

    #Accessors and observers for the class
    def set_ID(self, id ):
        #set ID
        self.ID = id
    def set_multiplicity(self, m ):
        #set multiplicity
        multiplicity = m


class Trial:
    def __init__(self):
        self.GRBs = [] #array of GRB objects
        self.GRBs_with_Opt = [] #GRBs with optical data

        #define a vector which has number of elements corresponding to the total number
        #of different GRB IDs in the X-Ray file and elements corresponding to the
        #total number of entries per individual GRB ID
        #Ex: GRBs_with_Opt.size() = N; GRBs_With_Opt = {1,5,...,n
        #This would correspond to the X-Ray file having N different GRB IDs,
        #and if the first two and last GRBs had IDs 050509B, 050709, XXXXXX, then
        #GRB 050509B would have 1 entry, 050709 would have 5 entries, GRB XXXXXX
        #would have N entries
        self.XRay_entries = []

        #define a vector identical to XRay_entries but for the optical file
        self.optical_entries = []

        #define a vector of GRBs filled with IDs that exist in optical data but not in X-Ray
        # data
        self.IDs_in_Opt_not_X = []

    # Loads X-Ray data file into the vector of GRBs
    def loadX_RayData(self, filename):
        #initialize variable for number of GRBs loaded
        counter = 0
        #initialize string variable for old ID used in GRB ID multiplicity determination
        old_ID = None
        #initialize counter used for determining number of data points for each unique ID
        entries_per_ID = 0

        # <insert exception code here>

        xrayData = pd.read_csv(str(filename), header=None)

        #display features
        print("\n****************** X-Ray GRB Data ******************")
        print(xrayData.rename(columns=dict(zip(range(5),["GRB ID","dt_X","Exp_X","F_x","sigma_X"]))))

        for line in xrayData.values:
            [ID, dtX, ExpX, Fx, sigmaX] = line
            
            grb = GRB(ID, dtX, ExpX, Fx, sigmaX)

            new_ID = ID

            #check if self is the first entry read from optical file
            if  counter == 0:
                #initialize old_ID to entry read from optical file if
                #self is the first entry from optical file
                old_ID = ID

            #check if ID has been reached
            if  new_ID == old_ID:
                #increment counter for number of entries per unique ID because
                #a multiplicity of the same idea has been found OR it is the first trial
                entries_per_ID += 1

            else:
                #construct Possibility object with corresponding data
                new_Possibility = Possibility(old_ID, entries_per_ID)

                #add in element to vector containing number of entires per unique
                #GRB ID read from optical file
                self.XRay_entries.append(new_Possibility)

                #reset counter used for determining number of data points for each unique ID
                #counter is reset to 1 because in order to reach self else clause, new_ID !=
                #old_ID, there is already a entry
                entries_per_ID = 1


            old_ID = ID
            #increment counter
            counter += 1

            #put the GRB object in the vector of GRBs
            self.GRBs.append(grb)

            #display loaded features
            #print(ID, dtX, ExpX, Fx, sigmaX)

        #add last pair to vector of multiplicities
        #construct Possibility object with corresponding data
        new_Possibility = Possibility(old_ID, entries_per_ID)

        #add in element to vector containing number of entires per unique
        #GRB ID read from optical file
        self.XRay_entries.append(new_Possibility)

        print("\nNumber of GRBs loaded:",counter)

        return counter

    # Loads Beta_X data file. Uses the GRB ID to locate the
    # GRB ID in the vector, then sets the GRB object to
    # also contain the correct Beta_X values
    def loadBeta_X(self, filename):

        #define variable for total loaded Beta X elements
        total_loaded = 0
        #initialize variable for number of successful pairs after Beta_X
        Beta_X_pairs = 0
        #initialize variable for rate of successful pairings
        pairing_rate = 0

        #<test to see if file opens successfully>
        
        BetaXData = pd.read_csv(filename, header=None)
        print("*** Beta_X Data ***")
        print(BetaXData)
        
        #read in data from file and assign them to corresponding variables
        #read until file no longer has any more data
        for line in BetaXData.values:
            [ID, Beta_X, Beta_X_upper_sigma, Beta_X_lower_sigma] = line

            #initialize success boolean to False
            success = False

            #increment counter for successful load
            total_loaded += 1

            # run through vector of GRBs
            for a in range(len(self.GRBs)):
                #match GRB IDs
                if  self.GRBs[a].GRB_ID == ID:
                    #assign attributes to appropriate GRB
                    self.GRBs[a].set_Beta_X(Beta_X)
                    self.GRBs[a].set_Beta_X_upper_sigma(Beta_X_upper_sigma)
                    self.GRBs[a].set_Beta_X_lower_sigma(Beta_X_lower_sigma)

                    #signify match
                    success = True
                    #increment success counter for successful pairing
                    Beta_X_pairs += 1

            #check if no match was found
            if  success == False :
                #notify user that no match was able to be found
                print("\nUnable to match GRB ID",ID,"with Beta_X",Beta_X)

        #calculate percent of loaded GRBs that are paired
        pairing_rate = (Beta_X_pairs /  len(self.GRBs))*100

        #notify user of loading statistics
        print("Number of loaded GRBs from Beta_X file:",total_loaded)
        print("Number of successful pairings with X-Ray Data:",Beta_X_pairs)
        print("Pairing Rate:",pairing_rate,"%")
        
        print("...Cleaning up X-Ray Entries...")
        self.clean_XRay_entries()

        return total_loaded

    # Loads optical data file. Uses the GRB ID to locate the
    # GRB ID in the vector, then sets the GRB object to
    # also contain the correct optical data values
    def loadOpticalData(self, filename):
        #define variables for GRB attributes

        #initialize variable for number of subjects loaded
        total_loaded = 0
        #initialize variable for rate of successful pairings
        pairing_rate = 0
        #initialize counter for no match found
        no_match_found_counter = 0

        optical_pairs = 0
        
        #initialize string variable for old ID used in GRB ID multiplicity determination
        old_ID = None
        #initialize counter used for determining number of data points for each unique ID
        entries_per_ID = 0
        #initialize variable for number of GRBs that are not in
        #BOTH the X-Ray and optical files
        disjoint = 0
        #initialize variable for a check on the number of GRBs that are not in
        #BOTH the X-Ray and optical files
        disjoint_confirm = 0

        #test to see if file opens successfully
        
        opticalData = pd.read_csv(filename, header=None)

        #display features
        print("*** Optical GRB Data ***")
        print(opticalData.rename(columns={0:"GRB ID",1:"dt_O [s]",2:"Telescope",3:"Instrument",4:"Filter",5:"Exposure Time [s]",6:"F_o [uJy]",7:"sigma_O [uJy]"}))

        #read in data from file and assign them to corresponding variables
        #read until file no longer has any more data
        for line in opticalData.values:
            [ID, dtO_hours, tel, inst, fil, ExpO, Fo, sigmaO] = line
            #transform optical dt measurement from hours into seconds
            dtO_seconds = 3600 * dtO_hours

            #display loaded features
            #print(line)

            #initialize variable for location in GRB vector
            location = 0
            #initialize counter for checking if any pairings were made at all
            check = 0

            #initialize new_ID to be the GRB ID read in from file
            new_ID = ID

            #check if self is the first entry read from optical file
            if  total_loaded == 0:
                #initialize old_ID to entry read from optical file if
                #self is the first entry from optical file
                old_ID = ID


            #check if ID has been reached
            if  new_ID == old_ID:
                #increment counter for number of entries per unique ID because
                #a multiplicity of the same idea has been found OR it is the first trial
                entries_per_ID += 1

            else:
                #construct Possibility object with corresponding data
                new_Possibility = Possibility(old_ID, entries_per_ID)

                #add in element to vector containing number of entires per unique
                #GRB ID read from optical file
                self.optical_entries.append(new_Possibility)

                #disjoint += check_ID(old_ID)

                #reset counter used for determining number of data points for each unique ID
                #counter is reset to 1 because in order to reach self else clause, new_ID !=
                #old_ID, there is already a entry
                entries_per_ID = 1

            old_ID = ID

            #increment counter for total loaded Optical elements
            total_loaded += 1

            while location < len(self.GRBs) and location != -1:
                location = self.matchGRB(ID, dtO_seconds, location)
                check += 1

                #print("\nBeginning checks")
                #print("\nLocation is:",location)

                #check if pairing is made
                if  location != -1:
                    #construct a GRB object that will be added
                    #into vector of GRBs with optical data
                    copy_grb = GRB(self.GRBs[location].GRB_ID, self.GRBs[location].dt_XRay, self.GRBs[location].ExpT_XRay, self.GRBs[location].F_x, self.GRBs[location].sigma_x)

                    #set corresponding GRB appropriate Beta_X parameters
                    copy_grb.set_Beta_X(self.GRBs[location].Beta_X)
                    copy_grb.set_Beta_X_lower_sigma(self.GRBs[location].Beta_X_lower_sigma)
                    copy_grb.set_Beta_X_upper_sigma(self.GRBs[location].Beta_X_upper_sigma)

                    #set corresponding GRB appropriate optical parameters
                    copy_grb.set_dt_Opt(dtO_seconds)
                    copy_grb.set_telescope(tel)
                    copy_grb.set_instrument(inst)
                    copy_grb.set_filter(fil)
                    copy_grb.set_Exp_Opt(ExpO)
                    copy_grb.set_F_o(Fo)
                    copy_grb.set_sigma_o(sigmaO)

                    #add newly-created GRB with optical data to
                    #vector of GRBs with optical data
                    self.GRBs_with_Opt.append(copy_grb)

                    #increment location
                    location += 1
                    #increment counter of successful pairing
                    optical_pairs += 1


            #check if no pairing was made
            if  check == 1:
                #notify user that no match was able to be found
                print("Unable to match GRB",ID,"with optical dt",dtO_hours,"[hr] =",dtO_seconds,"[s].")
                
                #increment counter for no match made
                no_match_found_counter += 1

        #add last pair to vector of multiplicities
        #construct Possibility object with corresponding data
        new_Possibility = Possibility(old_ID, entries_per_ID)

        #add in element to vector containing number of entires per unique
        #GRB ID read from optical file
        self.optical_entries.append(new_Possibility)

        self.total_possible_pairings = self.find_total_possible_pairings()

        #calculate percent of loaded GRBs that are paired
        pairing_rate = (optical_pairs /  self.total_possible_pairings )*100

        #notify user of loading statistics
        print("Number of loaded GRBs from optical data file:",int(total_loaded))
        print("Number of possible pairs:",int(self.total_possible_pairings))
        print("Number of successful pairs:",int(optical_pairs))
        print("Pairing Rate:",pairing_rate,"%")

        return total_loaded

    # Loads wavelength data file and pairs each GRB based on
    # telescope and filter information appropriate wavelength
    def loadWavelengthData(self, filename):
        #define variables for GRB attributes

        #initialize variable for number of subjects loaded
        loaded = 0
        #initialize variable for number of successful pairings loaded
        success_counter = 0
        #initialize variable for the counter to check pairing
        check_counter = 0

        #test to see if file opens successfully
        
        observationData = pd.read_csv(filename, header=None)

        #display features
        print("*** Wavelength Data ***")
        print(observationData.rename(columns=dict(zip(range(5),["Telescope","Instrument","Filter","Wavelength","Frequency"]))))

        #read in data from file and assign them to corresponding variables
        #read until file no longer has any more data
        for line in observationData.values:
            [telName, instrumentName, filterName, wavelength, frequency] = line
            #print(line)

            #call function to match GRB objects with appropriate wavelength data
            success_counter += self.matchFrequency(telName, instrumentName, filterName, wavelength, frequency)

            #increment counter for successful load
            loaded += 1

            #notify user that a match was able to be found
            #cout << "Able to match telescope " << telName << " with instrument " <<
            #instrumentName <<  " and with filter " << filterName << "." << endl


        print()
        #display those GRBs that couldn't get paired
        for a in range(len(self.GRBs_with_Opt)):
            #check to see if GRB has been paired with Beta_X, optical, not wavelength data
            if   self.GRBs_with_Opt[a].frequency_Opt == -1:
                print("\nGRB",self.GRBs_with_Opt[a].GRB_ID,"with optical dt",self.GRBs_with_Opt[a].dt_Opt,", telescope",self.GRBs_with_Opt[a].telescope,", instrument",self.GRBs_with_Opt[a].instrument,", with filter",self.GRBs_with_Opt[a].filter,"unpaired.")
                check_counter += 1

        #display loaded statistics
        print("Number of Wavelength Sets Loaded:",loaded)
        print("Number of Unsuccessfully Paired:",check_counter)
        print("Number of Successfully Paired:",success_counter)
        print("Overall Success Rate:",100 * (success_counter / self.total_possible_pairings ),"%")

        return loaded

    # A function used to calculate the optical to X-Ray spectral index
    # as well as the upper and lower bounds on the uncertainty
    def calculate_Beta_OX(self):
        #initialize variable for counter of successful Beta_OX calculations
        success_counter = 0
        #initialize counter for nan calculations of Beta_OX
        nan = 0

        #initialize variables necessary for calculation
        F_x = 0
        F_o = 0
        sigma_x = 0
        sigma_o = 0
        frequency_X = FREQUENCY_XRAY
        frequency_O = 0
        Beta_OX = 0
        sigma_OX_upper = 0
        sigma_OX_lower = 0

        #display features
        print("*** Beta_OX Data ***")
        print("GRB ID, F_x [uJy], sigma_X [uJy], F_o [uJy], sigma_o [uJy], Freq_X, Freq_O, Beta_OX, Upper sigma_OX, Lower sigma_OX")

        # run through vector of GRBs
        for a in range(len(self.GRBs_with_Opt)):
            #check to see if GRB has been fully paired
            if self.GRBs_with_Opt[a].frequency_Opt != -1:

                # Extract values from GRB objects in vector of GRBs
                F_x = self.GRBs_with_Opt[a].F_x
                F_o = self.GRBs_with_Opt[a].F_o
                frequency_O = self.GRBs_with_Opt[a].frequency_Opt
                sigma_x = self.GRBs_with_Opt[a].sigma_x
                sigma_o = self.GRBs_with_Opt[a].sigma_o

                try:
                    #calculate Beta_OX
                    Beta_OX = np.log(F_x / F_o) / np.log(frequency_X / frequency_O )
                except ZeroDivisionError:
                    #check to see if calculation returns NaN
                #if np.isnan(Beta_OX):
                    Beta_OX = 0
                    nan += 1

                try:
                    #calculate upper bound on uncertainty for Beta_OX
                    sigma_OX_upper = np.log( (1 + (sigma_x / F_x)) / (1 - (sigma_o / F_o))) / np.log( frequency_X / frequency_O )
                except ZeroDivisionError:
                    #check to see if calculation returns nan
                #if np.isnan(sigma_OX_upper):
                    sigma_OX_upper = 0

                try:
                    #calculate lower bound on uncertainty for Beta_OX
                    sigma_OX_lower = np.abs(np.log( (1 - (sigma_x / F_x)) / (1 + (sigma_o / F_o))) / np.log( frequency_X / frequency_O ))
                except ZeroDivisionError:
                    #check to see if calculation returns nan
                #if np.isnan(sigma_OX_lower):
                    sigma_OX_lower = 0

                #set calculated value into GRB Beta_OX attribute
                self.GRBs_with_Opt[a].set_Beta_OX(Beta_OX)
                self.GRBs_with_Opt[a].set_sigma_OX_upper(sigma_OX_upper)
                self.GRBs_with_Opt[a].set_sigma_OX_lower(sigma_OX_lower)

                print(self.GRBs_with_Opt[a].GRB_ID,self.GRBs_with_Opt[a].F_x,self.GRBs_with_Opt[a].sigma_x,self.GRBs_with_Opt[a].F_o,self.GRBs_with_Opt[a].sigma_o,self.GRBs_with_Opt[a].frequency_XRay,self.GRBs_with_Opt[a].frequency_Opt,self.GRBs_with_Opt[a].Beta_OX,self.GRBs_with_Opt[a].sigma_OX_upper,self.GRBs_with_Opt[a].sigma_OX_lower)

                #increment success counter for successful calculation
                success_counter += 1


        #display loaded statistics
        print("Number of Successful Beta_OX Calculations:",success_counter)
        print("Overall Success Rate:",100 * (success_counter / self.total_possible_pairings ),"%")

    # prints out the GRB vector by calling the report method
    # of each GRB. Used for debugging.
    def report(self):
        #run through vector of GRBs
        for grb in self.GRBs:
            grb.report()

    #writes paired GRB data to file
    def write_paired_data(self):
        #cast temporal percent difference to string
        dt = int(DT_PERCENT_DIF)
        percent_dif = str(dt)
        print("\nChoose a location to save the paired data tables.")
        savepath = diropenbox()
        #define variable for name of files
        filename_comprehensive = savepath + "\\Comprehensive_Paired_Data_Table_" + percent_dif + "%.csv"
        print("Saving to",filename_comprehensive)
        filename_terse = savepath + "\\GRB_Pairings-dt_" + percent_dif + "%.csv"

        #open file for comprehensive data file
        myFile_comprehensive = open(filename_comprehensive,"w+")

        #headers for file
        filecontent = ["GRB ID,X-Ray dt [hr],X-Ray Exposure Time [s],F_x [uJy],Sigma_x [uJy],Beta_X,Beta_X Upper Sigma,Lower Sigma,Optical dt [hr],Telescope,Instrument,Filter,Optical Exposure Time [s],F_o [uJy],Sigma_o [uJy],Wavelength_X [nm],Frequency_X [Hz],Wavelength_o [nm],Frequency_o [Hz],Beta_OX, Bound of Sigma_OX,Lower Bound of Sigma_OX,\n"]
        
        #run through vector of populated GRBs
        for a in range(len(self.GRBs_with_Opt)):
            #double check to see if GRB is fully populated
            if  self.GRBs_with_Opt[a].frequency_Opt != -1:
                filecontent.append(",".join(map(str, [self.GRBs_with_Opt[a].GRB_ID,self.GRBs_with_Opt[a].dt_XRay/3600,
                                            self.GRBs_with_Opt[a].ExpT_XRay,self.GRBs_with_Opt[a].F_x,self.GRBs_with_Opt[a].sigma_x,self.GRBs_with_Opt[a].Beta_X,
                                            self.GRBs_with_Opt[a].Beta_X_upper_sigma,self.GRBs_with_Opt[a].Beta_X_lower_sigma,self.GRBs_with_Opt[a].dt_Opt/3600,
                                            self.GRBs_with_Opt[a].telescope,self.GRBs_with_Opt[a].instrument,self.GRBs_with_Opt[a].filter,self.GRBs_with_Opt[a].ExpT_Opt,
                                            self.GRBs_with_Opt[a].F_o,self.GRBs_with_Opt[a].sigma_o,self.GRBs_with_Opt[a].frequency_XRay,self.GRBs_with_Opt[a].frequency_Opt,
                                            self.GRBs_with_Opt[a].Beta_OX,self.GRBs_with_Opt[a].sigma_OX_upper,self.GRBs_with_Opt[a].sigma_OX_lower])))
        myFile_comprehensive.write("\n".join(filecontent))

        #close file
        myFile_comprehensive.close()

        #open file for comprehensive data file
        myFile_terse = open(filename_terse,"w+")

        #print out headers for file
        filecontent = ["GRB ID,X-Ray dt [hr],Optical dt [hr],|dt_x - dt_o| [hr],Beta_X,Beta_X Upper Sigma, Lower Sigma,Beta_OX, Bound of Sigma_OX,Lower Bound of Sigma_OX"]

        #run through vector of populated GRBs
        for a in range(len(self.GRBs_with_Opt)):
            #double check to see if GRB is fully populated
            if  self.GRBs_with_Opt[a].frequency_Opt != -1:
                filecontent.append(",".join(map(str, [self.GRBs_with_Opt[a].GRB_ID,self.GRBs_with_Opt[a].dt_XRay/3600,self.GRBs_with_Opt[a].dt_Opt/3600,
                                            np.abs(self.GRBs_with_Opt[a].dt_Opt - self.GRBs_with_Opt[a].dt_XRay)/3600,self.GRBs_with_Opt[a].Beta_X,
                                            self.GRBs_with_Opt[a].Beta_X_upper_sigma,self.GRBs_with_Opt[a].Beta_X_lower_sigma,self.GRBs_with_Opt[a].Beta_OX,
                                            self.GRBs_with_Opt[a].sigma_OX_upper,self.GRBs_with_Opt[a].sigma_OX_lower])))
        myFile_terse.write("\n".join(filecontent))

        #close file
        myFile_terse.close()
    
    # A private function which is used to pair GRBs with their
    # appropriate optical wavelength based on telescope, instrument,
    # and filter information
    def matchFrequency(self, tel, inst, filt, wavelength, frequency):
        #initialize counter for successful pairing
        thatsapair = 0
        #initialize boolean for successful pairing
        success = False

        #run through vector of GRB
        for a in range(len(self.GRBs_with_Opt)):
            #test to see if GRB ID matches passed ID
            #also test to see if GRB has not yet been populated with wavelength parameters
            if  self.GRBs_with_Opt[a].telescope == tel and self.GRBs_with_Opt[a].instrument == inst and self.GRBs_with_Opt[a].filter == filt and self.GRBs_with_Opt[a].frequency_Opt == -1:
                #pair GRB with frequency data
                self.GRBs_with_Opt[a].set_frequency_Opt(frequency)

                #increment success counter
                thatsapair += 1
                #assign boolean to True to signify successful pairing
                success = True

        #return number of successful pairs
        return thatsapair

    # A private function which is used to determine if for a given
    # GRB ID in the optical data set, exists at least one GRB
    # with identical ID read in by the X-Ray Data set
    def check_ID(self, ID):
        #define boolean for match
        corresponding_ID_found = False
        #initialize counter for number of disjoint GRBs
        disjoint_counter = 0

        print("\n\nIn check_ID!")

        for kumquat in range(len(self.XRay_entries)):
            #check if passed optical IDs are identical to those in GRBs that have
            #already been paired with Beta_X data
            if  ID == self.XRay_entries[kumquat].ID:
                #set boolean to True if match found
                corresponding_ID_found = True


        #check if no match was found
        if  corresponding_ID_found == False:
            #add ID that is in optical data set but not X-Ray to
            #vector containing all such IDs
            self.IDs_in_Opt_not_X.append(ID)
            #increment counter for lack of match
            disjoint_counter += 1
            print("\n\nID in Optical but not in X-Ray :(")

        return disjoint_counter
    
    # A private function which is used to determine the total number of
    # possible pairings between optical and X-Ray data by, essence,
    # conducting matrix multiplication between vectors XRay_entries and
    # optical_entries
    def find_total_possible_pairings(self):
        #initialize counters for disjointed IDs
        in_Opt_not_XRay = 0
        in_XRay_not_Opt = 0

        #initialize variables for old and vector sizes
        old_XRay_entries_size = len(self.XRay_entries)
        new_XRay_entries_size = 0
        old_optical_entries_size = len(self.optical_entries)
        new_optical_entries_size = 0

        #initialize variable for total number of possibilities
        total_possibilities = 0

        #traverse through vector XRay_entries
        q = 0
        while q < len(self.XRay_entries):
            #initialize boolean signifying that ID is in XRay_entries
            #and also in optical_entries
            keeper = False

            #traverse through optical_entries vector
            n = 0
            while n < len(self.optical_entries):
                #check if the ID in XRay_entries and optical_entries are identical
                if  self.XRay_entries[q].ID == self.optical_entries[n].ID:
                    #assign boolean for keeper to True
                    keeper = True
                n += 1
            if keeper == False:
                #increment counter for number of IDs that exist in
                #XRay_entries but not in optical_entries
                in_XRay_not_Opt += 1

                #initialize iterator to beginning of XRay_entries vector
                it = 0
                #traverse through XRay_entries with iterator
                while it < len(self.XRay_entries):
                    # remove entity that had no pairing
                    if self.XRay_entries[it].ID == self.XRay_entries[q].ID:
                        # erase() invalidates the iterator, returned iterator
                        erased = self.XRay_entries.pop(it)

                        #check if index is at 0
                        if  q > 0 :
                            #shift index backwards one due to erased entity
                            q -= 1

                        else:
                            #increment iterator
                            #without self the iterator would remain stuck on
                            #index=0 and continuously delete entities
                            it += 1

                    else:
                        #increment iterator
                        it += 1
            q += 1

        #appropriately change variable for size of XRay_entries
        new_XRay_entries_size = len(self.XRay_entries)

        #traverse through vector optical_entries
        l = 0
        while l < len(self.optical_entries):
            #initialize boolean signifying that ID is in XRay_entries
            #and also in optical_entries
            keepme = False

            #traverse through optical_entries vector
            u = 0
            while u < len(self.XRay_entries):
                #check if the ID in XRay_entries and optical_entries are identical
                if  self.optical_entries[l].ID == self.XRay_entries[u].ID:
                    #assign boolean for keeper to True
                    keepme = True
                u += 1

            #check if no match was made
            if keepme == False:
                #increment counter for number of IDs that exist in
                #XRay_entries but not in optical_entries
                in_Opt_not_XRay += 1

                #initialize iterator to beginning of XRay_entries vector
                it_2 = 0
                #traverse through XRay_entries with iterator
                while it_2 < len(self.optical_entries):
                    # remove entity that had no pairing
                    if self.optical_entries[it_2].ID == self.optical_entries[l].ID:
                        # pop() deletes and returns the entry at index it_2
                        erased = self.optical_entries.pop(it_2)

                        #check if index is at 0
                        if  l > 0 :
                            #shift index backwards one due to erased entity
                            l -= 1

                        else:
                            #increment iterator
                            #without self the iterator would remain stuck on
                            #index=0 and continuously delete entities
                            it_2 += 1


                    else:
                        #increment iterator
                        it_2 += 1
            
            l += 1

        #appropriately change variable for size of optical_entries
        new_optical_entries_size = len(self.optical_entries)

        print("Final X-Ray Entries (left) and Optical Entries (right): ")
        print(pd.DataFrame([[self.XRay_entries[y].ID,self.XRay_entries[y].multiplicity,self.optical_entries[y].ID,self.optical_entries[y].multiplicity] for y in range(len(self.XRay_entries))]).rename(columns={0:"GRB ID",1:"Multiplicity (XRay)",2:"GRB ID",3:"Multiplicity (Opt)"}))

        print()

        #display appropriate statistics
        print("Original Number of Unique X-Ray GRBs with Beta_X Data:",old_XRay_entries_size)
        print("Number of Disjoint GRBs Removed from X-Rays:",in_XRay_not_Opt)
        print("Final Number of Unique X-Ray GRBs:",new_XRay_entries_size)
        print("Original Number of Unique Optical Entries:",old_optical_entries_size)
        print("Number of Disjoint GRBs Removed from Optical:",in_Opt_not_XRay)
        print("Final Number of Unique Optical GRBs:",new_optical_entries_size)

        #calculate the total number of possibilities that
        #could be obtained from perfect matching
        for g in range(len(self.optical_entries)):
            #multiply number of options in XRay_entries by
            #number of options in optical_entries per unique
            #GRB ID to determine total number of possibilities
            total_possibilities += self.XRay_entries[g].multiplicity * self.optical_entries[g].multiplicity

        #return number of erased
        return total_possibilities

    # A private function which is used to remove entities that haven't
    # been paired with Beta_X data from vector XRay_entries
    def clean_XRay_entries(self):
        #initialize counter for number we expect to keep
        should_not_keep = 0

        print("Original Number of Unique X-Ray GRBs:",len(self.XRay_entries))

        #traverse through vector XRay_entries
        q = 0
        while q < len(self.XRay_entries):
            keeper = False

            #traverse through vector GRBs that contain only some with BetaX data
            n = 0
            while n < len(self.GRBs):
                #check if the ID in XRay_entries and GRBs are identical and that
                #ID has been paired with BetaX data
                if self.XRay_entries[q].ID == self.GRBs[n].GRB_ID and self.GRBs[n].Beta_X != 31415926535:
                    #assign boolean for keeper to True
                    keeper = True
                n += 1
            #check if no match was made
            if keeper == False:
                #increment counter for number we should not keep
                should_not_keep += 1

                #initialize iterator to beginning of XRay_entries vector
                it = 0
                #traverse through XRay_entries with iterator
                while it < len(self.XRay_entries):
                    # remove entity that had no pairing
                    if self.XRay_entries[it].ID == self.XRay_entries[q].ID:
                        # pop() deletes and returns the entry at index it
                        erased = self.XRay_entries.pop(it)

                        #check if index is at 0
                        if  q > 0 :
                            #shift index backwards one due to erased entity
                            q -= 1

                        else:
                            #increment iterator
                            #without self the iterator would remain stuck on
                            #index=0 and continuously delete entities
                            it += 1

                    else:
                        #increment iterator
                        it += 1
            
            q += 1


        print("Final Number of Unique X-Ray GRBs:",len(self.XRay_entries))
        print("Number Removed for Lack of Beta_X Pairing:",should_not_keep)

    # A private function which is used to pair X-Ray data and
    # optical data to a particular GRB based on a small temporal
    # separation in measurement time
    def matchGRB(self, ID, dtO_s, location):
        #initialize bool variable to determine if match occurs
        paired = False

        #run through vector of GRB
        for a in range(location,len(self.GRBs)):
            #test to see if GRB ID matches passed ID
            #test to see if GRB has been populated with Beta_X data
            if  self.GRBs[a].GRB_ID == ID and self.GRBs[a].Beta_X != 31415926535:
                #test to see if temporal separation is small enough and GRB hasn't already
                #been populated with optical data
                if  (100 * np.abs(self.GRBs[a].dt_XRay - dtO_s) / dtO_s ) < DT_PERCENT_DIF :
                    #assign match to True to signify successful match
                    paired = True
                    #return location of self GRB
                    return a


        #test to see if no match occurred
        if  paired == False :
            #return -1 if no match occurred
            return -1

    # A private function used to return the location
    # in the  vector of GRBs of a GRB with a
    # particular GRB ID used in the loadOpticalData
    # to locate where a particular GRB is
    def findGRB(self, ID):
        #initialize bool variable to determine if match occurs
        match = False

        #run through vector of GRBs
        for a in range(len(GRBs)):
            #test to see if Subject ID matches passed ID
            if  GRBs[a].GRB_ID == ID:
                #assign match to True to signify successful match
                match = True
                #return location of self GRB
                return a

        #test to see if no match occurred
        if  match == False:
            #return -1 if no match occurred
            return -1


'''*********************************************************************
                      BEGIN MAIN FUNC CALLS
*********************************************************************'''

if __name__ == '__main__':
    #create a trial object
    t1 = Trial()

    #ask user for name of file they wish to load for X-Ray data
    print("\nPlease select the X-Ray data file.")
    XRayDataFile_name = fileopenbox()
    print("Selected",XRayDataFile_name)
    
    
    #call function for loading X-Ray data
    #pass name of file for X-Ray data
    t1.loadX_RayData(XRayDataFile_name)

    #formatting
    print("\n")

    #ask user for name of file they wish to load for Beta_X data
    print("Please select the Beta_X data file.")
    Beta_X_File_name = fileopenbox()
    print("Selected",Beta_X_File_name)

    #call function for loading X-Ray data
    #pass name of file for X-Ray data
    t1.loadBeta_X(Beta_X_File_name)

    #formatting
    print("\n")

    #ask user for desired temporal percent difference
    DT_PERCENT_DIF = float(input("Please enter the desired temporal percent difference (%): "))

    #ask user for name of file they wish to load for side effects
    print("Please select the optical data file.")
    OpticalData_name = fileopenbox()
    print("Selected",OpticalData_name)
    
    #call function for loading optical data
    #pass name of file for optical data
    t1.loadOpticalData(OpticalData_name)

    #formatting
    print("\n")

    #ask user for name of file they wish to load for side effects
    print("Please select the wavelength data file.")
    WavelengthData_name = fileopenbox()
    print("Selected",WavelengthData_name)

    #call function for loading wavelength data
    #pass name of file for wavelength data
    t1.loadWavelengthData(WavelengthData_name)

    print() #for formatting purposes
    
    input("Press any key to calculate Beta_OX.")

    #calculate X-Ray to optical spectral indices
    t1.calculate_Beta_OX()

    #write paired data to .csv file
    t1.write_paired_data()

    pass