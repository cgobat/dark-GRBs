 /**
  * Title: "Automating the Calculation of Beta_OX"
  *
  * Copyright (C) 2020 David Fitzpatrick
  *
  * From: "Analyzing Optically-Dark Short Gamma Ray Bursts"
  * (1) David Fitzpatrick, (2) Professor Alexander van der Horst, Ph.D.

  * 1. Georgetown University, Department of Physics,
  *    37 and O Streets NW, Washington D.C. 20057
  * 2. The George Washington University, Department of Physics,
  *    725 21 Street NW, Washington D.C. 20052
  *
  *  I hereby grant to Georgetown University and its agents the non-exclusive, worldwide
  *  right to reproduce, distribute, display and transmit my thesis in such tangible and
  *  electronic formats as may be in existence now or developed in the future. I retain all
  *  ownership rights to the copyright of the thesis including the right to use it in whole
  *  or in part in future works. I agree to allow the Georgetown University Department of
  *  Physics to serve as the institutional repository of my thesis and to make it available
  *  to the Georgetown University community through its website. I certify that the version
  *  that I have submitted is the same version that was approved by my senior research
  *  advisor.
  *
  * Description: In an effort to automate the calculation of the optical to X-Ray spectral
  * index (Beta_OX) of well-documented Gamma Ray Bursts (GRBs), this program loads in
  * multiple files containing disparate GRB characteristics (in order: X-ray flux data,
  * Beta_X data, optical flux data, and optical telescope filter data; fields correspond to
  * Tables in Fong et al. 2015) and pairs burst measurements based on ID number and
  * user-defined temporal separation between optical and X-Ray measurements.  The fully-
  * populated GRBs, with parameters which include a calculated value for Beta_OX, are then
  * written to .csv files for further analysis.
  *
  **/

//declare necessary preprocessor directives
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <ctime>

//initialize and declare necessary global constants

//set allowed temporal separation [hr]
double DT_PERCENT_DIF;
//set wavelength for X-Rays in nm
//this is based off of an energy of 1 keV and the fact that
// \lambda = hv
const double FREQUENCY_XRAY = 2.415e+17;

//initialize global variable for number of successful pairs after optical
//data is loaded
double optical_pairs = 0;

//initialize variable for total number of possible pairings; that is, if
//there was 100% accuracy for pairing, this would be how many pairs
double total_possible_pairings = 0;

using namespace std;

/**********************************************************************
                      BEGIN CLASS DEFINITIONS
***********************************************************************/

class GRB {

    public:
      // A constructor that creates a GRB beginning with
      // the GRB ID, X-Ray temporal separation, X-Ray
      // exposure time, and X-Ray uncertainty
      GRB(string, double, double, double, double);

      // Prints attributes of particular GRB neatly.  Items should
      // print in the order listed for private data members below
      void report();

      // Returns the GRB ID of the GRB
      string getGRB_ID();

      // Accessor and Observers for each field of the X-Ray data.
      double get_dt_XRay();
      double get_ExpT_XRay();
      double get_F_x();
      double get_sigma_x();
      string get_References_XRay();
      void set_Beta_X(double);
      double get_Beta_X();
      void set_Beta_X_upper_sigma(double);
      double get_Beta_X_upper_sigma();
      void set_Beta_X_lower_sigma(double);
      double get_Beta_X_lower_sigma();

      // Accessor and Observers for each field of the Optical data
      void set_dt_Opt(double);
      double get_dt_Opt();
      void set_telescope(string);
      string get_telescope();
      void set_instrument(string);
      string get_instrument();
      void set_filter(string);
      string get_filter();
      void set_Exp_Opt(double);
      double get_Exp_Opt();
      void set_F_o(double);
      double get_F_o();
      void set_sigma_o(double);
      double get_sigma_o();
      void set_References_Opt(string);
      string get_References_Opt();
      void set_frequency_XRay(double);
      double get_frequency_XRay();
      void set_frequency_Opt(double);
      double get_frequency_Opt();
      void set_Beta_OX(double);
      double get_Beta_OX();
      void set_sigma_OX_upper(double);
      double get_sigma_OX_upper();
      void set_sigma_OX_lower(double);
      double get_sigma_OX_lower();

    private:
        //17 attributes of the GRB class
        string GRB_ID;
        double dt_XRay;
        double ExpT_XRay;
        double F_x;
        double sigma_x;
        double Beta_X;
        double Beta_X_upper_sigma;
        double Beta_X_lower_sigma;
        string References_XRay;
        double dt_Opt;
        string telescope;
        string instrument;
        string filter;
        double ExpT_Opt;
        double F_o;
        double sigma_o;
        string References_Opt;
        double frequency_XRay;
        double frequency_Opt;
        double Beta_OX;
        double sigma_OX_upper;
        double sigma_OX_lower;

};

GRB::GRB(string id, double dt_X, double ExpT_X, double Fx, double sig_x)
{
      //create constructor by the following assignment statements
      GRB_ID = id;
      dt_XRay = dt_X;
      ExpT_XRay = F_x;
      F_x = Fx;
      sigma_x = sig_x;

      //initialize factors loaded after construction to 0 or NULL
      Beta_X = 31415926535;
      Beta_X_lower_sigma = 31415926535;
      Beta_X_upper_sigma = 31415926535;
      dt_Opt = 0;
      telescope = "NULL";
      instrument = "NULL";
      filter = "NULL";
      ExpT_Opt = 0;
      F_o = 0;
      sigma_o = 0;
      References_Opt = "NULL";
      References_XRay = "NULL";
      frequency_XRay = FREQUENCY_XRAY;
      Beta_OX = 0;
      sigma_OX_lower = 0;
      sigma_OX_upper = 0;
      //initialize optical frequency to -1
      //used to check for full population for the Trial print and write function
      frequency_Opt = -1;

}

void GRB::report()
{
    //set the numeric output formatting
    cout << fixed << showpoint << setprecision( 2 );
    //print out Subject's attributes from the first-loaded data file
    cout << GRB_ID << " " << dt_XRay << " " << ExpT_XRay << " " << F_x << " "
         << sigma_x << " " << Beta_X << " " << Beta_X_upper_sigma << " "
         << Beta_X_lower_sigma << " " << dt_Opt << " " << telescope << " " << instrument
         << " " << filter << " " << ExpT_Opt << " " << F_o << " " << sigma_o
         << " " << " " << frequency_XRay << " " << " " << frequency_Opt << " " << Beta_OX
         << " " << sigma_OX_upper << " " <<  sigma_OX_lower  <<endl;
}

string GRB::getGRB_ID()
{
    //return GRB ID of GRB
    return GRB_ID;
}

double GRB::get_dt_XRay()
{
    //return X-Ray dt of GRB
    return dt_XRay;
}
double GRB::get_ExpT_XRay()
{
    //return X-Ray Exposure time
    return ExpT_XRay;
}
double GRB::get_F_x()
{
    //return the X-Ray flux density
    return F_x;
}
double GRB::get_sigma_x()
{
    //return the X-Ray sigma value
    return sigma_x;
}
string GRB::get_References_XRay()
{
    //return X-Ray references
    return References_XRay;
}
void GRB::set_Beta_X(double b)
{
    //set X-Ray spectral flux density
    Beta_X = b;
}
double GRB::get_Beta_X()
{
    //get X-Ray spectral flux density
    return Beta_X;
}
void GRB::set_Beta_X_upper_sigma(double u)
{
    //set upper bound of X-Ray spectral flux density
    Beta_X_upper_sigma = u;
}
double GRB::get_Beta_X_upper_sigma()
{
    //return upper bound of X-Ray spectral flux density
    return Beta_X_upper_sigma;
}
void GRB::set_Beta_X_lower_sigma(double l)
{
    //set lower bound of X-Ray spectral flux density
    Beta_X_lower_sigma = l;
}
double GRB::get_Beta_X_lower_sigma()
{
    //return lower bound of X-Ray spectral flux density
    return Beta_X_lower_sigma;
}
void GRB::set_dt_Opt(double dt)
{
    //set optical dt of GRB
    dt_Opt= dt;
}

double GRB::get_dt_Opt()
{
    //return optical dt of GRB
    return dt_Opt;
}

void GRB::set_telescope(string tel)
{
    //set optical telescope name
    telescope = tel;
}

string GRB::get_telescope()
{
    //return optical telescope name
    return telescope;
}

void GRB::set_instrument(string i)
{
    //set optical instrument name
    instrument = i;
}

string GRB::get_instrument()
{
    //return optical instrument name
    return instrument;
}

void GRB::set_filter(string f)
{
    //set filter name
    filter = f;
}

string GRB::get_filter()
{
    //get filter name
    return filter;
}

void GRB::set_Exp_Opt(double e)
{
    //set optical exposure time
    ExpT_Opt = e;
}

double GRB::get_Exp_Opt()
{
    //return optical exposure time
    return ExpT_Opt;
}

void GRB::set_F_o(double f)
{
    //set optical flux density
    F_o = f;
}

double GRB::get_F_o()
{
    //get optical flux density
    return F_o;
}

void GRB::set_sigma_o(double s)
{
    //set optical flux density standard deviation
    sigma_o = s;
}

double GRB::get_sigma_o()
{
    //get optical flux density standard deviation
    return sigma_o;

}

void GRB::set_References_Opt(string r)
{
    //set optical references
    References_Opt = r;
}

string GRB::get_References_Opt()
{
    //get optical references
    return References_Opt;
}

void GRB::set_frequency_XRay(double f)
{
    //set X-Ray frequency
    frequency_XRay = f;
}

double GRB::get_frequency_XRay()
{
    //get X-Ray frequency
    return frequency_XRay;
}

void GRB::set_frequency_Opt(double wa)
{
    //set optical frequency
    frequency_Opt = wa;
}

double GRB::get_frequency_Opt()
{
    //get optical frequency
    return frequency_Opt;
}

void GRB::set_Beta_OX(double b)
{
    //set Beta_OX
    Beta_OX = b;
}

double GRB::get_Beta_OX()
{
    //get Beta_OX
    return Beta_OX;
}

void GRB::set_sigma_OX_upper(double s)
{
    //set upper bound on Beta_OX uncertainty
    sigma_OX_upper = s;
}

double GRB::get_sigma_OX_upper()
{
    //get upper bound on Beta_OX uncertainty
    return sigma_OX_upper;
}

void GRB::set_sigma_OX_lower(double s)
{
    //set lower bound on Beta_OX uncertainty
    sigma_OX_lower = s;
}

double GRB::get_sigma_OX_lower()
{
    //get lower bound on Beta_OX uncertainty
    return sigma_OX_lower;
}

// A class used for determining total number of possible outcomes
// between X-Ray and optical data
class Possibility {

    public:

      // A constructor consisting of the GRB ID number and its
      // corresponding multiplicity
      Possibility( string, int );

      //Accessors and observers for the class
      string get_ID();
      void set_ID(string);
      int get_multiplicity();
      void set_multiplicity(int);

    private:

      //variable for GRB ID
      string ID;
      //variable for multiplicity
      int multiplicity;
};

Possibility::Possibility( string id, int mult )
{
      //create constructor by the following assignment statements
      ID = id;
      multiplicity = mult;
}

string Possibility::get_ID()
{
    //return ID
    return ID;
}

void Possibility::set_ID( string id )
{
    //set ID
    ID = id;
}

int Possibility::get_multiplicity()
{
    //return multiplicity
    return multiplicity;
}

void Possibility::set_multiplicity( int m )
{
    //set multiplicity
    multiplicity = m;
}

class Trial {

    public:

      // Simple constructor
      Trial();

      // Loads X-Ray data file into the vector of GRBs
      int loadX_RayData(string);

      // Loads Beta_X data file. Uses the GRB ID to locate the
      // GRB ID in the vector, and then sets the GRB object to
      // also contain the correct Beta_X values
      int loadBeta_X(string);

      // Loads optical data file. Uses the GRB ID to locate the
      // GRB ID in the vector, and then sets the GRB object to
      // also contain the correct optical data values
      int loadOpticalData(string);

      // Loads wavelength data file and pairs each GRB based on
      // telescope and filter information appropriate wavelength
      int loadWavelengthData(string);

      // A function used to calculate the optical to X-Ray spectral index
      // as well as the upper and lower bounds on the uncertainty
      void calculate_Beta_OX();

      // prints out the GRB vector by calling the report method
      // of each GRB. Used for debugging.
      void report();

      //writes paired GRB data to file
      void write_paired_data();

    private:

      // A private function used to return the location
      // in the  vector of GRBs of a GRB with a
      // particular GRB ID used in the loadOpticalData
      // to locate where a particular GRB is
      int findGRB(string);

      // A private function which is used to pair X-Ray data and
      // optical data to a particular GRB based on a small temporal
      // separation in measurement time
      int matchGRB( string, double, int);

      // A private function which is used to determine if for a given
      // GRB ID in the optical data set, there exists at least one GRB
      // with identical ID read in by the X-Ray Data set
      int check_ID( string );

      // A private function which is used to remove entities that haven't
      // been paired with Beta_X data from vector XRay_entries
      int clean_XRay_entries();

      // A private function which is used to pair GRBs with their
      // appropriate optical wavelength based on telescope, instrument,
      // and filter information
      int matchFrequency( string, string, string, double, double );

      // A private function which is used to determine the total number of
      // possible pairings between optical and X-Ray data by, in essence,
      // conducting matrix multiplication between vectors XRay_entries and
      // optical_entries
      int find_total_possible_pairings();

      //define a vector of GRBs
      vector<GRB> GRBs;

      //define a vector of GRBs filled with optical data
      vector<GRB> GRBs_with_Opt;

      //define a vector which has number of elements corresponding to the total number
      //of different GRB IDs in the X-Ray file and elements corresponding to the
      //total number of entries per individual GRB ID
      //Ex: GRBs_with_Opt.size() = N; GRBs_With_Opt = {1,5,...,n}
      //This would correspond to the X-Ray file having N different GRB IDs,
      //and if the first two and last GRBs had IDs 050509B, 050709, and XXXXXX, then
      //GRB 050509B would have 1 entry, GRB 050709 would have 5 entries, and GRB XXXXXX
      //would have N entries
      vector<Possibility> XRay_entries;

      //define a vector identical to XRay_entries but for the optical file
      vector<Possibility> optical_entries;

      //define a vector of GRBs filled with IDs that exist in optical data but not in X-Ray
      // data
      vector<string> IDs_in_Opt_not_X;
};

Trial::Trial()
{
    //constructor is empty because it does not need any attributes
}


int Trial::loadX_RayData( string filename )
{
    //define variables for GRB attributes

    //define variable for GRB ID Number
    string ID;
    //define variable for X-Ray dt in seconds
    double dtX;
    //define variable for X-Ray exposure time
    double ExpX;
    //define variable for optical flux density
    double Fx;
    //define variable for std dev of optical flux density
    double sigmaX;

    //initialize variable for number of GRBs loaded
    int counter = 0;
    //initialize string variable for old ID used in GRB ID multiplicity determination
    string old_ID = "NULL";
    //initialize counter used for determining number of data points for each unique ID
    int entries_per_ID = 0;

    //define ifstream variable for file
    ifstream file;

    //open file specified by user
    file.open( filename.c_str() );

    //test to see if file opens successfully
    while ( !file )
    {
        //notify user of error and ask to input file name again
        cout << "Could not open file: " << filename << "." << endl;
        cout << "Please re-enter filename (or control-C to exit): ";
        cin >> filename; //assign user input to variable fileName
        file.open(filename.c_str()); //try again to open file specified by user
    }

    //display features
     cout << endl << right << setw(65) << "************************ X-Ray GRB Data ********"
          << "****************" << endl << endl;
     cout << right << setw(10) << "GRB ID" << right << setw(10) << "dt_X [s]"
          << right << setw(20) << "Exposure Time [s]" << right << setw(10)
          << "F_x [uJy]" << right << setw(15) << "sigma_X [uJy]" << endl << endl;
    //read in data from file and assign them to corresponding variables
    //read until file no longer has any more data
    while ( file >> ID >> dtX >> ExpX >> Fx >> sigmaX)
    {
       //create a GRB g and construct it using variables read from file
       GRB grb( ID, dtX, ExpX, Fx, sigmaX);

       //initialize new_ID to be the GRB ID read in from file
        string new_ID = ID;

        //check if this is the first entry read from optical file
        if ( counter == 0 )
        {
            //initialize old_ID to entry read from optical file if
            //this is the first entry from optical file
            old_ID = ID;
        }

        //check if new ID has been reached
        if ( new_ID.compare(old_ID) == 0 )
        {
            //increment counter for number of entries per unique ID because
            //a multiplicity of the same idea has been found OR it is the first trial
            entries_per_ID++;
        }
        else
        {
            //construct Possibility object with corresponding data
            Possibility new_Possibility( old_ID, entries_per_ID );

            //add in new element to vector containing number of entires per unique
            //GRB ID read from optical file
            XRay_entries.push_back(new_Possibility);

            //reset counter used for determining number of data points for each unique ID
            //counter is reset to 1 because in order to reach this else clause, new_ID !=
            //old_ID, meaning there is already a new entry
            entries_per_ID = 1;
        }

        old_ID = ID;
        //increment counter
        counter++;

       //put the GRB object in the vector of GRBs
       GRBs.push_back( grb );

       //display loaded features
       cout << right << setw(10) << ID << right << setw(10) << dtX  << right << setw(20)
            << ExpX << right << setw(10) <<  Fx << right << setw(15) << sigmaX <<endl;
    }

    //add last pair to vector of multiplicities
    //construct Possibility object with corresponding data
    Possibility new_Possibility( old_ID, entries_per_ID );

    //add in new element to vector containing number of entires per unique
    //GRB ID read from optical file
    XRay_entries.push_back(new_Possibility);

    cout << endl << "Number of GRBs loaded: " << counter;

    //close file
    file.close();

   return counter;
}

int Trial::loadBeta_X(string filename)
{
    //define variable for GRB ID
    string ID;
    //define variable for Beta_X
    double Beta_X;
    //define variable for Beta_X upper uncertainty
    double Beta_X_upper_sigma;
    //define variable for Beta_X lower uncertainty
    double Beta_X_lower_sigma;

    //define variable for total loaded Beta X elements
    int total_loaded = 0;
    //initialize variable for number of successful pairs after Beta_X
    double Beta_X_pairs = 0;
    //initialize variable for rate of successful pairings
    double pairing_rate = 0;

    //define ifstream variable for file
    ifstream file;

    //open file specified by user
    file.open( filename.c_str() );

    //test to see if file opens successfully
    while ( !file )
    {
        //notify user of error and ask to input file name again
        cout << "Could not open file: " << filename << "." << endl;
        cout << "Please re-enter filename (or control-C to exit): ";
        cin >> filename; //assign user input to variable fileName
        file.open(filename.c_str()); //try again to open file specified by user
    }

    //display features
     cout << endl << right << setw(60) << "************************Beta_X Data*************"
          << "*************"
          << endl << endl;
     cout << right << setw(10) << "GRB ID" << right << setw(10) << "Beta_X"
          << right << setw(20) << "Upper Uncertainty" << right << setw(20)
          << "Lower Uncertainty"<< endl << endl;
    //read in data from file and assign them to corresponding variables
    //read until file no longer has any more data

    while ( file >> ID >> Beta_X >> Beta_X_upper_sigma >> Beta_X_lower_sigma )
    {
        //display loaded features
        cout << right << setw(10) << ID << right << setw(10) << Beta_X
             << right << setw(20) << Beta_X_upper_sigma << right << setw(20)
             << Beta_X_lower_sigma << endl;

        //initialize success boolean to false
        bool success = false;

        //increment counter for successful load
        total_loaded++;

        // run through vector of GRBs
        for ( int a = 0; a < GRBs.size(); a++ )
        {
            //match GRB IDs
            if ( GRBs[a].getGRB_ID() == ID )
            {
                //assign attributes to appropriate GRB
                GRBs[a].set_Beta_X(Beta_X);
                GRBs[a].set_Beta_X_upper_sigma(Beta_X_upper_sigma);
                GRBs[a].set_Beta_X_lower_sigma(Beta_X_lower_sigma);

                //signify match
                success = true;
                //increment success counter for successful pairing
                Beta_X_pairs++;
            }
        }

        //check if no match was found
        if ( success == false )
        {
            //notify user that no match was able to be found
            cout << endl << "Unable to match GRB ID " << ID << " with Beta_X " << Beta_X
                 << endl << endl;
        }
    }

    //calculate percent of loaded GRBs that are paired
    pairing_rate = (Beta_X_pairs /  GRBs.size())*100;

    //notify user of loading statistics
    cout << endl << right << setw(58) << "Number of loaded GRBs from Beta_X file: "
         << right << setw(2) << total_loaded;
    cout << endl << right << setw(57) << "Number of successful pairings with X-Ray Data: "
         << right << setw(3) << Beta_X_pairs;
    cout << endl << right << setw(55) << fixed << showpoint << setprecision(1)
         << "Pairing Rate: " << right << setw(4) << pairing_rate << right << setw(1)
         << "%" << endl;

    cout << endl << right << setw(60) << "...Cleaning up X-Ray Entries..." << endl;
    clean_XRay_entries();

    //close file
    file.close();

    return total_loaded;
}

int Trial::loadOpticalData( string filename )
{
    //define variables for GRB attributes

    //define variable for GRB ID Number
    string ID;
    //define variable for optical dt in hours
    double dtO_hours;
    //define variable for optical dt in seconds
    double dtO_seconds;
    //define variable for optical telescope name
    string tel;
    //define variable for instrument name
    string inst;
    //define variable for optical filter name
    string fil;
    //define variable for optical exposure time in seconds
    double ExpO;
    //define variable for optical flux density in uJy
    double Fo;
    //define variable for optical flux density std dev in uJy
    double sigmaO;

    //initialize variable for number of subjects loaded
    double total_loaded = 0;
    //initialize variable for rate of successful pairings
    double pairing_rate = 0;
    //initialize counter for no match found
    int no_match_found_counter = 0;

    //initialize string variable for old ID used in GRB ID multiplicity determination
    string old_ID = "NULL";
    //initialize counter used for determining number of data points for each unique ID
    int entries_per_ID = 0;
    //initialize variable for number of GRBs that are not in
    //BOTH the X-Ray and optical files
    int disjoint = 0;
    //initialize variable for a check on the number of GRBs that are not in
    //BOTH the X-Ray and optical files
    int disjoint_confirm = 0;

    //define ifstream variable for file
    ifstream file;

    //open file specified by user
    file.open( filename.c_str() );

    //test to see if file opens successfully
    while ( !file )
    {
        //notify user of error and ask to input file name again
        cout << "Could not open file: " << filename << "." << endl;
        cout << "Please re-enter filename (or control-C to exit): ";
        cin >> filename; //assign user input to variable fileName
        file.open(filename.c_str()); //try again to open file specified by user
    }
    //display features
     cout << endl << right << setw(110) << "**********************************************"
          <<" Optical GRB Data **********************************************" << endl
          << endl;
     cout << right << setw(10) << "GRB ID" << right << setw(15) << "dt_O [s]"
          << right << setw(15) << "Telescope" << right << setw(15)
          << "Instrument" << right << setw(10) << "Filter" << right << setw(20)
          << "Exposure Time [s]" << right << setw(10) << "F_o [uJy]" << right
          << setw(15) << "sigma_O [uJy]" << endl << endl;

    //read in data from file and assign them to corresponding variables
    //read until file no longer has any more data
    while (file >> ID >> dtO_hours >> tel >> inst >> fil >> ExpO >> Fo >> sigmaO)
    {
        //transform optical dt measurement from hours into seconds
        dtO_seconds = 3600 * dtO_hours;

        //display loaded features
        cout << right << setw(10) << ID << right << setw(15) << dtO_seconds
             << right << setw(15) << tel << right << setw(15)
             << inst << right << setw(10) << fil << right << setw(20)
             << ExpO << right << setw(10) << Fo << right
             << setw(15) << sigmaO << endl;

        //initialize variable for location in GRB vector
        int location = 0;
        //initialize counter for checking if any pairings were made at all
        int check = 0;

        //initialize new_ID to be the GRB ID read in from file
        string new_ID = ID;

        //check if this is the first entry read from optical file
        if ( total_loaded == 0 )
        {
            //initialize old_ID to entry read from optical file if
            //this is the first entry from optical file
            old_ID = ID;
        }

        //check if new ID has been reached
        if ( new_ID == old_ID )
        {
            //increment counter for number of entries per unique ID because
            //a multiplicity of the same idea has been found OR it is the first trial
            entries_per_ID++;
        }
        else
        {
            //construct Possibility object with corresponding data
            Possibility new_Possibility( old_ID, entries_per_ID );

            //add in new element to vector containing number of entires per unique
            //GRB ID read from optical file
            optical_entries.push_back(new_Possibility);

            //disjoint += check_ID(old_ID);

            //reset counter used for determining number of data points for each unique ID
            //counter is reset to 1 because in order to reach this else clause, new_ID !=
            //old_ID, meaning there is already a new entry
            entries_per_ID = 1;
        }

        old_ID = ID;

        //increment counter for total loaded Optical elements
        total_loaded++;

        while ( location < GRBs.size() && location != -1)
        {
            location = matchGRB(ID, dtO_seconds, location);
            check++;

            //cout << "\nBeginning checks\n";
            //cout << "Location is: " << location << endl;

            //check if pairing is made
            if( location != -1 )
            {
                //construct a GRB object that will be added
                //into vector of GRBs with optical data
                GRB copy_grb( GRBs[location].getGRB_ID(), GRBs[location].get_dt_XRay(),
                             GRBs[location].get_ExpT_XRay(), GRBs[location].get_F_x(),
                             GRBs[location].get_sigma_x());

                //set corresponding GRB appropriate Beta_X parameters
                copy_grb.set_Beta_X(GRBs[location].get_Beta_X());
                copy_grb.set_Beta_X_lower_sigma(GRBs[location].get_Beta_X_lower_sigma());
                copy_grb.set_Beta_X_upper_sigma(GRBs[location].get_Beta_X_upper_sigma());

                //set corresponding GRB appropriate optical parameters
                copy_grb.set_dt_Opt(dtO_seconds);
                copy_grb.set_telescope(tel);
                copy_grb.set_instrument(inst);
                copy_grb.set_filter(fil);
                copy_grb.set_Exp_Opt(ExpO);
                copy_grb.set_F_o(Fo);
                copy_grb.set_sigma_o(sigmaO);

                //add newly-created GRB with optical data to
                //vector of GRBs with optical data
                GRBs_with_Opt.push_back(copy_grb);

                //increment location
                location ++;
                //increment counter of successful pairing
                optical_pairs++;
            }

        }
        //check if no pairing was made
        if ( check == 1)
        {
//            //notify user that no match was able to be found
//            cout << "Unable to match GRB " << ID << " with optical dt " << dtO_hours
//                 << " [hr] = " << dtO_seconds << " [s]." << endl << endl;

           //increment counter for no match made
           no_match_found_counter++;
        }
    }

    //add last pair to vector of multiplicities
    //construct Possibility object with corresponding data
    Possibility new_Possibility( old_ID, entries_per_ID );

    //add in new element to vector containing number of entires per unique
    //GRB ID read from optical file
    optical_entries.push_back(new_Possibility);

    total_possible_pairings = find_total_possible_pairings();

    //calculate percent of loaded GRBs that are paired
    pairing_rate = (optical_pairs /  total_possible_pairings )*100;

    //notify user of loading statistics
    cout << endl << right << setw(57) << "Number of loaded GRBs from optical data file: "
         << right << setw(3) << int(total_loaded);
    cout << endl << right << setw(57) << fixed << showpoint << setprecision(0)
         << "Number of possible pairs: " << right << setw(3)
         << int(total_possible_pairings);
    cout << endl << right << setw(57) << fixed << showpoint << setprecision(0)
         << "Number of successful pairs: " << right << setw(3) << int(optical_pairs);
    cout << endl << right << setw(55) << fixed << showpoint << setprecision(1)
         << "Pairing Rate: " << right << setw(5) << pairing_rate << "%";

    //close file
    file.close();

    return total_loaded;
}

int Trial::loadWavelengthData(string filename)
{
    //define variables for GRB attributes

    //define variable for telescope name
    string telName;
    //define variable for instrument name
    string instrumentName;
    //define variable for filter name
    string filterName;
    //define variable for wavelength in [nm]
    double wavelength;
    //define variable for frequency in [Hz]
    double frequency;

    //initialize variable for number of subjects loaded
    int loaded = 0;
    //initialize variable for number of successful pairings loaded
    int success_counter = 0;
     //initialize variable for the counter to check pairing
    int check_counter = 0;

    //define ifstream variable for file
    ifstream file;

    //open file specified by user
    file.open( filename.c_str() );

     //test to see if file opens successfully
    while ( !file )
    {
        //notify user of error and ask to input file name again
        cout << "Could not open file: " << filename << "." << endl;
        cout << "Please re-enter filename (or control-C to exit): ";
        cin >> filename; //assign user input to variable fileName
        file.open(filename.c_str()); //try again to open file specified by user
    }

    //display features
     cout << endl << right << setw(80) << "******************************* Wavelength Data"
          <<" *******************************"
          << endl << endl;
     cout << right << setw(20) << "Telescope" << right << setw(15) << "Instrument"
          << right << setw(10) << "Filter" << right << setw(15)
          << "Wavelength" << right << setw(20) << "Frequency" << endl << endl;

    //read in data from file and assign them to corresponding variables
    //read until file no longer has any more data
    while ( file >> telName >> instrumentName >> filterName >> wavelength >> frequency )
    {
        //display loaded features
        cout << right << setw(20) << telName << right << setw(15) << instrumentName
          << right << setw(10) << filterName << right << setw(15)
          << wavelength << right << setw(20) << frequency << endl;

        //call function to match GRB objects with appropriate wavelength data
        success_counter += matchFrequency( telName, instrumentName, filterName, wavelength,
                                           frequency );

        //increment counter for successful load
        loaded++;

        //notify user that a match was able to be found
        //cout << "Able to match telescope " << telName << " with instrument " <<
        //instrumentName <<  " and with filter " << filterName << "." << endl;
    }

    cout << endl;
    //display those GRBs that couldn't get paired
    for ( int a = 0; a < GRBs_with_Opt.size(); a++ )
    {
        //check to see if GRB has been paired with Beta_X, optical, but not wavelength data
        if (  GRBs_with_Opt[a].get_frequency_Opt() == -1 )
        {
            cout << "\nGRB " << GRBs_with_Opt[a].getGRB_ID() << " with optical dt "
                 << GRBs_with_Opt[a].get_dt_Opt() << ", telescope "
                 << GRBs_with_Opt[a].get_telescope() << ", instrument "
                 << GRBs_with_Opt[a].get_instrument() <<  ", and with filter "
                 << GRBs_with_Opt[a].get_filter() <<" unpaired." << endl;
            check_counter++;
        }
    }

    //display loaded statistics
    cout << endl << endl << right << setw(57) << "Number of Wavelength Sets Loaded: "
         << right << setw(3) << loaded << endl;
    cout << right << setw(57) << "Number of Unsuccessfully Paired: " << right << setw(3)
         << check_counter << endl;
    cout << right << setw(57) << "Number of Successfully Paired: " << right << setw(3)
         << success_counter << endl;
    cout << right << setw(55) << "Overall Success Rate: " << right << setw(4) << fixed
         << showpoint << setprecision(1)
         << 100 * (success_counter / total_possible_pairings ) << right << setw(1) << "%"
         << endl;

    //close file
    file.close();

    return loaded;
}

void Trial::calculate_Beta_OX()
{
    //initialize variable for counter of successful Beta_OX calculations
    int success_counter = 0;
    //initialize counter for nan calculations of Beta_OX
    int nan = 0;

    //initialize variables necessary for calculation
    double F_x = 0;
    double F_o = 0;
    double sigma_x = 0;
    double sigma_o = 0;
    double frequency_X = FREQUENCY_XRAY;
    double frequency_O = 0;
    double Beta_OX = 0;
    double sigma_OX_upper = 0;
    double sigma_OX_lower = 0;

    //display features
     cout << endl << right << setw(160) << "*********************************************"
          <<"********************************* Beta_OX Data ******************************"
          <<"*************************************************"
          << endl << endl;
     cout << right << setw(10) << "GRB ID" << right << setw(15)
          << "F_x [uJy]" << right << setw(15) << "sigma_X [uJy]" << right << setw(15)
          << "F_o [uJy]" << right << setw(15) << "sigma_o [uJy]" << right << setw(25)
          << "Freq_X" << right  << setw(25) << "Freq_O" << right << setw(10)
          << "Beta_OX" << right << setw(20) << "Upper sigma_OX" << right
          << setw(20) << "Lower sigma_OX" << endl;

    // run through vector of GRBs
    for ( int a = 0; a < GRBs_with_Opt.size(); a++ )
    {
        //check to see if GRB has been fully paired
        if (GRBs_with_Opt[a].get_frequency_Opt() != -1)
        {

            // Extract values from GRB objects in vector of GRBs
            F_x = GRBs_with_Opt[a].get_F_x();
            F_o = GRBs_with_Opt[a].get_F_o();
            frequency_O = GRBs_with_Opt[a].get_frequency_Opt();
            sigma_x = GRBs_with_Opt[a].get_sigma_x();
            sigma_o = GRBs_with_Opt[a].get_sigma_o();

            //calculate Beta_OX
            Beta_OX = log(F_x / F_o) / log(frequency_X / frequency_O );

            //check to see if calculation returns nan
            if (Beta_OX != Beta_OX )
            {
                Beta_OX = 0;
                nan++;
            }

            //calculate upper bound on uncertainty for Beta_OX
            sigma_OX_upper = log( (1 + (sigma_x / F_x)) / (1 - (sigma_o / F_o))) /
            log( frequency_X / frequency_O );
            //calculate lower bound on uncertainty for Beta_OX
            sigma_OX_lower = abs( log( (1 - (sigma_x / F_x)) / (1 + (sigma_o / F_o))) /
                                  log( frequency_X / frequency_O ));

            //check to see if calculation returns nan
            if (sigma_OX_upper != sigma_OX_upper)
            {
                sigma_OX_upper = 0;
            }
            //check to see if calculation returns nan
            if (sigma_OX_lower != sigma_OX_lower )
            {
                sigma_OX_lower = 0;
            }
            //set calculated value into GRB Beta_OX attribute
            GRBs_with_Opt[a].set_Beta_OX(Beta_OX);
            GRBs_with_Opt[a].set_sigma_OX_upper(sigma_OX_upper);
            GRBs_with_Opt[a].set_sigma_OX_lower(sigma_OX_lower);

            cout << right << setw(10) << GRBs_with_Opt[a].getGRB_ID() << right << setw(15)
                  << GRBs_with_Opt[a].get_F_x() << right << setw(15)
                  << GRBs_with_Opt[a].get_sigma_x() << right << setw(15)
                  << GRBs_with_Opt[a].get_F_o() << right << setw(15)
                  << GRBs_with_Opt[a].get_sigma_o() << right << setw(25)
                  << GRBs_with_Opt[a].get_frequency_XRay() << right  << setw(25)
                  << GRBs_with_Opt[a].get_frequency_Opt() << right << setw(10)
                  << GRBs_with_Opt[a].get_Beta_OX() << right << setw(20)
                  << GRBs_with_Opt[a].get_sigma_OX_upper() << right << setw(20)
                  << GRBs_with_Opt[a].get_sigma_OX_lower() << endl;

            //increment success counter for successful calculation
            success_counter++;
        }
    }

    //display loaded statistics
    cout << endl << endl << right << setw(57) << "Number of Successful Beta_OX"
         <<"Calculations: " << right << setw(3) << success_counter << endl;
    cout << right << setw(55) << "Overall Success Rate: " << right << setw(4) << fixed
         << showpoint << setprecision(1) << 100 * (success_counter /
         total_possible_pairings ) << right << setw(1) << "%" << endl;
}

void Trial::report()
{
    //run through vector of GRB
    for ( int a = 0; a < GRBs.size(); a++)
    {
        GRBs[a].report();
    }
}

void Trial::write_paired_data()
{
    //cast temporal percent difference to string
    int dt = static_cast<int>(DT_PERCENT_DIF);
    string percent_dif = to_string(dt);
    //define variable for name of files
    string filename_comprehensive = "./Written_Files/Comprehensive_Paired_Data_Table_"
    + percent_dif + "%.csv";
    string filename_terse = "./Written_Files/GRB_Pairings-dt_" + percent_dif + "%.csv";

    //define ofstream variables for the files to be written
    ofstream myFile_comprehensive;
    ofstream myFile_terse;

    //open file for comprehensive data file
    myFile_comprehensive.open(filename_comprehensive);

    //print out headers for file
    myFile_comprehensive << "GRB ID,X-Ray dt [hr],X-Ray Exposure Time [s],"
           << "F_x [uJy],Sigma_x [uJy],Beta_X,"
           << "Beta_X Upper Sigma,Beta_X Lower Sigma,"
           << "Optical dt [hr],Telescope,Instrument,Filter,"
           << "Optical Exposure Time [s],F_o [uJy], Sigma_o [uJy],"
           << "Wavelength_X [nm],Frequency_X [Hz],Wavelength_o [nm],"
           << "Frequency_o [Hz],Beta_OX,Upper Bound of Sigma_OX,"
           << "Lower Bound of Sigma_OX,\n";

    //run through vector of populated GRBs
    for ( int a = 0; a < GRBs_with_Opt.size(); a++)
    {
        //double check to see if GRB is fully populated
        if ( GRBs_with_Opt[a].get_frequency_Opt() != -1)
        {
        myFile_comprehensive << fixed << setprecision(2) << GRBs_with_Opt[a].getGRB_ID()
                             << "," << fixed << setprecision(2)
                             <<  GRBs_with_Opt[a].get_dt_XRay()/3600 << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_ExpT_XRay() << ","
                             << fixed << setprecision(2) << GRBs_with_Opt[a].get_F_x()
                             << "," << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_sigma_x() << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_Beta_X()
                             << "," << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_Beta_X_upper_sigma() << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_Beta_X_lower_sigma()
                             << "," << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_dt_Opt()/3600 << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_telescope()
                             << "," << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_instrument() << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_filter()
                             << "," << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_Exp_Opt() << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_F_o()
                             << "," << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_sigma_o() << "," << fixed
                             << GRBs_with_Opt[a].get_frequency_XRay() << "," << fixed
                             << setprecision(2)
                             << GRBs_with_Opt[a].get_frequency_Opt() << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_Beta_OX() << ","
                             << fixed << setprecision(2)
                             << GRBs_with_Opt[a].get_sigma_OX_upper() << "," << fixed
                             << setprecision(2) << GRBs_with_Opt[a].get_sigma_OX_lower()
                             << "\n";
        }
    }
    //close file
    myFile_comprehensive.close();

    //open file for comprehensive data file
    myFile_terse.open(filename_terse);

    //print out headers for file
    myFile_terse << "GRB ID,X-Ray dt [hr],"
           << "Optical dt [hr],"
           << "|dt_x - dt_o| [hr],"
           << "Beta_X,"
           << "Beta_X Upper Sigma,Beta_X Lower Sigma,"
           << "Beta_OX,Upper Bound of Sigma_OX,"
           << "Lower Bound of Sigma_OX,\n";

    for ( int a = 0; a < GRBs_with_Opt.size(); a++)
    {
        //check to see if GRB is fully populated
        if ( GRBs_with_Opt[a].get_frequency_Opt() != -1)
        {
        myFile_terse << GRBs_with_Opt[a].getGRB_ID() << "," << fixed << setprecision(2)
                     << GRBs_with_Opt[a].get_dt_XRay()/3600
               << "," << fixed << setprecision(2) << GRBs_with_Opt[a].get_dt_Opt()/3600
               << "," << fixed << setprecision(2)
               << abs(GRBs_with_Opt[a].get_dt_Opt() - GRBs_with_Opt[a].get_dt_XRay())/3600
               << "," << fixed << setprecision(2) << GRBs_with_Opt[a].get_Beta_X()
               << "," << fixed << setprecision(2)
               << GRBs_with_Opt[a].get_Beta_X_upper_sigma() << "," << fixed
               << setprecision(2) << GRBs_with_Opt[a].get_Beta_X_lower_sigma()
               << "," << fixed << setprecision(2) << GRBs_with_Opt[a].get_Beta_OX() << ","
               << fixed << setprecision(2) << GRBs_with_Opt[a].get_sigma_OX_upper()
               << "," << fixed << setprecision(2) << GRBs_with_Opt[a].get_sigma_OX_lower()
               << "\n";
        }
    }
    //close file
    myFile_terse.close();

}

int Trial::matchFrequency( string tel, string inst, string filt, double wavelength,
                          double frequency )
{
    //initialize counter for successful pairing
    int thatsapair = 0;
    //initialize boolean for successful pairing
    bool success = false;

    //run through vector of GRB
    for ( int a = 0; a < GRBs_with_Opt.size(); a++)
    {
        //test to see if GRB ID matches passed ID
        //also test to see if GRB has not yet been populated with wavelength parameters
        if ( GRBs_with_Opt[a].get_telescope() == tel && GRBs_with_Opt[a].get_instrument()
            == inst && GRBs_with_Opt[a].get_filter() == filt &&
             GRBs_with_Opt[a].get_frequency_Opt() == -1)
        {
            //pair GRB with frequency data
            GRBs_with_Opt[a].set_frequency_Opt(frequency);

            //increment success counter
            thatsapair++;
            //assign boolean to true to signify successful pairing
            success = true;
        }
    }

    //return number of successful pairs
    return thatsapair;

}

int Trial::check_ID( string ID )
{
    //define boolean for match
    bool corresponding_ID_found = false;
    //initialize counter for number of disjoint GRBs
    int disjoint_counter = 0;

    cout << "\n\nIn check_ID!" << endl << endl;

    for ( int kumquat = 0; kumquat < XRay_entries.size(); kumquat++ )
    {
        //check if passed optical IDs are identical to those in GRBs that have
        //already been paired with Beta_X data
        if ( ID.compare( XRay_entries[kumquat].get_ID() ) == 0 )
        {
            //set boolean to true if match found
            corresponding_ID_found = true;
        }
    }
    //check if no match was found
    if ( corresponding_ID_found == false )
    {
        //add ID that is in optical data set but not X-Ray to
        //vector containing all such IDs
        IDs_in_Opt_not_X.push_back(ID);
        //increment counter for lack of match
        disjoint_counter++;
        cout << "\n\nID in Optical but not in X-Ray :(" << endl << endl;
    }

    return disjoint_counter;

}

int Trial::clean_XRay_entries()
{
    //initialize counter for number we expect to keep
    int should_not_keep = 0;

    cout << endl << right << setw(58) << "Original Number of Unique X-Ray GRBs: "
         << right << setw(2) << XRay_entries.size();

    //traverse through vector XRay_entries
    for (int q = 0; q < XRay_entries.size(); q++ )
    {
        //initialize boolean signifying that XRay_entries[q] should be
        //kept to false
        bool keeper = false;

        //traverse through vector GRBs that contain only some with BetaX data
        for (int n = 0; n < GRBs.size(); n++ )
        {
            //check if the ID in XRay_entries and GRBs are identical and that
            //ID has been paired with BetaX data
            if ( XRay_entries[q].get_ID().compare(GRBs[n].getGRB_ID()) == 0 &&
                 GRBs[n].get_Beta_X() != 31415926535 )
            {
                //assign boolean for keeper to true
                keeper = true;
            }
        }
        //check if no match was made
        if ( keeper == false )
        {
            //increment counter for number we should not keep
            should_not_keep++;

            //initialize iterator to beginning of XRay_entries vector
            auto it = XRay_entries.begin();
            //traverse through XRay_entries with iterator
            while (it != XRay_entries.end())
            {
                // remove entity that had no pairing
                if ((*it).get_ID() == XRay_entries[q].get_ID() )
                {
                    // erase() invalidates the iterator, use returned iterator
                    it = XRay_entries.erase(it);

                    //check if index is at 0
                    if ( q > 0 )
                    {
                        //shift index backwards one due to erased entity
                        q -= 1;
                    }
                    else
                    {
                        //increment iterator
                        //without this the iterator would remain stuck on
                        //index=0 and continuously delete entities
                        ++it;
                    }
                }
                else
                {
                    //increment iterator
                    ++it;
                }
            }
        }
    }

    cout << endl << right << setw(58) << "Final Number of Unique X-Ray GRBs: " << right
         << setw(2) << XRay_entries.size() << endl;
    cout << right << setw(58) << "Number Removed for Lack of Beta_X Pairing: " << right
         << setw(2) << should_not_keep << endl;
}

int Trial::find_total_possible_pairings()
{
    //initialize counters for disjointed IDs
    int in_Opt_not_XRay = 0;
    int in_XRay_not_Opt = 0;

    //initialize variables for old and new vector sizes
    int old_XRay_entries_size = XRay_entries.size();
    int new_XRay_entries_size = 0;
    int old_optical_entries_size = optical_entries.size();
    int new_optical_entries_size = 0;

    //initialize variable for total number of possibilities
    int total_possibilities = 0;

    //traverse through vector XRay_entries
    for (int q = 0; q < XRay_entries.size(); q++ )
    {
        //initialize boolean signifying that ID is in XRay_entries
        //and also in optical_entries
        bool keeper = false;

        //traverse through optical_entries vector
        for (int n = 0; n < optical_entries.size(); n++ )
        {
            //check if the ID in XRay_entries and optical_entries are identical
            if ( XRay_entries[q].get_ID().compare(optical_entries[n].get_ID()) == 0 )
            {
                //assign boolean for keeper to true
                keeper = true;
            }
        }
        //check if no match was made
        if ( keeper == false )
        {
            //increment counter for number of IDs that exist in
            //XRay_entries but not in optical_entries
            in_XRay_not_Opt++;

            //initialize iterator to beginning of XRay_entries vector
            auto it = XRay_entries.begin();
            //traverse through XRay_entries with iterator
            while (it != XRay_entries.end())
            {
                // remove entity that had no pairing
                if ((*it).get_ID() == XRay_entries[q].get_ID() )
                {
                    // erase() invalidates the iterator, use returned iterator
                    it = XRay_entries.erase(it);

                    //check if index is at 0
                    if ( q > 0 )
                    {
                        //shift index backwards one due to erased entity
                        q -= 1;
                    }
                    else
                    {
                        //increment iterator
                        //without this the iterator would remain stuck on
                        //index=0 and continuously delete entities
                        ++it;
                    }
                }
                else
                {
                    //increment iteratorDecember 15
                    ++it;
                }
            }
        }
    }
    //appropriately change variable for size of XRay_entries
    new_XRay_entries_size = XRay_entries.size();

    //traverse through vector optical_entries
    for (int l = 0; l < optical_entries.size(); l++ )
    {
        //initialize boolean signifying that ID is in XRay_entries
        //and also in optical_entries
        bool keepme = false;

        //traverse through optical_entries vector
        for (int u = 0; u < XRay_entries.size(); u++ )
        {
            //check if the ID in XRay_entries and optical_entries are identical
            if ( optical_entries[l].get_ID().compare(XRay_entries[u].get_ID()) == 0 )
            {
                //assign boolean for keeper to true
                keepme = true;
            }
        }
        //check if no match was made
        if ( keepme == false )
        {
            //increment counter for number of IDs that exist in
            //XRay_entries but not in optical_entries
            in_Opt_not_XRay++;

            //initialize iterator to beginning of XRay_entries vector
            auto it_2 = optical_entries.begin();
            //traverse through XRay_entries with iterator
            while (it_2 != optical_entries.end())
            {
                // remove entity that had no pairing
                if ((*it_2).get_ID() == optical_entries[l].get_ID() )
                {
                    // erase() invalidates the iterator, use returned iterator
                    it_2 = optical_entries.erase(it_2);

                    //check if index is at 0
                    if ( l > 0 )
                    {
                        //shift index backwards one due to erased entity
                        l -= 1;
                    }
                    else
                    {
                        //increment iterator
                        //without this the iterator would remain stuck on
                        //index=0 and continuously delete entities
                        ++it_2;
                    }
                }
                else
                {
                    //increment iterator
                    ++it_2;
                }

            }
        }
    }

    //appropriately change variable for new size of optical_entries
    new_optical_entries_size = optical_entries.size();

    cout << endl << endl << right << setw(60) << "Final X-Ray Entries (left) and Optical"
         << "Entries ( right): " << endl;
    cout << endl << right << setw(10) << "GRB ID" << right << setw(15) << "Multiplicity"
         << right << setw(25) << "GRB ID" << right << setw(15) << "Multiplicity" << endl
         << endl;

    for ( int y = 0; y < optical_entries.size(); y++ )
    {
       cout << right << setw(10) << XRay_entries[y].get_ID() << right << setw(15)
            << XRay_entries[y].get_multiplicity() << right << setw(25)
            << optical_entries[y].get_ID() << right << setw(15)
            << optical_entries[y].get_multiplicity() << endl;
    }

    cout << endl;

    //display appropriate statistics
    cout << endl << right << setw(58) << "Original Number of Unique X-Ray GRBs with Beta_X"
         <<" Data: " << right << setw(2) << old_XRay_entries_size << endl;
    cout << right << setw(58) << "Number of Disjoint GRBs Removed from X-Rays: " << right
         << setw(2) << in_XRay_not_Opt << endl;
    cout << right << setw(58) << "Final Number of Unique X-Ray GRBs: " << right << setw(2)
         << new_XRay_entries_size << endl;
    cout << endl << right << setw(58) << "Original Number of Unique Optical Entries: "
         << right << setw(2) << old_optical_entries_size << endl;
    cout << right << setw(58) << "Number of Disjoint GRBs Removed from Optical: " << right
         << setw(2) << in_Opt_not_XRay << endl;
    cout << right << setw(58) << "Final Number of Unique Optical GRBs: " << right
         << setw(2) << new_optical_entries_size << endl;

    //calculate the total number of possibilities that
    //could be obtained from perfect matching
    for ( int g = 0; g < optical_entries.size(); g++ )
    {
        //multiply number of options in XRay_entries by
        //number of options in optical_entries per unique
        //GRB ID to determine total number of possibilities
        total_possibilities += XRay_entries[g].get_multiplicity() *
        optical_entries[g].get_multiplicity();
    }
    //return number of erased
    return total_possibilities;
}

int Trial::matchGRB( string ID, double dtO_s, int location )
{
    //initialize bool variable to determine if match occurs
    bool paired = false;

    //run through vector of GRB
    for ( int a = location; a < GRBs.size(); a++)
    {
        //test to see if GRB ID matches passed ID
        //test to see if GRB has been populated with Beta_X data
        if ( GRBs[a].getGRB_ID() == ID && GRBs[a].get_Beta_X() != 31415926535 )
        {
            //test to see if temporal separation is small enough and GRB hasn't already
            //been populated with optical data
            if ( ( 100 * abs( GRBs[a].get_dt_XRay() - dtO_s ) / dtO_s ) < DT_PERCENT_DIF )
            {
                //return location of this GRB
                return a;
                //assign match to true to signify successful match
                paired = true;
            }
        }
    }
    //test to see if no match occurred
    if ( paired == false )
    {
        //return -1 if no match occurred
        return -1;
    }
}

int Trial::findGRB(string ID)
{
    //initialize bool variable to determine if match occurs
    bool match = false;

    //run through vector of GRBs
    for ( int a = 0; a < GRBs.size(); a++)
    {
        //test to see if Subject ID matches passed ID
        if ( GRBs[a].getGRB_ID() == ID )
        {
            //return location of this GRB
            return a;
            //assign match to true to signify successful match
            match = true;
        }
    }
    //test to see if no match occurred
    if ( match == false )
    {
        //return -1 if no match occurred
        return -1;
    }
}

/**********************************************************************
                             BEGIN MAIN
***********************************************************************/

int main()
{
    //define variable for string type of desired temporal % difference
    string percent_dif;
    //define variable for file name of X-Ray data file
    string XRayDataFile_name;
    //define variable for file name of Beta_X data file
    string Beta_X_File_name;
    //define variable for file name of optical data file
    string OpticalData_name;
     //define variable for file name of wavelength data file
    string WavelengthData_name;
    //create a trial object
    Trial t1;


     //ask user for name of file they wish to load for X-Ray data
    cout << "Please enter the name of the X-Ray data file: ";
    //assign user input to variable for file name of subject data
    cin >> XRayDataFile_name;

    //call function for loading X-Ray data
    //pass name of file for X-Ray data
    t1.loadX_RayData(XRayDataFile_name);

    //formatting
    cout << endl << endl;

    //ask user for name of file they wish to load for Beta_X data
    cout << "Please enter the name of the Beta_X data file: ";
    //assign user input to variable for file name of subject data
    cin >> Beta_X_File_name;

    //call function for loading X-Ray data
    //pass name of file for X-Ray data
    t1.loadBeta_X(Beta_X_File_name);

    //formatting
    cout << endl << endl;

    //ask user for desired temporal percent difference
    cout << "Please enter the desired temporal percent difference [%]: ";
    //assign user input to constant for temporal percent difference
    cin >> DT_PERCENT_DIF;

    //ask user for name of file they wish to load for side effects
    cout << "Please enter the name of the optical data file: ";
    //assign user input to variable for file name for side effects
    cin >> OpticalData_name;

    //call function for loading optical data
    //pass name of file for optical data
    t1.loadOpticalData(OpticalData_name);

    //formatting
    cout << endl << endl;

    //ask user for name of file they wish to load for side effects
    cout << "Please enter the name of the wavelength data file: ";
    //assign user input to variable for file name for side effects
    cin >> WavelengthData_name;

    //call function for loading wavelength data
    //pass name of file for wavelength data
    t1.loadWavelengthData(WavelengthData_name);


    cout << endl; //for formatting purposes
    //set numeric output formatting
    cout << fixed << showpoint << setprecision(2);

    cout << endl << endl << "Press any key to calculate Beta_OX." << endl << endl;
    cin.ignore();
    cin.get();

    //calculate X-Ray to optical spectral indices
    t1.calculate_Beta_OX();

    //write paired data to .csv file
    t1.write_paired_data();

    return 0;
}
