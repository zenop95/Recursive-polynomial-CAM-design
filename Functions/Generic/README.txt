******************************************************************************

MATLAB Solar System Ephemeris Toolbox                             29 June 2001
                                             C.F. Walker, Kennedy Space Center

******************************************************************************

INTRODUCTION
------------

This README, along with the "help" information for each function, describes
how to install and use the MATLAB Solar System Ephemeris Toolbox.  The
functions in this toolbox enable the user to read and interpolate JPL
ephemerides in the MATLAB computing environment. 

This README file is liberally plagiarized from E.M. Standish's excellent
README file on the JPL-provided ephemerides and associated FORTRAN programs.
That README file contains information on available ephemerides, ephemeris
sizes, and references which describe the ephemerides which are absent from
this README.  It is suggested that the user read Dr. Standish's file for
background on the ephemerides.

The final goal of the installation process is the successful exectution of
"testeph.m".  This function exercises the functions which read and interpolate
planetary and lunar coordinates from a binary direct-read ephemeris file and
compares these results against corresponding numbers produced at JPL.
"testeph.m" uses the functirons which are of eventual interest to the user.

It is strongly suggested that a potential user first read through this README
in its entirety.  This will provide an idea of what is involved in both the
installation and in the usage of the toolbox.

This README contains the following sections:

INTRODUCTION
LOCATION OF THE FILES
WHAT TO DO
FILES TO BE RETRIEVED
BRIEF ITEM DESCRIPTION
CREATING THE BINARY EPHEMERIS FILE
VERIFYING PROPER INSTALLATION
USING THE TOOLBOX
ASSISTANCE

******************************************************************************

LOCATION OF THE FILES
---------------------

The JPL ASCII ephemeris header and data files are available over the Internet
from the anonymous ftp server:

                 ftp://ssd.jpl.nasa.gov

When connected, go to the "pub/eph/export" directory.

******************************************************************************

WHAT TO DO
----------

Installation consists of three major steps:

1. Retrieve the items, listed in the next section, from the anonymous ftp
   site.

2. Create the binary direct read ephemeris files.

3. Run "testeph.m" to verify proper installation.

******************************************************************************

FILES TO BE RETRIEVED
---------------------

Retrieve the following items from the anonymous ftp site:

	/pub/eph/export/ascii/ascSYYYY.XXX  (multiple files)
	/pub/eph/export/ascii/header.XXX
	/pub/eph/export/test-data/testpo.XXX

where "SYYYY" is the starting year of each particular ephemeris block ("S" is
the sign : "p" for + and "m" for -) and "XXX" is the ephemeris number: 200,
403, 405, or 406.  See choices below for different ephemerides.

The ASCII files come in 20-year blocks (about 3.5 Mbytes each). The blocks are
converted into binary blocks on the user's computer, using "asc2bin_eph.m";
then, contiguous binary blocks may be merged together to form a single
ephemeris file using "binmerge.m" or "mergemall.m".

******************************************************************************

BRIEF ITEM DESCRIPTION
----------------------

  "Contents.m"   : MATLAB-formatted contents file.

  Ephemeris Reading/Interpolation
  -------------------------------

  "pleph.m"      : This function gets the position and velocity of one body
                   relative to another body or gives nutations or librations.

  "state_eph.m"  : This function reads and interpolates the ephemeris at a
                   given julian ephemeris date and gives the heliocentric or
                   solar system barycentric state of the selected bodies
                   and/or nutations and librations.

  "init_eph.m"   : This function opens and reads the header of the selected
                   binary ephemeris file.

  "eph_global.m" : This function declares global variables used by the
                   ephemeris reading/interpolation functions.  You can have
                   access to these variables by typing "eph_global" at the
                   MATLAB command prompt.

  "eph_header_size.m": This function defines the dimensions of parameters in
                   the header section of a binary ephemeris file.

  "chebyval.m"   : This function evaluates Chebyshev polynomials.

  "chebyder.m"   : This function evaluates Chebyshev polynomial derivatives.
                      

  Ephemeris Toolbox Testing/Verification
  --------------------------------------

  "testeph.m"    : This function compares interpolation results with similar                                                         
                   runs made at JPL in order to ensure that the ephemeris is
                   installed and being read correctly.

                   This function exercises all of the ephemeris interpolation
                   functions.                   

  "testpo"       : Test results computed at JPL; these are input by the 
                   function "testeph" and are used for testing the ephemeris 
                   installation.  There is a different "testpo" for each 
                   different ephemeris; they must match or the test will not 
                   work correctly.

  "read_testpo.m": This function reads and parses "testpo".


  Binary Ephemeris File Creation
  ------------------------------

  "binmerge.m"   : Function to merge two adjoining binary ephemerides.

  "mergemall.m"  : Function to merge several adjoining binary ephemerides.

  "asc2bin_eph.m": A one-time conversion function which converts the ephemeris 
                   from ASCII format into binary form. 
                           
  "ascSYYYY.XXX" : ASCII ephemeris files from JPL Ephemeris DEXXX, each 
                   covering 20 years, starting in the year SYYYY ("p/m" for 
                   "+/-").  The 20-year blocks may be converted separately 
                   into binary ephemeris files using "asc2bin_eph". 
                   Subsequently, separate binary files may be merged into a
                   single ephemeris file using "binmerge.m" and "mergemall.m".

  "header.XXX"   : Header information for ephemeris deXXX, needed by
                   "asc2bin_eph.m".

It is suggested that you put the ".m" files in a separate directory (e.g.
"ephem") and make that directory part of your MATLAB path.

******************************************************************************

CREATING THE BINARY EPHEMERIS FILE
----------------------------------

The functions in this toolbox read and interpolate binary direct-read
ephemeris files.  These files must be created from the ASCII ephemeris files
described in the previous section.  This binary format is NOT the same as the
UNIX binary files available from the anonymous ftp site, even if you are using
a UNIX computer.

The format of the binary ephemeris file is NOT MATLAB's ".mat" format, so it
cannot be read using the "load" command.  It is meant to be directly accessed
by the functions in the toolbox.  Also, although the functions in the toolbox
should work in MATLAB on any type of computer, the binary ephemeris file
format is not the same on each type of computer.  That is you can use the
toolbox functions on any type of computer, but you must create separate binary
ephemeris files for each type of computer.

This process need only be performed one time.  Once the binary ephemeris files
are created, you may delete the ASCII ephemeris files.

The binary ephemeris file creation process may be completed using ephtool, a 
GUI-based tool for creating and using binary ephemeris files, or it may be
accomplished from the MATLAB command prompt as described below.

EPHTOOL METHOD

Type "ephtool" at the MATLAB command prompt to open the GUI.  Then click on the
"Create New File" button.  This will open a new window labeled "Create Binary
Ephemeris File."  If you have not already downloaded the ASCII ephemeris and
header files, click on "Download ASCII Ephemeris Files" to open a web browser
and download the files.  Next click on "Choose ASCII Ephemeris Files" which will
open a listing of files in the current directory.  Choose the files you wish to
use and click "OK."  Then click "Choose Ephemeris Header File" and do the same.
Then enter the binary ephemeris file name and click "Create Binary Ephemeris
File."

COMMAND PROMPT METHOD

To convert ASCII blocks into binary ephemeris files, use the function
"asc2bin_eph.m" as follows:

>> asc2bin_eph(HFILE,EFILES,BFILES)

where HFILE is the name of the ASCII header information file, "header.XXX";
EFILES is a cell array of ASCII ephemeris file names ("ascSYYYY.XXX"), and
BFILES is a cell array of desired binary output file names.

WARNING:  The ASCII ephemeris files use "D" as the exponent delimiter for
floating point numbers.  MATLAB does not recognize "D" as such, so the file
must be read line-by-line and parsed using the MATLAB function "fgetl.m"
rather than the faster "fscanf.m".  This is transparent to the user, but
unfortunately makes asc2bin_eph execute very slowly.  Fortunately, this step
is needed only one time.

The binary files created cover specific intervals of time.  They may be used
separately, or they may be merged together into a single large file.  To merge
two adjacent binary files into one file, use the function "binmerge.m" as
follows:

>> binmerge(FIN1,FIN2,FOUT)

where FIN1, and FIN2 are the names of the two binary ephemeris files you wish
to merge, and FOUT is the name for the merged binary ephemeris file.

Several adjacent files may be merged using the function "mergemall.m".  This
function calls "binmerge.m" successively to create a merged binary ephemeris
file.  Use "mergemall.m" as follows:

>> mergemall(FLIST,FOUT)

where FLIST is a cell array of the names of the binary ephemeris files you
wish to merge, and FOUT is the name desired for the merged binary ephemeris
file.

Once you have created merged binary ephemeris files, you may delete the
constituent binary ephemeris files if you wish.

******************************************************************************

VERIFYING PROPER INSTALLATION
-----------------------------

Once you have created the binary ephemeris file(s), verify proper creation of
these files and operation of the reading/interpolation functions using
"testeph.m" as follows:

>> testeph(EFILE,TFILE)

where EFILE is the name of the binary ephemeris file you wish to use and TFILE
is the name of the "testpo.XXX" file of JPL-generated results.

"testeph.m" computes ephemeris positions, velocities, nutations, and
librations and compares the results with those from "testpo.XXX".  Comparisons
will span the time span covered by the binary ephemeris file.  If any
comparison yields a difference larger than 1e-13 (in units of au or au/day),
an error message will be printed out.

Testing can also be accomplished using ephtool by clicking "Create New File"
and then "Test Binary Ephemeris File" and following the prompts.

WARNING:  Because some libration angles contain a large number of significant
digits, the JPL-generated value may be stored incorrectly by MATLAB, and
"testeph.m" may generate an error message.  If you see error messages only
involving ntarg=15, compare the "Computed Result" column against the actual
JPL-generated value in the "testpo.XXX" file.

******************************************************************************

USING THE TOOLBOX
-----------------

Detailed help for each function in the toolbox is available using the MATLAB
"help" utility:

>> help FUNCTION

where FUNCTION is the name of the function.

To use the reading/interpolation functions, "init_eph" must first be called to
open the ephemeris file, read its header, and define global variables:

>> init_eph(FNAME)

where FNAME is the name of the binary ephemeris file.

The functions of primary interest to the user are "pleph.m" and "state_eph.m".
Their help files are shown below:

PLEPH  Read JPL planetary ephemeris and give body states
   [R,V]=PLEPH(ET,NTARG,NCENT,KMFLAG) returns the position R
   and velocity V of body NTARG relative to body NCENT at Julian
   date ET.  The numbering convention for NTARG and NCENT is:

    1 -> Mercury          8 -> Neptune
    2 -> Venus            9 -> Pluto
    3 -> Earth           10 -> Moon
    4 -> Mars            11 -> Sun
    5 -> Jupiter         12 -> Solar system barycenter
    6 -> Saturn          13 -> Earth-Moon barycenter
    7 -> Uranus          14 -> Nutations (longitude and obliquity)
              15 -> Librations, if on file

   If nutations or librations are desired, set NCENT = 0.

   Units for position and velocity are AU and AU/day or, if KMFLAG
   is nonzero, km and km/s.  Units for nutations and librations
   are radians and radians/day.

------------------------------------------------------------------------------

STATE_EPH Read and interpolate JPL Planetary Ephemeris File
    STATE_EPH(JD,LIST,BARYFLAG,KMFLAG,EPHFILE) reads and
    interpolates a JPL planetary ephemeris file

    INPUTS: JD      Julian Date at desired interpolation epoch

            LIST    1-by-12 array indicating which bodies and what
                    type of interpolation is desired:
                    LIST(i) = 0 -> No interpolation for body i
                            = 1 -> Position only
                            = 2 -> Position and velocity

                    Designation of astronomical bodies
                     i = 1 -> Mercury
                       = 2 -> Venus
                       = 3 -> Earth-Moon barycenter
                       = 4 -> Mars
                       = 5 -> Jupiter
                       = 6 -> Saturn
                       = 7 -> Uranus
                       = 8 -> Neptune
                       = 9 -> Pluto
                       = 10-> Geocentric Moon
                       = 11-> Nutations in Earth longitude and 
                              obliquity
                       = 12-> Lunar librations (if on file)

           BARYFLAG  = 0 -> Planetary states are heliocentric
                    ~= 0 -> Planetary states are solar-system 
                            barycentric

           KMFLAG    = 0 -> Units are AU and AU/day
                    ~= 0 -> Units are km and km/s

           EPHFILE  Name of ephemeris file

    OUTPUTS: PV       3-by-13 vector of positions
             
             VV       3-by-13 vector of velocities

    States are relative to Earth mean equator and equinox of J2000
    if the DE number is >= 200; of B1950 if the DE number is < 200.

The GUI-based ephtool can also be used to get state data for solar system
bodies by pressing "Use Existing File" and following the prompts.

******************************************************************************

ASSISTANCE
----------

If you are really stuck, direct your questions to

Chuck Walker
VA-F3
Kennedy Space Center, FL  32899
321-476-3672 (Voice)
321-853-5528 (Fax)
mailto:charles.f.walker@nasa.gov

I shall try to answer your questions when I'm free from my normal obligations.
Please realize that I cannot provide customized service to each individual
user.

Please include your name, address, phone number and e-mail address.
