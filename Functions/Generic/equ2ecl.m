%Programmed by:         Dario Izzo
%                       Advanced Concepts Team
%Date:                  28/05/2004
%Revision:              1
%Tested by:             ----------
%
%
%Converts a column vector given in the equatorial coordinate frame to the
%earth ecliptic coordinate frame
%
%Usage:     recl=equ2ecl(requ)
%
%Inputs:    requ=coordinates in the equatorial J2000 frame
%
%
%Output:    recl=coordinates in the ecliptic frame
%
%X-ref:  none

function recl=equ2ecl(requ)

global incl %Earth axis inclination used for the 'ecl' options


    RR=[1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
    recl=RR*requ;


