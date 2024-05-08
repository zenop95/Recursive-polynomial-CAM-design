%Programmed by:         Dario Izzo
%                       Advanced Concepts Team
%Date:                  28/05/2004
%Revision:              1
%Tested by:             ----------
%
%
%Converts a column vector given in the ecliptic coordinate frame to the
%earth equatorial coordinate frame
%
%Usage:     requ=ecl2equ(recl)
%
%Inputs:    recl=coordinates in the ecliptic J2000 frame
%
%
%Output:    requ=coordinates in the equatorial frame
%
%X-ref:  none

function requ=ecl2equ(recl)

global incl %Earth axis inclination used for the 'ecl' options


    RR=[1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
    requ=RR'*recl;


