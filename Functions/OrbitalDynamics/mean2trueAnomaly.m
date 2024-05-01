function TA = mean2trueAnomaly(n,e,M,t)
% mean2trueAnomaly computes the true anomaly for a given mean anomaly in  
% the elliptic orbit
% INPUT: n [rad/s] = Mean motion of the orbit
%        e [-]     = eccentricity of the orbit
%        M [rad]   = mean anomaly
% 
% OUTPUT: TA [rad]  = True anomaly of the given time from periapsis

%Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
T  = 2*pi/n;                              %[s] orbital period
E  = fzero(@(E) M-E+e*sin(E), 0);         %[rad] eccentric anomaly;
TA = floor(t/T)*2*pi + 2*atan(sqrt((1+e)/(1-e))*tan(E/2));  %[rad] true anomaly;
end

