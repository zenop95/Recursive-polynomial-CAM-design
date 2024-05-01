function pp = generateInitSerra(orbit)

mu     = 398600.4418;    % [km^3/s^2]
HBR    = .006;           % [km]
n      = 1.4591e-4;
e      = 0.741;          %[-] eccentricity
theta0 = -3.071;         %[rad] initial true anomaly
omega  = 0;              %[rad] argument of periapsis
RAAN   = 0;              %[rad] right ascension node longitude
inc    = 0;              %[rad] inclination
a      = (mu/n^2)^(1/3); %[-] semimajor axis
T      = 2*pi/n;         %[s] orbital period

%% Relative initial conditions in Hill reference frame
r0 = [.04683;     0;  -.002986];
v0 = [.000000643; 0; .000001922];
c12 = 0; c13 = 2.5144;  c14 = 2.607e-3;  c15 = 0;         c16 = 8.0097e-3;
         c23 = 0;       c24 = 0;         c25 = -.0024e-3; c26 = 0;
                        c34 = 0.1574e-3; c35 = 0;         c36 = 0.4794e-3;
                                         c45 = 0;         c46 = 0.4928e-6;  
                                                          c56 = 0;
C0 = [42.9181        c12        c13         c14      c15          c16;
          c12     0.0649        c23         c24      c25          c26;
          c13        c23     0.1549         c34      c35          c36;
          c14        c24        c34   0.1613e-6      c45          c46;
          c15        c25        c35         c45    1e-10          c56;
          c16        c26        c36         c46      c56    1.5083e-6]/1e6;

%% Inertial state of the two objects
%inertial state of the satellite
x_p     = kepler2cartesian(a,e,RAAN,inc,omega,theta0,mu);
r_p     = x_p(1:3);
v_p     = x_p(4:6);
[l2ep,wp] = lvlh2eci(r_p,v_p);
dr      = l2ep*r0;
dv      = l2ep*v0;
r_s     = r_p - dr;
v_s     = v_p - dv;
x_s     = [r_s; v_s];
    
[l2es,ws] = lvlh2eci(r_s,v_s);
Dp        = [l2ep zeros(3); skew(wp) l2ep];
Ds        = [l2es zeros(3); skew(ws) l2es];
C0p = Dp'*C0*Dp/2;
C0s = Ds'*C0*Ds/2;

pp = struct( ...
            'orbit', orbit, ...
            'mu',    mu, ...
            'HBR',   HBR, ...
            'cart_s',   x_p, ...
            'cart_d',   x_s, ...
            'C0s',    C0p, ...
            'C0p',    C0s, ...
            'T',     T, ...
            'n',     n ...
            );
end