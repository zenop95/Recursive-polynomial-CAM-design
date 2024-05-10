function pp = generateInitCislunar(conj)
% Test cases from De Maria's thesis

%% Primary
Lsc     = 384405;                                                        % [km]   (1,1) Distance scaling constant
Tsc     = 375677;                                                        % [s]    (1,1) Time scaling constant
Vsc     = Lsc/Tsc;
scale = [Lsc*ones(3,1); Vsc*ones(3,1)];
x0p = [0.877951855; 0; -0.192527194; 0; 0.227333284; 0].*scale;
primary.x0 = x0p;
% c12 = 0.56662; c13 = 0.6678;  c14 = 5.9149e-6;  c15 = 6.7702e-7;  c16 = 3.6882e-6;
%                c23 = 0.50654; c24 = 4.1301e-6;  c25 = 6.5091e-7;  c26 = 2.5187e-6;
%                               c34 = 5.3721e-6;  c35 = 7.0304e-7;  c36 = 3.5061e-6;
%                                                 c45 = 5.0560e-12; c46 = 3.1068e-11;  
%                                                                   c56 = 1.9561e-12;
% primary.C0 = [0.76741        c12        c13         c14      c15          c16;
%                   c12     0.46737       c23         c24      c25          c26;
%                   c13        c23     0.66121        c34      c35          c36;
%                   c14        c24        c34   4.8742e-11     c45          c46;
%                   c15        c25        c35         c45    5.8238e-12     c56;
%                   c16        c26        c36         c46      c56    2.0591e-11];
primary.C0 = [[5.5e-2 2.5e-3 1.6e-3;
               2.5e-3 2.8e-2 1.4e-4;
               1.6e-3 1.4e-4 2.8e-2] zeros(3,3); zeros(3,6)];

primary.HBR    = 0.01;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0 = primary.x0;

%% Secondary
if strcmpi(conj,'parallel')
    x0s = [337790.22030754; 0; -74074.495159392; 
       0.00130390551665015; 0.0370542702726203; 0.000488964568743807];          % parallel case
else
    x0s = [337489.082821275;
                         0;
           -74008.45829412;
       0.00818586179084693;
         0.232625087737658;
        0.0030696981715676];        % perpendicular case
end
secondary.tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary.x0       = x0s;         
secondary.relState = [];                          % [km] [km/s] Relative cartesian state at TCA
% c12 = 0.69571; c13 = 0.70523;  c14 = 9.0699e-6;  c15 = -4.279e-7;   c16 = 4.883e-6;
%                c23 = 0.39546;  c24 = 4.7575e-6;  c25 = -2.8035e-7;  c26 = 2.5689e-6;
%                                c34 = 4.8979e-6;  c35 = -3.9782e-7;  c36 = 3.2715e-6;
%                                                  c45 = -2.8578e-12; c46 = 3.4489e-11;  
%                                                                     c56 = -7.0324e-13;
% secondary.C0 = [1.3049        c12        c13         c14      c15          c16;
%                     c12     0.40134      c23         c24      c25          c26;
%                     c13        c23     0.49791       c34      c35          c36;
%                     c14        c24        c34   6.3409e-11    c45          c46;
%                     c15        c25        c35         c45    2.4865e-12    c56;
%                     c16        c26        c36         c46     c56    2.6072e-11];
secondary.C0 = [[3.7e-2 -5e-3 1.2e-3;
               -5e-3 2.8e-2 -6.2e-4;
               1.2e-3 -6.2e-4 4.3e-2] zeros(3,3); zeros(3,6)];


secondary.HBR        = 0.01 + primary.HBR;         % [km]
secondary.mass       = 100;          % [kg] mass
secondary.A_drag     = 1;            % [m^2] drag surface area
secondary.Cd         = 2.2;          % [-] shape coefficient for drag
secondary.A_srp      = 1;            % [m^2] SRP surface area
secondary.Cr         = 1.31;         % [-] shape coefficient for SRP
secondary.x          = [];         % [-] 
secondary.covariance = [];         % [-] 
secondary.w          = 1;        
secondary.cdm        = true;        
secondary.ang        = false;        
%% Relative initial covariance in ECI reference frame

pp = struct( ...
            'mu',        403504.497, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'n_conj',    1, ...
            'tca_sep',   0, ...
            'T',         1 ...
            );
end