function pp = generateInitHalo()
% Test cases from De Maria's thesis

%% Primary
Lsc     = 384405;                                                        % [km]   (1,1) Distance scaling constant
Tsc     = 375677;                                                        % [s]    (1,1) Time scaling constant
Vsc     = Lsc/Tsc;
scale = [Lsc*ones(3,1); Vsc*ones(3,1)];
data =  table2array(readtable('./data/Halo_Jpl.csv'));
mu  = data(end);
x0p = data(2:7)'.*scale;
% x(:,1) = x0p;
% for j = 2:280
%     x(:,j) = propCr3bp(x(:,j-1),zeros(3,1),[0,0.01],mu);
% end
% plot3(x(1,:),x(2,:),x(3,:),'.')
primary.x0 = x0p;
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
primary.T     = data(9);
%% Secondary
relState = [0; 0; 0; -0.2; 0.1; 0].*scale;
x0s = x0p + relState;        

secondary.tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary.x0       = x0s;         
secondary.relState = [];                          % [km] [km/s] Relative cartesian state at TCA
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
            'mu',        mu, ...
            'Lsc',       Lsc, ...
            'Vsc',       Vsc, ...
            'Tsc',       Tsc, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'n_conj',    1, ...
            'tca_sep',   0, ...
            'T',         primary.T ...
            );
end