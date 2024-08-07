function pp = generateInitDro()

%% Primary
Lsc     = 384748;                                                               % [km]   (1,1) Distance scaling constant
Tsc     = 375700.4314474463;                                                    % [s]    (1,1) Time scaling constant
Vsc     = Lsc/Tsc;
scale = [Lsc*ones(3,1); Vsc*ones(3,1)];
x0p   = [0.839663412542214;   0.012238694287822; 0; 
         0.025170783226701;   0.485027083812405; 0].*scale;
x(:,1) = x0p([1,2,4,5])./scale([1,2,4,5]);
for j = 2:300
    x(:,j) = propCr3bp2dNoCtrl(x(:,j-1),[0,0.01], 0.012153731135914043);
end
plot(x(1,:),x(2,:),'.')
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
primary.cart0  = primary.x0;
T              = 2.522*Tsc;
%% Secondary
x0s = [0.839663412542214;   0.012238694287822; 1.1e-07;
      -0.0013895;          -0.18394916;       -0.003].*scale;      
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
            'mu',        0.012153731135914043, ...
            'Lsc',       Lsc, ...
            'Vsc',       Vsc, ...
            'Tsc',       Tsc, ...
            'primary',   primary, ...
            'secondary', secondary, ...
            'n_conj',    1, ...
            'tca_sep',   0, ...
            'T',         T ...
            );
end