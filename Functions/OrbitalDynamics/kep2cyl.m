function cyl = kep2cyl(kep, mu)
    % Extract Keplerian elements
    a = kep.a;  % Semi-major axis
    ecc = kep.ecc;  % Eccentricity
    inc = kep.inc;  % Inclination
    w = kep.w; % Argument of periapsis
    RAAN = kep.RAAN; % Argument of periapsis
    f = kep.theta; % True anomaly
    
    % Compute intermediate variables
    p = a * (1 - ecc^2);
    
    % Compute cylindrical coordinates
    cyl.th = f + w;
    cyl.Th = sqrt(mu * p);
    cyl.nu = RAAN;
    cyl.Nu = cyl.Th * cos(inc); 
    cyl.r = p/(1+ecc*cos(f));
    cyl.R = cyl.Th/p*ecc*sin(f);
end
