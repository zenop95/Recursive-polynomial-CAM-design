function kep = cyl2kep(cyl, mu)
    % Extract cylindrical coordinates
    r = cyl.r;
    th = cyl.th;
    nu = cyl.nu;
    R = cyl.R;
    Th = cyl.Th;
    Nu = cyl.Nu;
    
    % Compute intermediate variables
    i = acos(Nu / Th);
    cs = (-1 + Th^2 / (mu * r)) * cos(th) + (R * Th * sin(th)) / mu;
    ss = - (R * Th * cos(th)) / mu + (-1 + Th^2 / (mu * r)) * sin(th);
    e = sqrt(cs^2 + ss^2);
    p = Th^2 / mu;
    costrue = 1 / e * (p / r - 1);
    f = acos(costrue);
    
    if R < 0
        f = 2 * pi - f;
    end
    
    % Compute semi-major axis
    a = p / (1 - e^2);
    
    
    % Compute Keplerian elements
    kep.a = a;     % Semi-major axis
    kep.ecc = e;     % Eccentricity
    kep.inc = i;     % Inclination
    kep.w = th - f; % True anomaly
    if kep.w > 2*pi; kep.w = kep.w - 2*pi; end % argument of perigee
    kep.RAAN = nu;    % Argument of periapsis
    kep.theta = f;     % True anomaly
    kep.ex = e*cos(kep.w);     % Ecc vector x
    kep.ey = e*sin(kep.w);     % Ecc vector y
end
