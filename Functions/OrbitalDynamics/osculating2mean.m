function meanKep = osculating2mean(kep, mu, Lsc)

    % Constants
    J2 = 1.08262668e-3; % J2 coefficient
    J4 = 0;%-1.620e-6;
    rE = 6371 / Lsc; % Radius of Earth normalized

    % Convert Keplerian elements to cylindrical coordinates
    cyl = kep2cyl(kep, mu);

    % Extract cylindrical elements
    r  = cyl.r;
    th = cyl.th;
    nu = cyl.nu;
    R  = cyl.R;
    Th = cyl.Th;
    Nu = cyl.Nu;

    % Compute intermediate variables
    p = Th^2 / mu;
    ci = Nu / Th;
    si = sqrt(1 - ci^2);
    cs = (p/r - 1.0) * cos(th) + (R * Th * sin(th)) / mu;
    ss = -((R * Th * cos(th)) / mu) + (p/r - 1.0) * sin(th);
    e = sqrt(cs^2 + ss^2);
    eta = sqrt(1 - e^2);

    beta = 1 / (1 + eta);

    costrue = 1 / e * (p / r - 1);

    f = acos(costrue);

    if R < 0
        f = 2 * pi - f;
    end

    M = true2meanAnomaly(f, e);

    phi = f - M;

    % Calculate mean elements
    rMean = r + ( (rE^2 * beta * J2) / (2 * r) - (3 * rE^2 * beta * J2 * si^2) / (4 * r) + ...
                  (rE^2 * eta * J2 * mu^2 * r) / (Th^4) - (3 * rE^2 * eta * J2 * mu^2 * r * si^2) / (2 * Th^4) + ...
                  (rE^2 * J2 * mu) / (2 * Th^2) - (rE^2 * beta * J2 * mu) / (2 * Th^2) - ...
                  (3 * rE^2 * J2 * mu * si^2) / (4 * Th^2) + (3 * rE^2 * beta * J2 * mu * si^2) / (4 * Th^2) - ...
                  (rE^2 * J2 * mu * si^2 * cos(2 * th)) / (4 * Th^2) ) + (15 * (rE^4) * J4 / (8 * r^5)) * ((7 * (si^4)) - (6 * (si^2)) + 1) * cos(4 * th);


    thMean = th + ( (-3 * rE^2 * J2 * mu^2 * phi) / (Th^4) + (15 * rE^2 * J2 * mu^2 * phi * si^2) / (4 * Th^4) - ...
                    (5 * rE^2 * J2 * mu * R) / (2 * Th^3) - (rE^2 * beta * J2 * mu * R) / (2 * Th^3) + ...
                    (3 * rE^2 * J2 * mu * R * si^2) / (Th^3) + (3 * rE^2 * beta * J2 * mu * R * si^2) / (4 * Th^3) - ...
                    (rE^2 * beta * J2 * R) / (2 * r * Th) + (3 * rE^2 * beta * J2 * R * si^2) / (4 * r * Th) + ...
                    ( -(rE^2 * J2 * mu * R) / (2 * Th^3) + (rE^2 * J2 * mu * R * si^2) / (Th^3) ) * cos(2 * th) + ...
                    ( -(rE^2 * J2 * mu^2) / (4 * Th^4) + (5 * rE^2 * J2 * mu^2 * si^2) / (8 * Th^4) + ...
                      (rE^2 * J2 * mu) / (r * Th^2) - (3 * rE^2 * J2 * mu * si^2) / (2 * r * Th^2) ) * sin(2 * th) ) + ...
                      (5 * (rE^4) * J4 * mu * R / (8 * Th^3 * r^2)) * (7 * (si^4) - 6 * (si^2) + 1) * sin(4 * th);

    nuMean = nu + ( (3 * rE^2 * ci * J2 * mu^2 * phi) / (2 * Th^4) + (3 * rE^2 * ci * J2 * mu * R) / (2 * Th^3) + ...
                     (rE^2 * ci * J2 * mu * R * cos(2 * th)) / (2 * Th^3) + ...
                     ( (rE^2 * ci * J2 * mu^2) / (4 * Th^4) - (rE^2 * ci * J2 * mu) / (r * Th^2) ) * sin(2 * th) ) + ...
                     (5 * (rE^4) * J4 * ci * mu * R / (8 * Th^3 * r^2)) * (7 * (si^4) - 6 * (si^2) + 1) * sin(4 * th);

    RMean = R + ( -(rE^2 * beta * J2 * R) / (2 * r^2) + (3 * rE^2 * beta * J2 * R * si^2) / (4 * r^2) - ...
                   (rE^2 * eta * J2 * mu^2 * R) / (2 * Th^4) + (3 * rE^2 * eta * J2 * mu^2 * R * si^2) / (4 * Th^4) + ...
                   (rE^2 * J2 * mu * si^2 * sin(2 * th)) / (2 * r^2 * Th) );

    ThMean = Th + ( (rE^2 * J2 * mu^2 * si^2) / (4 * Th^3) - (rE^2 * J2 * mu * si^2) / (r * Th) ) * cos(2 * th) - ...
             ( rE^2 * J2 * mu * R * si^2 * sin(2 * th) ) / (2 * Th^2);

    NuMean = Nu; % No perturbations for normal velocity

    % Convert mean cylindrical coordinates to mean Keplerian elements
    meanCyl.r  = rMean;
    meanCyl.th = thMean;
    meanCyl.nu = nuMean;
    meanCyl.R  = RMean;
    meanCyl.Th = ThMean;
    meanCyl.Nu = NuMean;

    meanKep = cyl2kep(meanCyl, mu);
end


