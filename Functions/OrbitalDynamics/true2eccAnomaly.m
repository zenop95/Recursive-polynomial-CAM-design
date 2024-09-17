function E = true2eccAnomaly(theta, e)
    E = 2.0 * atan2(sqrt(1.0 - e^2) * sin(theta / 2.0), sqrt(1.0 + e^2) * cos(theta / 2.0));
end
