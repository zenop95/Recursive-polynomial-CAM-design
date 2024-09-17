function M = true2meanAnomaly(theta, e)
    E = true2eccAnomaly(theta, e);
    M = E - e * sin(E);
end
