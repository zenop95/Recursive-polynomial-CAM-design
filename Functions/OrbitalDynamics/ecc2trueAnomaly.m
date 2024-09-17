function theta = ecc2trueAnomaly(E, e)
    theta = 2.0 * atan2(sqrt(1.0 + e^2) * sin(E / 2.0), sqrt(1.0 - e^2) * cos(E / 2.0));
end
