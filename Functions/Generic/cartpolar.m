function [Rtransf] = cartpolar(Rold,coordtransf)

% File for coordinate transformation polar <-> cartesian
% cartesian -> polar
if (coordtransf == 1)
[Rtransf(2),Rtransf(3),Rtransf(1)] = cart2sph (Rold(1),Rold(2),Rold(3));
Rtransf(2) = qck(Rtransf(2));
elseif (coordtransf == 2)
% polar -> cartesian
[Rtransf(1),Rtransf(2),Rtransf(3)] = sph2cart (Rold(2),Rold(3),Rold(1));
else
    error('ERROR: not an available option for coordtransf in "cartpolar"!')
end