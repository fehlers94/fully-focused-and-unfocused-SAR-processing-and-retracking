T = 0.879

R=6371000;
phi0 = 0;
phi1 = V*T./(R+H);


MA = R*[sin(phi0) cos(phi0)]
MB = (R+H)*[sin(phi1) cos(phi1)]
AB = MB-MA

% normalize
ABn = AB./sqrt(sum(AB.^2))
MBn = MB./sqrt(sum(MB.^2))

% now the look angle
theta = acos(sum(ABn.*MBn))

antenna_pattern = exp(-gamma_x*theta.^2)