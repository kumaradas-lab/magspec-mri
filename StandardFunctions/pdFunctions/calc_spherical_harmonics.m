function orders = calc_spherical_harmonics(XX, YY, ZZ, r, lmax)

%% calculate RR, TTHETA and PPHI

RR = sqrt (XX.^2 + YY.^2 + ZZ.^2);
TTHETA = acos (ZZ ./ RR); % [0..pi]
PPHI = atan2 (YY, XX); % [-pi..pi[
RR = RR/r;

%% l=1
l = 1;
if lmax >= 1
    P1m = legendre (1, cos (TTHETA));
    
    % m = -1
    orders{l}{1} = RR .* squeeze (P1m(2,:,:,:)) .* sin (PPHI); % YY
    % m = 0
    orders{l}{2} = RR .* squeeze (P1m(1,:,:,:)); % ZZ
    % m = 1
    orders{l}{3} = RR .* squeeze (P1m(2,:,:,:)) .* cos (PPHI); % XX
end

%% l=2
l = 2;
if lmax >= l
    Plm = legendre (l, cos (TTHETA));
    
    % m = -2
    orders{l}{1} = RR.^l .* squeeze (Plm(3,:,:,:)) .* sin (2*PPHI);
    % m = -1
    orders{l}{2} = RR.^l .* squeeze (Plm(2,:,:,:)) .* sin (PPHI);
    % m = 0
    orders{l}{3} = RR.^l .* squeeze (Plm(1,:,:,:));
    % m = 1
    orders{l}{4} = RR.^l .* squeeze (Plm(2,:,:,:)) .* cos (PPHI);
    % m = 2
    orders{l}{5} = RR.^l .* squeeze (Plm(3,:,:,:)) .* cos (2*PPHI);
end

%% l=3
l = 3;
if lmax >= l
    Plm = legendre (l, cos (TTHETA));
    
    % m = -3
    orders{l}{1} = RR.^l .* squeeze (Plm(4,:,:,:)) .* sin (3*PPHI);
    % m = -2
    orders{l}{2} = RR.^l .* squeeze (Plm(3,:,:,:)) .* sin (2*PPHI);
    % m = -1
    orders{l}{3} = RR.^l .* squeeze (Plm(2,:,:,:)) .* sin (PPHI);
    % m = 0
    orders{l}{4} = RR.^l .* squeeze (Plm(1,:,:,:));
    % m = 1
    orders{l}{5} = RR.^l .* squeeze (Plm(2,:,:,:)) .* cos (PPHI);
    % m = 2
    orders{l}{6} = RR.^l .* squeeze (Plm(3,:,:,:)) .* cos (2*PPHI);
    % m = 3
    orders{l}{7} = RR.^l .* squeeze (Plm(4,:,:,:)) .* cos (3*PPHI);
end

%% l=4
l = 4;
if lmax >= l
    Plm = legendre (l, cos (TTHETA));
    
    % m = -4
    orders{l}{1} = RR.^l .* squeeze (Plm(5,:,:,:)) .* sin (4*PPHI);
    % m = -3
    orders{l}{2} = RR.^l .* squeeze (Plm(4,:,:,:)) .* sin (3*PPHI);
    % m = -2
    orders{l}{3} = RR.^l .* squeeze (Plm(3,:,:,:)) .* sin (2*PPHI);
    % m = -1
    orders{l}{4} = RR.^l .* squeeze (Plm(2,:,:,:)) .* sin (PPHI);
    % m = 0
    orders{l}{5} = RR.^l .* squeeze (Plm(1,:,:,:));
    % m = 1
    orders{l}{6} = RR.^l .* squeeze (Plm(2,:,:,:)) .* cos (PPHI);
    % m = 2
    orders{l}{7} = RR.^l .* squeeze (Plm(3,:,:,:)) .* cos (2*PPHI);
    % m = 3
    orders{l}{8} = RR.^l .* squeeze (Plm(4,:,:,:)) .* cos (3*PPHI);
    % m = 4
    orders{l}{9} = RR.^l .* squeeze (Plm(5,:,:,:)) .* cos (4*PPHI);
end

%% l=5
l = 5;
if lmax >= l
    Plm = legendre (l, cos (TTHETA));
    
    % m = -5
    orders{l}{1} = RR.^l .* squeeze (Plm(6,:,:,:)) .* sin (5*PPHI);
    % m = -4
    orders{l}{2} = RR.^l .* squeeze (Plm(5,:,:,:)) .* sin (4*PPHI);
    % m = -3
    orders{l}{3} = RR.^l .* squeeze (Plm(4,:,:,:)) .* sin (3*PPHI);
    % m = -2
    orders{l}{4} = RR.^l .* squeeze (Plm(3,:,:,:)) .* sin (2*PPHI);
    % m = -1
    orders{l}{5} = RR.^l .* squeeze (Plm(2,:,:,:)) .* sin (PPHI);
    % m = 0
    orders{l}{6} = RR.^l .* squeeze (Plm(1,:,:,:));
    % m = 1
    orders{l}{7} = RR.^l .* squeeze (Plm(2,:,:,:)) .* cos (PPHI);
    % m = 2
    orders{l}{8} = RR.^l .* squeeze (Plm(3,:,:,:)) .* cos (2*PPHI);
    % m = 3
    orders{l}{9} = RR.^l .* squeeze (Plm(4,:,:,:)) .* cos (3*PPHI);
    % m = 4
    orders{l}{10} = RR.^l .* squeeze (Plm(5,:,:,:)) .* cos (4*PPHI);
    % m = 5
    orders{l}{11} = RR.^l .* squeeze (Plm(6,:,:,:)) .* cos (5*PPHI);
end

end
