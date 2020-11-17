function[f_atm] = drag_force(dc, h, v_sate)
% input:
% 1. dc is the tabulated minimum and maximum density value (column 2 and column 3 (g/km^2))
%    at a given altitude h_i (column 1 (km))
% 2. h is the altitude of the satellite (km)
% 3. v_sate is the velocity of the satellite (3*1 vector) m/s`2
% output:
% f_atm is the atmospheric drag

% the atmosphere is radial symmetrical: rho(r,t) = rho(h)

% some constants in this job
dc(:,1) = dc(:,1) * 1000; % km to m
dc(:,2:3) = dc(:,2:3) / 1000^4; % g/km^3 to kg/m^3
AdM = 0.001;  % m`2/kg
C_d = 1; % no unit
v_atm = [0;0;0];

% find the precise altitude intervall
idx = find(h>dc(:,1));
if (h>1000000)
    rho_min = dc(50,2);
    rho_max = dc(50,3);
else
    if (dc(max(idx)) ~= h)
        idx = max(idx);
        % According to the Harris-Priester Density Model
        H_min = (dc(idx,1) - dc(idx+1,1)) / (log(dc(idx+1,2) / dc(idx,2)));
        H_max = (dc(idx,1) - dc(idx+1,1)) / (log(dc(idx+1,3) / dc(idx,3)));

        rho_min = dc(idx,2) * exp((dc(idx,1) - h) / H_min);
        rho_max = dc(idx,3) * exp((dc(idx,1) - h) / H_max);
    else
        idx = max(idx)+1;
        rho_min = dc(idx,2);
        rho_max = dc(idx,3);
    end
end


% A simplified estimation for atmospheric density
rho_atm = (rho_min + rho_max) / 2;

% atmospheric drag
f_atm = -1/2 * C_d * rho_atm * AdM * (v_sate - v_atm) * norm((v_sate - v_atm)); % m/s^2
end