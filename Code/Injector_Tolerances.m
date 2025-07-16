% ---------------------------------- %
% -------INJECTOR DESIGN CODE------- %
% ---------------------------------- %

%{
Copyright (c) 2025 Sam Owen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
%}

clc;
clear;
close all;

% Initial Geometry
rn_init = 2.3311e-3;
rv_init = 8.4094e-3;
rt_init = 1.4161e-3;
n = 6; % Number of tangential inlets
mdot_target = 0.25; % Target mass flow rate (kg/s)
rho = 998; % Fluid density (kg/m^3)
xi = 0.5; % Hydraulic loss coefficient

% Tolerances (Update Based on Requirements)
deltap_tolerance = 0.5; % bar
alphae_tolerance = 5; % degrees
eta_min = 2.99; % Minimum 
eta_max = 3; % Maximum

[~, deltap_target, alphae_target, eta_target] = injector_performance(rn_init, rv_init, rt_init, n, mdot_target, rho, xi);

fprintf('Target Pressure Drop: %.4f bar\n', deltap_target);
fprintf('Target Cone Angle: %.2f degrees\n', alphae_target);
fprintf('Target Eta: %.2f\n\n', eta_target);

% Geometry Range (Update Based on Requirements)
rn_range = (1:0.05:5) * 1e-3;
rv_range = (4:0.05:10) * 1e-3;
rt_range = (1:0.05:3) * 1e-3;

fprintf('Searching for valid geometries...\n');
found = false;

for rn = rn_range

    for rv = rv_range

        for rt = rt_range

            if rv - rt <= 0 || rt <= 0

                continue;

            end

            try

                [mu, deltap, alphae, eta] = injector_performance(rn, rv, rt, n, mdot_target, rho, xi);

                if abs(deltap - deltap_target) <= deltap_tolerance && ...
                   abs(alphae - alphae_target) <= alphae_tolerance && ...
                   eta >= eta_min && eta <= eta_max

                    fprintf('-----Valid Geometry Found-----\n');
                    fprintf('rn = %.2f mm, rv = %.2f mm, rt = %.2f mm\n', rn*1e3, rv*1e3, rt*1e3);
                    fprintf('Pressure Drop = %.4f bar\n', deltap);
                    fprintf('Spray Cone Angle = %.2f degrees\n', alphae);
                    fprintf('Eta = %.2f\n\n', eta);
                    found = true;

                end

            catch

                continue;

            end

        end

    end

end

if ~found

    fprintf('No matching geometries found within given tolerances.\n');
    
end


function [mu, deltap, alphae, eta] = injector_performance(rn, rv, rt, n, mdot, rho, xi)

    rin = rv - rt;
    eta = rin / rn;
    A = (rn * rin) / (n * rt^2);
    phi_guess = 0.8;
    phi = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A, phi_guess);
    mu = sqrt(((phi^3) * (eta^2)) / ((2 - phi) * (eta^2) + (xi * A^2 * phi^3)));
    mu0 = sqrt((phi^3) / (2 - phi));
    deltap = (1 / (2 * rho)) * (mdot / (pi * rn^2 * mu))^2 * 1e-5; 
    alphae = atan((mu0 * A) / sqrt(1 - ((mu0^2) * (A^2)))) * (360 / pi); 

end
