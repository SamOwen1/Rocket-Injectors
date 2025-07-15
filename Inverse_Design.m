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

% Inputs

rn = 2.3311e-3; % Nozzle Radius (m)
rv = 8.4094e-3; % Vortex Chamber Radius (m)
rt = 1.4161e-3; % Tangential Inlet Radius (m)
n = 6; % Number of Tangential Inlets
mdot = 0.25; % Mass Flow Rate (kg/s)
rho = 998; % Fluid density (kg/m^3)
xi = 0.5; % Hydraulic Loss Coefficient

% Outputs

rin = rv - rt; 
eta = rin / rn;
A = (rn * rin) / (n * rt ^ 2);
phi_guess = 0.8;
phi = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A, phi_guess);
mu = sqrt(((phi ^ 3) * (eta ^ 2)) / ((2 - phi) * (eta ^ 2) + (xi * A ^ 2 * phi ^ 3)));
mu0 = sqrt((phi ^ 3) / (2 - phi));
deltap = (1 / (2 * rho)) * (mdot / (pi * rn ^ 2 * mu)) ^ 2;
alphae = atan((mu0 * A) / sqrt(1 - ((mu0 ^ 2) * (A ^ 2)))) * (360 / pi);

fprintf('Pressure Drop: %4f bar\n', deltap * 1e-5);
fprintf('Spray Cone Angle: %4f degrees\n', alphae);
fprintf('Coefficient of Nozzle Opening: %4f \n', eta);
fprintf('Geometric Constant: %4f \n', A);
fprintf('Filling Efficiency: %4f \n', phi);
fprintf('Mass Flow Coefficient: %4f \n', mu);