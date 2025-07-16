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

% Total Flow Rates

Mdotf = 1.0; % Fuel Mass Flow Rate (kg/s)
Mdoto = 2.5; % Oxidiser Mass Flow Rate (kg/s)

%{
Here Tilt Angle is Defined As the Angle of the Injector Stream From the Vertical
%}

tot = 60; % Collision Angle (degrees)
f_tilt = 37.7784; % Fuel Tilt Angle (degrees)
o_tilt = tot - f_tilt; % Oxidiser Tilt Angle (degrees)
D = 5; % Impinging Distance Ratio to Average Orifice Diameter (3-7)

%{ 
Lower D Values Tend to Increase Performance at the Risk of Melting the Injector 
%}

% Number of Elements

nf = 20; 
no = nf;

% Fluid Parameters

rhof = 998; % Fuel Density (kg/m^3)
rhoo = 998; % Oxidiser Density (kg/m^3)
deltapf = 8e5; % Fuel Pressure Drop (Pa)
deltapo = 5e5; % Oxidiser Pressure Drop (Pa)
cdf = 0.7; % Fuel Discharge Coefficient
cdo = 0.7; % Oxidiser Discharge Coefficient

mdotf = Mdotf / nf;
mdoto = Mdoto / no;

rnf = sqrt(mdotf / (pi * cdf * sqrt(2 * rhof * deltapf)));
rno = sqrt(mdoto / (pi * cdo * sqrt(2 * rhoo * deltapo)));

fprintf('-----UNLIKE DOUBLETS-----\n\n')
fprintf('Fuel Element Radius: %.4f mm\n', rnf * 1e3);
fprintf('Fuel Element Number: %.4f \n\n', nf);
fprintf('Oxidiser Element Radius: %.4f mm\n', rno * 1e3);
fprintf('Oxidiser Element Number: %.4f \n\n', no);

vf = cdf * sqrt(2 * deltapf / rhof);
vo = cdo * sqrt(2 * deltapo / rhoo);

beta = atan((mdotf*vf*sin(f_tilt*pi/180)-mdoto*vo*sin(o_tilt*pi/180))/(mdotf*vf*cos(f_tilt*pi/180)+mdoto*vo*cos(o_tilt*pi/180)));
s = 0;

fprintf('Resultant Angle: %.4f degrees\n\n', beta * 180 / pi);

if beta < -1e-5

    angle = @(alpha) (mdotf*vf*sind(alpha)-mdoto*vo*sind(tot-alpha))/(mdotf*vf*cosd(alpha)+mdoto*vo*cosd(tot-alpha));
    sol = fzero(angle, f_tilt);
    fprintf('Resultant Angle Travels Towards FUEL Side\n')
    fprintf('Increase FUEL Tilt to %.4f degrees\n\n', sol)
    s = 1;

elseif beta > 1e-5

    angle = @(alpha) (mdotf*vf*sind(alpha)-mdoto*vo*sind(tot-alpha))/(mdotf*vf*cosd(alpha)+mdoto*vo*cosd(tot-alpha));
    sol = fzero(angle, f_tilt);
    fprintf('Resultant Angle Travels Towards OXIDISER Side\n')
    fprintf('Decrease FUEL Tilt to %.4f degrees\n\n', sol)
    s = 1;

end

d = D * (rnf + rno);
do = d * (tan(f_tilt*pi/180) + tan(o_tilt*pi/180));
fprintf('Impingement Distance: %.4f mm\n', d * 1e3);
fprintf('Orifice Seperation Distance: %.4f mm\n', do * 1e3);
fprintf('Fuel Tilt Angle: %.4f degrees\n', f_tilt);
fprintf('Oxidiser Tilt Angle: %.4f degrees\n\n', o_tilt);

a_rat = (rnf / rno) ^ 2;
a_rat_opt = ((rhoo / rhof) * (mdotf / mdoto) ^ 2) ^ 0.7;
fprintf('F-O Area Ratio: %.4f \n', a_rat);
fprintf('OPTIMUM F-O Area Ratio: %.4f \n\n', a_rat_opt);

if abs(a_rat_opt - a_rat) > (0.1 * a_rat_opt)

    fprintf('-----DESIGN INVALID-----\n\n')
    fprintf('Area Ratio Not Optimal Within 10 Percent\n')
    fprintf('Try Triplets or Change Pressure Drop\n')

else

    if s == 0

    fprintf('-----DESIGN VALID-----\n')

    end

end
