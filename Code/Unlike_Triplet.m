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
Mdoto = 2.0; % Oxidiser Mass Flow Rate (kg/s)

% Configuration: 0 for Oxidiser Centred, 1 for Fuel Centred

config = 1;

%{
Here Tilt Angle is Defined As the Angle of the Injector Stream From the Vertical
%}

tot = 60; % Collision Angle (degrees)
D = 5; % Impinging Distance Ratio to Average Orifice Diameter (3-7)

%{ 
Lower D Values Tend to Increase Performance at the Risk of Melting the Injector 
%}

% Number of Oxidiser Elements (Even if config = 1)
 
no = 20;

% Fluid Parameters

rhof = 998; % Fuel Density (kg/m^3)
rhoo = 998; % Oxidiser Density (kg/m^3)
deltapf = 8e5; % Fuel Pressure Drop (Pa)
deltapo = 8e5; % Oxidiser Pressure Drop (Pa)
cdf = 0.7; % Fuel Discharge Coefficient
cdo = 0.7; % Oxidiser Discharge Coefficient

if config == 0

    f_tilt = tot / 2; % Fuel Tilt Angle (degrees)
    o_tilt = 0; % Oxidiser Tilt Angle (degrees)
    nf = no * 2;
    fprintf('-----OXIDISER CENTRED-----\n\n')

elseif config == 1 

    f_tilt = 0; % Fuel Tilt Angle (degrees)
    o_tilt = tot / 2; % Oxidiser Tilt Angle (degrees)
    nf = no / 2;
    fprintf('-----FUEL CENTRED-----\n\n')

end

mdotf = Mdotf / nf;
mdoto = Mdoto / no;

vf = cdf * sqrt(2 * deltapf / rhof);
vo = cdo * sqrt(2 * deltapo / rhoo);

rnf = sqrt(mdotf / (pi * cdf * sqrt(2 * rhof * deltapf)));
rno = sqrt(mdoto / (pi * cdo * sqrt(2 * rhoo * deltapo)));

fprintf('-----UNLIKE TRIPLETS-----\n\n')
fprintf('Fuel Element Radius: %.4f mm\n', rnf * 1e3);
fprintf('Fuel Element Number: %.4f \n\n', nf);
fprintf('Oxidiser Element Radius: %.4f mm\n', rno * 1e3);
fprintf('Oxidiser Element Number: %.4f \n\n', no);

if config == 0

    d = D * (2 / 3) * (rnf + rnf + rno);
    do = d * (tan(f_tilt*pi/180) + tan(f_tilt*pi/180));
    fprintf('Impingement Distance: %.4f mm\n', d * 1e3);
    fprintf('F-F Orifice Seperation Distance: %.4f mm\n', do * 1e3);
    fprintf('Fuel Tilt Angle: %.4f degrees\n', f_tilt);
    fprintf('Oxidiser Tilt Angle: %.4f degrees\n\n', o_tilt);

    a_rat = (rno / rnf) ^ 2;
    MR = (rhoo * vo ^ 2 * rno ^ 2) / (rhof * vf ^ 2 * rnf ^ 2);
    fprintf('O-F Area Ratio: %.4f \n\n', a_rat);
    fprintf('O-F Momentum Ratio: %.4f \n\n', MR);
    fprintf('Fuel Velocity: %.4f m/s\n', vf);
    fprintf('Oxidiser Velocity: %.4f m/s\n\n', vo);

elseif config == 1

    d = D * (2 / 3) * (rnf + rno + rno);
    do = d * (tan(o_tilt*pi/180) + tan(o_tilt*pi/180));
    fprintf('Impingement Distance: %.4f mm\n', d * 1e3);
    fprintf('O-O Orifice Seperation Distance: %.4f mm\n', do * 1e3);
    fprintf('Fuel Tilt Angle: %.4f degrees\n', f_tilt);
    fprintf('Oxidiser Tilt Angle: %.4f degrees\n\n', o_tilt);

    a_rat = (rno / rnf) ^ 2;
    MR = (rhoo * vo ^ 2 * rno ^ 2) / (rhof * vf ^ 2 * rnf ^ 2);
    fprintf('O-F Area Ratio: %.4f \n\n', a_rat);
    fprintf('O-F Momentum Ratio: %.4f \n\n', MR);
    fprintf('Fuel Velocity: %.4f m/s\n', vf);
    fprintf('Oxidiser Velocity: %.4f m/s\n\n', vo);

end