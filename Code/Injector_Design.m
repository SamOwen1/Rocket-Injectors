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

% Choose Configuration: 0 for Open, 1 for Closed, 2 for Jet, 3 for Absent

op_cl1 = 1; % Inner Element Type
op_cl2 = 1; % Outer Element Type 

% Choose Propellent Configuration: 0 for Oxidiser Centred, 1 for Fuel Centred

prop = 0;

% Choose Coefficient of Nozzle Opening: <1 for Open, >1 for Closed

eta1 = 3; % Inner Element Coefficient of Nozzle Opening
eta2 = 2; % Outer Element Coefficient of Nozzle Opening

% Choose Mixing Type: 0 for Internal, 1 for External, 2 for Impinging

mix = 0; % Mixing Type

% Mixing Time ~0.2ms for Non-Hypergolic, ~0.1ms for Hypergolic

tao = 0.1e-3; % Valid for 0.2 < (mdot1 + mdot2) < 1.0 (kg/s)

% Input Fluid Parameters

mdot1 = 0.25; % Inner Element Mass Flow Rate (kg/s)
mdot2 = 0.10; % Outer Element Mass Flow Rate (kg/s)
deltap1 = 8e5; % Inner Element Pressure Drop (Pa)
deltap2 = 8e5; % Outer Element Pressure Drop (Pa)
rho1 = 998; % Inner Element Fluid Density (kg/m^3)
rho2 = 998; % Outer Element Fluid Density (kg/m^3)
n1 = 6; % Inner Element Number of Tangential Inlets
n2 = 4; % Outer Element Number of Tangential Inlets
xi1 = 0.5; % Inner Element Hydraulic Loss Coefficient
xi2 = 0.5; % Outer Element Hydraulic Loss Coefficient
wall = 1e-3; % Inner Element Wall Thickness (m)

%{
Choose Inner Element Spray Angle: 60-80 for Coaxial Internally Mixed
Larger for Simplex Inner Swirler or for External Mixing as alpha2<alpha1
%}

alpha1 = 60; % Inner Element Spray Angle
d = 7.5; % (alpha1 - alpha2) for mix == 1

cd1 = 0.65; % Jet Discharge Coefficient (0.60-0.85)
cd2 = 0.60; % Jet Discharge Coefficient (0.60-0.85)

% Calculations

if op_cl1 == 1

    phi_eq1 = @(phi1) ((sqrt(2) * (1 - phi1)) / sqrt(2 - phi1)) - sin(alpha1 * pi / 360);
    phi1 = fzero(phi_eq1, 0.8); 
    A1 = (sqrt(2) * (1 - phi1)) / (phi1 * sqrt(phi1)); 
    mu1 = sqrt(((phi1 ^ 3) * (eta1 ^ 2)) / ((2 - phi1) * (eta1 ^ 2) + (xi1 * A1 ^ 2 * phi1 ^ 3)));  
    rn1 = sqrt(mdot1 / (pi * mu1 * sqrt(2 * rho1 * deltap1)));
    rin1 = rn1 * eta1;
    rt1 = sqrt((rn1 * rin1) / (n1 * A1));
    rv1 = rin1 + rt1;
    rn1o = rn1 + wall;

elseif op_cl1 == 0

    A1 = eta1 / (n1 * (1 - eta1) ^ 2); 
    phi1 = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A1, 0.8); 
    mu1 = sqrt(((phi1 ^ 3) * (eta1 ^ 2)) / ((2 - phi1) * (eta1 ^ 2) + (xi1 * A1 ^ 2 * phi1 ^ 3)));
    alpha1 = (360 / pi) * asin(((sqrt(2) * (1 - phi1)) / sqrt((2 - phi1))));
    rn1 = sqrt(mdot1 / (pi * mu1 * sqrt(2 * rho1 * deltap1)));
    rin1 = rn1 * eta1;
    rt1 = sqrt((rn1 * rin1) / (n1 * A1));
    rv1 = rin1 + rt1;
    rn1o = rn1 + wall;

elseif op_cl1 == 2

    rn1 = sqrt(mdot1 / (pi * cd1 * sqrt(2 * rho1 * deltap1)));
    rn1o = rn1 + wall;
    alpha1 = 0;

end

% Spray Angle for Internal Mixing is 30-40 Degrees Less Than Isolated alpha2

if op_cl2 == 1 && (mix == 0 || mix == 2)

    alpha2 = 125; % Change for Desired Total Angle
    alphat = alpha2 - 35;
    
    phi_eq2 = @(phi2) ((sqrt(2) * (1 - phi2)) / sqrt(2 - phi2)) - sin(alpha2 * pi / 360);
    phi2 = fzero(phi_eq2, 0.8);
    A2 = (sqrt(2) * (1 - phi2)) / (phi2 * sqrt(phi2));
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3))); 
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

elseif op_cl2 == 0 && (mix == 0 || mix == 2)

    A2 = eta2 / (n2 * (1 - eta2) ^ 2); 
    phi2 = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A2, 0.8); 
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3)));
    alpha2 = (360 / pi) * asin(((sqrt(2) * (1 - phi2)) / sqrt((2 - phi2))));
    alphat = alpha2 - 35;
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

elseif op_cl2 == 0 && (mix == 1)

    A2 = eta2 / (n2 * (1 - eta2) ^ 2); 
    phi2 = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A2, 0.8); 
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3)));
    alpha2 = (360 / pi) * asin(((sqrt(2) * (1 - phi2)) / sqrt((2 - phi2))));
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

elseif op_cl2 == 1 && mix == 1 && op_cl1 ~= 2

    alpha2 = alpha1 - d;
    phi_eq2 = @(phi2) ((sqrt(2) * (1 - phi2)) / sqrt(2 - phi2)) - sin(alpha2 * pi / 360);
    phi2 = fzero(phi_eq2, 0.8);
    A2 = (sqrt(2) * (1 - phi2)) / (phi2 * sqrt(phi2));
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3))); 
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

elseif (op_cl2 == 1 || op_cl2 == 0) && mix == 1 && op_cl1 == 2

    alpha2 = 90; % Change for Desired Total Angle

    phi_eq2 = @(phi2) ((sqrt(2) * (1 - phi2)) / sqrt(2 - phi2)) - sin(alpha2 * pi / 360);
    phi2 = fzero(phi_eq2, 0.8);
    A2 = (sqrt(2) * (1 - phi2)) / (phi2 * sqrt(phi2));
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3))); 
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

elseif op_cl2 == 2

    if op_cl1 == 3

        rn1o = 0;

    end

    rn2 = sqrt((mdot2 / (pi * cd2 * sqrt(2 * rho2 * deltap2))) + (rn1o ^ 2));
    alpha2 = 0;

end

% Check for Errors

valid = 1;

if mix == 1 && op_cl2 ~= 3 && op_cl1 ~= 3

    if alpha1 < alpha2

        fprintf('No External Mixing\n');
        fprintf('-----DESIGN INVALID-----\n');
        valid = 0;

    end

end

if op_cl1 == 0 && abs(rv1 - rn1) > 1e-8

    fprintf('Incorrect eta1 Configuration\n');
    fprintf('-----DESIGN INVALID-----\n');
    valid = 0;

end

if op_cl2 == 0 && abs(rv2 - rn2) > 1e-8

    fprintf('Incorrect eta2 Configuration\n');
    fprintf('-----DESIGN INVALID-----\n');
    valid = 0;

end

if op_cl2 ~= 2 && op_cl2 ~= 3 && op_cl1 ~= 3

    rn_l2 = rn2 * sqrt(1 - phi2);  
    phiv0 = ((3 * phi2) - (2 * phi2 ^ 2)) / (2 - phi2);
    rv0_l2 = rn2 * sqrt(1 - phiv0); 
    rv0_l2_norm = sqrt(1 - phiv0);
    rv2_norm = rv2 / rn2;
    rv_l2 = rn2 * fzero(@(rvl) ((rv0_l2_norm^2 * (rv2_norm^2 - rvl^2)) / (rv2_norm^2 - rvl^2 - mu2^2)) - (rvl^2), 0.5); % Vortex Chamber

    if rn_l2 - rn1o < 0

        gvn = rn_l2 - rn1o;
        fprintf('Inner Element Not Accommodated in Nozzle Gas vortex: %.4f mm\n', gvn * 1e3);
        fprintf('-----DESIGN INVALID-----\n');
        valid = 0;

    elseif rv0_l2 - rn1o < 0

        gvv0 = rv0_l2 - rn1o;
        fprintf('Inner Element Not Accommodated in Head End Gas vortex: %.4f mm\n', gvv0 * 1e3);
        fprintf('-----DESIGN INVALID-----\n');
        valid = 0;

    elseif rv_l2 - rn1o < 0

        gvv = rv_l2 - rn1o;
        fprintf('Inner Element Not Accommodated in Vortex Chamber Gas Vortex: %.4f mm\n', gvv * 1e3);
        fprintf('-----DESIGN INVALID-----\n');
        valid = 0;

    elseif valid == 1

        gvn = rn_l2 - rn1o;
        fprintf('Nozzle Gas Vortex: %.4f mm\n', gvn * 1e3);
        gvv0 = rv0_l2 - rn1o;
        fprintf('Head End Gas Vortex: %.4f mm\n', gvv0 * 1e3);
        gvv = rv_l2 - rn1o;
        fprintf('Vortex Chamber Gas Vortex: %.4f mm\n\n', gvv * 1e3);
        
        if gvn < 3e-4

            fprintf('Nozzle Gas Vortex May be Invalid - Less Than 0.3mm\n\n');
            
        end

    end

end

% Find Recess

if mdot1 + mdot2 < 0.2 || mdot1 + mdot2 > 1.0

    fprintf('Mixing Time May be Invalid For Total Mass Flow Rate\n\n');

end


if prop == 0 && op_cl1 ~= 3 && op_cl2 ~= 3

    if (op_cl1 == 1 || op_cl1 == 0) && op_cl2 ~= 2

        Kmo = mdot1 / (mdot2 + mdot1);
        Kmf = mdot2 / (mdot2 + mdot1);
        lmix = sqrt(2)*tao*((Kmf*(mu2/phi2)*sqrt(deltap2/rho2)) + (Kmo*(mu1/phi1)*sqrt(deltap1/rho1)));
        del = ((rn2 - rn1) / tan(alpha1 * pi /360));

    end

elseif prop == 1 && op_cl1 ~= 3 && op_cl2 ~= 3

    if (op_cl1 == 1 || op_cl1 == 0) && op_cl2 ~= 2

        Kmf = mdot1 / (mdot2 + mdot1);
        Kmo = mdot2 / (mdot2 + mdot1);
        lmix = sqrt(2)*tao*((Kmo*(mu2/phi2)*sqrt(deltap2/rho2)) + (Kmf*(mu1/phi1)*sqrt(deltap1/rho1)));
        del = ((rn2 - rn1) / tan(alpha1 * pi /360));

    end

end

if mix == 0 && op_cl1 ~= 3 && op_cl2 ~= 3

    if (op_cl1 == 1 || op_cl1 == 0) && (op_cl2 == 1 || op_cl2 == 0)

        recess = lmix + del;

    end

elseif (mix == 1 || mix == 2) && op_cl1 ~= 3 && op_cl2 ~= 3

    if (op_cl1 == 1 || op_cl1 == 0) && (op_cl2 == 1 || op_cl2 == 0) 

        recess = del;

    end

end

if valid == 1

    fprintf('-----DESIGN VALID-----\n\n');

    if op_cl1 ~= 3

    fprintf('-----INNER ELEMENT-----\n\n')
    fprintf('Nozzle Radius: %.4f mm\n', rn1 * 1e3);
    fprintf('Outer Nozzle Radius: %.4f mm\n', rn1o * 1e3);

        if op_cl1 ~= 2

            fprintf('Vortex Chamber Radius: %.4f mm\n', rv1 * 1e3);
            fprintf('Tangential Inlet Radius: %.4f mm\n', rt1 * 1e3);
            fprintf('Number of Tangential Inlet: %.4f \n', n1);
            fprintf('Spray Angle: %.4f degrees\n', alpha1);
            fprintf('Mass Flow Coefficient: %.4f\n', mu1);

            fprintf('Nozzle Length > %.4f mm\n', 2e3 * rn1);
            fprintf('Vortex Chamber Length > %.4f mm\n', 2e3 * (rv1 - rt1));
            fprintf('Tangential Inlet Length > %.4f mm\n', 3e3 * rt1);

            if op_cl2 ~= 3 && op_cl2 ~= 2 && mix ~= 1

                fprintf('Recess: %.4f mm\n', recess * 1e3);
                fprintf('Recess Ratio: %.4f\n', recess / del);

            elseif op_cl2 == 2 && mix ~= 1

                del = ((rn2 - rn1) / tan(alpha1 * pi /360));
                fprintf('Recess > %.4f mm\n', del * 1e3);
                fprintf('Recess Ratio > %.4f\n', del / del);                

            end

        end

    end

    if op_cl2 ~= 3

    fprintf('\n-----OUTER ELEMENT-----\n\n')
    fprintf('Nozzle Radius: %.4f mm\n', rn2 * 1e3);

        if op_cl2 ~= 2

            fprintf('Vortex Chamber Radius: %.4f mm\n', rv2 * 1e3);
            fprintf('Tangential Inlet Radius: %.4f mm\n', rt2 * 1e3);
            fprintf('Number of Tangential Inlet: %.4f \n', n2);

            if (mix == 0 || mix == 2) && (op_cl1 ~= 3 && op_cl2 ~= 3) && (op_cl1 ~= 2)

                fprintf('Estimated Total Spray Angle: %.4f degrees\n', alphat);

            elseif mix == 1 || op_cl1 == 3

                fprintf('Spray Angle: %.4f degrees\n', alpha2);

            end

            fprintf('Mass Flow Coefficient: %.4f\n', mu2);

        end

    end

end

if op_cl1 == 2 && op_cl2 == 2 && valid == 1

    MR = (rho1 * (cd1 * sqrt(2 * deltap1 / rho1)) ^ 2) / (rho2 * (cd2 * sqrt(2 * deltap2 / rho2)) ^ 2);
    fprintf('\n1-2 Momentum Flux Ratio: %.4f\n', MR);

end
