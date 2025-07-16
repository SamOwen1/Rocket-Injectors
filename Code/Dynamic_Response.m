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

% Define Injector Parameters 

rn = 2.3311e-3; % Nozzle Radius (m)
rv = 8.4094e-3; % Vortex Radius (m)
rt = 1.4161e-3; % Tangential Inlet Radius (m)
lv = 70e-3; % Uniform Vortex Chamber Length (m)
ln = 20e-3; % Nozzle Length (m)
theta = 45 * (pi / 180); % Convergence Angle
lc = (rv - rn) / tan(theta); % Non-Uniform Vortex Chamber Length (m)
lt = 3 * rt; % Tangential Inlet Length (m)
rho = 998; % Fluid density (kg/m^3)
nu = 0.1; % Artificial Viscosity Coefficient
n = 6; % Number of tangential inlets

% Define Frequency Range (in Hz)

freq_range = linspace(0, 3000, 3000);
omega_range = freq_range .* (2 * pi);

% Define Throttle Values

throttle_values = [0.4, 0.6, 0.8, 1];
deltapnom = 8e5; % Nominal Pressure Drop (Pa)

% Colors for plotting

colors = lines(length(throttle_values));

figure;
hold on;

for p = 1:length(throttle_values)

    Psi = throttle_values(p);

    rin = rv - rt; 
    eta = rin / rn;
    A = (rn * rin) / (n * rt ^ 2);
    phi_guess = 0.8;
    phi = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A, phi_guess);

    mu = sqrt((phi ^ 3) / (2 - phi));
    mdot = mu * pi * (rn ^ 2) * sqrt(2 * rho * deltapnom);
    phiv0 = ((3 * phi) - (2 * phi ^ 2)) / (2 - phi);

    mdot = Psi * mdot;
    deltap = (1 / (2 * rho)) * (mdot / (mu * pi * rn ^ 2)) ^ 2;

    usum = sqrt((2 * deltap) / rho);
    rn_l_norm = sqrt(1 - phi);
    rv0_l_norm = sqrt(1 - phiv0);
    rn_norm = 1;
    rin_norm = rin / rn;
    rv_norm = rv / rn;

    rv_l_norm = fzero(@(r) ((rv0_l_norm^2 * (rv_norm^2 - r^2)) / (rv_norm^2 - r^2 - mu^2)) - (r^2), phi_guess);

    uvz = usum * sqrt(1 - (rv0_l_norm ^ 2) / (rv_l_norm ^ 2));
    unz = usum * sqrt(1 - (rv0_l_norm ^ 2) / (rn_l_norm ^ 2));
    uin = usum * rv0_l_norm / eta;
    C = usum * rv0_l_norm * rn;

    rt_norm = rt / rn;
    ln_norm = ln / rn;
    lv_norm = lv / rn;
    lt_norm = lt / rn;
    lc_norm = lc / rn;
    C_norm = C / (rn * uin);
    uvz_norm = uvz / uin;
    unz_norm = unz / uin;
    uin_norm = 1;

    cv_norm = uvz_norm + sqrt((C_norm ^ 2) * ((rv_norm ^ 2) - (rv_l_norm ^ 2)) / (2 * rv_l_norm ^ 4));
    cn_norm = unz_norm + sqrt((C_norm ^ 2) * ((rn_norm ^ 2) - (rn_l_norm ^ 2)) / (2 * rn_l_norm ^ 4));

    % Preallocate Response Magnitude

    response_magnitude = zeros(size(freq_range));

    % Loop Over Frequencies

    for idx = 1:length(omega_range)

        omega = omega_range(idx);
        omega_norm = (rn * omega) / uin;

        phiv = omega_norm * ((lv_norm + lc_norm) / cv_norm);
        phin = omega_norm * (ln_norm / cn_norm);

        Pi_refl = 1 - 2 * (sqrt(phi) / sqrt((rv_norm ^ 2) - (rv0_l_norm ^ 2)));

        V2_constant = 1 / (A * sqrt(2 * (rv_norm ^ 2 - rv0_l_norm ^ 2)));
        sum_terms_V2 = 1 / (1 - Pi_refl * exp(-2 * phiv * (1i + nu)));
        sum_terms_VN = exp(-phiv * (1i + nu)) / (1 - Pi_refl * exp(-2 * phiv * (1i + nu)));

        Pi_V2 = V2_constant * sum_terms_V2;
        Pi_VN = sum_terms_VN;

        Sh_V = omega_norm * rv_norm / uvz_norm;
        Sh_T = omega_norm * lt_norm / uin_norm;

        f_x = @(x) x .* Sh_V .* tan(pi * x / 2);
        Kdrop = 1 - (rv0_l_norm / rv_norm);
        Re_Pi_V3 = 2 * integral(@(x) (cos(f_x(x)) ./ (1 - Kdrop * x).^3) .* exp(-nu * f_x(x)), 0, 1, 'ArrayValued', true);
        Im_Pi_V3 = -2 * integral(@(x) (sin(f_x(x)) ./ (1 - Kdrop * x).^3) .* exp(-nu * f_x(x)), 0, 1, 'ArrayValued', true);
        Pi_V3 = Re_Pi_V3 + 1i * Im_Pi_V3;

        Pi_N = (1 - Pi_refl) * exp(-1i * phin);
        Pi_T = 0.5 * (1 - 1i * Sh_T) / (1 + Sh_T^2);

        Pi_Sigma = (rv_norm / rv0_l_norm)^2 * (Pi_T * Pi_VN * Pi_N) / (1 + 2 * Pi_T * (Pi_V2 + Pi_V3));

        response_magnitude(idx) = abs(Pi_Sigma);

    end

    % Detect Peaks

    [peakVals, peakLocs] = findpeaks(response_magnitude, freq_range, ...
        'MinPeakProminence', 0.01, 'MinPeakDistance', 100); % Adjust if Needed

    % Number of Resonant Modes

    n_modes = 5;

    % Bulk Mode 

    bulk_freq = peakLocs(1);
    bulk_mag = peakVals(1);

    % Print Bulk Mode Info

    fprintf('\nThrottle: %.1f\n', Psi);
    fprintf('Mass Flow Rate: %.4f kg/s\n', mdot);
    fprintf('Pressure Drop: %.4f bar\n', deltap * 1e-5);
    fprintf('  Bulk Mode: Frequency = %.2f Hz, Magnitude = %.3f\n', bulk_freq, bulk_mag);

    % Loop Over Next n_modes Peaks

    for k = 1:n_modes

        idx = k + 1; 

        if length(peakLocs) >= idx

            freq = peakLocs(idx);
            mag = peakVals(idx);
            fprintf('  Resonant Mode %d: Frequency = %.2f Hz, Magnitude = %.3f\n', k, freq, mag);

            if freq <= 400

                fprintf('    WARNING - Risk of Chugging - f < 400Hz\n');

            end

        else

            fprintf('  Resonant Mode %d: Not Detected\n', k);

        end

    end

    plot(freq_range, response_magnitude, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('$\\delta = %.1f$', Psi), ...
        'Color', colors(p, :));

end

xlabel('$f(Hz)$', 'Interpreter', 'latex');
ylabel('$\Pi_{\Sigma}$', 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on;
xlim([0 3000]);
