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

% Define Constant for Engine

a = 1087.88; % Local Speed of Sound (m/s)
R = 91e-3 / 2; % Radius of Chamber (m)
L = 90.406e-3; % Length (m)
threshold = 15000; % Frequency Threshold (Hz)

% Bessel Function Roots

alpha = [
    0.000, 1.220, 2.333, 3.238, 4.241;  % n=0
    0.586, 1.697, 2.717, 3.726, 4.731;  % n=1
    0.972, 2.135, 3.173, 4.192, 5.204;  % n=2
    1.337, 2.551, 3.612, 4.643, 5.662;  % n=3
    1.693, 2.995, 4.037, 5.082, 6.110   % n=4
];

q_values = 0:3;

% Initialize Array
f = zeros(5, 5, length(q_values)); % 5x5 Matrix

% Calculations

for m = 0:3

    for n = 0:3

        for i = 1:length(q_values)

            q = q_values(i);  

            alpha_mn = alpha(n+1, m+1);  
            
            f(m+1, n+1, i) = (a/2) * sqrt((alpha_mn / R)^2 + (q / L)^2);

        end

    end

end

fprintf('-----RESONANT FREQUENCIES-----\n\n')
fprintf('Range = %.4f - %.4f Hz\n\n', 0, threshold)

for m = 0:3

    for n = 0:3

        for i = 1:length(q_values)

            q = q_values(i);

            if f(m+1, n+1, i) < threshold

                fprintf('R = %d, T = %d, L = %d, f = %.4f Hz\n', m, n, q, f(m+1, n+1, i));

            end

        end

    end

end
