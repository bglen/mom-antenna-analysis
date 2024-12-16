%{

antenna_analysis.m
Brian Glen
12/4/24

Description:
This program analyzes a dipole antenna using both Pocklington's and
Hallen's equations, solved with method of moments. It sweeps both antenna
radius and dipole gap, and calculates the radiated power and radiations
patterns for each combination.

It implements the Electromagnetic Waves and Antennas MATLAB toolbox by Dr.
Sophocles J Orfanidis, Rutgers University

%}

clc;
clear all;

% add the EWA toolbox to the path
directory = fileparts(which(mfilename)); 
% Add  directory and all sub-directories to the path
addpath(genpath(directory));

%% Simulation Parameters
% runtime parameters
use_pocklington = 0;
use_hallens = 1;

check_convergence = 0;
calc_current_distribution = 0;
calc_input_impedance = 1;
sweep_antenna_radius = 0;
sweep_antenna_gap = 0;

% constant antenna parameters
lambda = 1.625; % wavelength [m]
antenna_length = 1.3; % length of antenna [m]
gap_voltage = 120; % time harmonic gap voltage [V]

% other parameters
eta_0 = 120 * pi; % intrinsic impedance of free space
k = 2 * pi / lambda; % wave number

%% Pocklington's Equation
if use_pocklington

    %% 1.a Check number of elements required for convergence of solution
    if check_convergence
        a = 0.1 / 1000; % radius [m]
        max_error = 0.001;
        N = check_convergence_pocklington(a, antenna_length, lambda, gap_voltage, eta_0, max_error);
    end
    
    %% 1.b Plot current distribution of antenna
    if calc_current_distribution

        a = 0.1 / 1000; % radius [m]
        N = 55; % number of elements
        
        [~, ~, z] = mom_parameters(N, antenna_length, lambda);

        % Magnetic Frill Source
        % Calculate input impedence of an ideal 0.8*lambda antenna
        r_rad = ideal_impedance(eta_0, (0.8*lambda), lambda); % input impedance of an ideal antenna is just radiation resistance
        
        % calculate required outer radius b for annular aperture from input impedence
        b = a * exp((2 * pi * r_rad) / eta_0); % [m]
        
        % Z-component of incident electric field
        E_z = zeros(N, 0);
        for n = 1:N
            E_z(n) = frill_e_field_simple(z(n), a, b, gap_voltage, k);
        end
        
        % calculate current distribution, plot result
        [~, ~, ~, ~] = current_distribution(N, "pocklington", a, antenna_length, lambda, E_z, gap_voltage, eta_0, 1);
        
    end
    
    %% 1.c Calculate input impedance
    if calc_input_impedance

        a = 0.1 / 1000; % radius [m]
        delta = 0.5 / 1000; % antenna gap [m]
        N = 55; % number of elements
        
        [~, ~, z] = mom_parameters(N, antenna_length, lambda);

        %% Magnetic Frill Source
        % Calculate input impedence of an ideal 0.8*lambda antenna
        r_rad = ideal_impedance(eta_0, (0.8*lambda), lambda); % input impedance of an ideal antenna is just radiation resistance
        
        % calculate required outer radius b for annular aperture from input impedence
        b = a * exp((2 * pi * r_rad) / eta_0); % [m]
        
        % Z-component of incident electric field
        E_z = zeros(N, 0);
        for n = 1:N
            E_z(n) = frill_e_field_simple(z(n), a, b, gap_voltage, k);
        end
        
        % calculate current distribution, no plotting
        [I, ~, ~, ~] = current_distribution(N, "pocklington", a, antenna_length, lambda, E_z, gap_voltage, eta_0, 0);

        % calculate input impedance
        Z_in = input_impedance(delta, z, I, gap_voltage);
    
        % ideal input impedance
        Z_in_ideal = ideal_impedance(eta_0, antenna_length, lambda); % input impedance of an ideal antenna is just radiation resistance
        fprintf("\nCalculated input impedance is %.2f [ohms]\n", Z_in);
        fprintf("Ideal input impedance is %.2f [ohms]\n", Z_in_ideal);
    end
    
    %% 1.d Sweep radius of antenna
    if sweep_antenna_radius
        a = [0.01, 0.05, 0.1, 0.5, 1, 4] / 1000; % [m]
        delta = 0.5 / 1000; % antenna gap [m]
        N = 55; % number of elements

        % Magnetic Frill Source
        % Calculate input impedence of an ideal 0.8*lambda antenna
        Z_source = ideal_impedance(eta_0, (0.8*lambda), lambda);

        sweep_radius(N, "pocklington", a, delta, antenna_length, lambda, gap_voltage, Z_source, eta_0);
    end
    
    %% 1.f Sweep antenna gap
    if sweep_antenna_gap
        a = 0.01 / 1000; % [m]
        delta = [0.01, 0.05, 0.1, 0.5, 1, 4] / 1000; % [m]
        N = 3255;
        
        % Magnetic Frill Source
        % Calculate input impedence of an ideal 0.8*lambda antenna
        Z_source = ideal_impedance(eta_0, (0.8*lambda), lambda);

        % sweep antenna
        sweep_gap(N, "pocklington", a, delta, antenna_length, lambda, gap_voltage, Z_source, eta_0);

    end
    
end

%% 2. Hallen's Equations
if use_hallens

    if check_convergence
        a = 0.1 / 1000; % radius [m]
        max_error = 0.001;
        N = check_convergence_hallen(a, antenna_length, lambda, gap_voltage, max_error);
    end
    
    if calc_current_distribution   
        a = 0.1 / 1000; % [m]
        delta = 0.5 / 1000; % [m]
        N = 213;

        [I, ~, ~, ~] = current_distribution(N, "hallen", a, antenna_length, lambda, 0, gap_voltage, eta_0, 1);
    end

    %% Calculate input impedence
    if calc_input_impedance
        a = 0.1 / 1000; % [m]
        delta = 0.5 / 1000; % [m]
        N = 213;

        % calculate current distribution
        [I, ~, ~, z] = current_distribution(N, "hallen", a, antenna_length, lambda, 0, gap_voltage, eta_0, 0);

        % calculated input impedance from current distribution
        Z_in = input_impedance(delta, z, I, gap_voltage);
    
        % ideal input impedance
        Z_in_ideal = ideal_impedance(eta_0, antenna_length, lambda); % input impedance of an ideal antenna is ~= radiation resistance

        fprintf("\nCalculated input impedance is %.2f [ohms]\n", Z_in);
        fprintf("Ideal input impedance is %.2f [ohms]\n", Z_in_ideal);
    end
    
    %% Sweep antenna radius
    if sweep_antenna_radius
        a = [0.01, 0.05, 0.1, 0.5, 1, 4] / 1000; % [m]
        delta = 1 / 1000; % [m]
        N = 213;

        sweep_radius(N, "hallen", a, delta, antenna_length, lambda, gap_voltage, 0, eta_0)
    end
    
    %% Sweep antenna gap
    if sweep_antenna_gap
        a = 0.01 / 1000; % [m]
        delta = [0.01, 0.05, 0.1, 0.5, 1, 4] / 1000; % [m]
        N = 213;

        % sweep antenna
        sweep_gap(N, "hallen", a, delta, antenna_length, lambda, gap_voltage, 0, eta_0);
    end


end

%% 3. Calculate Radiated Power

%% Wrapper Functions

% Increases the element size N until the required convergence error is
% acheived. Plots convergemce and returns N
function N = check_convergence_pocklington(a, antenna_length, lambda, gap_voltage, eta , max_error)
    N = 5; % initial number of elements, will always be odd
    beta = (2*pi) / lambda;
    
    converged = false;
    max_iterations = 200;
    i = 1;

    I_sum = zeros(1, max_iterations); % holds results
    error = zeros(1, max_iterations); % holds relative error between iterations
    elements = zeros(1, max_iterations); % holds the N values tested

    while converged == false && i <= max_iterations
        
        elements(i) = N;

        % for each iteration, calculate the sum of current on the antenna
        [~, ~, z] = mom_parameters(N, antenna_length, lambda);

        % Magnetic Frill Source
        r_rad = ideal_impedance(eta, (0.8*lambda), lambda); % input impedance of an ideal 0.8*lambda antenna
        
        % calculate required outer radius b for annular aperture from input impedence
        b = a * exp((2 * pi * r_rad) / eta); % [m]
        
        % Z-component of source electric field
        E_z = zeros(N, 0);
        for n = 1:N
            E_z(n) = frill_e_field_simple(z(n), a, b, gap_voltage, beta);
        end
        
        % calculate current distribution
        [I, ~, ~, ~] = current_distribution(N, "pocklington", a, antenna_length, lambda, E_z, gap_voltage, eta, 0);

        I_sum(i) = sum(abs(I));

        % Calculate a relative error between iterations
        if i >= 2 % only calculate error on 2nd iteration and beyond
            error(i) = abs((I_sum(i-1) - I_sum(i)) / N);

            fprintf("\nN = %d, Relative error = %f", N, error(i));

            if (error(i) <= max_error)
                converged = true; % end loop
                fprintf("\nSolution converged!");
            end
        end
        
        i = i + 1;

        if (i > max_iterations)
            fprintf("\nConvergence hit max iterations!");
        else
            N = N + 2; % go to next odd element size
        end
    end
    
    % trim arrays for plotting
    elements = elements(2:(i - 1));
    error = error(2:(i - 1));

    % Plot convergence of relative error
    figure;
    plot(elements, (error*100))
    xlabel("Number of elements")
    ylabel("Relative error")
    ytickformat('percentage')
    title(sprintf("Convergence of Solution - Pocklington Method, error <= %.3f", max_error))
end

% Increases the element size N until the required convergence error is
% acheived. Plots convergemce and returns N
function N = check_convergence_hallen(a, antenna_length, lambda, gap_voltage, eta, max_error)
    N = 5; % initial number of elements, will always be odd
    
    converged = false;
    max_iterations = 200;
    i = 1;

    I_sum = zeros(1, max_iterations); % holds results
    error = zeros(1, max_iterations); % holds relative error between iterations
    elements = zeros(1, max_iterations); % holds the N values tested

    while converged == false && i <= max_iterations
        
        elements(i) = N;
        
        % calculate current distribution
        [I, ~, ~, ~] = current_distribution(N, "hallen", a, antenna_length, lambda, 0, gap_voltage, eta, 0);

        I_sum(i) = sum(abs(I));

        % Calculate a relative error between iterations
        if i >= 2 % only calculate error on 2nd iteration and beyond
            error(i) = abs((I_sum(i-1) - I_sum(i)) / N);

            fprintf("\nN = %d, Relative error = %f", N, error(i));

            if (error(i) <= max_error)
                converged = true; % end loop
                fprintf("\nSolution converged!");
            end
        end
        
        i = i + 1;

        if (i > max_iterations)
            fprintf("\nConvergence hit max iterations!");
        else
            N = N + 2; % go to next odd element size
        end
    end
    
    % trim arrays for plotting
    elements = elements(2:(i - 1));
    error = error(2:(i - 1));

    % Plot convergence of relative error
    figure;
    plot(elements, (error*100))
    xlabel("Number of elements")
    ylabel("Relative error")
    ytickformat('percentage')
    title(sprintf("Convergence of Solution - Hallen Method, error <= %.3f", max_error))
end

%{
Calculate & plot current distribution for pocklington and hallen equations
 - E_z required for pocklington
%}
function [I, D, M, z] = current_distribution(N, method, a, antenna_length, lambda, E_z, gap_voltage, eta, plot_results)
    % Simulation setup
    [D, M, z] = mom_parameters(N, antenna_length, lambda);
    
    %% Solve for the unknown current distribution I(z)
    if method == "pocklington"
        [I, ~, ~] = pfield((antenna_length / lambda), (a / lambda), E_z, 'e', 't'); % using triangular basis function, might be better for impedance estimation
    elseif method == "hallen"
        [I, ~, ~] = hdelta((antenna_length / lambda), (a / lambda), M); 

        % Scale the normalized current distribution
        I = gap_voltage * I;
    else
        I = 0; % problem
        fprintf("Problem using current_distribution function");
    end

    % transpose
    I = I';
    
    % Plot Results
    if plot_results
        % plot current distribution
        figure;
        subplot(2, 1, 1);
        plot(z, abs(I))
        xlim([-antenna_length/2 antenna_length/2]);
        xlabel("z (m)")
        ylabel("Magnituide of Current (A)")
        title("Current Distribution Magnituide")
    
        subplot(2, 1, 2);
        plot(z, angle(I))
        xlim([-antenna_length/2 antenna_length/2]);
        xlabel("z (m)")
        ylabel("Phase Angle of Current (rad)")
        title("Current Distribution Phase Angle")
        
        % plot normalized radiation pattern
        figure;
        plot_radiation_pattern(100, z, I, lambda, eta, "Normalized Radiation Pattern");
    end
end

%{
Sweeps antenna calculations given antenna radius a array
 - E_z, r_rad, eta required for pocklington
%}
function sweep_radius(N, method, a, delta, antenna_length, lambda, gap_voltage, Z_source, eta)

    % Simulation setup
    [~, ~, z] = mom_parameters(N, antenna_length, lambda);
    I = zeros(1, length(N));
    Z_in = zeros(1, length(a)); % setup array for impedance results
    
    % setup figure for current distribution & rad pattern
    figure; 
    subplot(2, 1, 1);
    hold on;
    xlabel("z (m)")
    ylabel("Magnituide of current (A)")
    title("Current distribution throughout antenna")
    xlim([-antenna_length/2 antenna_length/2]);

    subplot(2, 1, 2);
    hold on;
    xlabel("z (m)");
    ylabel("Phase angle of current (rad)");
    title("Phase angle of current throughout antenna");
    xlim([-antenna_length / 2 antenna_length / 2]);
    
    for r = 1:length(a)
        if method == "pocklington"

            % Setup source for each radius a
            b = a(r) * exp((2 * pi * Z_source) / eta); % [m], calculate required outer radius b for annular aperture from input impedence
    
            % Z-component of incident electric field
            E_z = zeros(N, 0);
            for n = 1:N
                E_z(n) = frill_e_field_simple(z(n), a(r), b, gap_voltage, ((2 * pi) / lambda));
            end
    
            % Calculate current distribution
            [I, ~, ~, ~] = current_distribution(N, "pocklington", a(r), antenna_length, lambda, E_z, gap_voltage, eta, 0);

            % Calculate radiated power
            P = radiated_power(z, I, lambda, delta, eta);
            fprintf("\nRadiated Power is %.3f [W] for a = %0.2f mm", P, a(r)*1000);

        elseif method == "hallen"
    
            % Calculate current distribution
            [I, ~, ~, ~] = current_distribution(N, "hallen", a(r), antenna_length, lambda, 0, gap_voltage, eta, 0); 

            % Calculate radiated power
            P = radiated_power(z, I, lambda, delta, eta);
            fprintf("\nRadiated Power is %.3f [W] for a = %0.2f mm", P, a(r)*1000);

        else

            fprintf("problem sweeping radius!");

        end

        % Calculate impedance
        Z_in(r) = input_impedance(delta, z, I, gap_voltage); % result is negative because it is recieving power

        % Plot current distribution
        
        subplot(2, 1, 1);
        plot(z, abs(I), 'DisplayName', sprintf('a = %.4f m', a(r)));
    
        subplot(2, 1, 2);
        plot(z, angle(I), 'DisplayName', sprintf('a = %.4f m', a(r)));
        

        %figure;
        %plot_radiation_pattern(100, z, I, lambda, eta, sprintf("Rad Pattern, a = %0.4f m", a(r)));
    end

    
    % Add legends to current distribution subplots
    subplot(2, 1, 1);
    legend show;

    subplot(2, 1, 2);
    legend show;

    subplot(2, 1, 1);
    hold off;

    subplot(2, 1, 2);
    hold off;
    

    %% Plot imput impedance as function of antenna thickness
    figure;
    plot((a * 1000), Z_in);
    xlabel("Antenna radius [mm]");
    ylabel("Input impedance [ohms]");
    title("Antenna impedance as a function of radius");
end

%{
Sweeps antenna calculations given antenna gap delta array
 - E_z, r_rad, eta required for pocklington
%}
function sweep_gap(N, method, a, delta, antenna_length, lambda, gap_voltage, Z_source, eta)

    % Simulation setup
    [~, ~, z] = mom_parameters(N, antenna_length, lambda);
    I = zeros(1, length(N));
    Z_in = zeros(1, length(delta)); % setup array for impedance results
    
    % setup figure for current distribution
    figure;
    subplot(2, 1, 1);
    hold on;
    xlabel("z (m)")
    ylabel("Magnituide of current (A)")
    title("Current distribution throughout antenna")
    xlim([-antenna_length/2 antenna_length/2]);

    subplot(2, 1, 2);
    hold on; % Hold for current phase plot
    xlabel("z (m)");
    ylabel("Phase angle of current (rad)");
    title("Phase angle of current throughout antenna");
    xlim([-antenna_length / 2 antenna_length / 2]);

    for r = 1:length(delta)
        if method == "pocklington"

            % Setup source for each radius a
            b = a * exp((2 * pi * Z_source) / eta); % [m], calculate required outer radius b for annular aperture from input impedence
    
            % Z-component of incident electric field
            E_z = zeros(N, 0);
            for n = 1:N
                E_z(n) = frill_e_field_simple(z(n), a, b, gap_voltage, ((2*pi) / lambda));
            end
    
            % Calculate current distribution
            [I, ~, ~, ~] = current_distribution(N, "pocklington", a, antenna_length, lambda, E_z, gap_voltage, eta, 0);

            % Calculate radiated power
            P = radiated_power(z, I, lambda, delta(r), eta);
            fprintf("\nRadiated Power is %.3f [W] for delta = %0.2f mm", P, delta(r)*1000);

        elseif method == "hallen"
    
            % Calculate current distribution
            [I, ~, ~, ~] = current_distribution(N, "hallen", a, antenna_length, lambda, 0, gap_voltage, eta, 0);

            % Calculate radiated power
            P = radiated_power(z, I, lambda, delta(r), eta);
            fprintf("\nRadiated Power is %.3f [W] for delta = %0.2f mm", P, delta(r)*1000);


        else

            fprintf("problem sweeping gap!");

        end

        % Calculate impedance
        Z_in(r) = input_impedance(delta(r), z, I, gap_voltage);

        % Plot current distribution
        subplot(2, 1, 1);
        plot(z, abs(I), 'DisplayName', sprintf('delta = %.4f m', delta(r)));
    
        subplot(2, 1, 2);
        plot(z, angle(I), 'DisplayName', sprintf('delta = %.4f m', delta(r)));
    end

    % Add legends to current distribution subplots
    subplot(2, 1, 1);
    legend show;

    subplot(2, 1, 2);
    legend show;

    % Release holds
    subplot(2, 1, 1);
    hold off;

    subplot(2, 1, 2);
    hold off;

    %% 1.e Plot imput impedance as function of antenna gap
    figure;
    plot((delta * 1000), Z_in);
    xlabel("Antenna gap [mm]");
    ylabel("Input impedance [ohms]");
    title("Antenna impedance as a function of antenna gap");
end

%% Helper Functions

%{
    Plots the radiation pattern as a superposition of Hertzian dipoles for
    a given current distribution
%}
function plot_radiation_pattern(angle_samples, z, I, lambda, eta, plot_title)
    
    beta = (2 * pi) / lambda;
    R = 100; % arbitrary distance
    
    theta = linspace(0, 2*pi, angle_samples); % sample from 0 to 2*pi radians
    E_field = zeros(1, angle_samples);
    
    % Compute the electric field at each angle
    for i = 1:length(theta)
        for n = 1:length(z)
            R_n = R - z(n) * cos(theta(i)); % length to point at z(n)
            E_point = 1j * (beta * eta * I(n) / (4 * pi * R)) * ...
                      exp(-1j * beta * R) * exp(1j * beta * z(n)) * sin(theta(i)); % individual contribution to field
    
            E_field(i) = E_field(i) + E_point; % sum the contribution of each element
        end
    end
    
    % Normalize the magnitude of the pattern
    E_field = normalize(abs(E_field), "range");
    
    % Plot the pattern
    polarplot(theta, E_field);
    title(plot_title);
end

%{
    Calculates radiated power for antenna. Power is calculated as the
    superposition of radiated hertzian dipoles for each current
    distribution element.
%}
function P = radiated_power(z, I, lambda, delta, eta)
    P = 0;
    
    for n = 1:length(z)
        if (z(n) > (delta / 2) || z(n) < -(delta / 2)) % ensure no current distribution within the gap is counted
            P_z = real(I(n))^2 * (eta * pi / 3) * (z(n) * 2 / lambda)^2; % [W]
            P = P + P_z;
        end
    end
end

% Calculate parameters for method of moments
function [D, M, z] = mom_parameters(N, antenna_length, lambda)
    D = antenna_length / N; % distance between sample points

    if (D >= lambda / 2)
        print("Not enough samples");
    end

    M = (N - 1) / 2; % extent of samples from the origin of the antenna
    z = linspace(-M*D, M*D, N); % this array holds the sample distances along the z axis
end

% Calculates the impedance of an ideal dipole antenna
function r = ideal_impedance(eta, antenna_length, wavelength)
    % Source antenna with input impedence of an ideal 0.68*lambda antenna
    fun = @(theta, L, wavelength) (cos(L .* pi ./ wavelength .* cos(theta)) - cos(L .* pi ./ wavelength)).^2 ./ sin(theta); % define math function inside of integral
    r = eta / (2*pi) * integral(@(theta) fun(theta, antenna_length, wavelength), 0, pi); % evaluate integral of rad resistance from 0 to pi
end

%{
    Calculates electric field on antenna surface from frill generator, Implementation of Balanis 8-33
    where -l/2 <= z <= l/2
%}
function E_z = frill_e_field_simple(z, a, b, source_voltage, k)
    r_1 = sqrt(z^2 + a^2);
    r_2 = sqrt(z^2 + b^2);

    % source_voltage is voltage present at the input connection to the
    % antenna
    E_z = (-source_voltage / (log(b / a))) * ((exp(-j*k*r_1) / r_1) - (exp(-j*k*r_2) / r_2));
end

% Calculates input impedance, by finding the feed current at the edge of
% the gap
function Z_in = input_impedance(antenna_gap, z, I, gap_voltage)

    % Find the index corresponding to the edge of the gap on one side
    i = find(z < -antenna_gap / 2, 1, 'last'); % Last element before the gap

    I_in = abs(I(i));

    Z_in = gap_voltage / I_in;
end

