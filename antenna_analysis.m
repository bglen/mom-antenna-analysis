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

%% add the EWA toolbox to the path
directory = fileparts(which(mfilename)); 
% Add  directory and all sub-directories to the path
addpath(genpath(directory));

%% Define Unknown Antenna Parameters
a = 0.1 / 1000; % antenna radius [m]
lambda = 0.6; % wavelength [m]
length = 1.8 * lambda; % length of antenna
delta = 1 / 1000; % dipole gap [m]
gap_voltage = 125; % time harmonic gap voltage [V]

%% Method of Moments Parameters
D = lambda / 2; % desired distance between samples
N = round((2 * length / D) + 1); % number of elements. Calculate minimum N to cover length of antenna
if mod(N, 2) == 0
    N = N+1; % ensure N is odd
end

M = (N - 1) / 2; % extent of samples from the origin of the antenna

%% Magnetic Frill Source
% Source antenna with input impedence of an ideal 0.68*lambda antenna

% calculate radius b from characteristic impedence

% magnetic current density

% Z-component of incidient electric field
E_i_z = zeros(2*M);

% add electric field at gap location


%{
%% Pocklington's Equation and Analysis

% Sweep element density to determine what is needed for solution to
% converge

% Solve for the unknown antenna current distribution I(z)
pfield();

% Calculate imput impedance
pmat();

% Sweep antenna radius & plot results

% Sweep antenna gap with constant radius of 0.01 mm

%% Delta Gap Source
% Electric field in the gap
E_delta_gap = gap_voltage / delta;

% magnetic current density

% Z-component of incidient electric field
E_i_z = zeros(10);

%% Hallen's Equation and Analysis

% Change source to a delta gap generator

% Sweep element density to determine what is needed for solution to
% converge

% Solve for the unknown antenna current distribution I(z)

% Calculate imput impedance 

% Sweep antenna radius & plot results

% Sweep antenna gap with constant radius of 0.01 mm

%% Calculate Radiated Power

%% Radiation Pattern
%}
