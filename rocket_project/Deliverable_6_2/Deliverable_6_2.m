addpath(fullfile('..', 'src'));
addpath ('C:\Users\lolon\OneDrive\EPFL\RO-1ST YEAR 2023-2024\MPC\Casadi')
close all
clear all
clc

%% TODO: This file should produce all the plots for the deliverable

Ts = 1/40; % Higher sampling rate for this part!

rocket = Rocket(Ts);
H = 5; % Horizon length in seconds
nmpc = NmpcControl(rocket, H);

x0 = zeros(12, 1);
ref = [0.5, 0, 1, deg2rad(65)]';
Tf = 2.5;
rocket.mass = 1.75;
rocket.delay = 1; % 0 if not specified
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
rocket.anim_rate = 10;
ph = rocket.plotvis(T, X, U, Ref);