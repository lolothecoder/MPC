addpath ('C:\Users\lolon\OneDrive\EPFL\RO-1ST YEAR 2023-2024\MPC\Casadi')
addpath(fullfile('..', 'src'));

close all
clear all
clc

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/20;
rocket = Rocket(Ts);
H = 5; % Horizon length in seconds
nmpc = NmpcControl(rocket, H);

% MPC reference with default maximum roll = 15 deg
%ref = @(t_, x_) ref_TVC(t_);

% MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(50);
ref = @(t_, x_) ref_TVC(t_, roll_max);

% Evaluate once and plot optimal open−loop trajectory,
% pad last input to get consistent size with time and state
x0 = zeros(12,1); %Initial state
ref0 = [2 2 2 deg2rad(40)]'; % Reference Tracking for open loop simulaton

%x0 = [deg2rad([2 -2 0, -2 2 0]), 0 0 0, 0 0 0]';

[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref0);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref0); 

Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
rocket.anim_rate = 10;
ph = rocket.plotvis(T, X, U, Ref);
