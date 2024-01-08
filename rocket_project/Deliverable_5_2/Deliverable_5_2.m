addpath(fullfile('..', 'src'));
%% TODO: This file should produce all the plots for the deliverable

%close all
%clear all
%clc

Ts = 1/20; % Sample time 
Tf = 8; % Simulation time
rocket = Rocket(Ts);
H = 4; % Horizon
[xs, us] = rocket.trim(); % Finding trim point
sys = rocket.linearize(xs, us); % Linearization
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us); 

%% Design MPC controller
mpc_x = MpcControl_x(sys_x, Ts, H);
mpc_y = MpcControl_y(sys_y, Ts, H);
mpc_z = MpcControl_z(sys_z, Ts, H);
mpc_roll = MpcControl_roll(sys_roll, Ts, H);

% Merge four sub−system controllers into one full−system controller
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

%% Simulation

ref = [1.2, 0, 3, 0]';
x0 = [zeros(1, 9), 1 0 3]';
%x0 = zeros(12,1);
%ref = @(t_, x_)ref_TVC(t_);
rocket.mass=2.13;
rocket.mass_rate=-0.27;
[T, X, U, Ref, Zhat] = rocket.simulate_est_z(x0, Tf, @mpc.get_u, ref, mpc_z, sys_z);

% Plot pose
rocket.anim_rate = 15; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'With mass rate'; % Set a figure title
