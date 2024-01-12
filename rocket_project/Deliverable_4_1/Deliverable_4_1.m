addpath(fullfile('..', 'src'));


%close all
%clear all
%clc

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/20;
rocket = Rocket(Ts);
[xs,us] = rocket.trim();
sys = rocket.linearize(xs,us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys,xs,us);

%Design MPC controller
H = 8;
mpc_x = MpcControl_x(sys_x,Ts,H);
mpc_y = MpcControl_y(sys_y,Ts,H);
mpc_z = MpcControl_z(sys_z,Ts,H);
mpc_roll = MpcControl_roll(sys_roll,Ts,H);

% Merge four sub−system controllers into one full−system controller
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);

%Evaluate once and plot optimal open-loop trajectory,
%pad last input to get consistent size with time and state
x0 = zeros(12,1);
ref4 = [2, 2, 2, deg2rad(40)]';
[u, T_opt, X_opt, U_opt] = mpc.get_u(x0, ref4);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt,X_opt,U_opt,ref4);

%Setup reference function
ref = @(t_, x_) ref_TVC(t_);

%Simulate
Tf = 30;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);

% Visualize
rocket.anim_rate = 10; % Increase this to make the animation faster
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title

%% For Part 6.1, uncomment if you would like to try the code
% Ts = 1/20; % Sample time
% rocket = Rocket(Ts);
% H = 4; % Horizon length in seconds
% 
% [xs, us] = rocket.trim();
% sys = rocket.linearize(xs, us);
% [sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);
% 
% mpc_x = MpcControl_x(sys_x, Ts, H);
% mpc_y = MpcControl_y(sys_y, Ts, H);
% mpc_z = MpcControl_z(sys_z, Ts, H);
% mpc_roll = MpcControl_roll(sys_roll, Ts, H);
% 
% % Merge four sub−system controllers into one full−system controller
% mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
% 
% 
% % Setup reference function
% Tf = 30;
% roll_max = deg2rad(50);
% ref = @(t_, x_) ref_TVC(t_, roll_max);
% x0 = zeros(12,1);
% 
% % Evaluate once and plot optimal open−loop trajectory,
% % pad last input to get consistent size with time and state
% [u, T_opt, X_opt, U_opt] = mpc.get_u(x0, ref);
% U_opt(:,end+1) = nan;
% ph = rocket.plotvis(T_opt, X_opt, U_opt, ref);
% 
% [T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
% 
% 
% % Plot pose
% rocket.anim_rate = 10; % Increase this to make the animation faster
% ph = rocket.plotvis(T, X, U, Ref);
% ph.fig.Name = 'Merged lin. MPC in nonlinear simulation with 50 deg roll_max'; % Set a figure title


