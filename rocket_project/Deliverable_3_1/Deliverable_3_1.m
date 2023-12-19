addpath(fullfile('..', 'src'));

%close all
%clear all
%clc

%% TODO: This file should produce all the plots for the deliverable

Ts = 1/20;
Tf = 7;
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

%u_x = mpc_x.get_u(x_x);
% 
% [u, T_opt, X_opt, U_opt] = mpc_x.get_u(x);
% U_opt(:,end+1) = NaN;
% 
% % X_opt = ...
% % U_opt = ...
% ph = rocket.plotvis_sub(T_opt,X_opt,U_opt,sys_x,xs,us);



%% X 
% x0 = [0,0,0,3]'; 
% [T, X_sub, U_sub] = rocket.simulate_f(sys_x, x0, Tf, @mpc_x.get_u, 0);
% ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us);


%% Y
% x0 = [0,0,0,3]'; 
% [T, X_sub, U_sub] = rocket.simulate_f(sys_y, x0, Tf, @mpc_y.get_u, 0);
% ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_y, xs, us);

%% Z 
x0 = [0,3]';
[T, X_sub, U_sub] = rocket.simulate_f(sys_z, x0, Tf, @mpc_z.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_z, xs, us);

%% Roll
% x0 = [0, deg2rad(40)]'; 
% [T, X_sub, U_sub] = rocket.simulate_f(sys_roll, x0, Tf, @mpc_roll.get_u, 0);
% ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_roll, xs, us);

