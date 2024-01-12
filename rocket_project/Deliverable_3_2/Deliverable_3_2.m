addpath(fullfile('..', 'src'));

%close all
%clear all
%clc

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/20;
Tf = 10;
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


x_ref_open_roll = deg2rad(35);
x_ref_open = -4;

%% X open loop
[u, T_opt, X_opt, U_opt] = mpc_x.get_u([0,0,0,0]',x_ref_open);
U_opt(:,end+1) = NaN;

xs_x = repmat(xs([1:4],1), 1, length(X_opt));

X_opt = X_opt + xs_x;
U_opt = U_opt + us(1);
ph = rocket.plotvis_sub(T_opt,X_opt,U_opt,sys_x,xs,us,x_ref_open);

%% Y open loop
[u, T_opt, X_opt, U_opt] = mpc_y.get_u([0,0,0,0]',x_ref_open);
U_opt(:,end+1) = NaN;

xs_x = repmat(xs([5:8],1), 1, length(X_opt));

X_opt = X_opt + xs_x;
U_opt = U_opt + us(2);
ph = rocket.plotvis_sub(T_opt,X_opt,U_opt,sys_y,xs,us,x_ref_open);

%% Z open loop
[u, T_opt, X_opt, U_opt] = mpc_z.get_u([0,0]',x_ref_open);
U_opt(:,end+1) = NaN;

xs_x = repmat(xs([9:10],1), 1, length(X_opt));

X_opt = X_opt + xs_x;
U_opt = U_opt + us(3);
ph = rocket.plotvis_sub(T_opt,X_opt,U_opt,sys_z,xs,us,x_ref_open);

%% Roll open loop
[u, T_opt, X_opt, U_opt] = mpc_roll.get_u([0,0]',x_ref_open_roll);
U_opt(:,end+1) = NaN;

xs_x = repmat(xs([11:12],1), 1, length(X_opt));

X_opt = X_opt + xs_x;
U_opt = U_opt + us(4);
ph = rocket.plotvis_sub(T_opt,X_opt,U_opt,sys_roll,xs,us,x_ref_open_roll);

%% X closed loop
x0 = [0,0,0,0]'; 
x_ref = -4;
[T, X_sub, U_sub] = rocket.simulate_f(sys_x, x0, Tf, @mpc_x.get_u, x_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us, x_ref);


%% Y closed loop
x0 = [0,0,0,0]'; 
y_ref = -4; 
[T, X_sub, U_sub] = rocket.simulate_f(sys_y, x0, Tf, @mpc_y.get_u, y_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_y, xs, us, y_ref);

%% Z closed loop
x0 = [0,0]'; 
z_ref = -4;
[T, X_sub, U_sub] = rocket.simulate_f(sys_z, x0, Tf, @mpc_z.get_u, z_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_z, xs, us, z_ref);

%% Roll closed loop
x0 = [0,0]'; 
roll_ref = deg2rad(35); 
[T, X_sub, U_sub] = rocket.simulate_f(sys_roll, x0, Tf, @mpc_roll.get_u, roll_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sys_roll, xs, us, roll_ref);

