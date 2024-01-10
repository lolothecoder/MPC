classdef NmpcControl < handle
    
    properties
        solver
        nx, nu, N
        nlp_x0
        nlp_lbx, nlp_ubx
        nlp_lbg, nlp_ubg
        nlp_p
        
        T_opt
        sol
        idx
        
        % Delay compensation
        rocket
        expected_delay
        mem_u
        
        % Warmstart
        nlp_lam_x0
        nlp_lam_g0
    end
    
    methods
        function obj = NmpcControl(rocket, tf, expected_delay)
            
            if nargin < 3, expected_delay = 0; end
            
            import casadi.*
            
            N_segs = ceil(tf/rocket.Ts); % MPC horizon
            nx = 12; % Number of states
            nu = 4;  % Number of inputs
            
            % Decision variables (symbolic)
            N = N_segs + 1; % Index of last point
            X_sym = SX.sym('X_sym', nx, N); % state trajectory
            U_sym = SX.sym('U_sym', nu, N-1); % control trajectory)
            
            % Parameters (symbolic)
            x0_sym  = SX.sym('x0_sym', nx, 1);  % initial state
            ref_sym = SX.sym('ref_sym', 4, 1);  % target position
            
            % Default state and input constraints
            ubx = inf(nx, 1);
            lbx = -inf(nx, 1);
            ubu = inf(nu, 1);
            lbu = -inf(nu, 1);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            new_inf = 500;
            
            %Overwriting Bound Values
            ubx = zeros(nx, 1) + new_inf;
            lbx = zeros(nx, 1) - new_inf;
            ubu = zeros(nu, 1) + new_inf;
            lbu = zeros(nu, 1) - new_inf;
            
            
            f_discrete = @(x,u) RK4(x,u,rocket.Ts,@rocket.f);
            % Cost
            cost = 0;
            
            % Equality constraints (Casadi SX), each entry == 0
            eq_constr = [ ; ];
            
            % Inequality constraints (Casadi SX), each entry <= 0
            ineq_constr = [ ; ];

            % For box constraints on state and input, overwrite entries of
            % lbx, ubx, lbu, ubu defined above
            
            %Bounds on the input
            ubu(1,1) = 0.26;
            lbu(1,1) = -0.26;
            
            ubu(2,1) = 0.26;
            lbu(2,1) = -0.26;
            
            ubu(3,1) = 80;
            lbu(3,1) = 50;
            
            ubu(4,1) = 20;
            lbu(4,1) = -20;
            
            %Bounds on the states
            ubx(5,1) = deg2rad(75);
            lbx(5,1) = -deg2rad(75);
            
            wx = 5;
            wy = 5;
            wz = 1;
            alpha = 5;
            beta = 5;
            gamma = 200;
            vx = 1;
            vy = 1;
            vz = 1;
            x = 500;
            y = 500;
            z = 800;
            
            delta_1 = 10;
            delta_2 = 10;
            P_avg = 0.0001;
            P_diff = 0.001;
            
            Q = diag([wx wy wz alpha beta gamma vx vy vz x y z]);
            R = diag([delta_1 delta_2 P_avg P_diff]);
            
            gamma_ref = ref_sym(4,1);
            x_ref = ref_sym(1,1);
            y_ref = ref_sym(2,1);
            z_ref = ref_sym(3,1);
            reference = [0 0 0 0 0 gamma_ref 0 0 0 x_ref y_ref z_ref]';
            
            cost_init = (x0_sym-reference)'*Q*(x0_sym-reference);
            
            eq_constr = [X_sym(:,1)- f_discrete(x0_sym, zeros(nu,1))]; 
            cost_iterations = 0;
            for i = 1:N-1
                eq_constr = [eq_constr; X_sym(:,i+1) - f_discrete(X_sym(:,i), U_sym(:,i))];
                ineq_constr = [ineq_constr; lbu - U_sym(:,i); U_sym(:,i) - ubu];
                ineq_constr = [ineq_constr; lbx - X_sym(:,i); X_sym(:,i) - ubx];

                %Cost calculation
                cost_iterations= cost_iterations +(X_sym(:,i)-reference)'*Q*(X_sym(:,i)-reference) + (U_sym(:,i))'*R*(U_sym(:,i));
            end

            ineq_constr = [ineq_constr; lbx - (X_sym(:,N)); X_sym(:,N) - ubx];

            [xs, us] = rocket.trim();
            sys = rocket.linearize(xs,us);
            sys = c2d(sys, rocket.Ts);
            
            sys_linear = LTISystem('A', sys.A, 'B', sys.B);
            
            for i = 1:nu
                sys_linear.u.max(i) = ubu(i,1);
                sys_linear.u.min(i) = lbu(i,1);
            end
            sys_linear.x.max(5) = ubx(5,1);
            sys_linear.x.min(5) = lbx(5,1);
            
            sys_linear.x.penalty = QuadFunction(Q); 
            sys_linear.u.penalty = QuadFunction(R);
            
            Qf_mpt = sys_linear.LQRPenalty.weight;
            cost_terminal = (X_sym(:,N)-reference)'*Qf_mpt*(X_sym(:,N)-reference);

            cost = cost_init + cost_iterations + cost_terminal;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ---- Assemble NLP ------
            nlp_x = [X_sym(:); U_sym(:)];
            nlp_p = [x0_sym; ref_sym];
            nlp_f = cost;
            nlp_g = [eq_constr; ineq_constr];
            
            nlp = struct('x', nlp_x, 'p', nlp_p, 'f', nlp_f, 'g', nlp_g);
            
            % ---- Setup solver ------
            opts = struct('ipopt', struct('print_level', 0), 'print_time', false);
            obj.solver = nlpsol('solver', 'ipopt', nlp, opts);
            
            % ---- Assemble NLP bounds ----
            obj.nlp_x0  = zeros(size(nlp_x));
            
            obj.nlp_ubx = [repmat(ubx, N, 1); repmat(ubu, (N-1), 1)];
            obj.nlp_lbx = [repmat(lbx, N, 1); repmat(lbu, (N-1), 1)];
            
            obj.nlp_ubg = [zeros(size(eq_constr)); zeros(size(ineq_constr))];
            obj.nlp_lbg = [zeros(size(eq_constr)); -inf(size(ineq_constr))];
            
            obj.nlp_p = [zeros(size(x0_sym)); zeros(size(ref_sym))];
            
            obj.nlp_lam_x0 = [];
            obj.nlp_lam_g0 = [];
            
            obj.nx = nx;
            obj.nu = nu;
            obj.N = N;
            obj.T_opt = linspace(0, N * rocket.Ts, N);
            
            obj.idx.X = [1, obj.N * obj.nx];
            obj.idx.U = obj.idx.X(2) + [1, (obj.N-1) * obj.nu];
            obj.idx.u0 = obj.idx.U(1) + [0, obj.nu-1];
            
            % Members for delay compensation
            obj.rocket = rocket;
            obj.expected_delay = expected_delay;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % Could optimize using optimal u in last step
            % solve linearize problem for initial guess
            u_init = zeros(4, 1); % Replace this by a better initialization
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.mem_u = repmat(u_init, 1, expected_delay);
        end
        
        function [u, T_opt, X_opt, U_opt] = get_u(obj, x0, ref)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            delay = obj.expected_delay;
            mem_u = obj.mem_u;
            
            % Delay compensation: Predict x0 delay timesteps later.
            % Simulate x_ for 'delay' timesteps
            % Euler integration to go in the future
            x_ = x0;
            % ...
       
            x0 = x_;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute solution from x0
            obj.solve(x0, ref);
            
            % Evaluate u0
            nlp_x = obj.sol.x;
            id = obj.idx.u0;
            u = full( nlp_x(id(1):id(2)) );      
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % Delay compensation: Save current u
            if obj.expected_delay > 0
               % obj.mem_u = ...
            end
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargout > 1, T_opt = obj.get_T_opt(); end
            if nargout > 2, X_opt = obj.get_X_opt(); end
            if nargout > 3, U_opt = obj.get_U_opt(); end
            return
            
            % Additional evaluation
            % Complete trajectory
            % % X_opt = full(reshape(nlp_x(idx_X(1):idx_X(2)), obj.nx, obj.N));
            % % U_opt = full(reshape(nlp_x(idx_U(1):idx_U(2)), obj.nu, obj.N - 1));
            % %
            % % cost_opt = full(sol.f);
            % % constr_opt = full(sol.g);
            % %
            % % stats = obj.solver.stats;
        end
        
        function solve(obj, x0, ref)
            
            % ---- Set the initial state and reference ----
            obj.nlp_p = [x0; ref];     % Initial condition
            obj.nlp_x0(1:obj.nx) = x0; % Initial guess consistent
            
            % ---- Solve the optimization problem ----
            args = {'x0', obj.nlp_x0, ...
                'lbg', obj.nlp_lbg, ...
                'ubg', obj.nlp_ubg, ...
                'lbx', obj.nlp_lbx, ...
                'ubx', obj.nlp_ubx, ...
                'p', obj.nlp_p, ...
                %                 'lam_x0', obj.nlp_lam_x0, ...
                %                 'lam_g0', obj.nlp_lam_g0
                };
            
            obj.sol = obj.solver(args{:});
            if obj.solver.stats.success ~= true
                solve_status_str = obj.solver.stats.return_status;
                fprintf([' [' class(obj) ': ' solve_status_str '] ']);
                obj.sol.x(obj.idx.u0) = nan;
            end
            
            % Use the current solution to speed up the next optimization
            obj.nlp_x0 = obj.sol.x;
            obj.nlp_lam_x0 = obj.sol.lam_x;
            obj.nlp_lam_g0 = obj.sol.lam_g;
        end
        function T_opt = get_T_opt(obj)
            T_opt = obj.T_opt;
        end
        function X_opt = get_X_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.X;
            X_opt = full(reshape(nlp_x(id(1):id(2)), obj.nx, obj.N));
        end
        function U_opt = get_U_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.U;
            U_opt = full(reshape(nlp_x(id(1):id(2)), obj.nu, obj.N - 1));
        end
    end
end

