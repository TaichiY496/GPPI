classdef stuart_landau
    % This class defines the stuart_landau oscillator.
    % If you want to try other oscillators, 
    % make another class similar to this
    % and write the equation and parameters there.
    % Author:
    % Taichi Yamamoto
    % yamamoto-taichi913@g.ecc.u-tokyo.ac.jp
    properties
        % parameters of dynamics
        alpha = 2;
        beta = 1;
        % parameters for plot
        x_lim = [-1.6,1.6];
        y_lim = [-1.6,1.6];
        delta = 0.02; % mesh size
        plot_ind = [1,2]; % plot (x(1), x(2))
        x_tick = [-1.6,-0.8,0,0.8,1.6];
        y_tick = [-1.6,-0.8,0,0.8,1.6];
        x_name = "$x$"; y_name = "$y$";
        % parameters for data generation
        time = 2.5;
        n = 500;
        x_lc_0 = [1;0]; % origin of phase function is approximately here (not accurate)
        varsigma_phase = [-0.2,0,0.2];
        name = "sl";
        dt = 0.005;
        % threshold for period measuring:
        % (noisy timescale) << th_M*dt << (period) should hold.
        % see utils.period_noise()
        th_M = 100;
        % parameters estimated numerically (just define the variables and put zeros).
        T = 0; omega = 0;
        scale = [];
        x_lc = [];
    end
    methods
        function f = func(obj,x)
            % dynamics ODE
            % x: 2 * N, N state vectors
            % f: 2 * N, dx/dt values for each state vector
            f1 = x(1,:) - obj.alpha * x(2,:) - (x(1,:) - obj.beta * x(2,:)) .* (x(1,:).^2 + x(2,:).^2);
            f2 = obj.alpha * x(1,:) + x(2,:) - (obj.beta * x(1,:) + x(2,:)) .* (x(1,:).^2 + x(2,:).^2);
            f = [f1;f2];
        end
        
        function r = cond(~, x1, x2)
            % define Poincare surface for measuring period
            % should return 1 if and only if the orbit (x1->x2) passes the surface
            % surface: y = 0 & x > 0
            r = x1(1) > 0 && (x1(2) * x2(2) < 0);
        end

        function theta = phase_calc(obj, x)
            % Analytical phase function
            % You don't need to define this function to try GPPI.
            % x : 2 * N, N state vectors
            % theta : 1 * N, phase values for each state vector
            r2 = sum(x.^2, 1); % r^2
            theta = atan2(x(2,:), x(1,:)) - log(r2) * obj.beta/2;
            theta = rem(theta, 2*pi);
            index = theta < 0;
            theta(index) = theta(index) + 2*pi;
        end
    end
end