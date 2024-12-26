classdef learn
    % This class consists of functions for DPR method.
    % Author:
    % Norihisa Namura (2022) 
    % namura.n.aa@m.titech.ac.jp
    %
    % If you use this code, please cite the following paper:
    % Estimating asymptotic phase and amplitude functions of limit-cycle oscillators from time series data
    % Norihisa Namura, Shohei Takata, Katsunori Yamaguchi, Ryota Kobayashi, and Hiroya Nakao
    % Phys. Rev. E 106, 014204 – Published 7 July 2022
    % https://doi.org/10.1103/PhysRevE.106.014204
    methods (Static)
        function dxdt = second_order_diff(x,dt)
            dxdt = zeros(size(x));
            dxdt(:,2:end-1) = (x(:,3:end) - x(:,1:end-2)) / (2*dt);
            dxdt(:,1) = (-3*x(:,1) + 4*x(:,2) - x(:,3)) / (2*dt);
            dxdt(:,end) = (3*x(:,end) - 4*x(:,end-1) + x(:,end-2)) / (2*dt);
        end

        function U = polynomial(x,p)
            % This function calculates vector of polynomials U(X) for M time steps.
            % input: 
            % - x: States (matrix: 2*M)
            % - p: Order of polynomials
            % output:
            % - U: Polynomials (matrix: 2*M)

            order = 1:p+1;

            N = length(x(1,:)); % dimension
            U = zeros(N,sum(order));
            for i = 1:p+1
                for j = 1:order(i)
                    U(:, sum(order(1:i)) -j+1 ) = x(1,:).^(j-1) .* x(2,:).^(order(i)-j); % powered n(i)-1
                end
            end
        end
        
        function [dUdx] = polynomial_diff(x,p)
            % This function calculates vector of partial derivative of polynomials for M time steps.
            % input: 
            % - x: States (matrix: 2*M)
            % - p: Order of polynomials
            % output:
            % - dUdx: Polynomials (matrix: M*P)

            order = 1:p+1;

            N = length(x(1,:)); % dimension
            dUdx = zeros(2,N,sum(order));
            for i = 1:p+1
                for j = 1:order(i)
                    if j < 2
                        dUdx(1,:, sum(order(1:i)) -j+1 ) = 0;
                    else
                        dUdx(1,:, sum(order(1:i)) -j+1 ) = (j-1) * x(1,:).^(j-2) .* x(2,:).^(order(i)-j); % powered n(i)-1
                    end
                    if order(i) - j - 1 < 0
                        dUdx(2,:, sum(order(1:i)) -j+1 ) = 0;
                    else
                        dUdx(2,:, sum(order(1:i)) -j+1 ) = x(1,:).^(j-1) .* ((order(i) - j) * x(2,:).^(order(i)-j-1)); % powered n(i)-1
                    end
                end
            end
        end

        function [H,f,l] = eval_phase(x,dxdt,mu,sig,omega,gamma,p)
            % This function calculates the objective function for estimating phase.
            % input: 
            % - x: Data used for calculation (matrix: 2*length)
            % - dxdt: Time derivative of data (matrix: 2*length)
            % - mu: Mean
            % - sig: Standard deviation
            % - omega: Natural frequency
            % - gamma: Weight for regularization
            % - p: Order of polynomials
            % output:
            % - x_lc: States of limit cycle
            % - theta_lc: Phases of limit cycle corresponding to "x_lc"
            % - x_lc_noise: States of limit cycle with noise

            U_x = learn.polynomial(x,p);
            [k,l] = size(U_x);
            
            pd = learn.polynomial_diff(x,p);
            dUdt = zeros(size(U_x));
            for i = 1:k
                dUdt(i,:) = dxdt(1,i) * pd(1,i,:) + dxdt(2,i) * pd(2,i,:);
            end
            
            % standardization
            U_x(:,2:end) = (U_x(:,2:end) - mu) ./ sig;
            dUdt(:,2:end) = dUdt(:,2:end) ./ sig;
            
            H_c = [dUdt,U_x*omega];
            H_s = [-U_x*omega,dUdt];
            C_Trans_C = H_c.'*H_c + H_s.'*H_s;
            
            H = 2*(C_Trans_C + gamma * eye(size(C_Trans_C)));
            f = zeros(length(H),1);
        end
        
        function theta = phase_model(X,z1,z2,p,mu,sig)
            % Estimate phase on input points using learned parameters
            U = learn.polynomial(X,p);
            U(:,2:end) = (U(:,2:end) - mu) ./ sig; % standardization
            
            fc = sum(U.*z1',2);
            fs = sum(U.*z2',2);
            theta = atan2(fs,fc) + 2*pi;
            theta = rem(theta,2*pi);
        end
        
        function [z1,z2,mu,sig] = learn_phase(x,dxdt,phase_abs,omega,p,gamma)
            % This function estimates the coefficient vectors for the phase function. 
            % input: 
            % - x: Data used for calculation (matrix: 2*length)
            % - dxdt: Time derivative of data (matrix: 2*length)
            % - phase_abs: State where phase is zero
            % - omega: Natural frequency
            % - p: Order of polynomials
            % - gamma: Weight for regularization (default: 0)
            % output:
            % - z1: 
            % - omega: Natural frequency

            if ~exist("gamma",'var')
                gamma = 0;
            end

            U_x = learn.polynomial(x,p);
            [~,mu,sig] = zscore(U_x(:,2:end));
            [H,f,l] = learn.eval_phase(x,dxdt,mu,sig,omega,gamma,p);
            
            Aeq_c = squeeze(learn.polynomial(phase_abs,p));
            Aeq_s = squeeze(learn.polynomial(phase_abs,p));
            
            Aeq_c(2:end) = (Aeq_c(2:end) - mu) ./ sig;
            Aeq_s(2:end) = (Aeq_s(2:end) - mu) ./ sig;
            Aeq = [Aeq_c,zeros(1,length(Aeq_c));zeros(1,length(Aeq_c)),Aeq_s];
            beq = [1,0].';
            
            options = optimoptions('quadprog','Display','iter');
            z = quadprog(H,f,[],[],Aeq,beq,[],[],[],options);
            %options = optimoptions('lsqlin','Display','iter');
            %z = lsqlin(2*H,f,[],[],Aeq,beq,[],[],[],options);
            z1 = z(1:l);
            z2 = z(l+1:2*l);
        end
    end
end
