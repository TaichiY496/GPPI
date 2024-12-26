classdef funcs
    % This class consists of functions for:
    % - simulation of dynamics
    % - numerical phase function calculation
    % - estimation by proposed method GPPI
    % - evaluation of results
    % Author:
    % Taichi Yamamoto
    % yamamoto-taichi913@g.ecc.u-tokyo.ac.jp
    methods (Static)
        function R2 = coefficient_determination(y,y_est)
            % This function calculates the coefficient of determination.
            % input: 
            % - y: True signals (vector)
            % - y_est: Estimated signals (vector with same length as y)
            % output:
            % - R2: Coefficient of determination

            y_mean = mean(y);
            SSR = sum((y - y_est).^2);
            SST = sum((y - y_mean).^2);
            R2 = 1 - SSR / SST;
        end

        function theta = theta_adjust(theta)
            % limit theta to [-pi, pi]
            cond = theta > pi;
            theta(cond) = theta(cond) - 2*pi;
            cond = theta < -pi;
            theta(cond) = theta(cond) + 2*pi;
        end
        
        function fx = runge_kutta_4(x,dt,cla)
            % RK4 method.
            % input: 
            % - x: State
            % - dt: Time interval
            % - cla: Class of a limit-cycle oscillator
            % output:
            % - fx: State at next time
            ka = cla.func(x);
            kb = cla.func(x + dt*ka/2);
            kc = cla.func(x + dt*kb/2);
            kd = cla.func(x + dt*kc);
            fx = x + dt*(ka + 2*kb + 2*kc + kd)/6;
        end

        function [T,omega] = period_noise(eta,windowsize,cla,s,n)
            % estimate the period and natural frequency from noisy data.
            % input: 
            % - eta: Strength of noise
            % - windowsize: Windowsize of filtering
            % - cla: Class of limit-cycle oscillators
            % - s : Rotate s times before measure
            % - n : Rotate n times for measure
            % output:
            % - T: Period
            % - omega: Natural frequency

            dt = cla.dt;
            M = 100000000; % max iteration
            ini = cla.x_lc_0; % initial point
            d = size(ini,1); % dim of dynamics
            xx = zeros(d,windowsize+1); % clean states
            xx_filtered = zeros(d,3); % filtered states

            % filter parameters
            a = 1;
            b = (1/windowsize)*ones(1,windowsize);
            
            xx(:,1) = ini;
            for i = 1:windowsize
                xx(:,i+1) = funcs.runge_kutta_4(xx(:,i),dt,cla);
            end
            xx_noise = utils.add_noise(xx,eta); % noisy states

            % filter states
            tmp = filter(b,a,xx_noise,[],2);
            xx_filtered(:,2:3) = tmp(:,end-1:end); 
            
            xx = xx(:,end);
            xx_noise = xx_noise(:,2:end);
            
            count = 0;
            icount = 0;
            for i = 1:M-1
                xx_noise(:,1:end-1) = xx_noise(:,2:end);
                xx_filtered(:,1:end-1) = xx_filtered(:,2:end);
                xx = funcs.runge_kutta_4(xx,dt,cla); % time evolution of clean states
                xx_noise(:,end) = utils.add_noise(xx,eta); % add noise
                tmp = filter(b,a,xx_noise,[],2);
                xx_filtered(:,end) = tmp(:,end);
                cond_var1 = cla.cond(xx_filtered(:,1),xx_filtered(:,2));
                % Avoid too much lap-counts due to noise
                if i > icount + cla.th_M && cond_var1
                    count = count + 1;
                    icount = i;
                    if count == s
                        i_start = i;
                    elseif count == s+n
                        T = (i - i_start) * dt / n;
                    	break;
                    end
                end
            end
            omega = 2*pi/T;
        end

        function [T,omega] = period(cla,s,n)
            % the period and natural frequency by exact data.
            % input: 
            % - cla: Class of limit-cycle oscillators
            % - s : Rotate s times before measure
            % - n : Rotate n times for measure
            % output:
            % - T: Period
            % - omega: Natural frequency

            dt = cla.dt;
            th_M = cla.th_M;
            M = 100000000; % max iteration
            x2 = cla.x_lc_0; % initial point
            d = size(x2,1); % dim of dynamics
            x1 = zeros(d,1);

            count = 0;
            icount = 0;
            for i = 1:M-1
                x1(:) = x2(:);
                x2(:) = funcs.runge_kutta_4(x1,dt,cla); % time evolution of clean states
                cond_var1 = cla.cond(x1,x2);
                if i > icount + th_M && cond_var1
                    count = count + 1;
                    icount = i;
                    if count == s
                        i_start = i;
                    elseif count == (s+n)
                        T = (i - i_start) * dt / n;
                    	break;
                    end
                end
            end
            omega = 2*pi/T;
        end
        
        function [x_lc,theta_lc,x_lc_noise,theta_lc_noise] = phase_map(eta,windowsize,cla,omega_noise)
            % calculates the limit cycle and its phase.
            % input: 
            % - eta: Strength of noise
            % - windowsize: Windowsize of filtering
            % - cla: Class of a limit-cycle oscillator
            % - omega_noise: the frequency estimated by data
            % 
            % output:
            % - x_lc: States of limit cycle
            % - theta_lc: Phases of limit cycle corresponding to "x_lc"
            % - x_lc_noise: States of limit cycle with noise
            % - theta_lc_noise: Phases of limit cycle estimated by omega_noise

            dt = cla.dt; T = cla.T;
            initial = cla.x_lc_0;
            d = size(initial,1);
            n = floor(T/dt);
            % run 2 periods
            xx = zeros(d,n*2+windowsize); xx(:,1) = cla.x_lc_0;
            for i=1:size(xx,2)-1
                xx(:,i+1) = funcs.runge_kutta_4(xx(:,i),dt,cla);
            end
            xx_noise = utils.add_noise(xx,eta);
            xx_filtered = filter(ones(1,windowsize)/windowsize, 1, xx_noise, [], 2);
            l = (windowsize-1)/2;
            x_lc = xx(:,n+1:n*2);
            x_lc_noise = xx_filtered(:,n+l+1:n*2+l);
            theta_lc = cla.omega * dt * (0:n-1);
            theta_lc_noise = omega_noise * dt * (0:n-1);
        end

        function [x_valid, theta_data, valid] = phase_data(x,cla,x_lc,theta_lc,scale,precise,th)
            % estimates the phase of datapoints (GPPI step2)
            dt = cla.dt; omega = cla.omega;
            m = size(x, 2);
            x_last = squeeze(x(:,m,:));
            [theta_last, invalid] = funcs.phase_cycle(x_last,x_lc,theta_lc,scale,precise,th);
            valid = ~invalid;
            x_valid = x(:,:,valid);
            n_valid = size(x_valid,3);
            theta_data = zeros(m,n_valid);
            theta_data(m,:) = theta_last(valid);
            for i = 1:m-1
                theta_data(i,:) = theta_data(m,:) - omega*(m-i)*dt;
                theta_data(i,:) = rem(theta_data(i,:)+10*pi, 2*pi);%fix to [0,2pi]
            end
        end

        function [theta_map, invalid] = phase_cycle(X,x_lc,theta_lc,scale,precise,th)
            % estimates phase on limit cycle.
            % normalize by scale
            x_lc_n = x_lc ./ scale;  xx = X ./ scale;
            index = dsearchn(x_lc_n.',xx.');
            theta_map = theta_lc(index);
            l = sum((x_lc_n(:, index) - xx).^2, 1);
            invalid = l > th;
            if precise
                % linear interpolation
                L = size(x_lc_n,2);
                ind2 = mod(index,L)+1;
                ind3 = mod(index-2,L)+1;
                cond = sum((x_lc_n(:, ind2) - xx).^2, 1) > sum((x_lc_n(:, ind3) - xx).^2, 1);
                ind2(cond) = ind3(cond);
                a = xx - x_lc_n(:,index);
                b = x_lc_n(:,ind2) - x_lc_n(:,index);
                k = sum(a.*b, 1) ./ sum(b.^2, 1);
                dtheta = funcs.theta_adjust(theta_lc(ind2) - theta_lc(index));
                theta_map = theta_map + k .* dtheta;
                %theta_map = rem(theta_map+2*pi,2*pi);
                invalid = sum((a - k.*b).^2,1) > th;
            end
            
        end

        function [theta_map, invalid] = phase_numerical(X,x_lc,theta_lc,cla,M,precise)
            % This function calculates the phase function numerically.
            % input: 
            % - X: d * N, States to obain phase
            % - x_lc: d * L, States on limit cycle
            % - theta_lc: L * 1, Phases on limit cycle
            % - cla: Class of a limit-cycle oscillator
            % - M: number of RK4 steps to simulate
            % - precise: if 1, the phase on limit cycle is linear interpolated
            % output:
            % - theta_map: N * 1, Phase function numerically
            % - invalid: N * 1 (logic), theta_map(i) might be invalid
            %   when invalid(i) == 1

            dt = cla.dt; scale = cla.scale;
            xx = X;
            for i = 1:M
                xx = funcs.runge_kutta_4(xx,dt,cla);
            end
            [theta_map, invalid] = funcs.phase_cycle(xx,x_lc,theta_lc,scale,precise,0.0001);
            theta_map = theta_map +2*pi - rem(cla.omega*M*dt,2*pi);
            theta_map = rem(theta_map, 2*pi);
        end

        function [sin_model,cos_model] = learn_GP(X,theta,scale,GPoptions)
            X = X./scale;
            sin_model = fitrgp(X', sin(theta), GPoptions{:});
            cos_model = fitrgp(X', cos(theta), GPoptions{:});
        end

        function [theta, fs, fc, sd] = phase_GP(X,sin_model,cos_model,scale)
            X = X ./ scale;
            [fs, fssd] = predict(sin_model, X');
            [fc, fcsd] = predict(cos_model, X');
            theta = rem(atan2(fs,fc) + 2*pi, 2*pi);
            sd = sqrt(cos(theta).^2 .* fssd + sin(theta).^2 .* fcsd);
        end

        function ind_V_lc = phase_golf(target,dtheta)
            % finds an index close to target with less cost.
            % dtheta : n_V * n_lc
            if target > max(dtheta(1,:))
                % can't reach to target(>0) with cost 1 -> assume max
                [ps,argps] = max(dtheta,[],2);
                [tmp,ind_V] = max(ps);
                if target < tmp
                    [~,ind_V] = min(abs(ps-target));
                end
                ind_lc = argps(ind_V);
            elseif target > min(dtheta(1,:))
                % can reach to target with cost 1
                ind_V = 1;
                [~,ind_lc] = min(abs(dtheta(1,:)-target));
            else
                % can't reach to target(<0) with cost 1 -> assume min
                [ps,argps] = min(dtheta,[],2);
                [tmp,ind_V] = min(ps);
                if target > tmp
                    [~,ind_V] = min(abs(ps-target));
                end
                ind_lc = argps(ind_V);
            end
            ind_V_lc = [ind_V, ind_lc];
        end
    end
end