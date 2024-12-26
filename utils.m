classdef utils
    % This class mainly consists of functions for
    % - data sampling
    % - data generation
    % - visualization
    % - other utils
    % Author:
    % Taichi Yamamoto
    % yamamoto-taichi913@g.ecc.u-tokyo.ac.jp
    methods (Static)
        function initials = gen_initials(cla, n)
            % This function generates initial points for time evolution.
            initials = [cla.x_lim(2)-cla.x_lim(1),0;0,cla.y_lim(2)-cla.x_lim(1)] * rand([2,n]);
            initials(1,:) = initials(1,:) + cla.x_lim(1);
            initials(2,:) = initials(2,:) + cla.y_lim(1);
        end
        
        function [x,dxdt,x_data] = gen_data(eta,M,initials,windowsize,dt,cla)
            % This function generates simulation data.
            % input: 
            % - eta: Strength of noise
            % - M: Number of data points per one time series
            % - initials: Initial points for time evolution
            % - windowsize: Windowsize of filtering
            % - dt: Time interval
            % - cla: Class of a limit-cycle oscillator
            % 
            % output:
            % - x: Data used for calculation (matrix: 2*(M-a)*n)
            % - dxdt: Time derivative of data (matrix: 2*(M-a)*n)
            % - x_data: Whole data for used (matrix: 2*Mn)

            dn = size(initials);
            x_tmp = zeros(dn(1),M,dn(2)); % n: number of initial point 
            x_tmp(:,1,:) = initials; % set initial states
            
            % time evolution
            for i = 1:M-1
                x_tmp(:,i+1,:) = funcs.runge_kutta_4(squeeze(x_tmp(:,i,:)),dt,cla);
            end
            x_tmp = utils.add_noise(x_tmp,eta); % add noise
            x_data = x_tmp; % log used data
            
            if eta > 0 % if noise exists
                l = fix(windowsize/2) + 1;
                x = x_tmp(:,l:end-l+1,:);
                [dxdt] = polynomial_interpolation(x_tmp,M,dt,windowsize);
            else % no noise
                x = x_tmp(:,2:M-1,:);
                dxdt = (x_tmp(:,3:M,:) - x_tmp(:,1:M-2,:)) / (2*dt);
            end
        end
        
        function x_noise = add_noise(x,eta)
            x_noise = x + randn(size(x)) * eta;
        end
        
        function [X,X1,X2,area_size] = mesh_grid(cla)
            % M * N data grid is determined by cla
            % X : 2 * MN
            % X1, X2 : M * N
            % area_size = [M, N]
            x1 = cla.x_lim(1):cla.delta:cla.x_lim(2);
            x2 = cla.y_lim(1):cla.delta:cla.y_lim(2);
            [X1,X2] = meshgrid(x1,x2);
            area_size = size(X1);
            X = horzcat(X1(:),X2(:));
            X = X.';
        end

        function [ind, fast, slow] = sample_data_fast(x,cla,fast_time,slow_time,t_h)
            % depends on cla.fast_region(x)
            % (use only for HH model data)
            dt = cla.dt;
            m = round(t_h/dt); %size(x,2);
            fs = round(fast_time/dt); ss = round(slow_time/dt);
            cond = cla.fast_region(x);
            slow = zeros(size(cond)); slow(1:fs:m,:) = 1;
            fast = zeros(size(cond)); fast(1:ss:m,:) = 1; 
            cond = cond(:); slow = ~cond & slow(:); fast = cond & fast(:);
            ind = slow | fast;
        end

        % functions for showing figures
        function fig_axis(cla,titletxt)
            % Setup a panel
            if cla.name == "sl"
                axis equal
            end
            xticks(cla.x_tick); xlim(cla.x_lim);
            yticks(cla.y_tick); ylim(cla.y_lim);
            xticklabels(compose("%.1f", cla.x_tick));
            yticklabels(compose("%.1f", cla.y_tick));
            t = title(titletxt); t.FontSize = 20;
            ax = gca; ax.FontSize = 20;
            ax.TickLabelInterpreter = "latex";
            ax.TickDir = "out";
            ax.Layer = "top";
            ax.Box = "off"; ax.LineWidth = 1;
            t = xlabel(cla.x_name, Interpreter="latex"); t.FontSize = 30;
            t = ylabel(cla.y_name, Interpreter="latex"); t.FontSize = 30;
        end
        
        function fig = fig_data(cla,x,sz,color,titletxt,show_cycle)
            % Scatter plot of data points
            % color can be array of colors, resulting in colorful plot
            x_lc = cla.x_lc;
            ind = cla.plot_ind;
            fig = figure(); fig.Position(3:4) = [400,300];
            hold on
            scatter(x(ind(1),:),x(ind(2),:),sz,color, "filled");
            if show_cycle
                plot(x_lc(ind(1),:),x_lc(ind(2),:),"k","LineWidth",2);
            end
            hold off
            utils.fig_axis(cla,titletxt);
            set(gcf,"visible","on");
        end

        function x_show = reduce_data(x,x_lc,scale,th1,th2,iter)
            % reduces too many data points for scatter plot
            x_n = x ./ scale; x_lc_n = x_lc ./ scale;
            index = dsearchn(x_lc_n.',x_n.');
            l = sum((x_lc_n(:, index) - x_n).^2, 1);
            index = l > th1;
            x_show = x(:,index);
            for i=1:iter
                x1 = x_show(:,1:2:end);
                x2 = x_show(:,2:2:end); 
                x1_n = x1 ./scale; x2_n = x2 ./ scale;
                index = dsearchn(x1_n.',x2_n.');
                l = sum((x1_n(:, index) - x2_n).^2, 1);
                index = l > th2;
                x_show = [x1, x2(:,index)];
            end
            fprintf("data reduced: %d -> %d", size(x,2), size(x_show,2))
        end
        
        function fig = fig_phase_heat(cla,theta,titletxt)
            x_lc = cla.x_lc;
            fig = figure(); fig.Position(3:4) = [400,300];
            hold on
            colormap("hsv")
            imagesc("XData",cla.x_lim, "YData",cla.y_lim, "CData",theta)
            ind = cla.plot_ind;
            plot(x_lc(ind(1),:),x_lc(ind(2),:),"k","LineWidth",2);
            colorbar('Ticks',[0,pi/2,pi,3*pi/2,2*pi],...
                'TickLabels',["$0$","$\frac{1}{2}\pi$","$\pi$","$\frac{3}{2}\pi$","$2\pi$"],...
                TickLabelInterpreter='latex');
            clim([0,2*pi]);
            hold off
            utils.fig_axis(cla,titletxt);
            set(gcf,"visible","on");
        end

        function fig = fig_phase_error(cla,diff,colors,titletxt)
            x_lc = cla.x_lc;
            ind = cla.plot_ind;
            fig = figure(); hold on
            tmp = diff(:);
            tmp = funcs.theta_adjust(tmp);
            tmp = reshape(abs(tmp), size(diff));
            colormap(colors)
            imagesc("XData",cla.x_lim, "YData",cla.y_lim, "CData",tmp)
            plot(x_lc(ind(1),:),x_lc(ind(2),:),"k","LineWidth",2);
            colorbar('Ticks',[0,pi/4,pi/2],...
                'TickLabels',["$0$","$\frac{1}{4}\pi$","$\frac{1}{2}\pi$"],...
                TickLabelInterpreter='latex');
            clim([0,pi/2]);
            hold off
            utils.fig_axis(cla,titletxt);
            set(gcf,"visible","on");
        end
    end
end