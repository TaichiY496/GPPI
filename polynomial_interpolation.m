% Noisy time series interpolation for DPR method.
% Author:
% Norihisa Namura (2022)
% namura.n.aa@m.titech.ac.jp
function [dxdt] = polynomial_interpolation(x,M,dt,windowsize)
    % This function calculates the time derivative of the data.
    % input: 
    % - x: Data (tensor: 2*M*length)
    % - M: Number of data points per one time series
    % - dt: Time interval
    % - windowsize: Windowsize of filtering
    % output: 
    % - dxdt: Time derivative of data (matrix: 2*length)
    t = 0:dt:M*dt-dt;
    t = t.';
    
    n = numel(x)/2/M;
    x1 = squeeze(x(1,:,:));
    x2 = squeeze(x(2,:,:));
    td = zeros(M-windowsize+1,windowsize);
    xd1 = zeros(M-windowsize+1,windowsize,n);
    xd2 = zeros(M-windowsize+1,windowsize,n);
    for i = 1:windowsize
        td(:,i) = t(i:end-windowsize+i);
        xd1(:,i,:) = x1(i:end-windowsize+i,:);
        xd2(:,i,:) = x2(i:end-windowsize+i,:);
    end

    p = 1; % order of regression
    s = size(xd1); % dim: M-windowsize+1,windowsize,n
    c1 = zeros(s(1),s(3),p+1);
    c2 = zeros(s(1),s(3),p+1);
    for j = 1:n
        for i = 1:M-windowsize+1
            A1 = polycat(td(i,:).',p);
            A2 = polycat(td(i,:).',p);
            b1 = squeeze(xd1(i,:,j)).';
            b2 = squeeze(xd2(i,:,j)).';
            c1(i,j,:) = pinv(A1) * b1;
            c2(i,j,:) = pinv(A2) * b2;
        end
    end
    
    dx1dt = zeros(s(1),s(3));
    dx2dt = zeros(s(1),s(3));
    
    l = fix(windowsize/2) + 1;
    t = t(l:end-l+1,:); 
    for j = 1:n
        for i = 1:p
            dx1dt(:,j) = dx1dt(:,j) + i * squeeze(c1(:,j,i+1)) .* power(t,i-1);
            dx2dt(:,j) = dx2dt(:,j) + i * squeeze(c2(:,j,i+1)) .* power(t,i-1);
        end
    end
    dxdt = cat(3,dx1dt,dx2dt);
    dxdt = permute(dxdt,[3,1,2]);
    %dxdt = dxdt(:,:);
end

function polycat = polycat(x,p)
    polycat = zeros(length(x),p+1);
    for i = 1:p+1
        polycat(:,i) = power(x,i-1);
    end
end