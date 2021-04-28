% Solving IVP of first order ODEs with Euler's Method.
% 
%     u' = f(u,t)
%     u(t_0) = u_0
%     
% NOTE: under the assumption t_0 = 0
% 
% FUNCTION:   Eulers
% 
% INPUT:
%     - func: ODE; f(u,t)
%     - x0: initial value; u_0
%     - h: time step
%     - n: number of steps
%     
% OUTPUT:
%     - t: history of time
%     - u: history of numerical approximation of u(t_n)
% 
% EXAMPLE:
%         g = @(u, t) -(10^6)*sinh(u-cos(t));
%         x0 = 0;
%         T = 2^(-10); % total time
%         h = 2^(-18);
%         n = T/h;
%         
%         [times, solns] = Eulers(g, x0, h, n);

function [t, u] = Eulers(func, x0, h, n)
    
    u = zeros(1,n);
    t = zeros(1,n);
    
    % Define u(1) by manually. 
    u(1) = x0;
    
    % Discretization.
    % t_n = t_0 + nh
    % Under the assumption t_0 = 0
    for k = 1: n
        t(k) = k*h;     
    end
    
    % Euler's Method
    for k = 1: n-1
        u(k+1) = u(k) + h*func(u(k),t(k));
    end 
    
end
