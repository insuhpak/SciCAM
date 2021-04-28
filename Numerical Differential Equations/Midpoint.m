% Solving IVP of first order ODEs with two-step midpoint method.
% 
%     u' = f(u,t)
%     u(t_0) = u_0
%     
% 
% FUNCTION:   Midpoint
% 
% INPUT:
%     - func: ODE; f(u,t)
%     - t0: initial value
%     - x0: initial value; u_0
%     - h: time step
%     - n: number of steps
%     
% OUTPUT:
%     - t: history of time
%     - u: history of numerical approximation of u(t_n)
% 
% EXAMPLE:
%         g = @(u, t) -1*u;
%         t0= 0;
%         x0 = 1;
%         T = 1=20; % total time
%         h = 0.2;
%         n = T/h;
%         
%         [times, solns] = Midpoint(g, t0, x0, h, n);

function [t, u] = Midpoint(func, t0, x0, h, n)

    u = zeros(1,n);
    t = zeros(1,n);
    
    % Define u(1) by hand. 
    % Define time by hand.
    u(1) = x0;
    for k = 1: n
        t(k) = t0 + k*h;     %t_n = t_0 + nh
    end
    u(2) = exp(-h);
    
    % Euler's Method
    for k = 2: n-1
        u(k+1) = u(k-1) + 2*h*func(u(k),t(k));
    end 

end