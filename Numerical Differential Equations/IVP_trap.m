% Solving IVP of first order ODEs with TRAPEZOIDAL method.
% 
%     u' = f(u,t)
%     u(t_0) = u_0
%     
% NOTE: this function uses newton_solver.m authored by Hongyun Wang
% 
% FUNCTION:   IVP_trap()
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
%         g = @(u, t) -(10^6)*sinh(u-cos(t));
%         t0= 0;
%         x0 = 0;
%         T = 10; % total time
%         h = 0.1;
%         n = T/h;
%         
%         [times, solns] = IVP_trap(g, t0, x0, h, n);

function [t, u] = IVP_trap(func, t0, x0, h, n)

    tol = 1e-16;

    % Store u and t-values in arrays.
    u = zeros(1,n);
    t = zeros(1,n);
    
    % Define initial condition by hand. 
    % Define time by hand.
    u(1) = x0;
    t(1) = t0;
    for k = 2:n
        t(k) = (k-1)*h;     %t_n = t_0 + nh
    end
    
    % Function used in Newton's method.
    Newton_func = @(x, c1, c2, c3, h) c1 + (h/2)*func(c1, c2) + (h/2)*func(x,c3) - x;
    
    % Trapezoidal Method
    for k = 1: n-1
        
        % Define the constant used in Newton_g.
        c1 = u(k);
        c2 = t(k);
        c3 = t(k+1);
        
        % Call Newton's method to solve for u(k+1).
        [xr, iter, flag] = newton_solver( @(x) Newton_func(x, c1, c2, c3, h), x0, tol);
        
        % Trapezoidal Method
        u(k+1) = u(k) + (h/2)*( func(u(k), t(k)) + func(xr, t(k+1)) );
    end 
    
end 