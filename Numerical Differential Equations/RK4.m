% Solving IVP with classic 4th-order Runge-Kutta method.
% 
% FUNCTION:   IVP_trap()
% 
% INPUT:
%     - system: 1st order system of the IVP
%     - t0: initial value
%     - init_vals: initial values of IVP
%     - h: time step
%     - n: number of steps
%     
% OUTPUT:
%     - t: history of time
%     - num_vals: history of numerical approximation
% 
% EXAMPLE:
%         Given IVP: (change into 1st order system)
%         y''- mu(1-y^2)y' + y = 0
%         y(0)=y0, y'(0)=v0
%         
%         system = @(w, mu, t) [w(2) , mu * (1-w(1)^2) * w(2) - w(1)] ;
%         mu = 0.5;
%         t0= 0;
%         y0 = 0.5;
%         v0 = 1;
%         init_vals = [y0, v0];
%         T = 30; % total time
%         h = 0.025;
%         n = T/h;
%         
%         [times, solns] = RK4( @(w,t) system(w, mu, t), t0, init_vals, h, n);

function [t, num_vals] = RK4(sys, t0, init_vals, h, n )

    step = 4;
    
    A = [0 0 0 0; 
        0.5 0 0 0;
        0 0.5 0 0;
        0 0 1 0];
    
    c = [0 0.5 0.5 1];
    
    b = [(1/6) (1/3) (1/3) (1/6)];
    
    sys_dim = size(init_vals, 2);
    
    k = zeros(step, sys_dim);
    
    num_vals = zeros(n+1, sys_dim);
    num_vals(1, :) = init_vals;
    
    t = zeros(n+1, 1);
    t(1) = t0;
    
    for i = 1 : n
        
        for j = 1: step
            k(j, :) = h * sys( init_vals + A(j, 1:(j-1)) * k(1:(j-1), :), t0+c(j)*h);
        end
        
        init_vals = init_vals + b*k;
        num_vals(i+1, :) = init_vals;
        t(i+1) = t(i) + h;
        
    end
    
    
end