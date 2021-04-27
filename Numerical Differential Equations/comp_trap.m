% Approximate an integral using the COMPOSITE TRAPEZOIDAL rule on a uniform grid.
% 
% FUNCTION: comp_trap
% 
% INPUT: 
%     - func:  function you want to integrate
%     - a: lower bound of integral
%     - b: upper bound of integral
%     - N: number of grids
%     
% OUTPUT:
%     - num = approximate solution to integral
%     
%     EXAMPLE: func = @(x) sqrt( 1 + (cos(x))^2 ) * exp( cos(x) );
%              a = 0;
%              b = 0;
%              N = 4;
%              
%             [soln] = comp_trap(func, a, b, N);


function [num] = comp_trap(func, a, b, N)
    
        % Step size
        h = (b-a) / N;
        
        % Define each point based on step.
        x = a:h:b;
        dimx = size(x);
        n = dimx(2);
        
        % Composite trapeziodal method.
        sum = func(x(1)) + func(x(n));
        
        for j = 2: n-1
            sum = sum + 2*func(x(j));
        end 
        
        num = sum * h * 0.5;
        
end 
