%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y0 = ProbFy0(x, params)
%
% Input arguments
%   x             ... spatial coordinates
%   params        ... [struct] as provided by ProbGetParams().
% Output arguments: None
% 
% Returns the initial values of all species in (x(1)).  Can be made
% space-dependent (to have inhomogeneous initial conditions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y0 = ProbFy0(x, params)

switch params.icswitch    
    case 1
        sig = 0.1;
        mass = 10.0;
        y0 = [(mass/(sig*sqrt(2.0*pi)))*exp(-0.5*((x-2.5)^2.0)/(sig^2.0));(mass/(sig*sqrt(2.0*pi)))*exp(-0.5*((x-1.5)^2.0)/(sig^2.0))];
    case 2
        y0 = [1.0+0.2*(rand-0.5);1.0+0.2*(rand-0.5)];
    case 3
        sig = 0.2;
        mass = 10.0;
        y0 = [(1.0+0.2*(rand-0.5))*(mass/(sig*sqrt(2.0*pi)))*exp(-0.5*((x-2.5)^2.0)/(sig^2.0));(1.0+0.2*(rand-0.5))*(mass/(sig*sqrt(2.0*pi)))*exp(-0.5*((x-1.5)^2.0)/(sig^2.0))];
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
