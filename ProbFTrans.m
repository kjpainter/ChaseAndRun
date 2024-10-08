%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pij = ProbFTrans(i, j, params, varargin)
%
% Input arguments
%   i             ... variable who's transport phenomena are being
%                     described
%   j             ... variable who's concentration/density is determining
%                     the transport of variable i
%   params        ... [struct] as provided by ProbGetParams().
%   varargin      ... {1} contains the results for the variables
%                          -> varargin{1}(:,:,c_m) is a matrix containing
%                          the concentration of MSC in all x and y 
%                     {2} is the patchId 
%                     {3} is either x or y
% Output arguments: None
% 
% Return function values of parameter function p_{i,j}(...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pij, terms] = ProbFTrans(i, j, params, varargin)


if ((i == params.eq.n1) & (j == params.eq.n1))
    n1 = varargin{1}(:,params.eq.n1); 
    n2 = varargin{1}(:,params.eq.n2);
    pij = params.n1_D*ones(size(n1));  
elseif ((i == params.eq.n2) & (j == params.eq.n2))
    n1 = varargin{1}(:,params.eq.n1); 
    n2 = varargin{1}(:,params.eq.n2);
    pij = params.n2_D*ones(size(n2));    
else
  error(['ProbFTrans(' num2str(i) ',' num2str(j) ...
	 ') should not have been called.']);
end

return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
