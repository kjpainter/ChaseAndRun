function [G] = ProbFNonLocal(i,params,varargin)


switch i
 case params.eq.n1
  n1 = varargin{1}(:,params.eq.n1);
  n2 = varargin{1}(:,params.eq.n2);
  switch params.hswitch
      case 1
          G1 = params.S_n1n1*n1;
          G2 = params.S_n1n2*n2;
  end
 case params.eq.n2
  n1 = varargin{1}(:,params.eq.n1);
  n2 = varargin{1}(:,params.eq.n2);
  switch params.hswitch
      case 1
          G1 = params.S_n2n1*n1;
          G2 = params.S_n2n2*n2;
  end
 
 otherwise
  error('Subroutine should not have been called.');
end
G = {G1, G2};
return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
