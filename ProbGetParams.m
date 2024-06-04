%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [grd, tdr, params] = ProbGetParams(varargin)
%
% Input arguments
%   varargin      ... 
% Output arguments
%   grd           ... information on grid composition and size
%   tdr           ... information on the tdr system
%   params        ... parameter values of the problem
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grd, tdr, params] = ProbGetParams(varargin)

grd    = struct();
tdr    = struct();
params = struct();

%
% The model consists of 1 PDEs implemented in the following order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indices of PDEs                %%%%
params.eq.n1  = 1;
params.eq.n2  = 2;
%%% End of Indices of PDEs         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% process input parameters
switch length(varargin)
  case 0
    inputStruct.filled = false;
  case 1
    inputStruct = varargin{1};
    inputStruct.filled = true;
  otherwise
    error('Too many arguments in call to ProbGetParams().');
end


%
% initialise struct params
%
params.domainlength = setValue('domainlength', 400, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.domainlength = ' ...
      num2str(params.domainlength) '.']);
params.gridCells = setValue('gridCells', 25, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.gridCells = ' ...
      num2str(params.gridCells) ' grid cells per unit length.']);
params.BCs = setValue('BCs', 'pp', inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.BCs = ' ...
      num2str(params.BCs) '.']);
params.alpha_n1n1 = setValue('alpha_n1n1', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.alpha_n1n1 = '...
      num2str(params.alpha_n1n1) '.']);
params.alpha_n1n2 = setValue('alpha_n1n2', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.alpha_n1n2 = '...
      num2str(params.alpha_n1n2) '.']);
params.alpha_n2n1 = setValue('alpha_n2n1', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.alpha_n2n1 = '...
      num2str(params.alpha_n2n1) '.']);
params.alpha_n2n2 = setValue('alpha_n2n2', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.alpha_n2n2 = '...
      num2str(params.alpha_n2n2) '.']);  
params.xi_n1n1 = setValue('xi_n1n1', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.xi_n1n1 = '...
      num2str(params.xi_n1n1) '.']);
params.xi_n1n2 = setValue('xi_n1n2', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.xi_n1n2 = '...
      num2str(params.xi_n1n2) '.']);
params.xi_n2n1 = setValue('xi_n2n1', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.xi_n2n1 = '...
      num2str(params.xi_n2n1) '.']);
params.xi_n2n2 = setValue('xi_n2n2', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.xi_n2n2 = '...
      num2str(params.xi_n2n2) '.']);  
params.r = setValue('r', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.r = '...
      num2str(params.r) '.']);
params.delta = setValue('delta', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.delta = '...
      num2str(params.delta) '.']);
params.fswitch = setValue('fswitch', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.fswitch = '...
      num2str(params.fswitch) '.']);
params.hswitch = setValue('hswitch', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.hswitch = '...
      num2str(params.hswitch) '.']);
params.oswitch = setValue('oswitch', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.oswitch = '...
      num2str(params.oswitch) '.']);
params.icswitch = setValue('icswitch', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.icswitch = '...
      num2str(params.icswitch) '.']);
params.index = setValue('index', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.index = '...
      num2str(params.index) '.']);
  
%
% 1. initialise struct tdr
%
tdr.versionString = '1.6pre-AdhPack-1.2';
tdr.size   = 2; % number of species (PDEs)
%tdr.tvec   = [[0:0.25:5],[6:1:20],[25:5:40]]; % initial, output and final time.
tdr.tvec   = [0:0.025:50.0];
%%% CODES for describing a function f %%%%
Zero     = 0;              % f \equiv 0
Const    = 1;              % f \equiv const
DependsT = 2;              % f = f(t)
DependsS = 4;              % f = f(x,y)
DependsU = 8;              % f = f(u)
DependsTS = DependsT + DependsS; %  f = f(t, x, y)
DependsTU = DependsT + DependsU; %  f = f(t, u)
DependsSU = DependsS + DependsU; %  f = f(x, y, u)
DependsTSU = DependsT + DependsS + DependsU; %  f = f(t, x, y, u)

% describe reaction function ProbFReac.m 
% if only Zero entries then ProbFReac.m is never called
tdr.FReac.depends = Zero * ones(tdr.size, 1); 
tdr.FReac.depends(params.eq.n1) = [ DependsU ];
tdr.FReac.vec = true;  % true if ProbFReac.m can be called vectorized

% describe transport function ProbFTrans.m
% if only Zero entries then ProbFTrans.m is never called
tdr.FTrans.depends = Zero * ones(tdr.size, tdr.size)
tdr.FTrans.depends(params.eq.n1, params.eq.n1) = DependsU; % const. diff. of n1
tdr.FTrans.depends(params.eq.n2, params.eq.n2) = DependsU; % const. diff. of n2
%

tdr.FTrans.vec = true; % true if ProbFTrans.m can be called vectorized

% describe nonlocal function ProbFNonLocal.m
% if only Zero entries then ProbFNonLocal.m is never called
% (ProbFNonLocal.m encodes the function g under the integral)
tdr.FNonLocal.depends(params.eq.n1) = DependsU;
tdr.FNonLocal.depends(params.eq.n2) = DependsU;

% describe nonlocal function ProbFNonLocal2.m
% if only Zero entries then ProbFNonLocal2.m is never called and
% value 1 is assumed.
% (ProbFNonLocal2.m encodes the function p multiplying the integral)
tdr.FNonLocal.depends2 = Zero * ones(tdr.size, 1);

switch params.oswitch
    case 1
        params.S_n1n1=params.alpha_n1n1/(2.0*params.xi_n1n1);
        params.S_n1n2=params.alpha_n1n2/(2.0*params.xi_n1n2);
        params.S_n2n1=params.alpha_n2n1/(2.0*params.xi_n2n1);
        params.S_n2n2=params.alpha_n2n2/(2.0*params.xi_n2n2);
        
        tdr.FNonLocal.R          = zeros(tdr.size); 
        tdr.FNonLocal.R(params.eq.n1,params.eq.n1) = params.xi_n1n1;
        tdr.FNonLocal.R(params.eq.n1,params.eq.n2) = params.xi_n1n2;
        tdr.FNonLocal.R(params.eq.n2,params.eq.n1) = params.xi_n2n1;
        tdr.FNonLocal.R(params.eq.n2,params.eq.n2) = params.xi_n2n2;
        maxrange = max(max(tdr.FNonLocal.R));
% describe Omega function (the same for all)
        tdr.FNonLocal.OmegaStr = cell(tdr.size);
        tdr.FNonLocal.OmegaStr{params.eq.n1,params.eq.n1} = ['@(r)(' num2str(params.xi_n1n1,'%25.16d') '* (r <= ' num2str(tdr.FNonLocal.R(params.eq.n1,params.eq.n1),'%25.16d') ') )'];
        tdr.FNonLocal.OmegaStr{params.eq.n1,params.eq.n2} = ['@(r)(' num2str(params.xi_n1n2,'%25.16d') '* (r <= ' num2str(tdr.FNonLocal.R(params.eq.n1,params.eq.n2),'%25.16d') ') )'];
        tdr.FNonLocal.OmegaStr{params.eq.n2,params.eq.n1} = ['@(r)(' num2str(params.xi_n2n1,'%25.16d') '* (r <= ' num2str(tdr.FNonLocal.R(params.eq.n2,params.eq.n1),'%25.16d') ') )'];
        tdr.FNonLocal.OmegaStr{params.eq.n2,params.eq.n2} = ['@(r)(' num2str(params.xi_n2n2,'%25.16d') '* (r <= ' num2str(tdr.FNonLocal.R(params.eq.n2,params.eq.n2),'%25.16d') ') )'];
end


% specify boundary condition (only periodic supported)
tdr.FNonLocal.BCs = params.BCs;
tdr.FNonLocal

%
% 2. initialise struct params
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters of the initial conditions    %%%
params.IC.n1      = 1.0; % initial PC density
params.IC.n2      = 1.0; % initial NC density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters of the differential equation %%%
% reaction terms
switch params.fswitch
    case 1
        params.selectReaction = 'none';
    case 2 
        params.selectReaction = 'simple';
end

% cell random motility
params.n1_D = 1.0;
params.n2_D = 1.0;
%
% params struct completed
%


%
% 3. define the spatial domain via its patches and store data in grd struct
%

% We use a 1-patch
grd.nop = 1;             % Number Of Patches
grd.isAxiSymmetric = false;
grd.is1D = true;
grd.y0  = 0.0;
grd.ny  = round(params.domainlength*params.gridCells);
grd.dy  = (params.domainlength) / grd.ny;
switch params.BCs
 case 'pp'
  grd.ngb = [ 1 1 ];
 otherwise
  error('Unsupported BCs');
end
grd
%
% grd struct completed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = setValue(fieldName, defaultValue, inputStruct, subStructName)
  val = defaultValue;
  if isfield(inputStruct, subStructName)
    tmpStruct = getfield(inputStruct, subStructName);
    if isfield(tmpStruct, fieldName)
      val = getfield(tmpStruct, fieldName);
    end
  end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
