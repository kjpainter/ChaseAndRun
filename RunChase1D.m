
function [t,y]=NCPC1D(varargin)
% function [t,y]=CellAdhesionRepulsion_1D(varargin)
%
% This is the main function to run a simulation and can just be
% called with one argument (default values given below):
%    params.domainlength    = 1;
%    params.gridCells       = 400;
%    params.BCs             = 'pp';
%    inputStruct.params = params;
%    CellAdhesionRepulsion_1D(inputStruct);
% 
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : 
%* Date created  : 2008, Sep 30
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.0 initial version (Alf Gerisch)
%*
%*******************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2006-2008 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%******************************************************************************


clear functions
format compact;

oldpath=path();          % save current Matlab path
% include mTDR base directory in path and setup mTDR system
path('mTDR-1.6pre-AdhPack-1.2', path); 
tdrSetup;
% add rowmap to Matlab path
path('~/Matlab/rowmap-0.93',path); 


% print TDR, rowmap and problem directories used 
% as convenience for the user
which tdrInit            % is in tdr/tdr directory
which rowmap
which ProbGetParams      % is in problem directory

% Select time integration scheme
integrator = 'ode15s';
% Set tolerance for time integration
tol=1e-6;
% Set output function
clear function localOutputfun; % to clear persistent variables
outputfun=@localOutputfun; % This function is defined below.

% init the TDR problem (initialises the data in TDRP and returns 
% the initial value y0 and the vector of output times tspan)

clear function ProbFy0; % to reinitialize persistent variables
[y0, tspan] = tdrInit(varargin{:});

% normalise initial masses
tempsum(1) = sum(y0(1:2:length(y0)));
tempsum(2) = sum(y0(2:2:length(y0)));
y0(1:2:length(y0)) = 1*y0(1:2:length(y0))*0.5*length(y0)/tempsum(1);
y0(2:2:length(y0)) = 1*y0(2:2:length(y0))*0.5*length(y0)/tempsum(2);


timerReset();           % reset all timer
switch integrator
 case 'ode45'
  disp('Running ode45 ...');
  % install output function for odesolver
  clear options;
  options=[];
  options = odeset(options, 'OutputFcn', outputfun);
  options = odeset(options, 'AbsTol', tol, 'RelTol', tol);
  tic
  [t,y] = ode45(@tdrFdgl, tspan, y0, options);
  toc
 case 'ode15s'
  disp('Running ode15s ...');
  % install output function for odesolver
  clear options;
  options=[];
  options = odeset(options, 'OutputFcn', outputfun);
  options = odeset(options, 'AbsTol', tol, 'RelTol', tol);
  tic
  [t,y] = ode45(@tdrFdgl, tspan, y0, options);
  toc
  
 otherwise
  error('unknown integrator');
end
%%%%%%%%%%%%%%%%%%%%%%%%
timerPrint('main::');   % print current values of timer

% reset Matlab path
path(oldpath);
  
return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = localOutputfun(t, y, flag)
global TDRP;           % global struct for TDR Problem data

persistent saved;
saveEvery = 100;
doPlot=false;
doPlot=true;

res = 0;

v1 = TDRP.params.index;
v2 = TDRP.params.xi_n1n2;
v3 = TDRP.params.xi_n2n1;

if (v1<10)
datafile =['Simulation00' num2str(v1) '.mat'];    
end
if (v1>=10 & v1<100)
datafile =['Simulation0' num2str(v1) '.mat'];    
end
if (v1>=100)
datafile =['Simulation' num2str(v1) '.mat'];    
end

if (strcmp(flag,'init'))
  disp('Received flag ''init'' in output function.');
  saved.params = TDRP.params;
  saved.tdr    = TDRP.tdr;
  saved.grd    = TDRP.grd;
  saved.t      = t(1);
  
  if (size(y,2)~=[1])
    error('wrong size of y')
  end
  saved.y = y';
  save(datafile, 'saved');  
  disp(['Saved to ' datafile ' at time t = ' num2str(t(1))]);
  if (doPlot)
    plotfun(saved.t(end), saved.y(end,:));
    %disp('Hit button to continue...');
  end
  
  return
end

if (strcmp(flag,'done'))
  disp('Received flag ''done'' in output function.');
  %i = saved.params.S_n1n1;
  %datafile =['sr_o' num2str(saved.params.oswitch) '_f' num2str(saved.params.fswitch) '_h' num2str(saved.params.hswitch) '_ics' num2str(saved.params.icswitch) '_alpha' num2str(i) '.mat'];  
  save(datafile, 'saved');

  if (1 == 1)
      eval(['load ' datafile])
      x = cell2mat(saved.grd.yc);
      t = saved.t;
      u = saved.y(:,1:2:end);
      v = saved.y(:,2:2:end);    
      map(1,1:3) = [1 1 1];
      map(2:102,3) = 0.75;
      map(103:203,1) = 1;    
  end

  umax = max(max(u-0.1));
  vmax = max(max(v-0.1));

  figure(2)
  set(gcf,'paperpositionmode','auto')
  set(gcf,'position',[400 200 280 450])
  clf
  subplot('position',[0.1 0.045 0.35 0.935])
  p1 = surf(x,t,zeros(size(u)),'AlphaData',(1+max(-1,log10(u+v))).*min(20,max(u-v,0)).^0.75,...
    'FaceAlpha','flat',...
    'FaceColor','blue',...
    'edgecolor','none');
  hold on
  p2 = surf(x,t,zeros(size(v)),'AlphaData',(1+max(-1,log10(u+v))).*min(20,max(v-u,0)).^0.75,...
    'FaceAlpha','flat',...
    'FaceColor','red',...
    'edgecolor','none');
  axis([0 10 0 ceil(max(t))])
%  axis([0 10 0 100])
  
  view([0 -90]) 
  grid off
  box on
  set(gca,'fontsize',10)
%  times = ([0 0.5 1 2 5 10 15 20])
  times = 0.5*([0 5 10 20 40 60 80 100])  
  for profiles = 1:8
     ii = find(t>=times(profiles),1); 
     subplot('position',[0.55 0.985-0.1175*profiles 0.4 0.1125]) 
     area(x,u(ii,:),'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)
     hold on
     plot(x,u(ii,:),'color',[0 0 1],'linewidth',2)
     area(x,v(ii,:),'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3)
     plot(x,v(ii,:),'color',[1 0 0],'linewidth',2)
     ylimit = max(max(u(ii,:)),max(v(ii,:)));

     if (ylimit<20) yl = 25; end
     if (ylimit>=20 & ylimit< 40) yl = 50; end
     if (ylimit>=40 & ylimit< 80) yl = 100; end
     if (ylimit>=80 & ylimit< 160) yl = 200; end
     if (ylimit>=160 & ylimit< 320) yl = 400; end
     if (ylimit>=320) yl = 625; end

     axis([0 max(ceil(x)) 0 1.2*yl])
     text(0.7*max(ceil(x)),0.95*yl,['t =' num2str(t(ii))],'fontsize',10);
     if (profiles<8) set(gca,'xtick',[],'ytick',[ceil(0.8*yl)],'fontsize',10); else set(gca,'ytick',[ceil(0.8*yl)],'fontsize',10); end
  end

  eval(['print -djpeg100 -r300 ' datafile '.jpg'])
  return
end

disp(['Reached time t = ' num2str(t)]);

if (size(y,2)~=[1])
  error('wrong size of y')
end
if (size(t) ~= [1 1])
  error('wrong size of t');
end

itmp=length(saved.t)+1;
saved.t(itmp) = t;
saved.y(itmp,:) = y';

if (~mod(length(saved.t)-1, saveEvery))
  save(datafile, 'saved');
  disp(['Saved to ' datafile ' at time t = ' num2str(t(end))]);
end

if (doPlot)
    plotfun(saved.t(end), saved.y(end,:));
end

return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = plotfun(t, y)

global TDRP

%disp([t min(y) max(y), sum(y)]);
figure(1)
clf;
%set(gcf, 'Position', [0 0 600 800]);
%set(gcf,'renderer','zbuffer');
set(gcf,'DoubleBuffer','on');

%grid.x1 = TDRP.grd.xc{1};
grid.x2 = TDRP.grd.yc{1};

n1 = y(1:2:end);
n2 = y(2:2:end);
% solution plot
ntot = n1+n2;
[nmax xmax] = max(ntot)

plot(grid.x2,n1,'b-');
hold on
plot(grid.x2,n2,'r-');

title(['t=' num2str(t)])
axis([0 ceil(max(grid.x2)) 0 1.2*nmax])
%axis([grid.x2(xmax)-7 grid.x2(xmax)+3 0 1.2*nmax])
drawnow
pause(0.0001)

res = 0;


return;
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

