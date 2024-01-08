
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Begin = 1;
%fin = SYSInitial(Begin);
global activase; 
activase = 0.006; 
global options1;
global tglobal;
time = 6000;%tglobal;
global GLight;
GLight=1000;
global RuACT_EPS_com;
RuACT_EPS_com = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the initial condition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RuACT_Con = RuACT_Ini(Begin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Getting the Velocity from RuACT_Rate.m %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global RuACT_OLD_TIME;
global RuACT_TIME_N;
global RuACT_VEL;
global RuACT_CON;

RuACT_OLD_TIME = 0;
RuACT_TIME_N = 1;

RuACT_VEL = zeros(1,3);    % Clean memory
RuACT_CON = zeros(3,1);    % Clean memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Tt,d] = ode15s(@RuACT_MB,[0,time],RuACT_Con,options1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the output                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
    for m = 1:4
        subplot(1,4,m);p = plot(Tt,d(:,m),'.');ylabel('mM');xlabel(' second');
    end
    