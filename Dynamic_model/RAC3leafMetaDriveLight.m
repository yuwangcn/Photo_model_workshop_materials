
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AandW=RAC3leafMetaDriveLight(inputfile,CultivarNo,Pst,PRca)
%clear all;
 global lightM;
 global Vmaxdata;
STdata=importdata('DynamicParameter_Cassava.txt');
CasSTdata=STdata.data;
Vdata=importdata('EnzymeVmax.txt');
Vmaxdata=Vdata.data;
Envdata=importdata('Environment.txt');
EnvData=Envdata.data;
Lightd=importdata(inputfile);
Lightdata=Lightd.data;


lightM(:,1)=Lightdata(1:800,1);
lightM(:,2)=Lightdata(1:800,2);

simNo=CultivarNo;
dVcmax=CasSTdata(simNo,1);
dJmax=CasSTdata(simNo,2);
dKd=CasSTdata(simNo,3);
dKi=CasSTdata(simNo,4);
dSlop=CasSTdata(simNo,5);
dInter=CasSTdata(simNo,6);

Begin = 1;
global RAInteg;
RAInteg=PRca;% 1=consider Rac; 0=Rubisco always actived 
global activase; 
activase = 0.002; 
global RuACT_EPS_com;
RuACT_EPS_com = 1;

%%%%%%%%%%%%%RuACT_EPS_com
global BallBerryInterceptC3;
global BallBerrySlopeC3
BallBerryInterceptC3=1.6*dInter;%0.008;%WY201804%Ball 1988
BallBerrySlopeC3=1.6*dSlop*100;%WY201804  %9.29;%Ball 1988%10.5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ki_Gs;
global kd_Gs;
global GsResponse;
ki_Gs=1/dKi;%0.9*60;%6.9*60*exp(1);%2*60;0.9 %WY201804 cassava
kd_Gs=1/dKd;%6.5*60*exp(1);%1*60;4.1
GsResponse=Pst; %%if GsResponse=0 Ball berry model no time dependent Gs response ; %%if GsResponse=1 time dependent Gs response, using ki_Gs and kd_Gs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Para_mata;
global PhotosynthesisType;
Para_mata=1;%%if Para_mata=1, C3 Metabolic model and Gs model integrated  if Para_mata=0 Farquhar mdoel and gs model
PhotosynthesisType=1;% 1:C3 model 2:C4 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leaf model parameters

global Air_CO2;
global Air_O2;
global WeatherWind;
global Radiation_PAR;
global Vcmax25;
global Jmax25;
global WeatherTemperature;
global WeatherRH;
global Radiation_NIR;
global Radiation_LW;
global PhiLeaf

Convert=1E6/(2.35E5); 
WeatherTemperature=EnvData(1);
Air_CO2=EnvData(2);
Air_O2=EnvData(3);
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=2000/Convert*0.85*0.85;%10*i;
Radiation_NIR=0;
Radiation_LW=0;
Vcmax25=dVcmax;
Jmax25=dJmax;
PhiLeaf=-0;%Mpa Min=-2.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global O2_cond

PCO2=Air_CO2;
O2_cond=Air_O2/1000*1.26;
Plight=Radiation_PAR*Convert;
PTemp=WeatherTemperature;

global CI
CI=Air_CO2*0.7/(3 * 10^4);
global CO2_Env;
global CO2_cond;
global LI;
global Jmax;
global alfa;
global fc;
global Theta;
global beta;
CO2_Env=PCO2;
CO2_cond = CO2_Env/(3 * 10^4);%CO2_Env*0.7/(3 * 10^4);% for leaf model input is CI
LI=Plight/1000;
Jmax=dJmax/1000;
alfa=0.85;
fc=0.15;
Theta=0.7;
beta=0.7519;

global Tp;
Tp=PTemp;
Begin = 1;
global tglobal;     % The total running time
tglobal =1*3600;%;12*3600
% global options1 
options1 = odeset('RelTol',1e-1,'AbsTol',1e-2);
time = tglobal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variables used for obtaining flux and concentration data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     

global PS_OLD_TIME;
global PS_TIME_N;
global PS_VEL;
PS_OLD_TIME = 0;
PS_TIME_N= 0;
PS_VEL = zeros(1,1);

global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
PR_OLD_TIME = 0;
PR_TIME_N = 1;
PR_VEL = zeros(1,1);

global SUCS_OLD_TIME;
global SUCS_TIME_N;
global SUCS_VEL;
global SUCS_CON;

SUCS_OLD_TIME = 0;
SUCS_TIME_N = 1;
SUCS_VEL = zeros(1,3);    % Clean memory
SUCS_CON = zeros(3,1);    % Clean memory

global TIME_N;
global OLD_TIME;
global Gs_VEL;
TIME_N=0;
OLD_TIME=0;
Gs_VEL=zeros(1,10);
global RuACT_OLD_TIME;
global RuACT_TIME_N;
global RuACT_VEL;
global RuACT_CON;

RuACT_OLD_TIME = 0;
RuACT_TIME_N = 1;

RuACT_VEL = zeros(1,3);    % Clean memory
RuACT_CON = zeros(3,1);    % Clean memory

%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%

global GP; 
GP = 0;
CMs = RAC3leafMetaIni(Begin);

global PR_PS_com;             % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global ATPActive;
ATPActive = 0;

global PSPR_SUCS_com;        % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

global RedoxReg_RA_com;
RedoxReg_RA_com = 0;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%

CM_Param = 0;

SUCS_Param = zeros(2,1);
SUCS_Param(1) = 1;
SUCS_Param(2) = 1;

PS_PR_Param = 0;


[Tt,d] = ode15s(@RAC3leafMetaMB,[0,time],CMs,PS_PR_Param, SUCS_Param);
   
global d_plot;
d_plot=d;

global Tt_plot;
Tt_plot = Tt;
global Result;
Result =[Tt,d]; 
Eb=d(:,3);
gm=0.7;
Sc=3*10^4;
%vinf=gm*Sc*10^(-3)*(Ci/(3 * 10^4)-MC_CO2);
%vH2Ob=Gbw*(Eb-Ea)/(Pressure / 1000.0)*10^6.0;
global Gs_VEL;

Gsw=1.6*Result(:,5);
ESaturation = 0.611 * exp(17.502 * Result(:,6)./ (Result(:,6) + 240.97));
Eb=Result(:,4);
Pressure=101325.0;
figure;
subplot(2,3,1); plot(lightM(:,1)*60,lightM(:,2),'k');title('PAR');xlim([0,60]);ylim([0,1800]);
subplot(2,3,2); plot(Result(:,1)/60,Result(:,2),'k');title('Ci');ylim([0,Air_CO2*1.2]);
subplot(2,3,3); plot(Result(:,1)/60,Result(:,6),'k');title('Tleaf');
subplot(2,3,4); plot(Result(:,1)/60,Result(:,5),'k');title('Gs');
subplot(2,3,5); plot(Result(:,1)/60,gm*Sc*(Result(:,2)/(3 * 10^4)-Result(:,9)),'k');title('A');ylim([0,40])
subplot(2,3,6); plot(Result(:,1)/60,Gsw.*(ESaturation-Eb)./(Pressure / 1000.0)*10^6.0,'k');title('Transpiration');


figure;
subplot(3,4,1); plot(Tt/60,d(:,8+1),'k');title('RuBP');
subplot(3,4,2); plot(Tt/60,d(:,8+2),'k');title('PGA');
subplot(3,4,3); plot(Tt/60,d(:,8+3),'k');title('DPGA');
subplot(3,4,4); plot(Tt/60,d(:,8+4),'k');title('T3P');
subplot(3,4,5); plot(Tt/60,d(:,8+5),'k');title('FBP');
subplot(3,4,6); plot(Tt/60,d(:,8+6),'k');title('E4P');
subplot(3,4,7); plot(Tt/60,d(:,8+7),'k');title('S7P');
subplot(3,4,8); plot(Tt/60,d(:,8+8),'k');title('SBP');
subplot(3,4,9); plot(Tt/60,d(:,8+13),'k');title('HexP');
subplot(3,4,10); plot(Tt/60,d(:,8+14),'k');title('PenP');

% [row2,col2]=size(Result);
% AandW(simNo,1)=Result(row2,7)';
% AandW(simNo,2)=Result(row2,8)';

%AMeM=CM_OUT(Tt,d(:,9:43));
% 
% PSPR_SUCS_com = 0;
% IModelCom;
end
%end


