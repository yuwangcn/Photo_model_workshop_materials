%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright   2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


% CM_mb.m  This model includes the mass balance equations for the full model of the carbon metabolism inlcuding SUCS and PSPR
function CM_DYDT = CM_mb(t, CM_Con, PS_PR_Param, SUCS_Param)

 
global WeatherTemperature;


SUCS_Param(1) = light; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the rates    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CM_Vel = CM_Rate(t, CM_Con, SUCS_Param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the rate of different reactions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1	=	CM_Vel(1);
v2	=	CM_Vel(2);
v3	=	CM_Vel(3);
NONE=	CM_Vel(4);
v5	=	CM_Vel(5);
v6	=	CM_Vel(6);
v7	=	CM_Vel(7);
v8	=	CM_Vel(8);
v9	=	CM_Vel(9)	;
v10	=	CM_Vel(10);
v13	=	CM_Vel(11);
v16	=	CM_Vel(12);
v23	=	CM_Vel(13);
v31	=	CM_Vel(14);
v32	=	CM_Vel(15);
v33	=	CM_Vel(16);
v111 = CM_Vel(17)  ;
v112 = CM_Vel(18)  ;
v113 = CM_Vel(19)  ;
v121 = CM_Vel(20)  ;
v122 = CM_Vel(21)  ;
v123 = CM_Vel(22)  ;
v124 = CM_Vel(23)  ;
v131 = CM_Vel(24)  ;
v1in = CM_Vel(25) ;
v2out= CM_Vel(26) ;

v51	=	CM_Vel(27);  %	DHAP+GAP --FBP
v52	=	CM_Vel(28);  %	FBP --F6P + Pi
v55	=	CM_Vel(29);  %	G1P+UTP --OPOP+UDPG 
v56	=	CM_Vel(30);  %	UDPG+F6P--SUCP + UDP
v57	=	CM_Vel(31);  %	SUCP--Pi + SUC
v58	=	CM_Vel(32);  %	F26BP--F6P + Pi
v59	=	CM_Vel(33);  %	F6P + ATP --ADP + F26BP
v60	=	CM_Vel(34);  %	ATP+UDP --UTP + ADP
v62	=	CM_Vel(36);  %	SUC SINK 
vdhap_in=	CM_Vel(37);%	DHAP IN
vgap_in	=	CM_Vel(38);%	GAP Export from chloroplast
vpga_in	=	CM_Vel(39);%	PGA export from chloroplast
vpga_use=	CM_Vel(40);%	PGA utilisation in chloroplast
vatpf   =   CM_Vel(41);    %	ATP synthesis rate

NetAssimilation=   CM_Vel(42);
vCO2b=   CM_Vel(43);
vCO2s=   CM_Vel(44);
vH2Ob=   CM_Vel(45);
vH2Os=   CM_Vel(46);
EnergyBalanceResidual=   CM_Vel(47);
vCO2total=   CM_Vel(48);
vH2Ototal=   CM_Vel(49);
vgs=   CM_Vel(50);
vinf=   CM_Vel(51);

VolMC=0.01;
VolMchl=0.02;
Volper=0.00045;
Volmito=0.00045;

Cpwater=4.184;%Jg-1c-1
Vandp=140;%gm-2 leaf thickness 200um leaf density 0.7*10^3kg m-3
Vol_airspace=0.04;%L
Molar_Volume=22.4/273*(WeatherTemperature+273);
Delta_Ci=(vCO2s-NetAssimilation)/(Vol_airspace/Molar_Volume);
Delta_Cb=0;%vCO2b-vCO2s;
Delta_Eb=0;%vH2Os-vH2Ob;
Delta_Gs=vgs;
Delta_Tleaf=EnergyBalanceResidual/Cpwater/Vandp;
Delta_H2Oou=vH2Ob/10^6;%umol->mol
Delta_CO2in=vCO2b/10^6;
Delta_MC_CO2=(vinf-v1)/VolMchl;
LeafMB=zeros(8,1);
LeafMB(1)=Delta_Ci;
LeafMB(2)=Delta_Cb;
LeafMB(3)=Delta_Eb;
LeafMB(4)=Delta_Gs;
LeafMB(5)=Delta_Tleaf;
LeafMB(6)=Delta_H2Oou;
LeafMB(7)=Delta_CO2in;
LeafMB(8)=Delta_MC_CO2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the mass balance equation%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = zeros(23,1);
% RuBP
tmp(1) = (v13-v1-v111)/VolMchl;
% PGA;
tmp(2) = (2*v1-v2-v32 + v113 + v111)/VolMchl;
% DPGA
tmp(3) = (v2-v3)/VolMchl;
% T3P
tmp(4) = (v3 - 2* v5 -v7 -v8 - v10 - v31-v33)/VolMchl;
% FBP;
tmp(5) = (v5-v6)/VolMchl;
% E4P;
tmp(6) = (v7-v8)/VolMchl;
%S7P;
tmp(7) = (v9-v10)/VolMchl;
%SBP;
tmp(8) = (v8 - v9)/VolMchl;
%  ATP;
tmp(9) = (v16 - v2 - v23 - v13- v113)/VolMchl;    
% NADPH;
tmp(10) = 0;
% CO2;
tmp(11) = 0;
% O2;
tmp(12) = 0;
% HexP;
tmp(13) = (v6 - v7 - v23)/VolMchl;
%  PenP;
tmp(14) = (v7 + v10 * 2 - v13)/VolMchl;
% Gcea 
tmp(15) = (v1in  - v113)/VolMchl;
% Gca 
tmp(16) =  (v112 - v2out)/VolMchl;
% PgCa 
tmp(17) = (v111 - v112)/VolMchl;
% Gcac 
tmp(18) = (v2out - v121)/Volper;
% Goac
tmp(19) = (v121 - v122- v124)/Volper;
% Serc 
tmp(20) = (v131 - v122)/Volper;
% Glyc 
tmp(21) = (v122 + v124 - 2*v131)/Volper;
% Hprc 
tmp(22) = (v122 - v123)/Volper;
% Gceac 
tmp(23) = (v123 - v1in)/Volper;
PSPR_DYDT = tmp;

SUCS_mb	    =	    zeros(12,1)	;				
SUCS_mb(1)	=	(v31 + v33 -  2* v51)/VolMC		;   %	T3Pc
SUCS_mb(2)	=	(v51-v52)/VolMC		;   %	FBPc
SUCS_mb(3)	=	(v52 -v55 -v59 + v58	-v56)/VolMC	;   %	HexPc
SUCS_mb(4)	=	(v59 - v58)/VolMC		;   %	F26BPc
SUCS_mb(5)	=	(0)/VolMC;% vatpf - v59 - v60		        ;   %	ATPc
SUCS_mb(6)	=	(0)/VolMC		        ;   %	ADPc
SUCS_mb(8)	=	(v55 - v56)/VolMC		;   %	UDPGc
SUCS_mb(9)	=	(0)/VolMC;% v60 - v55		;   %	UTPc
SUCS_mb(10)	=	(v56 - v57)/VolMC		;   %	SUCP
SUCS_mb(11)	=	(v57 - v62)/VolMC		;   %	SUC
SUCS_mb(12)	=	(v32 - vpga_use)/VolMC ;   %	pgaC
SUCS_DYDT=SUCS_mb;

CM_DYDT = zeros(43,1);
for m = 1:8
    CM_DYDT(m) = LeafMB(m);
end
for m = 1:23
    CM_DYDT(m+8) = PSPR_DYDT(m);
end
for m = 1:12
    CM_DYDT(m+31) = SUCS_DYDT(m);
end


