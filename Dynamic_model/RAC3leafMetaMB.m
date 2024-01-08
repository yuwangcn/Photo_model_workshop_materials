
function CM_DYDT = CM_mb(t, CM_Con, PS_PR_Param, SUCS_Param)
VolMC=0.01;
VolMchl=0.02;
Volper=0.00045;
Volmito=0.00045;
global Radiation_PAR;
 
global WeatherTemperature;
Convert=1E6/(2.35E5); 

 global lightM;

 lightTime=(lightM(:,1))*60*60;

 [maxt,positionm]=max(lightTime);
 if t<=maxt
 Radiation_PAR=interp1(lightTime,lightM(:,2),t,'linear')/Convert*0.85*0.85;

 end

fini = cdn (t);
light = Radiation_PAR*Convert/1000;

SUCS_Param(1) = light; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RuACT_Con=CM_Con(44:47);
global PSPR_RA_O2;
global PSPR_RA_CO2;
PSPR_RA_CO2=CM_Con(8);
PSPR_RA_O2=CM_Con(20);
%global PSPR_RA_CA;
global PS2RA_ATP;
PS2RA_ATP=CM_Con(17);
RuACT_Vel = RuACT_Rate(t,RuACT_Con);

RAv1	=	RuACT_Vel	(	1	)	;	%	v1	The rate of ER activation due to Rubisco activase
RAvn1	=	RuACT_Vel	(	2	)	;	%	vn1	The rate of E inactiavtion due to binding of RuBP
RAv7	=	RuACT_Vel	(	3	)	;	%	v7	The rate of formation of ECMR from ECM by binding of RuBP
RAvn7	=	RuACT_Vel	(	4	)	;	%	vn7	The rate of actiavtion of ECMR by Rubisco activase
RAv6_1	=	RuACT_Vel	(	5	)	;	%	v6_1	The rate of RuBP carboxylation
RAv6_2	=	RuACT_Vel	(	6	)	;	%	v6_2	The rate of RuBP oxygenation


RuACT_mb	=	zeros(4,1)	;				
RuACT_mb	(	1	)	=	(RAvn1 - RAv1)/VolMchl	;	%	ER		
RuACT_mb	(	2	)	=	(RAv1- RAvn1 - RAv7 + RAvn7 + RAv6_1 + RAv6_2)/VolMchl	;	%	EAF		
RuACT_mb	(	3	)	=	(RAv7 - RAvn7 - RAv6_1 - RAv6_2)/VolMchl	;	%	ECMR	
RuACT_mb	(	4	)	=	(RAv6_1 + RAv6_2 + RAv1 - RAvn1 + RAvn7 - RAv7)/VolMchl;	%	RuBP	

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
Rd= CM_Vel(52);


Cpwater=4.184;%Jg-1c-1
Vandp=140;%gm-2 leaf thickness 200um leaf density 0.7*10^3kg m-3
%leaf volume leaf thickness 200um~0.2L
%bundary layer
%leaf air space 20%
Vol_airspace=0.04;%L

Molar_Volume=22.4/273*(WeatherTemperature+273);
Delta_Ci=(vCO2s-NetAssimilation)/(Vol_airspace/Molar_Volume);
Delta_Cb=0;%vCO2b-vCO2s;
Delta_Eb=0;%vH2Os-vH2Ob;
Delta_Gs=vgs;
Delta_Tleaf=EnergyBalanceResidual/Cpwater/Vandp;
Delta_H2Oou=vH2Ob/10^6;%umol->mol
Delta_CO2in=vCO2b/10^6;
Delta_MC_CO2=(vinf+Rd-v1)/VolMchl;

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
tmp(1) = (v13-v1-v111)/VolMchl;%WY201902
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
%tmp(11) = (vinf+Rd-v1)/VolMchl;
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



%%%%%%%%%%%%%%%%%%%%%
CM_DYDT = zeros(47,1);
for m = 1:8
    CM_DYDT(m) = LeafMB(m);
end
for m = 1:23
    CM_DYDT(m+8) = PSPR_DYDT(m);
end
for m = 1:12
    CM_DYDT(m+31) = SUCS_DYDT(m);
end
CM_DYDT(44:47,1)=RuACT_mb;

