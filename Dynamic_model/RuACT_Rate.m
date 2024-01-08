%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%   Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%   CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%   China Institute of Genomic Biology and Department of Plant Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%   University of Illinois at Urbana Champaign
%   Global Change and Photosynthesis Research Unit, USDA/ARS, 1406 Institute of Genomic Biology, Urbana, IL 61801, USA.
 
%   This file is part of e-photosynthesis.
 
%    e-photosynthesis is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; 
 
%    e-photosynthesis is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
 
%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function RuACT_Vel = RuACT_Rate(t,RuACT_Con)

global RuACT_RC;										
k1	=	RuACT_RC	(	1	)			;	%	The rate constant of the activation of the Rubisco bound with RuBP. This step is associated with the ARubisco activase content or activity; i.e. 	Lazar (1999), 0.25~1 *10^(9)			
kn1	=	RuACT_RC	(	2	)			;	%	The rate constant of E inactivation by binding of RuBP;	Lazar 1999, with a lifetime of 5 ns at closed reaction center			
km1	=	RuACT_RC	(	3	)			;	%	The michaelis menton constant for RuBP with E.	Reference needed, a guess			
Ke2	=	RuACT_RC	(	4	)			;	%	Data from Mate et al 1996. Unit: micormolar;	Reference needed, a guess			
Ke3	=	RuACT_RC	(	5	)			;	%	Data from Mate et al 1996. Unit: micormolar;				
k6	=	RuACT_RC	(	6	)			;	%	micromolar per meter square per second, transfered to unit				
kc	=	RuACT_RC	(	7	)			;	%	Michaelis menton constant for CO2				
ko	=	RuACT_RC	(	8	)			;	%	Michaelis menton constant for O2				
k7	=	RuACT_RC	(	9	)			;	%	The rate constant for ecm to ecmr			
kr	=	RuACT_RC	(	10	)			;	%	The apparaent michaelis menton constant for RuBP			
							
ER	=	RuACT_Con	(	1	) 		;	%	The concentration of inactive ER		
Eaf	=	RuACT_Con	(	2	)			;	%	The total concentration of  E, EC, AND ECM		
ECMR	=	RuACT_Con	(	3	)		;	%	The concentration of ECMR	
RT	=	RuACT_Con	(	4	)		;	%	The concentration of ECMR	%%WY2018 RuBP->Rt
							
global RuACT_Pool  ;		
ET = RuACT_Pool(1) ;	
Rac = RuACT_Pool(2);	
C = RuACT_Pool (3) ;	
O = RuACT_Pool (4) ;
MT = RuACT_Pool (5);

global activase;
% activase = 80  ; 

global RROEA_EPS_com;

if RROEA_EPS_com == 1
    
    global RROEA_Pool;
    Rac  = RROEA_Pool(10);      

    global RROEA2RuACT_RuAC;
    activase = RROEA2RuACT_RuAC * 14364;       
end

global RuACT_EPS_com;
global PSPR_RA_O2;
global PSPR_RA_CO2;
global FIBF_RA_Mg;
global PS_C_CA
global PS2RA_ATP;
global O2_cond

if RuACT_EPS_com == 0
    ATP = 1.45;%WY2018 1.45
    ADP = 1.5 - ATP;
    RatioDT = ADP/ATP;
else
     C = PSPR_RA_CO2;
     O = O2_cond;

 
    ATP = PS2RA_ATP;
    ADP = PS_C_CA - ATP;
    RatioDT = ADP/ATP*1/10;

end


CA = 1;
CB = Ke3 + Ke2 * Ke3/C + Eaf -RT; %%WY2018 MT->RT
CC = -RT * (Ke3 + Ke2 * Ke3/C);%%WY2018 MT->RT
RuBP = (-CB +( CB ^2 - 4 * CA * CC )^0.5)/(2 * CA);%% WY2018 M->RuBP
%%%%WY2018 EC calculation
EC=-(C*Ke3 - C*Eaf - (C^2*Eaf^2 + 2*C^2*Eaf*Ke3 - 2*C^2*Eaf*MT + C^2*Ke3^2 + 2*C^2*Ke3*MT + C^2*MT^2 + 2*C*Eaf*Ke2*Ke3 + 2*C*Ke2*Ke3^2 + 2*C*Ke2*Ke3*MT + Ke2^2*Ke3^2)^(1/2) + C*MT + Ke2*Ke3)/(2*(C + Ke2));												
E = EC /C * Ke2;												
ECM=EC*MT/(Ke3+EC);%%WY2018%ECM = EC * M/Ke3;->ECM=EC*MT/Ke3/(1+EC/Ke3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if activase < 10^(-6)
    RCA = 0;
else
%   WY2018 
RCA=activase;
end

if RatioDT > 1
    RatioDT = 1;
elseif RatioDT < 0.1
    RatioDT = 0;
end

FATP = 1 - RatioDT;


v1	=	k1* ER + RCA * ER * FATP;
if v1<0
    v1	=0;
end
factor_n7 = 1; 
kn7 = 0.5 * factor_n7; 

vn7	=	ECMR * kn7; 
if vn7<0
    vn7	=0;
end

if RuBP>0%%%%WY201902
    vn1	=	kn1 * E * RuBP;  
    v7	=	k7 * ECM * RuBP	;	 
    v6_1	=	ECMR * k6 * C/(C+kc * (1 + O/ko));		 
    v6_2	=	ECMR * k6/3 * O/(O+ko * (1 + C/kc));
    v1	=	k1* ER + RCA * ER * FATP;
    if v1<0
    v1	=0;
    end
    factor_n7 = 1; 
    kn7 = 0.5 * factor_n7; 

    vn7	=	ECMR * kn7; 
    if vn7<0
    vn7	=0;
    end
else
   vn1=0;
   v7=0;
   v6_1=0;
   v6_2=0;
   vn7=0;
   v1=0;
end
Percent = 1-ER/ET;
%Percent = (ECMR+ECM)/ET;


global RuACT_OLD_TIME;
global RuACT_TIME_N;
global RuACT_VEL;
global RuACT_CON;

if (RuACT_TIME_N ==0)
    RuACT_TIME_N = 1;
end

if (t > RuACT_OLD_TIME)
    RuACT_TIME_N = RuACT_TIME_N + 1;
    RuACT_OLD_TIME = t;
end

RuACT_VEL	(	RuACT_TIME_N	,	1	)	=	t;
RuACT_VEL	(	RuACT_TIME_N	,   2	)	=	v1	;	 							
RuACT_VEL	(	RuACT_TIME_N	,   3	)	=	vn1	;	 							
RuACT_VEL	(	RuACT_TIME_N	,   4	)	=	v7	; 							
RuACT_VEL	(	RuACT_TIME_N	,   5	)	=	vn7	;	 							
RuACT_VEL	(	RuACT_TIME_N	,   6	)	=	v6_1	;	 					
RuACT_VEL	(	RuACT_TIME_N	,   7	)	=	v6_2	;	 


RuACT_CON(RuACT_TIME_N,1) = t;
RuACT_CON(RuACT_TIME_N,2) = E;
RuACT_CON(RuACT_TIME_N,3) = EC;
RuACT_CON(RuACT_TIME_N,4) = ECM;
RuACT_CON(RuACT_TIME_N,5) = Percent;

RuACT_Vel	(	1	)	=	v1	;	 								
RuACT_Vel	(	2	)	=	vn1	;	 						
RuACT_Vel	(	3	)	=	v7	;	 						
RuACT_Vel	(	4	)	=	vn7	; 								
RuACT_Vel	(	5	)	=	v6_1;	 								
RuACT_Vel	(	6	)	=	v6_2;	 	

global RuACT2RA_v61;
global RuACT2RA_v62;
global RuACT2RA_v1;
global RuACT2RA_vn1;
global RuACT2RA_v7;
global RuACT2RA_vn7;

RuACT2RA_v61 = v6_1;
RuACT2RA_v62 = v6_2;
RuACT2RA_v1 = v1;
RuACT2RA_vn1 = vn1;
RuACT2RA_vn7 = vn7;
RuACT2RA_v7 = v7;




global RuACT2PS_Percent; 
RuACT2PS_Percent = Percent; 