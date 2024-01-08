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



function RuACT_Con = RuACT_Ini(begin)
		
k1	=	0.006/20/2;	%	The rate constant of the activation of the Rubisco bound with RuBP. This step is associated with the ARubisco activase content or activity;
kn1	=	1.6 * 10^(-3)	; 
km1	=	20 * 10 ^ ( -6)	;	 
Ke2	=	0.1	;	    	
Ke3	=	1.600	;	 	
k6	=	2.5*1.5  	;	 		
kc	=	0.016	;	%	Michaelis menton constant for CO2, 0.016 mM				
ko	=	0.448	;	%	Michaelis menton constant for O2, .448 mM				
k7	=	k6 * 10;	%	The rate constant for ecm to ecmr		mM-1 s-1   
kr	=	20 * 10^(-3)	;	

global RuACT_RC;
RuACT_RC = zeros(5,1);

% The rate constant used in the model											
RuACT_RC	(	1	)	=	k1	;	%	The rate constant of the activation of the Rubisco bound with RuBP. This step is associated with the ARubisco activase content or activity;
RuACT_RC	(	2	)	=	kn1	;	%	The rate constant of E inactivation by binding of RuBP;		
RuACT_RC	(	3	)	=	km1	;	%	The michaelis menton constant for RuBP with E.	
RuACT_RC	(	4	)	=	Ke2	;	%	Data from Mate et al 1996. Unit: micormolar;	
RuACT_RC	(	5	)	=	Ke3	;	%	Data from Mate et al 1996. Unit: micormolar;				
RuACT_RC	(	6	)	=	k6	;	%	micromolar per meter square per second, transfered to unit				
RuACT_RC	(	7	)	=	kc	;	%	Michaelis menton constant for CO2				
RuACT_RC	(	8	)	=	ko	;	%	Michaelis menton constant for O2				
RuACT_RC	(	9	)	=	k7	;	%	The rate constant for ecm to ecmr	
RuACT_RC	(	10	)	=	kr	;	%	The apparaent michaelis menton constant for RuBP	

factor = 1; 
ER	=	0.2 * 4 *factor	;	% 	The concentration of inactive ER			
Eaf	=	0.05 * 4 *factor	;	%	The total concentration of E, EC, AND ECM			
ECMR	=	0.05	* 4 *factor;	% 	The concentration of ECMR			
RuBP = 2;

% global RuACT_EPS_com;
% global PS2RA_RuBP_ini;
% 
% if RuACT_EPS_com ==1
%     RuBP = PS2RA_RuBP_ini;
% end
% 								
% Assign value to a variable that is transferred to the program				
RuACT_Con = zeros(3,1);

RuACT_Con	(	1	)	=	ER	;	%	The concentration of inactive ER
RuACT_Con	(	2	)	=	Eaf	;	%	The total concentration of E, EC, AND ECM
RuACT_Con	(	3	)	=	ECMR	;	%	The concentration of ECMR
RuACT_Con	(	4	)	=	RuBP	;	%	The concentration of ECMR

ET	=	0.3	  * 4   ;	    % 	The total concentraiton of E, ER, EC, ECM, ECMR	, mM	
Rac	=	0.0056;	    %	The concentration of the activase, mM	
C = 0.012;              %   mM
O = 0.260;              %   mM
M = 5;                   

global RuACT_Pool;		
RuACT_Pool = zeros(2,1);	
RuACT_Pool (1) = ET;	
RuACT_Pool (2) = Rac;	
RuACT_Pool (3) = C;
RuACT_Pool (4) = O;	
RuACT_Pool (5) = M;	