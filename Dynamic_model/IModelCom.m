%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright ©  2007

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

% This function set the default components used for the full photosynthesis model. 
function ModelComb = IModelCom;

global RuACT_EPS_com;
RuACT_EPS_com = 0;

global BF_FI_com;     % The combination of BF and FI model 
BF_FI_com = 0;

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 0;

global FIBF_PSPR_com; % 1 means that the overall EPS model is used. 0 means partial model of FIBF is used. 
FIBF_PSPR_com = 0;    

global ATPActive;
ATPActive = 0;

global RedoxReg_RA_com;
RedoxReg_RA_com = 0;

global XanCycle_BF_com;
XanCycle_BF_com = 0;

global RROEA_EPS_com;
RROEA_EPS_com = 0;

global StomCond_TrDynaPS_com;
StomCond_TrDynaPS_com = 0;

global PSPR_SUCS_com;
PSPR_SUCS_com = 0;  

global trDynaPS_SUCS_com;
trDynaPS_SUCS_com = 0;

global EPS_SUCS_com;
EPS_SUCS_com = 0;

ModelComb = 1;