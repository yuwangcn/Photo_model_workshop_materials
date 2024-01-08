% This function is used to store all the required environmental variables, such as light, CO2, O2, humidity as such. 
%   Copyright ? 2007


function fini = cdn (t)

global RUBISCOTOTAL;        % Total concentration of Rubisco.
global RUBISCOMETHOD;       % The method for calculation of Rubisco catalyzed reaction
RUBISCOTOTAL = 3;           
RUBISCOMETHOD = 2;          % 1: Use enzyme concentration for calculation
                            % 2: Use the michaelis menton and enzyme concentration together for calculation
global VolRatioStCyto;
VolRatioStCyto = 1;

global GP;
global gp2condition_RuCon;

if GP ==1
    RUBISCOTOTAL = gp2condition_RuCon;
end


% First get the generic environmental conditions
                            
%global CO2_cond;
global O2_cond;
global GLight;
% global V16;


%CO2Temp = 550*0.7;%375*0.7;%270;  % default is 270 ppm, which corresponds to an atmospheric CO2 concentration of 360 ppm. 
O2Temp = 0.21;  % default is 0.21, i.e. 21%. 

light = 1000; 
GLight = light;


fini = 1;
