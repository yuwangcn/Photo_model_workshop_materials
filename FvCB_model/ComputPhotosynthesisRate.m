function PhotosynthesisRate=ComputPhotosynthesisRate(Vcmax25,Jmax25,Rd25,LeafTemperature,PAR,Ci,Air_O2)

R=8.314472E-3;%Gas constant KJ mole^{-1} K^{-1}
PhotosynthesisTheta=0.76;
%Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute C3 photosynthesis
Rate_TPu = 23.0; %u moles/m2 leaf area/s
%Temporary variables
Rd = Rd25 * exp(18.72 - 46.39 / (R * (LeafTemperature + 273.15)));
LeafTemperatureKelvin = LeafTemperature + 273.15; %Leaf temperature in K
GammaStar = exp(19.02 - 37.83 / (R * LeafTemperatureKelvin));
Ko = exp(12.38 - 23.72 / (R * LeafTemperatureKelvin));
Kc = exp(35.98 - 80.99 / (R * LeafTemperatureKelvin));	
Vcmax = Vcmax25 * exp(26.36 - 65.33 / (R * LeafTemperatureKelvin));
PhiPS2 = 0.385 + 0.02166 * LeafTemperature - 3.37 * LeafTemperature^2.0 / 10000.0;% Match PS_FIT
%I = Convert * Radiation_PAR * PhiPS2 * 0.5;
I =  PAR * PhiPS2 * 0.5;
ThetaPS2 = PhotosynthesisTheta + 0.01713 * LeafTemperature - 3.75 * LeafTemperature^2.0 / 10000.0; % Match PS_FIT
Jmax = Jmax25 * exp(17.71 - 43.9 / (R * LeafTemperatureKelvin));
J = (I + Jmax - sqrt((I + Jmax)^2.0 - 4.0 * ThetaPS2 * I * Jmax)) / (2.0 * ThetaPS2);
LeafAc = (1.0 - GammaStar / Ci) * (Vcmax * Ci) /(Ci + Kc * (1.0 + Air_O2 / Ko));%Rubisco limited photosynthesis
LeafAj = (1.0 - GammaStar / Ci) * (J * Ci) /(4.5 * Ci + 10.5 * GammaStar); %Light limited photosynthesis
if (LeafAj < 0.0)
    LeafAj = 0.0;
end
LeafAp = (3.0 * Rate_TPu) / (1.0 - GammaStar /Ci); %TPU limited photosynthesis
if (LeafAp < 0.0)
    LeafAp=0.0;
end
GrossAssimilation=min(min(LeafAc,LeafAj),LeafAp);%Minimum of three limitations
NetAssimilation=GrossAssimilation-Rd;
  
    if (isinf(GrossAssimilation) || isnan(GrossAssimilation))
%         fprintf(LogOutputFile, "Error in ComputeC3Photosynthesis for Leaf ID: %f: GrossAssimilation=%f\n",
%                 Photosynthesis->LeafID, LeafMassFlux->GrossAssimilation);
        GrossAssimilation = 0.0;
        NetAssimilation =-Rd;

    end

PhotosynthesisRate(1)=NetAssimilation;
% PhotosynthesisRate(2)=GrossAssimilation;
% PhotosynthesisRate(3)=Rd;
% PhotosynthesisRate(4)=GammaStar;
end