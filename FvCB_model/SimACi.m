Vcmax25 =137;
Jmax25 =175;
Rd25 = 1;
LeafTemperature=25;
PAR=1800;%light intensity
Air_O2=210.0;%O2 concertation

%%ACI curve%%
SimCi=50:50:1000;
[Row, col]= size(SimCi);
for i= 1:col
    SimA(i,1)=SimCi(i);
    SimA(i,2)=ComputPhotosynthesisRate(Vcmax25,Jmax25,Rd25,LeafTemperature,PAR,SimCi(i),Air_O2);
end 

figure;
plot(SimA(:,1), SimA(:,2));
xlabel('Ci (\mumol mol^-^1)');
ylabel('A (\mumol m^-^2 s^-^1)');

