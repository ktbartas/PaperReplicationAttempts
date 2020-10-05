%10/4/20 
%trying to recreate fig 3a from Tyson JJ paper
clc;close all;

%amino acid concentration
aa=1; %this is a constant but idk what the value is supposed to be

%define rates - given in table 2
k1 = .015*.02/aa; %min^-1
%^^^IF YOU CHANGE THE INITIAL VALUES, change .02 above to whatever number is
%printed repeatedly below after you run 
k2=0;%min^-1
k3=200/.02;%min^-1
%^^^IF YOU CHANGE THE INITIAL VALUES, change .02 above to whatever number is
%printed repeatedly below after you run 
k4=180; %adjustable - will change later
kprime4=0.018;%min^-1
k5=0;%min^-1
k6=1; %adjustable - will change later
k7=0.6;%min^-1
k9=k6+10.;%min^-1
k8=k9+10.;%min^-1

%equations - given in table 1
eq1 =@(C2,CP,pM,M,Y,YP) k6*M - k8*C2 + k9*CP   ;
eq2 =@(C2,CP,pM,M,Y,YP) -k3*CP*Y + k8*C2 - k9*CP  ;
eq3 =@(C2,CP,pM,M,Y,YP) k3*CP*Y - pM*(kprime4 + k4*(M/(C2 + CP + pM + M))^2) + k5*M   ;
eq4 =@(C2,CP,pM,M,Y,YP) pM*(kprime4 + k4*(M/(C2 + CP + pM + M))^2) - k5*M - k6*M   ;
eq5 =@(C2,CP,pM,M,Y,YP) k1*aa - k2*Y - k3*CP*Y   ;
eq6 =@(C2,CP,pM,M,Y,YP) k6*M - k7*YP   ;

[time, ex] = ode45(@(a,b)[eq1(b(1),b(2),b(3),b(4),b(5),b(6));eq2(b(1),b(2),b(3),b(4),b(5),b(6));eq3(b(1),b(2),b(3),b(4),b(5),b(6));eq4(b(1),b(2),b(3),b(4),b(5),b(6));eq5(b(1),b(2),b(3),b(4),b(5),b(6));eq6(b(1),b(2),b(3),b(4),b(5),b(6))], [0,100], [0,0.01,0,0.01,0.01,0] );
%add pM, M, Y and YP together for Ytotal
Ytot=ex(:,3)+ex(:,4)+ex(:,5)+ex(:,6);
Ctot=ex(:,1)+ex(:,2)+ex(:,3)+ex(:,4)
%divide by CT
Ytotvdiv=Ytot./ Ctot;

Mfordiv =ex(:,4);
%divide M by CT
Mdiv=Mfordiv./ Ctot;

figure; hold on;
plot(time,Ytotvdiv,'-r') % red for Ytot
plot(time,Mdiv,'-', 'color', [0.5 0 1]) % purple/blue for M
legend('[Ytot]/[CT]','[M]/[CT]')
ylabel('Molecular concentration (micromolar? molar?)')
xlabel('Time (minutes)')
