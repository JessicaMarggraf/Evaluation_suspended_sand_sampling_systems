% 
%	test of the validity of the CL formula for suspension
%   impact on the transit time for the Delft bottle
% 

clear
close all

graphic = 1;
% 1: Error as a function of Z/h
% 2: same with data from Jessica Marggraf
% 3: comparison calculation Rouse number with data from Jessica Marggraf
% 4 :boxplot error as a function of z/h

%========================================================================== 
%
%  Variables
%

zh = linspace(0,1,100);
% relative water depth

PRM = [0.25 0.5 1 2];
DPRM = [(PRM(2)+PRM(1))/2 (PRM(3)+PRM(2))/2 (PRM(4)+PRM(3))/2];
% Rouse number

DT = 0.1;
% water depth

ErUp1 = 1/6/PRM(1)*DT*1./(1-zh).*(1-exp(-6*PRM(1)*(1-zh)));
ErUp2 = 1/6/PRM(2)*DT*1./(1-zh).*(1-exp(-6*PRM(2)*(1-zh)));
ErUp3 = 1/6/PRM(3)*DT*1./(1-zh).*(1-exp(-6*PRM(3)*(1-zh)));
ErUp4 = 1/6/PRM(4)*DT*1./(1-zh).*(1-exp(-6*PRM(4)*(1-zh)));

%========================================================================== 
%
%  Experimental data from Marggraf (2024)
%

load ./BD_Marggraf.dat
%load ./BD_Marggraf.txt

yM = BD_Marggraf(:,1);
% abscissa

zM = BD_Marggraf(:,2);
% sampling height

hM = BD_Marggraf(:,3);
% water depth

zhM = zM./hM;
selectzh1 = find(zhM>0.8);
selectzh2 = find((zhM>0.6).*(zhM<=0.8));
selectzh3 = find((zhM>0.4).*(zhM<=0.6));
selectzh4 = find((zhM>0.2).*(zhM<=0.4));
selectzh5 = find(zhM<=0.2);
% relative sampling depth

uM = BD_Marggraf(:,4);
% local horizontal velocity

UM = BD_Marggraf(:,5);
% depth averaged velocity

TsM = BD_Marggraf(:,6);
% sampling duration

TmM = BD_Marggraf(:,7);
DTM = TmM./TsM;
% sampling move duration

d50M = BD_Marggraf(:,8)*1e-6;
% median grain size of the suspension

RoM = BD_Marggraf(:,9);
% Rouse number (calcalated by Marggraf based on the fit on vertical concentrations)


kst = 0.05;
% roughness height

%--------------------------------------------------------------------------
%  constants
%

nu = 1e-06;
% kinematic viscosity of water

g = 9.81;
% acceleration of gravity

kappa = 0.41;
% Von Karman constant

rhos = 2650;
rho = 1000;
s = rhos/rho;
% sediment and water density

%--------------------------------------------------------------------------
%  sediment and hydraulic parameters 
%

D = (g*(s-1)/nu^2).^(1/3).*d50M;
% sedientologic diameter

thetacr = 0.24./D + 0.055*(1-exp(-0.02*D));
% critical Shields parameter for the inception of movement (Soulsby, Whitehouse)

ACam = 24.6; BCam = 0.96; nCam = 1.53;
RCam = ( sqrt( (ACam/BCam).^(2/nCam)/4+(4/3/BCam*D.^3).^(1/nCam) ) ...
         - (ACam/BCam).^(1/nCam)/2 ).^nCam;
Ws = RCam*nu./d50M;
% settling velocity (Camenen, 2007)

%--------------------------------------------------------------------------
%  Bed shear stress
%

z0 = kst/30;
% roughness length

fc = 2*(kappa./(1+log(z0./hM))).^2;
% dimensionless friction coefficient

theta = 0.5*fc.*UM.^2./((s-1).*d50M*g);
% current related skin Shields parameter

ust = sqrt(theta.*(s-1).*d50M*g);
% shear velocity

Rou = Ws./(kappa.*ust);
% Rouse number

selectR1 = find(Rou<=DPRM(1));
selectR2 = find((Rou<=DPRM(2)).*(Rou>DPRM(1)));
selectR3 = find((Rou<=DPRM(3)).*(Rou>DPRM(2)));
selectR4 = find(Rou>DPRM(3));

ErUpDTM = 1/6./Rou*1./(1-zhM).*(1-exp(-6*Rou.*(1-zhM)));
ErUpM = ErUpDTM.*DTM;
% sampling error

%==========================================================================
% 
%  graphics 
% 
 
%------------------------------------------------------------------------------
if graphic == 1
    
subplot(1,1,1) 
plot(zh,ErUp1./DT,'k-',zh,ErUp2./DT,'r--',zh,ErUp3./DT,'b-.',zh,ErUp4./DT,'g:','LineWidth',2)
xlabel('$z / h$ [-]','FontSize',16,'Interpreter','Latex'); 
ylabel('$E_{rM} / \Delta{T}$ [-]','FontSize',16,'Interpreter','Latex');
axis([0 1 0 1])
%set(gca,'XTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])
set(gca,'Fontsize',14)
legend('$P_R=0.25$','$P_R=0.5$','$P_R=1.0$','$P_R=2.0$','Location','NorthWest','Interpreter','Latex')
grid on

%------------------------------------------------------------------------------
elseif graphic == 2
    
subplot(1,1,1) 
plot(zh,ErUp1./DT,'k-',zh,ErUp2./DT,'r--',zh,ErUp3./DT,'b-.',zh,ErUp4./DT,'g:','LineWidth',2), hold on
%plot(zhM,ErUpDTM,'ko'), hold on
plot(zhM(selectR1),ErUpDTM(selectR1),'ko'), hold on
plot(zhM(selectR2),ErUpDTM(selectR2),'ro'), hold on
plot(zhM(selectR3),ErUpDTM(selectR3),'bo'), hold on
plot(zhM(selectR4),ErUpDTM(selectR4),'go'), hold on
xlabel('$z / h$ [-]','FontSize',16,'Interpreter','Latex'); 
ylabel('$E_{rM} / \Delta{T}$ [-]','FontSize',16,'Interpreter','Latex');
axis([0 1 0 1])
set(gca,'Fontsize',14)
legend('$P_R=0.25$','$P_R=0.5$','$P_R=1.0$','$P_R=2.0$','Location','NorthWest','Interpreter','Latex')
grid on

%------------------------------------------------------------------------------
elseif graphic == 3
    
subplot(1,1,1) 
plot(RoM,Rou,'ko'), hold on
plot([0 2],[0 2],'k-')
xlabel('$P_R$ (Margraff) [-]','FontSize',16,'Interpreter','Latex'); 
ylabel('$P_R$ (Camenen) [-]','FontSize',16,'Interpreter','Latex');
axis([0 2 0 2])
set(gca,'Fontsize',14)
grid on

%------------------------------------------------------------------------------
elseif graphic == 4  

A = ErUpM(selectzh1);    
B = ErUpM(selectzh2);    
C = ErUpM(selectzh3);    
D = ErUpM(selectzh4);    
E = ErUpM(selectzh5);    
group  = [0.9*ones(size(A)); 0.7*ones(size(B)); 0.5*ones(size(C)); ...
          0.3*ones(size(D)); 0.1*ones(size(E))];

subplot(1,1,1) 
boxplot([A;B;C;D;E],group,'Orientation','Horizontal'), hold on
xlabel('$E_{rM}$ [-]','FontSize',16,'Interpreter','Latex'); 
ylabel('$z/h$ [-]','FontSize',16,'Interpreter','Latex');
axis([0 0.6 0 6])
set(gca,'Fontsize',14)
grid on

%------------------------------------------------------------------------------
end    

