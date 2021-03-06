% Cable Equation
% Kiersten Scott
% BENG 260
% close all
% clear all


% Solve equation linear cable equation

L = 0.04; % [cm]
tau=.00064e-3;
x = linspace(0,L,1000);  %[cm]
t = linspace(0,20*tau,10000); %[s]

% V(x,t)= transmembrane potential [mV]
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
m = 0;
do=-7.7500e-04;
soldl1 = pdepe(m,@cable,@cableic,@cablebc,x,t, [],do);
do=-7.7500e-03;
soldl2 = pdepe(m,@cable,@cableic,@cablebc,x,t, [],do);
% --------- PLOT SOLUTION -------------------------------------
figure
surf(x,t,soldl1, 'EdgeColor','none')
hcb=colorbar()
xlabel('Time (ms)')
ylabel('x(cm)')
zlabel('V (mv)')
colorTitleHandle = get(hcb,'Title');
titleString = 'V (V)';
set(colorTitleHandle ,'String',titleString);
set(gca,'Fontsize',15)

figure
transparency = fliplr(linspace(0.2, 0.9, length(x)-10));
i=1
for k = 5:200:length(x)-10
    plot(t/tau,soldl1(:,k),'color',transparency(k).*[0 0 1],'LineWidth',2);
    hold on
    legendString(i)='x = '+string(round(x(k),2))+' cm';
    i=i+1;
end
xlabel('t/\tau')
ylabel('V')
set(gca,'Fontsize',15)
legend(legendString)
title('Linearly Decreasing Diameter')
% xlim([0 1])

figure
plot(t./tau, sol2(:,500),'LineWidth',2)
hold on
plot(t./tau,soldl1(:,500),':','LineWidth',2)
xlabel('t/\tau')

xlabel('t/\tau')
ylabel('V (V)')
set(gca,'Fontsize',15)
legend('Constant Diameter','Linear Diameter')



figure
hold on
plot(t.*10^3,soldl1(:,1),'LineWidth',2)
plot(t.*10^3,soldl1(:,end/10),'LineWidth',2)
plot(t.*10^3,soldl1(:,end/2),'LineWidth',2)
plot(t.*10^3,soldl1(:,end),'LineWidth',2)
ylabel('V(\mum)')
xlabel('Time (ms)')
legend (string(round(x(2),3))+' cm',string(round(x(end/10),2))+' cm',string(round(x(end/2),2))+' cm',string(x(end))+' cm')
set(gca,'Fontsize',15)


figure
lambdad=0.1;
d=-7.7500e-04*10*x+0.00037;
plot(x,d,'LineWidth',2);
ylabel('d(x)')
xlabel('X (\mum)')
set(gca,'Fontsize',15)
% legend ('t=0', 't=end')

% ---- FUNCTIONS ------------------

% Linear cable equation
% (tau/lmda2)dudt = d^2u/dX^2 - (u/lmda2)
function [c,f,s] = cable(x,t,u,dudx,do)
rl=0.1e3;
cm=1e-6;
% do=-7.7500e-04
d=do*x+0.00037;
dddx=do;
tau=0.05;

c=4*rl*cm/d;
f=dudx;


s=-u*4*rl*cm/(tau*d)+dudx*2/d*dddx;



%
% c = tau/lmda2;
% f = dudx;
% s = -(lmda2^-1)*u;

end

% Initial conditon
function u0 = cableic(x,do)

% u0 = 50;
L=0.04;
if x > 0.99*L/2 && x < 1.01*L/2
% if x == 0
    u0 = 50; %mV
else
    u0 = 0; %mV
end

%     l = sqrt(lmda2);
% %     I =  0; % injected current
% %     u0 = ri*I*sqrt(l)*exp(-x./l);
%     u0 = 50*exp(-x./l);
end

% Boundary equations
function [pl,ql,pr,qr] = cablebc(xl,ul,xr,ur,t,do)

%pl = 50*(1-exp(-(t./tau)));
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end