% Linear Cable Equation
% Kiersten Scott
% BENG 260
% close all
clear all
% constants (Basser 1993 table 1)

% Solve equation linear cable equation 
% tau*dVdt = lbda^2*d^2V/dx^2 - V
L = 0.04; % [cm]
tau=.00064e-3;
x = linspace(0,L,1000);  %[cm]
t = linspace(0,200*tau,10000); %[s]
% V(x,t)= transmembrane potential [mV]
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
m = 0;
do=0.00006
sol1 = pdepe(m,@cable,@cableic,@cablebc,x,t, [],do);
do=0.00037
sol2 = pdepe(m,@cable,@cableic,@cablebc,x,t, [],do);

% --------- PLOT SOLUTION -------------------------------------
figure
surf(x,t,sol1, 'EdgeColor','none')
colorbar()
xlabel('time (ms)')
ylabel('x(cm)')
zlabel('V (mv)')

figure
surf(x,t,sol1, 'EdgeColor','none')
colorbar()
xlabel('time (ms)')
ylabel('x(cm)')
zlabel('V (mv)')

figure
transparency = fliplr(linspace(0.2, 0.9, length(x)-10));
i=1
for k = 5:200:length(x)-10
    plot(t/tau,sol1(:,k),'color',transparency(k).*[0 0 1],'LineWidth',2);
    hold on
    legendString(i)='x = '+string(round(x(k),2))+' cm';
    i=i+1;
end
xlabel('t/tau')
ylabel('V')
set(gca,'Fontsize',15)
legend(legendString)
title('Linearly Decreasing Diameter')



figure
plot(t./tau, sol1(:,500),'LineWidth',2)
hold on
plot(t./tau,sol2(:,500),':','LineWidth',2)
xlabel('t/\tau')

xlabel('t/\tau')
ylabel('V (V)')
set(gca,'Fontsize',15)
legend('d=0.00006 cm','d=0.00037')





% figure
% hold on
% % plot(t/tau,sol(:,2)./max(sol(:,2)),'LineWidth',2)
% % plot(t/tau,sol(:,end/10)./max(sol(:,end/10)),'LineWidth',2)
% % plot(t/tau,sol(:,end/2)./max(sol(:,end/2)),'LineWidth',2)
% % plot(t/tau,sol(:,end)./max(sol(:,end)),'LineWidth',2)
% plot(t.*10^3,sol(:,2),'LineWidth',2)
% plot(t.*10^3,sol(:,end/10),'LineWidth',2)
% plot(t.*10^3,sol(:,end/2),'LineWidth',2)
% plot(t.*10^3,sol(:,end),'LineWidth',2)
% ylabel('V(\mum)') 
% xlabel('Time (ms)')
% legend (string(round(x(2),2))+' cm',string(round(x(end/10),2))+' cm',string(round(x(end/2),2))+' cm',string(x(end))+' cm')
% set(gca,'Fontsize',15)
% xlim([0 1])
figure
lambdad=0.1;
d=0.00037*(-x/0.00037+1)
plot(x,d,'LineWidth',2);
ylabel('d(x)/d_0')
xlabel('X (\mum)')
% legend ('t=0', 't=end')

% ---- FUNCTIONS ------------------

% Linear cable equation 
% (tau/lmda2)dudt = d^2u/dX^2 - (u/lmda2) 
function [c,f,s] = cable(x,t,u,dudx,do)
rl=0.1e3;
cm=1e-6;
d=do;
% d=0.00037*(-x/0.02+1)
tau=0.05;

c=4*rl*cm/d;
f=dudx;

s=-u*4*rl*cm/(tau*d);



% 
% c = tau/lmda2;
% f = dudx;
% s = -(lmda2^-1)*u;

end

% Initial conditon
function u0 = cableic(x,do)
L=0.04
% u0 = 50;
if x > 0.99*L/2 && x < 1.01*L/2
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