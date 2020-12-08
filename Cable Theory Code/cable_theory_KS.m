% Linear Cable Equation
% Stephanie Sincomb 
% BENG 260
% close all
% constants (Basser 1993 table 1)
do = 1.5e-3; % outer myelin diameter [cm]
di = 0.7*do; % inner diameter
Cm = 5e-6; % myelin capacitance per unit area [muF/ cm^2]
Rm = 100; % myelin resistance [kOhms-cm^2]]
Ri = 0.14; % axoplasm resistance [kOhm-cm]

% per unit length
ri = (4*Ri)/(pi*di^2); % [kOhm/cm]
rm = Rm/(pi*di) ; % [kOhms-cm]
cm = Cm*pi*di; % [muF/cm]

% time constant 
tau = rm*cm; %[s]

% space constant lambda^2
lmda2 = rm/ri; %[cm^2]

% Solve equation linear cable equation 
% tau*dVdt = lbda^2*d^2V/dx^2 - V
L = 100*do; % [cm]
x = linspace(0,4*L,100);  %[cm]
t = linspace(0,3*tau,1000); %[s]

% V(x,t)= transmembrane potential [mV]
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
m = 0;
sol = pdepe(m,@cable,@cableic,@cablebc,x,t, [], tau, lmda2, Rm,rm,ri);

% --------- PLOT SOLUTION -------------------------------------
surf(x,t,sol, 'EdgeColor','none')
colorbar()
xlabel('time (ms)')
ylabel('x(cm)')
zlabel('V (mv)')

figure
transparency = linspace(0.2, 0.9, length(x)-10);
i=1
for k = 5:20:length(x)-10
   plot(t/tau,sol(:,k)/max(sol(:,k)),'color',transparency(k).*[0 0 1],'LineWidth',2); 
   hold on
   legendString(i)='x/\lambda^2 = '+string(round(x(k)/sqrt(lmda2),2))+' cm';
   i=i+1;
end
xlabel('t/\tau')
ylabel('V/Vmax') 
set(gca,'Fontsize',15)
legend(legendString)
% xlim([0 1])


figure
for k = 5:100:length(t)
    plot(x./sqrt(lmda2), sol(k,:),'LineWidth',2); 
    hold on
end
xlabel('x')


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
d=(-tanh(5*(x./sqrt(lmda2)))+1)
plot(x,d,'LineWidth',2);
ylabel('d(x)/d_0')
xlabel('X (\mum)')
% legend ('t=0', 't=end')

% ---- FUNCTIONS ------------------

% Linear cable equation 
% (tau/lmda2)dudt = d^2u/dX^2 - (u/lmda2) 
function [c,f,s] = cable(x,t,u,dudx,tau, lmda2,Rm,rm,ri)

do = 1.5e-3; % outer myelin diameter [cm]

di=0.7*do*(-tanh(2*(x/sqrt(lmda2)))+1)
Cm = 5e-6; % myelin capacitance per unit area [muF/ cm^2]
Rm = 100; % myelin resistance [kOhms-cm^2]]
Ri = 0.14; % axoplasm resistance [kOhm-cm]

ri = (4*Ri)/(pi*di^2); % [kOhm/cm]
rm = Rm/(pi*di) ; % [kOhms-cm]
cm = Cm*pi*di; % [muF/cm]
% 
% % time constant 
tau = rm*cm; %[s]

% % per unit length
% 
% % d=sin(x*0.2/lambdad);
% 
c = di*4*(tau)*ri/rm;
f=di^2*dudx;
s=-(4*ri*di/(rm))*u;
% 
% c = tau/lmda2;
% f = dudx;
% s = -(lmda2^-1)*u;

end

% Initial conditon
function u0 = cableic(x,tau, lmda2, Rm,rm,ri)

% u0 = 50;
    if x == 0
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
function [pl,ql,pr,qr] = cablebc(xl,ul,xr,ur,t,tau, lmda2, Rm,rm,ri)

    %pl = 50*(1-exp(-(t./tau)));
    pl = 0;
    ql = 1;
    pr = 0;
    qr = 1; 
end