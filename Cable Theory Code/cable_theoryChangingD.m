% Linear Cable Equation
close all 
clear all

% global variables
global ri tau lmda2 do rm lambdad

% constants (Basser 1993 table 1)
do = 0.5e-3; % outer myelin diameter [cm]
di = 0.7*do; % inner diameter
Cm = 5e-3; % myelin capacitance per unit area [muF/ cm^2]
Rm = 100; % myelin resistance [kOhms-cm^2]]
Ri = 0.14; % axoplasm resistance [kOhm-cm]

% per unit length
ri = (4*Ri)/(pi*di^2); % [kOhm/cm]
rm = Rm/(pi*di) ; % [kOhms-cm]
cm = Cm*pi*di; % [muF/cm]
lambdad=0.5
% time constant 
tau = rm*cm; %[s]

% space constant lambda^2
lmda2 = rm/ri;

% Solve equation linear cable equation 
% tau*dVdt = lbda^2*d^2V/dx^2 - V
L = 0.2; % [cm]
x = linspace(0,L,1000);  %[ cm]
t = linspace(0,tau,10000); %[s]

% V(x,t)= transmembrane potential [mV]
opts = odeset('RelTol',1e-5,'AbsTol',1e-7)
m =0;
sol = pdepe(m,@cable,@cableic,@cablebc,x,t);

% Linear cable equation 
% (tau/lmda2)dudt = d^2u/dX^2 - (u/lmda2)

%% plotting
figure
surf(x,t,sol, 'EdgeColor','none')
h=colorbar()
xlabel(h, 'V (mV)')
xlabel('X(cm)')
ylabel('Time (ms)')
zlabel('V (mV)')
set(gca,'Fontsize',15)
figure
hold on
plot(t,sol(:,1),'LineWidth',2)
plot(t,sol(:,end/4),'LineWidth',2)
plot(t,sol(:,end/2),'LineWidth',2)
plot(t,sol(:,end),'LineWidth',2)
ylabel('V(\mum)') 
xlabel('Time (ms)')
legend ('0 cm','0.05 cm', '0.1 cm','0.2 cm')
set(gca,'Fontsize',15)

figure
hold on
plot(x,sol(1,:),'LineWidth',2)
plot(x,sol(end,:),'LineWidth',2)
ylabel('V')
xlabel('space (\mum)')
legend ('t=0', 't=end')

figure
plot(x,1./(1+exp(x.*0.2/lambdad)));
ylabel('d(x)')
xlabel('X (\mum)')
legend ('t=0', 't=end')
%%

function [c,f,s] = cable(x,t,u,dudx)
global tau lmda2 do rm ri lambdad

% 
c = (tau)/lmda2; f = dudx; s = -(lmda2^-1)*u;

% For the chagning diameter case
d=1/(1+exp(x*0.2/lambdad));
c = d*4*(tau)*ri/rm;
f=d^2*dudx;
s=-(4*ri*d/(rm))*u;

end

% Initial conditons
function u0 = cableic(x)
% 
    if x > 0.01 && x < 0.012
        u0 = 50; %mV
    else
        u0 = 0; %Mv
    end
% u0=50
    global ri lmda2
    l = sqrt(lmda2);
%     I = 50 % injected current
%     u0 = ri*I*sqrt(l)*exp(-x./l);
%     u0 = 50*exp(-x./l);
end

% Boundary equations
function [pl,ql,pr,qr] = cablebc(xl,ul,xr,ur,t)
%.50*(tanh(20*(t-1))+1)/2
%     pl =50*(tanh(20*(t-1))+1)/2; 
%     ql = 1;
    pl = ul; 
    ql = 0;
    pr = ur;
    qr = 0; 
end

