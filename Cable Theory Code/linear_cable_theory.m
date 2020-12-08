% Linear Cable Equation
% Stephanie Sincomb 
% BENG 260

for ct = 1%:6 % to vary the diameter

% constants (Basser 1993 table 1)
do = 1.5e-3; % outer myelin diameter [cm]
di = 0.7*do;
%di = (1-(0.1*ct))*do; % inner diameter
%di = 0.7*do; % inner diameter

% Solve equation linear cable equation 
% tau*dVdt = lbda^2*d^2V/dx^2 - V
L = 100*do; % [cm]
x = linspace(0,9*L,500);  %[cm]
dx = diff(x); dx = dx(1);
t = linspace(0,1.5,10000); %[ms]

% V(x,t)= transmembrane potential [mV]
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
m =0;
sol = pdepe(m,@cable,@cableic,@cablebc,x,t,[],di,L,dx);

% --------- PLOT SOLUTION -------------------------------------
%plot_cable_theory % plot IC 50 mV

%plot_cable_theory_2 % plot for IC 50 mV 

%plot_cable_theory_3 % plot for Injected current

ii = find(t<=0.5); ii=ii(end);
kk = find(x<=(1.1*L)); kk = kk(end);
[solmax, idx] =max(sol(:,kk));
yval(ct) = sol(ii,kk)./solmax;
xval(ct) = x(kk);

end

% ---- FUNCTIONS ------------------

% Linear cable equation 
% (tau/lmda2)dudt = d^2u/dX^2 - (u/lmda2) 
function [c,f,s] = cable(x,t,u,dudx,di,L,dx)

Cm = 5e-3; % myelin capacitance per unit area [muF/ cm^2]
Rm = R_m(x,L); % myelin resistance [kOhms-cm^2]]
Ri = 0.14; % axoplasm resistance [kOhm-cm]

% per unit length
ri = (4*Ri)/(pi*di^2); % [kOhm/cm]
rm = Rm/(pi*di) ; % [kOhms-cm]
cm = Cm*pi*di; % [muF/cm]

% time constant 
tau = rm*cm; %[ms]

% space constant lambda^2
lmda2 = rm/ri; %[cm^2]

c = tau/lmda2;
f = dudx;

s = -(lmda2^-1)*u + ((rm/10)/(lmda2))*ie(x,t,tau); %including 
%s = -(lmda2^-1)*u;

end

% Initial conditon
function u0 = cableic(x,di,L,dx)

u0 = 0;
%     if x == 0
%         u0 = 50; %mV
%     else
%         u0 = 0; %mV
%     end

end

% Boundary equations
function [pl,ql,pr,qr] = cablebc(xl,ul,xr,ur,t,di,L,dx)

    %pl = 50*(1-exp(-(t./tau)));
    pl = 0; %Amps
    ql = 1;
    pr = 0;
    qr = 1; 
end

% injected current function
function i= ie(x,t,tau)

if round(x*100) ==7  && t >= 0.5/4 && t <= 0.5/2
%     v0 = 50e-3; %V
%     do = 1.5e-3; % outer myelin diameter [cm]
%     di = 0.7*do; % inner diameter
%     Ri = 0.14; % axoplasm resistance [kOhm-cm]
%     Rm = pi*di*rm;
%     i = v0*pi*sqrt(di^3/(Ri*1e3*Rm*1e3));
    i = 1; %Amps
else 
    i = 0;
end

end

function Rm = R_m(x,L,dx)

if x<=L
    Rm = 100; % kOhm-cm^2
else
    %Rm = 100; %kOhm-cm^2
    Rm = 20; %kOhm-cm^2
end

end