% Author: Min Suk Lee
% Date: 10/25/2020
% Description: Modeling affects of demyelination on neuron signal
% propagation using HH-model
%----------------------------------------------------------------
clear;
clc;
close all;

%filename CapstoneProject
%numerical solution of the space-clamped Hodgkin-Huxley equations
global check;
global t1p t2p ip; %parameters for the function izero(t)
in_HH
in_mhnv

% pNodal = (load('NodalCurrentDensity.mat','ihh').ihh);
% 
% figure;
% time = (1:klokmax).*dt;
% plot(time,pNodal);
% ylabel('current density [muA/cm^2]');

figure;
for klok=1:klokmax
  t=klok*dt;                      %note time
  m=snew(m,alpham(v),betam(v),dt); %update m
  h=snew(h,alphah(v),betah(v),dt); %update h
  n=snew(n,alphan(v),betan(v),dt); %update n
  gNa=gNabar*(m^3)*h;    %sodium conductance
  gK=gKbar*(n^4);    %potassium conductance
  g=gNa+gK+gLbar;         %total conductance
  gE=gNa*ENa+gK*EK+gLbar*EL;         %gE=g*E  (current density)
  im(klok) = gm*v; %current density (muA/cm^2)
  %save old value of v for checking purposes:
  v_old=v;
  %update v:
  %v=(v+(dt/C)*(gE+pNodal(klok)))/(1+(dt/C)*g);
  v=(v+(dt/C)*(gE+izero(t)))/(1+(dt/C)*g);
  if(check)
    E=gE/g;
    chv=C*(v-v_old)/dt+g*(v-E)-izero(t)
  end
  %store results for future plotting:
  mhn_plot(:,klok)=[m h n]';
  v_plot(klok)=v;
  t_plot(klok)=t;  
end
subplot(2,1,1),plot(t_plot,v_plot)
legend('membrane potential')
ylabel('Membrane potential [mV]')
subplot(2,1,2),plot(t_plot,mhn_plot)
legend('m-gate','h-gate','n-gate')
xlabel('time [ms]')
ylabel('gating variable')


figure;
time = (1:klokmax).*dt;
plot(time,im);
ylabel('current density [muA/cm^2]');
%save('NodalCurrentDensity','ihh')
