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
in_HH;
in_mhnv;
%pNodal = (load('NodalCurrentDensity.mat','ihh').ihh)

%V = zeros(numnodes,klokmax);
ismyelin = [true, true, true, true true true, false, false, false, false, true, true, true true];
myelinlen = [L,L,L,L3,A,B,NL,NL,NL,NL,L,L,L,L];
for i = 1:length(ismyelin)
    
    for klok=1:klokmax
      t=klok*dt;                      %note time
      x=klok*dx;
      m=snew(m,alpham(v),betam(v),dt); %update m
      h=snew(h,alphah(v),betah(v),dt); %update h
      n=snew(n,alphan(v),betan(v),dt); %update n
      gNa=gNabar*(m^3)*h;    %sodium conductance
      gK=gKbar*(n^4);    %potassium conductance
      g=gNa+gK+gLbar;         %total conductance
      gE=gNa*ENa+gK*EK+gLbar*EL;         %gE=g*E  (current density muA/cm^2)
      
      if ismyelin(i) == false % demyelinated regions
        im(klok) = gm*v*(10^3)/myelinlen(i); %current density (muA/cm) (x10^3 for converting to mu)
        C = Cn*10^6/(pi*d);
      else
        ihh(klok)=gNa*(v-ENa)+gK*(v-EK)+gLbar*(v-EL); % muA/cm^2
        im(klok) = (ihh(klok)*pi*d*NL/dx+(dx-NL)/dx*v*gm*(10^3))/myelinlen(i); % muA/cm
        C = Cn*10^6/(pi*d);
        %C = Cstar*10^6/(pi*d); %(muF/cm^2)
        %im(klok) = ihh(klok);
      end
      
      %save old value of v for checking purposes:
      v_old=v;
      %update v:
      if i == 1
        v=(v+(dt/C)*(gE+izero(t)))/(1+(dt/C)*g);
      else
        v=(v+(dt/C)*(gE+pI(klok)))/(1+(dt/C)*g);
      end
      
      if(check)
        E=gE/g;
        chv=C*(v-v_old)/dt+g*(v-E)-izero(t);
      end
      
      %store results for future plotting:
      mhn_plot(:,klok)=[m h n]';
      v_plot(klok)=v;
      t_plot(klok)=t;  
    end
    vhold, vstart = v_plot(end);
    in_mhnv;
    pI = im;
    max(pI)
    MatIhh(i,:) = ihh;
    MatI(i,:) = im;
    MatV(i,:) = v_plot;
%     subplot(numnodes,1,i),plot(t_plot,pI);
%     ylabel('pI [muA/cm]')
%     legend;
%     hold on;
%     subplot(2,1,1),plot(t_plot,v_plot), hold on;
%     legend('membrane potential')
%     ylabel('Membrane potential [mV]')
%     subplot(2,1,2),plot(t_plot,mhn_plot), hold on;
%     legend('m-gate','h-gate','n-gate')
%     xlabel('time [ms]')
%     ylabel('gating variable')


end
figure;
subplot(3,1,1), plot(t_plot,MatV);
ylabel('[mV]');
title('Membrane potential');
legend('IN1','IN2','IN3','IN(Shortened)','INA','INB','D1','Others Show No Action Potential');
subplot(3,1,2), plot(t_plot,MatI);
ylabel('Current [muA]')
title('Ionic Membrane Current Density for Demyelinated Regions');
legend('IN1','IN2','IN3','IN(Shortened)','INA','INB','Others Show No Current');
subplot(3,1,3), plot(t_plot,MatIhh);
ylabel('Current [muA]')
xlabel('time [ms]')
title('Hodgkin-Huxley Current Density');
legend('IN1','IN2','IN3','IN(Shortened)','INA','INB','D1 (No Current)','D2 (No Current)','D3 (No Current)','D4 ((Right Side) Yellow Spike)','Others Show No Current');

