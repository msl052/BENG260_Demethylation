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

demyelinated = false;
numnodes = 7;
%V = zeros(numnodes,klokmax);
for i = 1:numnodes
    
    if i == 4 | i==5 | i==6
        demyelinated = true; 
    end
    
    for klok=1:klokmax
      t=klok*dt;                      %note time
      m=snew(m,alpham(v),betam(v),dt); %update m
      h=snew(h,alphah(v),betah(v),dt); %update h
      n=snew(n,alphan(v),betan(v),dt); %update n
      gNa=gNabar*(m^3)*h;    %sodium conductance
      gK=gKbar*(n^4);    %potassium conductance
      g=gNa+gK+gLbar;         %total conductance
      gE=gNa*ENa+gK*EK+gLbar*EL;         %gE=g*E  (current density muA/cm^2)
      
      if demyelinated == true
        im(klok) = gm*v/NL; %current density (muA/cm)
      else
        ihh(klok)=gNa*(v-ENa)+gK*(v-EK)+gLbar*(v-EL); % muA/cm^2
        im(klok) = (ihh(klok)*pi*d*NL/dx+(dx-NL)/dx*v*gm)/NL; % muA/cm
        %im(klok) = ihh(klok)
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
        chv=C*(v-v_old)/dt+g*(v-E)-izero(t)
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
subplot(2,1,1), plot(t_plot,MatV);
ylabel('Membrane potential [mV]')
subplot(2,1,2), plot(t_plot,MatI);
ylabel('Current [muA]')
xlabel('time [ms]')