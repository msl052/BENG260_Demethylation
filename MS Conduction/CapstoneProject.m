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
labels = ["2","3","A","B","D1","D2","D3","D4","4"];
ismyelin = [true,true,true,true,false,false,false,false,true];
myelinlen = [L,1600/10000,400/10000,400/10000,NL,NL,NL,NL,L];

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
      
      ihh(klok)=gNa*(v-ENa)+gK*(v-EK)+gLbar*(v-EL); % muA/cm^2
      if ismyelin(i) == false % demyelinated regions
        im(klok) = gm*(v+69.92)*(10^3)/myelinlen(i); %current density (muA/cm^2) (x10^3 for converting to mu)
        C = Cstar*10^6/(pi*d);
        %C = Cstar*10^6;
        %C = Cm*10^6/(pi*d);
      else
        im(klok) = (ihh(klok)*pi*d*NL/dx+(dx-NL)/dx*v*gm*(10^3))/myelinlen(i); % muA/cm^2
        C = Cm*10^6/(pi*d);
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
      
      %Matizero(klok) = izero(t);
      
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
    MatIhh(i,:) = ihh;
    MatI(i,:) = im;
    MatV(i,:) = v_plot;
    Matmhn(i,:,:) = mhn_plot;

%     figure;
%     hold on;
%     plot(t_plot,mhn_plot);
%     legend('m-gate','h-gate','n-gate')
%     xlabel('time [ms]')
%     ylabel('gating variable')


end

figure;
subplot(2,1,1), plot(t_plot,MatV);
[M,I] = max(MatV,[],2);
for i = 1:length(labels)
    if M(i) > -69
        text(t_plot(I(i)),MatV(i,I(i)),labels(i));
    end
end
ylabel('Membrane potential [mV]')
legend(labels);
subplot(2,1,2), plot(t_plot,MatI);
[M,I] = max(MatI,[],2);
for i = 1:length(labels)
    if M(i) > 1
        text(t_plot(I(i)),MatI(i,I(i)),labels(i));
    end
end
Ipks = islocalmax(MatI,2);
ylabel('Current [muA/cm^2]')
legend(labels);
%subplot(3,1,3), plot(t_plot,MatIhh);
%ylabel('Current density [muA/cm^2]')
xlabel('time [ms]')

% figure;
% hold on;
% plot(t_plot,Matizero);
% xlabel('time [ms]')
% ylabel('stimulus [mua/cm^2]')
% ylim([0,100]);

% figure;
% hold on;
% plot(t_plot,Matmhn);
% legend('m-gate','h-gate','n-gate')
% xlabel('time [ms]')
% ylabel('gating variable')

