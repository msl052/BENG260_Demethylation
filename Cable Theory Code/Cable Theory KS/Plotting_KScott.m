figure
plot(t./tau, sol1(:,500),'LineWidth',2)
hold on
plot(t./tau,sol2(:,500),'LineWidth',2)
plot(t./tau,soldl(:,500),':','LineWidth',2)
plot(t./tau,soldl2(:,500),':','LineWidth',2)
xlabel('t/\tau')

xlabel('t/\tau')
ylabel('V (V)')
set(gca,'Fontsize',15)
legend('d=0.00006 cm','d=0.00037','d=-7.75e-04*x+0.00037','d=-7.75e-03*x+0.00037')
title('Voltage For x=L/2')


%%

figure
plot(t./tau, sol1(:,500),'LineWidth',2)
hold on
plot(t./tau,sol2(:,500),'LineWidth',2)
plot(t./tau,solAtand(:,500),'color', [.5 0 .5],'LineWidth',2)

xlabel('t/\tau')

xlabel('t/\tau')
ylabel('V (V)')
set(gca,'Fontsize',15)
legend('d=0.00006 cm','d=0.00037','d=0.00037*1/pi*(atan(1000*(-x+L/2))+pi/2)')
title('Voltage For x=L/2')
%%
close all
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated3.gif';

for n=linspace(1,length(t),100)
    % Draw plot for y = x.^n
    plot(x, sol1(n,:),x, sol2(n,:),'LineWidth',2)

%     plot(x, sol2(n,:),'LineWidth',2)
    ylim([0 max(max(sol1(:,:)))])
    xlabel('x(cm)')
    ylabel('V (V)')
    legend('d=0.00006 cm','d=0.00037 cm')
    
    set(gca,'Fontsize',15)
   drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end
