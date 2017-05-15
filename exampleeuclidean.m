[azel,curvatures, poleclass]=PCAlenght;
figure
subplot(2,2,1)
plot(azel(944,:),'r')
hold on
plot(azel([531,527],:)','b')
legend('Pole','Air')
xlabel('Time [ms]')
ylabel('Angle [Degrees]')
title('[Azimuth,elevation]')
subplot(2,2,2)
plot(curvatures(944,:),'r')
hold on
plot(curvatures([531,527],:)','b')
title('[K horizontal,K coronal]')
legend('Pole','Air')
xlabel('Time [ms]')
ylabel('Curvature')
subplot(2,2,3)
plot(azel(1007,:),'r')
hold on
plot(azel([546,547,548,549,551,552],:)','b')
legend('Pole','Air')
xlabel('Time [ms]')
ylabel('Angle [Degrees]')
title('[Azimuth,elevation]')
subplot(2,2,4)
plot(curvatures(1007,:),'r')
hold on
plot(curvatures([546,547,548,549,551,552],:)','b')
title('[K horizontal,K coronal]')
legend('Pole','Air')
xlabel('Time [ms]')
ylabel('Curvature')