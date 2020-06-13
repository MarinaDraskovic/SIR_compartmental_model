clear; close all; clc;
addpath('Data Serbia');

[beginDate, total_cases, active_cases, new_cases, total_recovered, total_deaths]=getDataSerbia();
t=1:length(total_cases); %Number of days

figure(1)
subplot(211) 
hold on; grid on;
plot(t, total_cases,'b-o',t, active_cases,'r-o', 'LineWidth',1);
plot(t, total_recovered,'g-o',t, total_deaths,'k-o', 'LineWidth',1);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})
ylabel('Cases');
title({'Data Serbia'});
legend ('Total cases', 'Active cases', 'Total recovered', 'Total deaths');

subplot(212)
semilogy(t, total_cases,'b-o', 'LineWidth',1);
hold on; grid on;
semilogy(t, active_cases,'r-o', 'LineWidth',1);
semilogy(t,total_recovered,'g-o', 'LineWidth',1);
semilogy(t, total_deaths,'k-o', 'LineWidth',1);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})
ylabel('Cases');
title({'Data Serbia, Log scale'});
legend ('Total cases', 'Active cases', 'Total recovered', 'Total deaths');

figure(2)
bar(new_cases);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})
ylabel('Cases');
title({'Daily new cases'});
