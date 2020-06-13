%% Initialization
clear; close; clc;
format long;

addpath('Data China');
[beginDate, ~, active_cases, total_recovered] = getDataChina();   
% addpath('Data Serbia');
% [beginDate, ~, active_cases, ~, total_recovered, ~] = getDataSerbia();

adata=active_cases';
rdata =total_recovered';
times=1:1:length(adata); 

%% Model without control 
% Parameter estimation of basic SIR model
% paramests = [beta, gamma, k]
[paramests, test, yest] = SIR_model_paramest(adata, rdata, times);

% Models with control 
dy1 = SIR_model_control_personal_protection(paramests, adata);
dy2 = SIR_model_control_treatment(paramests, adata);
dy3 = SIR_model_control_persProtection_treatment(paramests, adata);
dy4 = SIR_model_control_vaccination(paramests, adata);

%% Plotting SIR model and paramest
figure(1)
hold on; grid on;
plot(1:length(adata), adata, 'ko', 'LineWidth', 2);
plot(test, yest(:,2)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]});
ylabel('Infected cases');
title('Param. ests: SIR model China with added k, Fit using PSO, Without control');
legend('Real data',...
    'Model with parameter estimates', 'Location', 'northeast', 'FontSize', 12);

figure(2)
hold on; grid on;
plot(1:length(rdata), rdata, 'ko', 'LineWidth', 2);
plot(test, yest(:,3)*paramests(3), 'y', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]});
ylabel('Recovered cases');
title('Param. ests: SIR model China with added k, Fit using PSO, Without control');
legend('Real data',...
    'Model with parameter estimates', 'Location', 'northwest', 'FontSize', 12);

figure(3)
hold on; grid on;
plot(test, yest*paramests(3) ,'LineWidth',2)
plot(1:length(total_recovered), total_recovered, 'oy', 'LineWidth',2);
plot(1:length(active_cases), active_cases, 'or', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Cases');
title('Param. ests: SIR model China with added k, Fit using PSO, Without control');
legend('Susceptible','Infected','Recovered', 'Real recovered cases','Real total cases',...
    'Location','northwest','FontSize',12)

%% Ploting models with and without control
figure(4) % Personal protection
sgtitle('China: Added personal protection as control parameter');
subplot(221);
hold on; grid on;
plot(dy1(1,:), dy1(2,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,1)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Suspectible Cases'); 
legend('S - with control', 'S - without control',...
    'Location','northeast','FontSize',12);

subplot(222);
hold on; grid on;
plot(dy1(1,:), dy1(3,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,2)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Infectious Cases');
legend('I - with control', 'I - without control',...
    'Location','northeast','FontSize',12);

subplot(223);
hold on; grid on;
plot(dy1(1,:), dy1(4,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,3)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Recovered Cases');
legend('R - with control', 'R - without control',...
    'Location','northwest','FontSize',12);

subplot(224); grid on;
plot(dy1(1,:), dy1(5,:), 'r', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('u1^*');

%_________________________________________________
figure(5) % Treatment
sgtitle('China: Added treatment as control parameter');
subplot(221);
hold on; grid on;
plot(dy2(1,:), dy2(2,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,1)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Suspectible Cases'); 
legend('S - with control', 'S - without control',...
    'Location','northeast','FontSize',12);

subplot(222);
hold on; grid on;
plot(dy2(1,:), dy2(3,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,2)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Infectious Cases');
legend('I - with control', 'I - without control',...
    'Location','northeast','FontSize',12);

subplot(223);
hold on; grid on;
plot(dy2(1,:), dy2(4,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,3)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Recovered Cases');
legend('R - with control', 'R - without control',...
    'Location','northwest','FontSize',12);

subplot(224); grid on;
plot(dy2(1,:), dy2(5,:), 'r', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('u2^*');

%_________________________________________________________
figure(6) % Personal protection and treatment
sgtitle('China: Added personal protection and treatment as control parameters');
subplot(231);
hold on; grid on;
plot(dy3(1,:), dy3(2,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,1)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Suspectible Cases'); 
legend('S - with control', 'S - without control',...
    'Location','northeast','FontSize',12);

subplot(232);
hold on; grid on;
plot(dy3(1,:), dy3(3,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,2)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Infectious Cases');
legend('I - with control', 'I - without control',...
    'Location','northeast','FontSize',12);

subplot(233);
hold on; grid on;
plot(dy3(1,:), dy3(4,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,3)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Recovered Cases');
legend('R - with control', 'R - without control',...
    'Location','northwest','FontSize',12);

subplot(2,3,4.5); grid on;
plot(dy3(1,:), dy3(5,:), 'r', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('u1^*');

subplot(2,3,5.5);
grid on;
plot(dy3(1,:), dy3(6,:), 'r', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('u2^*');

%_________________________________________________
figure(7) % Vaccination
sgtitle('China: Added vaccination as control parameter');
subplot(221);
hold on; grid on;
plot(dy4(1,:), dy4(2,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,1)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Suspectible Cases'); 
legend('S - with control', 'S - without control',...
    'Location','northeast','FontSize',12);

subplot(222);
hold on; grid on;
plot(dy4(1,:), dy4(3,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,2)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Infectious Cases');
legend('I - with control', 'I - without control',...
    'Location','northeast','FontSize',12);

subplot(223);
hold on; grid on;
plot(dy4(1,:), dy4(4,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,3)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Recovered Cases');
legend('R - with control', 'R - without control',...
    'Location','northwest','FontSize',12);

subplot(224); grid on;
plot(dy4(1,:), dy4(5,:), 'r', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('u3^*');