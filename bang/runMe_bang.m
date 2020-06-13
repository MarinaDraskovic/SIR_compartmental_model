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
dy1 = bang_personal_protection(paramests, adata);
dy2 = bang_treatment(paramests, adata);
dy3 = bang_vaccination(paramests, adata);

figure(1) % Personal protection
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

%__________________________________________________
figure(2) % Treatment
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
figure(3) % Vaccination
sgtitle('China: Added vaccination as control parameters');
subplot(221);
hold on; grid on;
plot(dy3(1,:), dy3(2,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,1)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Suspectible Cases'); 
legend('S - with control', 'S - without control',...
    'Location','northeast','FontSize',12);

subplot(222);
hold on; grid on;
plot(dy3(1,:), dy3(3,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,2)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Infectious Cases');
legend('I - with control', 'I - without control',...
    'Location','northeast','FontSize',12);

subplot(223);
hold on; grid on;
plot(dy3(1,:), dy3(4,:)*paramests(3), 'r', 'LineWidth',2);
plot(test, yest(:,3)*paramests(3), 'b', 'LineWidth', 2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('Recovered Cases');
legend('R - with control', 'R - without control',...
    'Location','northwest','FontSize',12);

subplot(224); grid on;
plot(dy3(1,:), dy3(5,:), 'r', 'LineWidth',2);
xlabel({'Time [days]'; ['Starting date: ' datestr(beginDate)]})  
ylabel('u3^*');
