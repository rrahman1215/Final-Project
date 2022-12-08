clc; clear all; close all;
%% Data Information
data_inc = readtable('time_series_df.csv','PreserveVariableNames',true); % Daily Reported Monkeypox Data loading data
format long        % sdypecifying higher precision
 
%% Data Information
%N=7.9*10^(9);              % World Population 
qdata = table2array(data_inc(:,3));      % define array with y−coordinates of the data
tdata =table2array(data_inc(:,1));       %1:1:length(qdata); % % define array with t−coordinates of the data
%tforward =tdata'; % 1:1:length(qdata);            % t mesh for the solution of the differential equation
% %0.94705776, 0.60745201
% beta = 0.94705776;
% gamma = 0.60745201;
beta = 1.75;
gamma = 0.650;
k= [beta gamma];
P=7.837*10^(9);          % World Population 
% Define the func11tions as vector 
f1 = @(t,Y) -beta*Y(1)*Y(2)/P;
f2 = @(t,Y) beta*Y(1)*Y(2)/P-gamma*Y(2);
f3 = @(t,Y) gamma*Y(2);
Fun = @(t,Y) [f1(t,Y); f2(t,Y); f3(t,Y)];

% Initial condition
Y0 = [P-3;3;0];                    % Initial value of S, I , and R 
a = 0;                           % Time start 
b = 33;                          % Time end
t = 1:1:32;                        % time observed vector
tspan=0:1:33;
N = 33;                         % Number of iterations  
[t1,Y1] = eulsys(Fun,a,b,Y0,N);  % Forward Euler scheme 
[t2,Y2] = eulmodsys(Fun,a,b,Y0,N)% Modified Euler's method
[t3,Y3] = rk4sys(Fun,a,b,Y0,N);  % 4th-order Runge-Kutta method 
[t4,Y4] = ode45(@(t,y)(model_1(y,k,P)),tspan,[P-3  3  0]);
figure(1)
plot(t1, Y1(2,:),'-ro', 'LineWidth', 2);
hold on 
legend('I', 'Location','northwest');
%xlabel('hours');
%print('tmp', '-dpdf');  print('tmp', '-dpng');

%figure(2)
plot(t2, Y2(2,:),'-go', 'LineWidth', 2);
%hold on 
legend('I','Location','northwest');
%xlabel('hours');
%print('tmp', '-dpdf');  print('tmp', '-dpng');

%figure(3)
plot(t3, Y3(2,:),'-bo', 'LineWidth', 2);
%hold on 
legend('I','Location','northwest');
%xlabel('hours');
%print('tmp', '-dpdf');  print('tmp', '-dpng');


%figure(4)
plot(t4, Y4(:,2),'-mo', 'LineWidth', 2);
%hold on 
legend('I','Location','northwest');
%xlabel('hours');
%print('tmp', '-dpdf');  print('tmp', '-dpng');
legend('Euler','Mod Euler','RK4','ode45')
set(gca,'Fontsize',15)
hold off


I= Y3(2,:);
C= 10^(-4)*I;

filename = 'test.mat';
save(filename,'I')