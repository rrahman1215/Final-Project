%https://basicmedicalkey.com/fitting-models-to-data/

function MonkepoxSirFittingv
% This function fits the Monkeypox data 2022 and SIR model 
clear alll; close all; clc
data_inc = readtable('time_series_df.csv','PreserveVariableNames',true); % Daily Reported Monkeypox Data loading data
format long        % specifying higher precision
 
%% Data Information
%beta = 1.75;
%gamma = 0.650;
%days = 33;                 % Simulate for D Weeks
N=7.837*10^(9);              % World Population 
%N = 7000*10^(9);              % World Population 
eta = 10^(-5);
%%
qdata = table2array(data_inc(:,3));      % define array with y−coordinates of the data
                                         % observed data
tdata =table2array(data_inc(:,1));       %1:1:length(qdata); % % define array with t−coordinates of the data

tforward = tdata' %1:1:length(qdata);            % t mesh for the solution of the differential equation

%tmeasure =1.0e+8*[0.000000000003000 0.000000000008964 0.000000000026782 0.000000000080019 0.000000000239083 0.000000000714342 0.000000002134338 0.000000006377054 0.000000019053600 0.000000056929060 0.000000170094773 0.000000508215488 0.000001518464806 0.000004536921747 0.000013555546170 0.000040501423148 0.000121008581905 0.000361526405454 0.001079936064991 0.003224478372744 0.009614668104520 0.028553783711477 0.083799460095178 0.237667480227084 0.614705951506542 1.288918054578621 1.919298274467157 2.002320865547142 1.632300231611211 1.160912122207342 0.768796635774286 0.490126085055269 0.305971199743162]'
                  % actual multiplication 1.0e+12*[0.000000000003000]


%B(1:33) %[1:101:3301]';  % selects the points in the solution
                          % corresponding to the t values of tdata

%initial values of parameters to be fitted
beta  = 1.75;  %7; %1.85; %1.75;       % beta 
gamma = 1.0 ; %65;  % 5.650;      % alpha 
k=[beta gamma]
 %------------------------------------------------------------%
function dy = model_1(t,y,k)                        % DE

           beta = k(1);                                % Assignes the parameters in the DE the current
           gamma = k(2);                                % value of the parameters 
           dy = zeros(2,1);                         % assigns zeros to dy 
           dy(1) = -beta*y(1)*y(2)/N;                  % RHS of first equation
           dy(2) =  beta*y(1)*y(2)/N- gamma*y(2);          % RHS of second equation
end 

 %------------------------------------------------------------%

function error_in_data = moder(k)                   % computing the error in the data
                 [T Y] = ode45(@(t,y)(model_1(t,y,k)),tforward,[N-3 3]);
                                                    % solves the DE; output is written in T and Y 
                     q = eta*Y(:,2);           % assignts the y−coordinates of the solution at
                       %q= Y(:,2);                   % at the t−coordinates of tdata
        error_in_data = sum((q-qdata).^2)           %computes SSE
end

% Finding the optimum parameters k= [beta gamma] by fitting qdata and model values ql

[k1,fval] =  fminsearch(@moder,k,optimset); % minimization routine; assigns the new

disp(k1);                          %print final values of beta and gamma

fprintf('the value of k % .9f\n',k1)
%k2=[1.001409504007677, 1]
%k3 =[0.94705776, 0.60745201]

[T Y] = ode45(@(t,y)(model_1(t,y,k1)),tforward,[N-3  3]);
                        % solving the DE with the final values of the
                        % parameters
yint = eta*Y(:,2); % computing the y−coordinates corresponding to the t data
%------------------------------------------------------------%

%subplot(1,2,2)
plot(tdata,qdata,'-ro','linewidth',2);
hold on
plot(tdata,yint,'-bo','linewidth',2);
legend('obs','est', 'Location', 'Best')
%xlabel('time(per week)');         % plotting final fit
%ylabel('Number of cases');
%axis([3 14 0 350]);
set(gca,'Fontsize',15)
end 
