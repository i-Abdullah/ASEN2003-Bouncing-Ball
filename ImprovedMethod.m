%% info

% this script is to re-estimate the Coefficient of restitution for the pong ball
% and also a tennis ball using an improved method, which depends on using audio waves to 
% estimate the time between bounces and total time it took the ball to stop

% Done by:
% 
% - Sarah Foley
% - Hugo Stetz
% - Alexander Lowry
% - Abdulla Al Ameri



%% read data

clear
clc
close all


Data = xlsread('Data/ImprovedMethodData.xlsx');

trial_pong = Data(:,1);
trial_golf = Data(:,3);
time_pong = Data(:,2); % time to stop.
time_golf = Data(:,4); % time to stop.

%% estimate errors:

% let the random error in each individual ball for time to stop be std because they
% in theory should have the same time.

random_error_pong = std(time_pong);
random_error_golf = std(time_golf);

% the insturment error really depends on when we decide to make the
% cut-off, this method is really really accurate, so we can take this error
% to be a very small fraction of a second!

inst_error_pong = 0.001 ;
inst_error_golf = 0.001 ;


time_total_error_pong = sqrt ( (random_error_pong)^2 + (inst_error_pong)^2 ) ;
time_total_error_golf = sqrt ( (random_error_golf)^2 + (inst_error_golf)^2 ) ;

h0_inches = 36 ; %inches.
g = 386.09 ; % gravity in inches/s^2

h0_error = 1/14 ; % the same for all trials


%% calculating e:


syms h0 ts 

e_stop_error(h0,ts) = (ts - sqrt((2*h0)/g))/(ts + sqrt((2*h0)/g));
%partial of e_stop with respect to h0;
Partial_stop_h0 = diff(e_stop_error,h0);
Partial_stop_ts = diff(e_stop_error,ts);


% for ping-pong ball

% do the math!

for i=1:length(time_pong)
    
e_pong(i) = (time_pong(i) - sqrt((2*h0_inches)/g))/(time_pong(i) + sqrt((2*h0_inches)/g));
total_error_pong(i) = double( sqrt ( (Partial_stop_h0(h0_inches,time_pong(i)) * h0_error) .^2 + ((Partial_stop_ts(h0_inches,time_pong(i)) * time_total_error_pong).^2)));

end

% repeat for golf ball.

for i=1:length(time_golf)
    
e_golf(i) = (time_golf(i) - sqrt((2*h0_inches)/g))/(time_golf(i) + sqrt((2*h0_inches)/g));
total_error_golf(i) = double( sqrt ( (Partial_stop_h0(h0_inches,time_golf(i)) * h0_error) .^2 + ((Partial_stop_ts(h0_inches,time_golf(i)) * time_total_error_golf).^2)));

end


%% plot results

figure(1)

errorbar(trial_golf,e_golf,total_error_golf,'-.*')
hold on
errorbar(trial_pong,e_pong,total_error_pong,'-.s')
grid minor
ylabel('Coefficient of restitution and error')
xlabel('Trial')
title('Coefficient of restitution for golf and ping-pong ball')
legend('e_g_o_l_f','e_p_o_n_g','Location','SouthEast')
xlim([ 0.5 10.5])

%% printout results

fprintf('e_golf: %0.3f',mean(e_golf) );
fprintf(' ± %0.3f \n', mean(total_error_golf) );

fprintf('e_pong: %0.3f',mean(e_pong) );
fprintf(' ± %0.3f \n', mean(total_error_pong) );
