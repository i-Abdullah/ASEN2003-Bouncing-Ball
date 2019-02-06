%% info


%% get data

clear
clc
close all

% read data
Data = xlsread('ExperimentalData.xlsx');

Trial = Data(:,1); % trial number
Bounce1_Time = Data(:,2); % time between the first two bounces
Bounce2_Time = Data(:,3);
TotalTime = Data(:,4);
Error_time_values = Data(:,5);
Height_firstbounce = Data(:,7); % inches, 
Height_secondboune = Data(:,6); % inches, 


h0_inches = 36 ; %inches.
h0_SI = 0.9144 ; % meters
g = 9.81 ; % gravity
g = 386.09 ; % gravity in inches/s^2


%% error analysis : all error sources

h0_error = 1/4 ; % the same for all trials
TotalTime_error = 1/30 ; % changes per trial
Height_error = 1/4 ; % 1/14 inche is our error from height, eyeballed from the video
Time_error_each_bounce = 1/60; %using the same thing

% alternative method for error:

Height_firstbounce_std = std(Height_firstbounce);
Height_secondbounce_std = std(Height_secondboune);
Error_time_values_std = std(Error_time_values)+0.5; %std + 0.5 time to reaction?
Bounce1_Time_std = std(Bounce1_Time);
Bounce2_Time_std = std(Bounce2_Time);
%% error analysis : functions


syms h0 ts hn tn tn_1 hn_1 tn n
% hn_1 tn_1 is previous height, and time

e_stop_error(h0,ts) = (ts - sqrt((2*h0)/g))/(ts + sqrt((2*h0)/g));
e_height_error(hn,h0,n) = (hn/h0)^(1/(2*n)) ;
e_bonces_error(tn,tn_1) = ( tn / tn_1 ) ; 


%% e : time to stop

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% error analysis : derivatives: e_Stop

%partial of e_stop with respect to h0;
Partial_stop_h0 = diff(e_stop_error,h0);
Partial_stop_ts = diff(e_stop_error,ts);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% do the math!

for i=1:length(Trial)
    
e_stop(i) = (TotalTime(i) - sqrt((2*h0_inches)/g))/(TotalTime(i) + sqrt((2*h0_inches)/g));
Error_stop(i) = double(sqrt ( ((Partial_stop_h0(h0_inches,TotalTime(i)) * h0_error ))^2 + ((Partial_stop_ts(h0_inches,TotalTime(i)) * Error_time_values_std)))^2 );

end

%% e : time of bounce

% error analysis : derivatives: e_bounces


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%partial of e_bouncs with respect to tn, tn_1;

Partial_bounces_tn = diff(e_bonces_error,tn);
Partial_bounces_tn_1 = diff(e_bonces_error,tn_1);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


for i=1:length(Trial)
    
e_bounces(i) = Bounce2_Time(i) / Bounce1_Time(i) ;
Error_bouncs(i) = double(sqrt ( ((Partial_bounces_tn(Bounce2_Time(i),Bounce1_Time(i)) * Bounce2_Time_std )).^2 + ((Partial_bounces_tn_1(Bounce2_Time(i),Bounce1_Time(i)) * Bounce1_Time_std)).^2) );

end

%% e : height of bounce

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% error analysis : derivatives: e_height

%partial of e_height with respect to tn, tn_1;
Partial_height_hn = diff(e_height_error,hn);
Partial_height_h0 = diff(e_height_error,h0);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for i=1:length(Trial)
    
    % the equation says * n but we only using the first bounce of each
    % trial so n = 1 ;
    
e_height(i) = ( Height_firstbounce(i) / h0_inches ) ^ ( 1 / ( 2 )) ;
e_height_2(i) = ( Height_secondboune(i) / Height_firstbounce(i) ) ^ ( 1 / ( 4 )) ;

Error_height(i) = double(sqrt ( ((Partial_height_hn(Height_firstbounce(i),h0_inches,1) * Height_error ))^2 + ((Partial_height_h0(Height_firstbounce(i),h0_inches,1) * h0_error))^2 ));
Error_height_2(i) = double(sqrt ( ((Partial_height_hn(Height_secondboune(i), Height_firstbounce(i),1) * Height_secondbounce_std ))^2 + ((Partial_height_h0(Height_secondboune(i), Height_firstbounce(i),1) * Height_firstbounce_std))^2 ));


end


%% plot error


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure(1)

errorbar([1:length(Trial)],e_height_2,Error_height_2,'-*')
hold on
errorbar([1:length(Trial)],e_bounces,Error_bouncs,'-o')
hold on
errorbar([1:length(Trial)],e_stop,Error_stop,'-^')
grid minor
ylabel('Coefficient of restitution (unitless)')
xlabel('Trial')
title('Coefficient of restitution and error')
legend('e_h_e_i_g_h_t','e_b_o_u_n_c_e_s','e_s_t_o_p','Location','SouthEast')


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

