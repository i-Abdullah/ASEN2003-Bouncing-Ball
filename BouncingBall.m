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
Error_time_values = Data(:,5)+0.05;
Height = Data(:,6); % inches, 


h0_inches = 36 ; %inches.
h0_SI = 0.9144 ; % meters
g = 9.81 ; % gravity
g = 386.09 ; % gravity in inches/s^2


%% error analysis : all error sources

h0_error = 2/7 ; % the same for all trials
TotalTime_error = Error_time_values ; % changes per trial
Height_error = 1/14 ; % 1/14 inche is our error from height, eyeballed from the video
Time_error_each_bounce = Error_time_values; %using the same thing

%% error analysis : functions


syms h0 ts hn tn tn_1 hn_1 tn n
% hn_1 tn_1 is previous height, and time

e_stop_error(h0,ts) = (ts - sqrt((2*h0)/g))/(ts + sqrt((2*h0)/g));

e_height_error(hn,h0,n) = (hn/h0)^(1/(2*n)) ;

e_bonces_error(tn,tn_1) = ( tn / tn_1 ) ; 


%% error analysis : derivatives: e_Stop


%partial of e_stop with respect to h0;
Partial_stop_h0 = diff(e_stop_error,h0);
Partial_stop_ts = diff(e_stop_error,ts);




%% e : time to stop


for i=1:length(Trial)
    
e_stop(i) = (TotalTime(i) - sqrt((2*h0_inches)/g))/(TotalTime(i) + sqrt((2*h0_inches)/g));
Error_stop(i) = double(sqrt ( ((Partial_stop_h0(h0_inches,TotalTime(i)) * h0_error ))^2 + ((Partial_stop_ts(h0_inches,TotalTime(i)) * Error_time_values(i))))^2 );

end

%% error analysis : derivatives: e_bounces


%partial of e_bouncs with respect to tn, tn_1;
Partial_bounces_tn = diff(e_bonces_error,tn);
Partial_bounces_tn_1 = diff(e_bonces_error,tn_1);




%% e : time of bounce


for i=1:length(Trial)
    
e_bounces(i) = Bounce2_Time(i) / Bounce1_Time(i) ;
Error_bouncs(i) = double(sqrt ( ((Partial_bounces_tn(Bounce2_Time(i),Bounce1_Time(i)) * Error_time_values(i) )).^2 + ((Partial_bounces_tn_1(Bounce2_Time(i),Bounce1_Time(i)) * Error_time_values(i))).^2) );

end

%% error analysis : derivatives: e_height


%partial of e_height with respect to tn, tn_1;
Partial_height_hn = diff(e_height_error,hn);
Partial_height_h0 = diff(e_height_error,h0);



%% e : height of bounce

for i=1:length(Trial)
    
    % the equation says * n but we only using the first bounce of each
    % trial so n = 1 ;
    
e_height(i) = ( Height(i) / h0_inches ) ^ ( 1 / ( 2 )) ;

Error_height(i) = double(sqrt ( ((Partial_height_hn(Height(i),h0_inches,i) * Height_error ))^2 + ((Partial_height_h0(Height(i),h0_inches,i) * h0_error))^2 ));


end

