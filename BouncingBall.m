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
Error_time = Data(:,5);
Height = Data(:,6); % inches, 


h0_inches = 36 ; %inches.
h0_SI = 0.9144 ; % meters
g = 9.81 ; % gravity
g = 386.09 ; % gravity in inches/s^2
%% e : time to stop


for i=1:length(Trial)
    
e_stop(i) = (TotalTime(i) - sqrt((2*h0_inches)/g))/(TotalTime(i) + sqrt((2*h0_inches)/g));

end

%% e : time of bounce


for i=1:length(Trial)
    
e_bounces(i) = Bounce2_Time(i) / Bounce1_Time(i) ;

end

%% e : height of bounce

for i=1:length(Trial)
    
    % the equation says * n but we only using the first bounce of each
    % trial so n = 1 ;
    
e_height(i) = ( Height(i) / h0_inches ) ^ ( 1 / ( 2 )) ;

end


%% error analysis : all error sources

h0_error = 2/7 ; % the same for all trials
TotalTime_error = Error_time ; % changes per trial
Height = 1/14 ; % 1/14 inche is our error from height, eyeballed from the video
Time_error_each_bounce = Error_time; %using the same thing

%% error analysis : e_stop


syms h0 ts 

e_stop_error(h0,ts) = (ts - sqrt((2*h0)/g))/(ts + sqrt((2*h0)/g));



