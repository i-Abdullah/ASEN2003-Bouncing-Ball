%% info:
% this script is to estimate the coefficient of restitution of a ball using
% an image trackign software.


%% import data:

clear
clc
close all;


Data1 = xlsread('Data/Trial_1_Data_Image_Tracking.xlsx');
Data7 = xlsread('Data/Trial_7_Data_Image_Tracking.xlsx');
Data8 = xlsread('Data/Trial_8_Data_Image_Tracking.xlsx');

Data1(1,:) = [];

Data1(:,2) = Data1(:,2)*39.3701; % convert from inches to meters
Data7(:,2) = Data7(:,2)*39.3701; % convert from inches to meters
Data8(:,2) = Data8(:,2)*39.3701; % convert from inches to meters

% correct outlayer in all data set by making them hit the ground at height
% 0

% Data8(2:end,2) = Data8((2:end),2)  ; 

Data8(:,2) = Data8(:,2) - min(Data8(:,2));
Data7(:,2) = Data7(:,2) - min(Data7(:,2));
Data1(:,2) = Data1(:,2) - min(Data1(:,2));

% zero the time

Data1(:,1) = Data1(:,1) - Data1(1,1) ;
Data7(:,1) = Data7(:,1) - Data7(1,1) ;
Data8(:,1) = Data8(:,1) - Data8(1,1) ;


% h0 is just the first value in each data

h0_trial_1 = Data1(1,2);
h0_trial_7 = Data7(1,2);
h0_trial_8 = Data8(1,2);


% height of second bounce, it happens after 0.4 seconds so we can start
% looking after that

% get time bigger than 0.4

time_index_1 = find(Data1(:,1)>0.4);
time_index_7 = find(Data7(:,1)>0.4);
time_index_8 = find(Data8(:,1)>0.6); %for this one it happens after 0.6 seconds


% hn is height of the second bounce.

hn_Data1 = max(Data1(time_index_1(1):end,2));
hn_Data7 = max(Data7(time_index_7(1):end,2));
hn_Data8 = max(Data8(time_index_8(1):end,2));

%% get e values:

e_height_1 = ( hn_Data1 / h0_trial_1 ) ^ ( 1 / ( 2 )) ;
e_height_7 = ( hn_Data7 / h0_trial_7 ) ^ ( 1 / ( 2 )) ;
e_height_8 = ( hn_Data8 / h0_trial_8 ) ^ ( 1 / ( 2 )) ;

%% error analysis : all error sources

% error from instrument
h0_error_insturment = 1/4 ; % the same for all trials
hn_error_insturment = 1/4 ; % 1/4 inche is our error from height, eyeballed from the video

% random error

% std of all the three h0's because in theory they should be the same
% this's the same concept for hn

h0_error_random = std([ h0_trial_1 h0_trial_7 h0_trial_8 ]);
hn_error_random = std([ hn_Data1 hn_Data7 hn_Data8 ]);

% total error:

h0_total_error = sqrt ( (h0_error_insturment).^2 + (h0_error_random).^2 );
hn_total_error = sqrt ( (hn_error_insturment).^2 + (hn_error_random).^2 );

% error:

syms hn h0

e_height_error(hn,h0) = (hn/h0)^(1/(2));
Partial_height_hn = diff(e_height_error,hn);
Partial_height_h0 = diff(e_height_error,h0);

total_error_Data1 = double(sqrt ( ((Partial_height_h0(h0_trial_1,hn_Data1)) * h0_total_error ).^2 + (((Partial_height_hn(h0_trial_1,hn_Data1)) * hn_total_error).^2) ));
total_error_Data7 = double(sqrt ( ((Partial_height_h0(h0_trial_7,hn_Data7)) * h0_total_error ).^2 + (((Partial_height_hn(h0_trial_7,hn_Data7)) * hn_total_error).^2) ));
total_error_Data8 = double(sqrt ( ((Partial_height_h0(h0_trial_8,hn_Data8)) * h0_total_error ).^2 + (((Partial_height_hn(h0_trial_8,hn_Data8)) * hn_total_error).^2) ));



%% plot results

figure(1)

plot(Data1(:,1),Data1(:,2),'-.*')
hold on
plot(Data7(:,1),Data7(:,2),'-.s')
hold on
plot(Data8(:,1),Data8(:,2),'-.o')
hold on
hold on
plot([0:0.1:1.3], zeros(1,length([0:0.1:1.3])),'LineWidth',1.3)
legend('Trial 1','Trial 8','Trial 8','Ground')
grid minor
xlabel('time (second)');
ylabel('Height (inches)');
ylim([-1 40])
xlim([0 1.05])
title('Trajectory of the ball from height h0')

figure(2)
errorbar(1,e_height_1,total_error_Data1,'s')
hold on
errorbar(2,e_height_7,total_error_Data7,'^')
hold on
errorbar(3,e_height_8,total_error_Data8,'o')
grid minor
title('Image tracking coefficient of restitution with error bars')
xticks([ 1 2 3 ])
xticklabels({'Trial 1','Trial 7','Trial 8'})
xlim([ 0.5 3.5])
ylabel('Coefficient of restitution (unitless)')


%% printout resutls

fprintf('e_height_1: %0.3f',e_height_1 );
fprintf(' ± %0.4f \n', total_error_Data1 );

fprintf('e_height_7: %0.3f',e_height_7 );
fprintf(' ± %0.4f \n', total_error_Data7 );

fprintf('e_height_8: %0.3f',e_height_8 );
fprintf(' ± %0.4f \n', total_error_Data8 );
