%% info:
% this script is to estimate the coefficient of restitution of a ball using
% an image trackign software.


%% import data:

clear
clc
close all;


Data1 = xlsread('ImageTracking/Trial_1_Data_Image_Tracking.xlsx');
Data7 = xlsread('ImageTracking/Trial_7_Data_Image_Tracking.xlsx');
Data8 = xlsread('ImageTracking/Trial_8_Data_Image_Tracking.xlsx');


Data1(:,2) = Data1(:,2)*39.3701; % convert from inches to meters
Data7(:,2) = Data7(:,2)*39.3701; % convert from inches to meters
Data8(:,2) = Data8(:,2)*39.3701; % convert from inches to meters


%% plot results

plot(Data1(:,1),Data1(:,2))
hold on
plot(Data7(:,1),Data7(:,2))
hold on
plot(Data8(:,1),Data8(:,2))
legend('1','2','3')