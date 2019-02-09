%% info

% this script is to re-estimate the Coefficient of restitution for the bong ball
% and also a tennis ball using an improved method, which depends on using audio waves to 
% estimate the time between bounces and total time it took the ball to stop


%% read data

clear
clc
close all


Data = xlsread('Data/ImprovedMethodData.xlsx');

