clear 
clc
close all

load target_output_FDDL
y1 = output;
t1 = target;

load target_output_HC
y2 = output(:,2);
t2 = target(:,2);

I = eye(2);

figure
plotconfusion(I(t1,:)',I(y1,:)')

figure
plotconfusion(I(t2,:)',I(y2,:)')

