function [ H ] = HIST9(sample,opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if opt == true
X = abs(sample);
else
X = sample;
end
    
[H,~]=histcounts(X,linspace(0,1,10));


end

