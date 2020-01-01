function [diffusionParam] = constantDiffusion(xdata, timeStamp)
%CONSTANTDIFFUSION Summary of this function goes here
%   Detailed explanation goes here
scaleVal = 1e15;
diffusionParam = scaleVal .* ones(size(xdata));
end
