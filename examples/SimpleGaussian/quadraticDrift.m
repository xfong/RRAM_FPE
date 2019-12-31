function [driftVector] = quadraticDrift(xdata, timeStamp)
%QUADRATICDRIFT Summary of this function goes here
%   Detailed explanation goes here
minPt = 50e-9;
scaleVal = 50;
driftVector = (-1.0 .* scaleVal) .* (xdata - minPt.*ones(size(xdata)));
end