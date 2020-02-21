close all;clear all;clc;
% This script generates the Weibull distribution parameters to fit to data
% in 
Resistances=[777873.70 530487.70 1255166.38 4354728.05 203746.85];
Percentiles=[0.5 0.3 0.7 0.98 0.02];
logR = log10(Resistances);
mu0 = logR(1);
s1 = (mu0 - logR(2)) ./ (log1p(-Percentiles(2)) - log(Percentiles(2)));
s2 = (mu0 - logR(3)) ./ (log1p(-Percentiles(3)) - log(Percentiles(3)));
s3 = (mu0 - logR(4)) ./ (log1p(-Percentiles(4)) - log(Percentiles(4)));
s4 = (mu0 - logR(5)) ./ (log1p(-Percentiles(5)) - log(Percentiles(5)));
s0 = 0.1942;

numPoints=1e5;
normDist=sort(random('Normal',0,1,numPoints,1));
logDist=sort(random('Logistic',mu0,s0,numPoints,1));
figure(1);
probplot(logDist,'noref');
xlim([4 8]);
ylim([-3 3]);