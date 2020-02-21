close all;clear all;clc;
% This script generates the Weibull distribution parameters to fit to data
% in 
Resistances=[777873.70 530487.70 1255166.38 4354728.05 203746.85];
Percentiles=[0.5 0.3 0.7 0.98 0.02];
logR = log10(Resistances);
k1 = log(log1p(-Percentiles(3)) ./ log1p(-Percentiles(1))) ./ (logR(3) - logR(1));
k2 = log(log1p(-Percentiles(2)) ./ log1p(-Percentiles(1))) ./ (logR(2) - logR(1));
k3 = log(log1p(-Percentiles(3)) ./ log1p(-Percentiles(2))) ./ (logR(3) - logR(2));
k4 = log(log1p(-Percentiles(4)) ./ log1p(-Percentiles(5))) ./ (logR(4) - logR(5));
k5 = log(log1p(-Percentiles(4)) ./ log1p(-Percentiles(1))) ./ (logR(4) - logR(1));
k6 = log(log1p(-Percentiles(5)) ./ log1p(-Percentiles(1))) ./ (logR(5) - logR(1));

%k_res = mean([k2 k3 k4]);
k_res = 3.0;
lamb1 = logR(1) ./ ((-log1p(-Percentiles(1))).^(1 ./ k_res));
lamb2 = logR(2) ./ ((-log1p(-Percentiles(2))).^(1 ./ k_res));
lamb3 = logR(3) ./ ((-log1p(-Percentiles(3))).^(1 ./ k_res));
lamb4 = logR(4) ./ ((-log1p(-Percentiles(4))).^(1 ./ k_res));
lamb5 = logR(5) ./ ((-log1p(-Percentiles(5))).^(1 ./ k_res));
lamb_equi = mean([lamb1 lamb2 lamb3 lamb4]);

numPoints=1e4;
normDist = sort(random('Normal',0,1,numPoints,1));
wblDist = sort(random('Weibull',7.0,3.5,numPoints,1));

ind_low=1;
ind_high=1;
for idx=1:length(wblDist)
    if (wblDist(idx) > 5)
        if (ind_low == 1)
            ind_low = idx;
        end
    end
    if (wblDist(idx) < 7)
        ind_high = ind_high + 1;
    end
end
wblArr = wblDist(ind_low:ind_high);
resMap = 10.^(wblArr);
figure(1)
ax1=axis();
ah1=probplot(wblArr,'noref');
xlim([4 8]);
ylim([-3 3]);