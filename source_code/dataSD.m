function [initSTD, initSTDsmooth, inittSTD, initMean, initMeanSmooth] = dataSD(window, dt, data, dataSmooth, excludedIndices, prevDataLength)

window = round(window/dt);
rows = floor(length(data)/window);
vert = repmat(window*(0:rows-1)',1,window);
horz = repmat((1:window),rows,1);
inds = vert + horz;
arrayData = data(inds);
arrayDataSmooth = dataSmooth(inds);

arrayExcl = logical(sum(excludedIndices(inds),2));
STDtemp = std(arrayData,1,2);
meanTemp = mean(arrayData,2);
initSTD = STDtemp(~arrayExcl);
initMean = meanTemp(~arrayExcl);
STDtemp = std(arrayDataSmooth,1,2);
meanTemp = mean(arrayDataSmooth,2);
initSTDsmooth = STDtemp(~arrayExcl);
initMeanSmooth = meanTemp(~arrayExcl);
tFile = (1:length(data))*dt - dt;
tFile = prevDataLength*dt - dt + tFile(inds(:,end));
inittSTD = tFile(~arrayExcl);