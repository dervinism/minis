function SD = stdMinis(timeWindow, dt, data, excludedTimes, detectionParameters)

excludedIndices = round(excludedTimes/dt) + 1;
exclIndLogical = zeros(1, length(data));
exclIndLogical(excludedIndices) = 1;
circ = detectionParameters.SWstart + round(detectionParameters.BLduration/dt);
response = ones(1, 1 + 2*circ);
excludedIndices = logical(conv(exclIndLogical, response, 'same'));
excludedIndices(1 : (detectionParameters.SWstart + detectionParameters.refractoryPeriod)/dt - 1) = true;

window = round(timeWindow/dt);
rows = floor(length(data)/window);
vert = repmat(window*(0:rows-1)',1,window);
horz = repmat((1:window),rows,1);
inds = vert + horz;
arrayData = data(inds);

arrayExcl = logical(sum(excludedIndices(inds),2));
SD = std(arrayData,1,2);
SD = mean(SD(~arrayExcl));