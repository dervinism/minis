function [V, simV, shapes, tailMinis, simulatedEvents1D, simulatedEvents2D, simulatedEvents1D_RT, SD, freq, noiseSpectrum] = resimulateMinis(tailMinis, shapes,...
    lengthRatio, simulationParameters, detectionParameters, classificationParameters, baseline, distributionType, noiseProperties, V, excludedTimes,...
    smoothWindow, amplitudeArraySim, filtN, waveform, costFuncStruct, parallelCores)

tailMinis = round(tailMinis/lengthRatio);
[simV, V, ~, shapesTail] = simulateMinis(baseline, simulationParameters, distributionType, tailMinis, noiseProperties.dt, noiseProperties.baseline, V,...
    excludedTimes, smoothWindow, parallelCores, classificationParameters.amplitudeArray, amplitudeArraySim, classificationParameters.riseTimeArray);
tailMinis = tailMinis*lengthRatio;

if costFuncStruct.costScale(end)
    options.SD = [0 1 0 0];
    [simulatedEvents, ~, noiseSpectrum] = detectMinis(V, excludedTimes, detectionParameters, filtN, waveform, parallelCores, options);
    SD = simulatedEvents(1,22);
else
    [simulatedEvents, ~, noiseSpectrum] = detectMinis(V, excludedTimes, detectionParameters, filtN, waveform, parallelCores);
    SD = 0;
end
freq = noiseSpectrum(4,:);
noiseSpectrum = noiseSpectrum(2,:);

[simulatedEvents1D, simulatedEvents1D_RT, simulatedEvents2D] = classifyMinis(simulatedEvents(:,4), simulatedEvents(:,12), classificationParameters);
simulatedEvents1D = lengthRatio*simulatedEvents1D;
simulatedEvents2D = lengthRatio*simulatedEvents2D;
simulatedEvents1D_RT = lengthRatio*simulatedEvents1D_RT;

shapes = [shapes; shapesTail];