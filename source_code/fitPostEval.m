function [fitness, fitnessNeg, fitness1D_LDNeg, fitness1_5D_LDNeg, fitness2D_LDNeg, fitness1D_LWDNeg, fitness1_5D_LWDNeg, fitness2D_LWDNeg, fitness1D_LRDNeg,...
        fitness1_5D_LRDNeg, fitness2D_LRDNeg, fitness1D_LSNeg, fitness1_5D_LSNeg, fitness2D_LSNeg, difference1DNeg, difference1D_RTNeg, difference2DNeg,...
        distances1DNeg, distances2DNeg, simulatedEvents1DNeg, simulatedEvents2DNeg, simulatedEvents1D_RTNeg, targetEvents1DNeg, targetEvents1D_RTNeg,...
        targetEvents2DNeg, thrFitness] = fitPostEval(fitness, fitnessFunc, threshold, fitLevel, fitAmp, fitRT, fitFact, dimension, V, noiseProperties,...
        targetEventsNeg, minis2D, meanMinisRT, medianMinisRT, maxDev, lengthRatio, searchParameters, classificationParameters, amplitudeArray,...
        amplitudeArraySim, riseTimeArray, waveform, parallelCores, currentNeg, thrFitness)

fitStep = 10*threshold.thrValue;

% Negative SAD:
if fitLevel(3) && isempty(currentNeg)
    [fitnessNeg, fitness1D_LDNeg, fitness1_5D_LDNeg, fitness2D_LDNeg, fitness1D_LWDNeg, fitness1_5D_LWDNeg, fitness2D_LWDNeg, fitness1D_LRDNeg,...
        fitness1_5D_LRDNeg, fitness2D_LRDNeg, fitness1D_LSNeg, fitness1_5D_LSNeg, fitness2D_LSNeg, difference1DNeg, difference1D_RTNeg, difference2DNeg,...
        distances1DNeg, distances2DNeg, simulatedEvents1DNeg, simulatedEvents2DNeg, simulatedEvents1D_RTNeg, targetEvents1DNeg, targetEvents1D_RTNeg,...
        targetEvents2DNeg] = negativeFitness(fitnessFunc, dimension, V, targetEventsNeg, noiseProperties, lengthRatio, searchParameters,...
        classificationParameters, waveform, parallelCores);
    fitness = fitness + fitnessNeg;
elseif ~isempty(currentNeg)
    simulatedEvents1DNeg = sum(currentNeg,1);
    simulatedEvents2DNeg = currentNeg;
    simulatedEvents1D_RTNeg = sum(currentNeg,2)';
    
    targetEvents1DNeg = sum(targetEventsNeg,1);
    targetEvents2DNeg = targetEventsNeg;
    targetEvents1D_RTNeg = sum(targetEventsNeg,2)';
    
    [fitnessNeg, fitness1D_LDNeg, fitness1_5D_LDNeg, fitness2D_LDNeg, fitness1D_LWDNeg, fitness1_5D_LWDNeg, fitness2D_LWDNeg, fitness1D_LRDNeg, fitness1_5D_LRDNeg,...
        fitness2D_LRDNeg, fitness1D_LSNeg, fitness1_5D_LSNeg, fitness2D_LSNeg, difference1DNeg, difference1D_RTNeg, difference2DNeg, distances1DNeg,...
        distances2DNeg] = fitnessEval(targetEvents1DNeg, simulatedEvents1DNeg, targetEvents1D_RTNeg, simulatedEvents1D_RTNeg, targetEvents2DNeg,...
        simulatedEvents2DNeg, fitnessFunc, dimension);
    fitness = fitness + fitnessNeg;
else
    fitnessNeg = [];
    fitness1D_LDNeg = [];
    fitness1_5D_LDNeg = [];
    fitness2D_LDNeg = [];
    fitness1D_LWDNeg = [];
    fitness1_5D_LWDNeg = [];
    fitness2D_LWDNeg = [];
    fitness1D_LRDNeg = [];
    fitness1_5D_LRDNeg = [];
    fitness2D_LRDNeg = [];
    fitness1D_LSNeg = [];
    fitness1_5D_LSNeg = [];
    fitness2D_LSNeg = [];
    difference1DNeg = [];
    difference1D_RTNeg = [];
    difference2DNeg = [];
    distances1DNeg = [];
    distances2DNeg = [];
    simulatedEvents1DNeg = [];
    simulatedEvents2DNeg = [];
    simulatedEvents1D_RTNeg = [];
    targetEvents1DNeg = [];
    targetEvents1D_RTNeg = [];
    targetEvents2DNeg = [];
end

% Positive SAD:
if fitness < threshold.thrValue
    sadFitness = fitFact;
else
    sadFitness = 1;
end
fitness = fitness*sadFitness + 2*fitStep;

% Maximum deviation of the positive SAD:
if sadFitness == fitFact && fitLevel(4) && maxDev < threshold.devValue
    devFitness = fitFact;
    fitness = fitness - 2*fitStep;
    fitness = fitness*devFitness + 2*fitStep;
else
    devFitness = 1;
end

% Runs test:
if sadFitness == fitFact && ((fitLevel(4) && devFitness == fitFact) || ~fitLevel(4)) && fitLevel(5) && fitAmp == 1 && fitRT == 1
    runsFit = fitFact;
    fitness = fitness - 2*fitStep;
    fitness = fitness*runsFit + 2*fitStep;
else
    runsFit = 1;
end

% Skewness:
skewness = meanMinisRT/medianMinisRT;
if fitLevel(1) && sadFitness == fitFact && (((fitLevel(4) && devFitness == fitFact) || ~fitLevel(4)) && ((fitLevel(5) && runsFit == fitFact) || ~fitLevel(5)))
    if skewness <= 1
        fitness = (skewness/riseTimeArray(end))*fitStep + fitStep;
        skew = true;
    else
        skew = false;
    end
elseif sadFitness == fitFact && (((fitLevel(4) && devFitness == fitFact) || ~fitLevel(4)) && ((fitLevel(5) && runsFit == fitFact) || ~fitLevel(5)))
    skew = false;
end

% Cliff:
if fitLevel(2) && exist('skew','var') && ~skew
    cliff = detectCliff(minis2D);
    if cliff
        fitness = cliff*fitStep;
    end
elseif exist('skew','var') && ~skew
    cliff = 'false';
end

% Threshold fitness:
if strcmpi(threshold.state,'on') && exist('skew','var') && exist('cliff','var') && ~skew && ~cliff && isempty(thrFitness)
    fitness = .25*threshold.thrValue;
    thrRange = amplitudeArray(1): amplitudeArray(2)-amplitudeArray(1) :threshold.upBound;
    weights = fliplr(thrRange)/thrRange(end);
    thrFitness = weights(amplitudeArray == amplitudeArraySim(1));
elseif isempty(thrFitness)
    thrFitness = 1;
end
fitness = fitness*thrFitness;
end