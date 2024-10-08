function [fitness, fitnessNeg, fitness1D_LDNeg, fitness1_5D_LDNeg, fitness2D_LDNeg, fitness1D_LWDNeg, fitness1_5D_LWDNeg, fitness2D_LWDNeg, fitness1D_LRDNeg,...
        fitness1_5D_LRDNeg, fitness2D_LRDNeg, fitness1D_LSNeg, fitness1_5D_LSNeg, fitness2D_LSNeg, difference1DNeg, difference1D_RTNeg, difference2DNeg,...
        distances1DNeg, distances2DNeg, simulatedEvents1DNeg, simulatedEvents2DNeg, simulatedEvents1D_RTNeg, targetEvents1DNeg, targetEvents1D_RTNeg,...
        targetEvents2DNeg] = fitPreEval(fitness, fitnessFunc, threshold, fitLevel, fitAmp, fitRT, fitFact, dimension, V, noiseProperties,...
        targetEventsNeg, minis2D, meanMinisRT, medianMinisRT, errAmp, errRT, maxDev, maxDevRT, maxDevAmp, lengthRatio, searchParameters,...
        classificationParameters, amplitudeArray, amplitudeArraySim, riseTimeArray, waveform, SD, parallelCores, currentNeg) %#ok<INUSL>

fitStep = .4*threshold.thrValue;

% SD:
if (~isempty(threshold.SDupbound) && SD > threshold.SDupbound) || (~isempty(threshold.SDlobound) && SD < threshold.SDlobound)
    fitness = fitness + 2*fitStep;
end

% Skewness:
skewness = meanMinisRT/medianMinisRT;
if fitLevel(1) && skewness <= 1
    fitness = fitness + skewness*fitStep;
end

% Cliff:
if fitLevel(2)
    cliff = detectCliff(minis2D);
    if cliff(1)
        fitness = fitness + cliff(1)*fitStep;
    end
    if cliff(2)
        fitness = fitness + cliff(2)*fitStep;
    end
end

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
    fitAdd = (errAmp - threshold.thrValueAmp)/threshold.thrValueAmp;
    if fitAdd < 0
        fitAdd = 0;
    end
    sadFitness = min([fitFact + fitAdd 1]);
    if errAmp < threshold.thrValueAmp
        fitAdd = (errRT - threshold.thrValueRT)/threshold.thrValueRT;
        if fitAdd < 0
            fitAdd = 0;
        end
        sadFitnessAmp = min([fitFact + fitAdd 1]);
        if errRT < threshold.thrValueRT
            sadFitnessRT = fitFact;
        else
            sadFitnessRT = 1;
        end
    else
        sadFitnessAmp = 1;
        sadFitnessRT = 1;
    end
else
    sadFitness = 1;
    sadFitnessAmp = 1;
    sadFitnessRT = 1;
end
fitness = fitness*sadFitness*sadFitnessRT*sadFitnessAmp;

% Maximum deviation of the positive SAD:
if sadFitness == fitFact && sadFitnessAmp == fitFact && sadFitnessRT == fitFact && fitLevel(4) && maxDev < threshold.devValue
    fitAdd = (maxDevAmp - threshold.devValueAmp)/threshold.devValueAmp;
    if fitAdd < 0
        fitAdd = 0;
    end
    devFitness = min([fitFact + fitAdd 1]);
    if maxDevAmp < threshold.devValueAmp
        fitAdd = (maxDevRT - threshold.devValueRT)/threshold.devValueRT;
        if fitAdd < 0
            fitAdd = 0;
        end
        devFitnessAmp = min([fitFact + fitAdd 1]);
        if maxDevRT < threshold.devValueRT
            devFitnessRT = fitFact;
        else
            devFitnessRT = 1;
        end
    else
        devFitnessAmp = 1;
        devFitnessRT = 1;
    end    
else
    devFitness = 1;
    devFitnessRT = 1;
    devFitnessAmp = 1;
end
fitness = fitness*devFitness*devFitnessRT*devFitnessAmp;

% Runs test:
if sadFitness == fitFact && sadFitnessAmp == fitFact && sadFitnessRT == fitFact &&...
        ((fitLevel(4) && devFitness == fitFact && devFitnessAmp == fitFact && devFitnessRT == fitFact) || ~fitLevel(4)) &&...
        fitLevel(5) && fitAmp == 1 && fitRT == 1
    runsFit = fitFact;
else
    runsFit = 1;
end
fitness = fitness*runsFit;
end