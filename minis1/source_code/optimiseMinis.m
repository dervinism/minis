function [distParameters, fitness] = optimiseMinis(wd, targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters,...
    simulationParameters, optimisationParameters, classificationParameters, filtering, firstCost, parallelCores)
% optimiseMinis is an automated search algorithm for finding the optimal
% shape of the distribution for generating miniature excitatory
% postsynaptic potentials (mEPSPs or, simply, minis). The user can choose
% between finding the most optimal uni-dimensional (amplitudes only or
% amplitudes and 10-90% rise times simultaneously) or two-dimensional
% (amplitudes combined with 10-90% rise times) distributions and their
% actual shapes.
%

clear fitnessMinis
rng('shuffle');



%% Process the target file:
% Load the target file:
targetProperties = loadABF(targetFilename);
dt = targetProperties.dt;

% Determine filtering mode:
filtT.state = filtering.state;
filtT.nSweeps = targetProperties.hd.lActualEpisodes;
filtT.excludedTimes = targetExcludedTimes;
if isfield(filtering, 'filtfs')
  filtT.filtfs = filtering.filtfs;
end

% Determine excluded times for the target file:
sweepDuration = (length(targetProperties.sweep)*dt - dt)/targetProperties.hd.lActualEpisodes;
targetExcludedTimes = calcExcludedTimes(sweepDuration, targetProperties.hd.lActualEpisodes, 1000*targetExcludedTimes.startPulse,...
    1000*targetExcludedTimes.endPulse, 1000*targetExcludedTimes.startGlitch, 1000*targetExcludedTimes.endGlitch, dt);



%% Process the noise file:
% Load the noise file:
noiseProperties = loadABF(noiseFilename);
dt = noiseProperties.dt;

% Determine filtering mode:
filtN.state = filtering.state;
filtN.nSweeps = noiseProperties.hd.lActualEpisodes;
filtN.excludedTimes = noiseExcludedTimes;
if isfield(filtering, 'filtfs')
  filtN.filtfs = filtering.filtfs;
end

% Determine excluded times for a noise file:
sweepDuration = (length(noiseProperties.sweep)*dt - dt)/noiseProperties.hd.lActualEpisodes;
noiseExcludedTimes = calcExcludedTimes(sweepDuration, noiseProperties.hd.lActualEpisodes, 1000*noiseExcludedTimes.startPulse,...
    1000*noiseExcludedTimes.endPulse, 1000*noiseExcludedTimes.startGlitch, 1000*noiseExcludedTimes.endGlitch, dt);



%% Initialise the optimisation parameters
detectionParameters.sampleInterval = dt;
detectionParameters.smoothWindow = round(detectionParameters.smoothWindow/dt);

distributionType = optimisationParameters.distType;
baseline = optimisationParameters.distBaseline;

noGenerations = optimisationParameters.options.nGenerations;
draw = optimisationParameters.options.figureDisplay;
fullParallel = optimisationParameters.options.fullParallel;
tauRange = optimisationParameters.options.tauRange;
cluster = optimisationParameters.options.cluster;
clusterProfile = optimisationParameters.options.clusterProfile;
cliff = optimisationParameters.options.cliff;



%% Detect events in a noise-alone file:
if strcmpi(filtN.state, 'on')
    [noiseProperties.sweep, ~, f2] = filterMinis(noiseProperties.sweep, noiseProperties.dt, filtN, true);
    close(f2);
end

SD = [0 0];
options.SD = [0 1 0 0];
if strcmpi(baseline, 'Subtracted') || ~isempty(optimisationParameters.SDlobound) || ~isempty(optimisationParameters.SDupbound)
    filtering.state = 'off';
    waveformN.estimate = false;
    [noiseEvents, filtNoise] = detectMinis(noiseProperties.sweep, noiseExcludedTimes, detectionParameters, filtering, waveformN, parallelCores, options);
    if strcmpi(baseline, 'Subtracted')
        noiseProperties.sweep = filtNoise;
        [~, ~, noiseDistribution] = classifyMinis(noiseEvents(:,4), noiseEvents(:,12), classificationParameters);
        noiseDataLength = length(noiseProperties.sweep) - length(noiseProperties.excludedTimes);
    end
    if ~isempty(optimisationParameters.SDlobound) || ~isempty(optimisationParameters.SDupbound)
        SD(2) = noiseEvents(1,22);
    end
end



%% Detect events in a target file:
if tauRange && isempty(simulationParameters.tau_PSPm)
    waveformT.estimate = 2;
    waveformT.tau_m = simulationParameters.tau_m;
    waveformT.riseTimeArrayExt = classificationParameters.riseTimeArrayExt;
else
    waveformT.estimate = false;
end

if strcmpi(filtT.state, 'on')
    [targetProperties.sweep, ~, f1] = filterMinis(targetProperties.sweep, targetProperties.dt, filtT, true);
    close(f1);
end

filtering.state = 'spectrum';
filtering.nSweeps = filtT.nSweeps;
filtering.excludedTimes = filtT.excludedTimes;
filtT.excludedTimes = targetExcludedTimes;
[targetEvents, ~, tSpectrum, waveformT] = detectMinis(targetProperties.sweep, targetExcludedTimes, detectionParameters, filtering, waveformT,...
    parallelCores, options);
tSpectrum = tSpectrum(2,:);
targetDataLength = length(targetProperties.sweep) - length(targetExcludedTimes);
[targetEvents1D, targetEvents1D_RT, targetEvents2D] = classifyMinis(targetEvents(:,4), targetEvents(:,12), classificationParameters);

if strcmpi(baseline, 'Subtracted')
    baselineDistribution = (noiseDataLength/targetDataLength)*targetEvents2D - noiseDistribution;
    baselineDistribution = [zeros(1, size(baselineDistribution,2)); baselineDistribution];
    baselineDistribution = [zeros(size(baselineDistribution,1), 1) baselineDistribution];
    baselineDistribution(baselineDistribution < 0) = 0;
    if strcmpi(baseline, 'Subtracted')
        noiseProperties.baseline = baselineDistribution;
    else
        noiseProperties.baseline = length(classificationParameters.amplitudeArray);
    end
else
    noiseProperties.baseline = length(classificationParameters.amplitudeArray);
end

if ~isempty(optimisationParameters.SDlobound) || ~isempty(optimisationParameters.SDupbound)
    SD(1) = targetEvents(1,22);
end



%% Determine optimisation parameter bounds and the initial range:
if strcmpi(distributionType, 'Normal') || strcmpi(distributionType, 'Gaussian')
    nVars = 6;
elseif strcmpi(distributionType, 'Bimodal Normal') || strcmpi(distributionType, 'Bimodal Gaussian')
    nVars = 12;
elseif strcmpi(distributionType, 'Trimodal Normal') || strcmpi(distributionType, 'Trimodal Gaussian')
    nVars = 18;
elseif strcmpi(distributionType, 'Skew-normal') || strcmpi(distributionType, 'Log Normal')
    nVars = 8;
elseif strcmpi(distributionType, 'Bimodal Skew-normal') || strcmpi(distributionType, 'Bimodal Log Normal')
    nVars = 16;
elseif strcmpi(distributionType, 'Trimodal Skew-normal') || strcmpi(distributionType, 'Trimodal Log Normal')...
        || strcmpi(distributionType, 'Quadrimodal Normal') || strcmpi(distributionType, 'Quadrimodal Gaussian')
    nVars = 24;
end

% bounds on minis distribution parameters
loBound = optimisationParameters.options.bounds(1,:);
upBound = optimisationParameters.options.bounds(2,:);
loBound(nVars+1:end) = [];
upBound(nVars+1:end) = [];

% bounds on tau_sy1 (synaptic rise constant)
nVars = nVars + 1;
loBound = [loBound .15];
upBound = [upBound 3];

% bounds on tau_m
if tauRange
    nVars = nVars + 1;
    loBound = [loBound 0.9*simulationParameters.tau_m];
    if ~isempty(waveformT)
        upBound = [upBound 1.1*waveformT.parameters.tau_m];
    else
        upBound = [upBound simulationParameters.tau_PSPm];
    end
end

initRange = [loBound; upBound];



%% Choose the population size (the number of individuals generated) in GA:
% Aim to have at least 40.
populationSize1 = max(nVars*10, 40);



%% Create the data folder:
if ~isdeployed
    cd(wd);
end
dataDir = 'data';
if ~exist([wd filesep dataDir],'dir')
    mkdir([wd filesep dataDir]);
end



%% Initialise a customised output function for GA:
outPath = strcat(wd, '/data');
outputFcn = @(options, state, flag)customOutput(options, state, flag, outPath);



%% Load data for additional visualisation:
if ~isfield(optimisationParameters.options, 'optimisationData') || strcmpi(optimisationParameters.options.optimisationData, '')
    if isfield(optimisationParameters.options, 'headless') && optimisationParameters.options.headless
        button = 'No';
    else
        button = questdlg('Do you want to supply data for optimisation?','Load File','Yes','No','Yes');
    end
else
    button = 'Yes';
end
if strcmpi(button, 'Yes')
    if ~isfield(optimisationParameters.options, 'optimisationData') || strcmpi(optimisationParameters.options.optimisationData, '')
        [dataFilename, dataPathname, filterIndex] = uigetfile({'*.mat','MAT files (*.mat)'},'Data file', 'data_for_optimisation.mat');
    else
      [dataPathname, dataFilename] = fileparts(optimisationParameters.options.optimisationData);
      filterIndex = true;
    end
    if filterIndex
        dataFilename = fullfile(dataPathname, dataFilename);
        load(dataFilename, 'optimData');
        
        % Estimate centiles and the range:
        filtfs(1) = filtT.nSweeps;
        filtfs(2) = filtN.nSweeps;
        [estimated, diffAmps, diffAmpsDev, diffAmpsNeg, diffAmpsNegDev, diffRTs, diffRTsDev, diffRTsNeg, diffRTsNegDev, diffTwoDs, diffTwoDsDev,...
            diffTwoDsNeg, diffTwoDsNegDev, minAmps, minAmpsDev, minAmpsNeg, minAmpsNegDev, minRTs, minRTsDev, minRTsNeg, minRTsNegDev, minTwoDs,...
            minTwoDsDev, minTwoDsNeg, minTwoDsNegDev, maxAmps, maxAmpsDev, maxAmpsNeg, maxAmpsNegDev, maxRTs, maxRTsDev, maxRTsNeg, maxRTsNegDev,...
            maxTwoDs, maxTwoDsDev, maxTwoDsNeg, maxTwoDsNegDev] = compareSweepsOptim(size(optimData.initAmps,1), filtfs, optimData.initAmps,...
            optimData.initAmpsNeg, optimData.initRTs, optimData.initRTsNeg, optimData.initTwoDs, optimData.initTwoDsNeg, optimData.prct, optimData.prctNeg);
        if estimated
            measuredUF = [diffAmps; diffAmpsDev; diffAmpsNeg; diffAmpsNegDev; diffRTs; diffRTsDev; diffRTsNeg; diffRTsNegDev; diffTwoDs;...
                diffTwoDsDev; diffTwoDsNeg; diffTwoDsNegDev; minAmps; minAmpsDev; minAmpsNeg; minAmpsNegDev; minRTs; minRTsDev; minRTsNeg;...
                minRTsNegDev; minTwoDs; minTwoDsDev; minTwoDsNeg; minTwoDsNegDev; maxAmps; maxAmpsDev; maxAmpsNeg; maxAmpsNegDev; maxRTs;...
                maxRTsDev; maxRTsNeg; maxRTsNegDev; maxTwoDs; maxTwoDsDev; maxTwoDsNeg; maxTwoDsNegDev];
            [bottom, ~, diffAmpsBottom, diffAmpsDevBottom, diffAmpsNegBottom, diffAmpsNegDevBottom, minAmpsBottom, minAmpsDevBottom, minAmpsNegBottom,...
                minAmpsNegDevBottom, maxAmpsBottom, maxAmpsDevBottom, maxAmpsNegBottom,...
                maxAmpsNegDevBottom] = compareSweepsAmpsOptim(size(optimData.initAmps,1), filtfs, optimData.initAmps, optimData.initAmpsNeg,...
                optimData.prct, optimData.prctNeg, .5);
            measuredUFBottom = [diffAmpsBottom, diffAmpsDevBottom, diffAmpsNegBottom, diffAmpsNegDevBottom, minAmpsBottom, minAmpsDevBottom,...
                minAmpsNegBottom, minAmpsNegDevBottom, maxAmpsBottom, maxAmpsDevBottom, maxAmpsNegBottom, maxAmpsNegDevBottom];
            [mid, ~, diffAmpsMid, diffAmpsDevMid, diffAmpsNegMid, diffAmpsNegDevMid, minAmpsMid, minAmpsDevMid, minAmpsNegMid, minAmpsNegDevMid,...
                maxAmpsMid, maxAmpsDevMid, maxAmpsNegMid, maxAmpsNegDevMid] = compareSweepsAmpsOptim(size(optimData.initAmps,1), filtfs,...
                optimData.initAmps, optimData.initAmpsNeg, optimData.prct, optimData.prctNeg, .9);
            measuredUFMid = [diffAmpsMid, diffAmpsDevMid, diffAmpsNegMid, diffAmpsNegDevMid, minAmpsMid, minAmpsDevMid, minAmpsNegMid, minAmpsNegDevMid,...
                maxAmpsMid, maxAmpsDevMid, maxAmpsNegMid, maxAmpsNegDevMid];
            [top, ~, diffAmpsTop, diffAmpsDevTop, diffAmpsNegTop, diffAmpsNegDevTop, minAmpsTop, minAmpsDevTop, minAmpsNegTop, minAmpsNegDevTop,...
                maxAmpsTop, maxAmpsDevTop, maxAmpsNegTop, maxAmpsNegDevTop] = compareSweepsAmpsOptim(size(optimData.initAmps,1), filtfs,...
                optimData.initAmps, optimData.initAmpsNeg, optimData.prct, optimData.prctNeg, .98);
            measuredUFTop = [diffAmpsTop, diffAmpsDevTop, diffAmpsNegTop, diffAmpsNegDevTop, minAmpsTop, minAmpsDevTop, minAmpsNegTop, minAmpsNegDevTop,...
                maxAmpsTop, maxAmpsDevTop, maxAmpsNegTop, maxAmpsNegDevTop];
            histoData.Amps = optimData.Amps(filtfs(2));
            histoData.Amps = histoData.Amps{1};
            histoData.RTs = optimData.RTs(filtfs(2));
            histoData.RTs = histoData.RTs{1};
            TwoDsTemp = optimData.TwoDs(filtfs(2));
            TwoDsTemp = TwoDsTemp{1};
            TwoDs = zeros(size(TwoDsTemp{1},1), size(TwoDsTemp{1},2), numel(TwoDsTemp));
            for iPlane = 1:numel(TwoDsTemp)
                TwoDs(:,:,iPlane) = TwoDsTemp{iPlane};
            end
            histoData.TwoDs = TwoDs;
            histoData.SDs = SD;
        else
            measuredUF = [];
            measuredUFBottom = [];
            bottom = [];
            measuredUFMid = [];
            mid = [];
            measuredUFTop = [];
            top = [];
            histoData = [];
        end
    else
        measuredUF = [];
        measuredUFBottom = [];
        bottom = [];
        measuredUFMid = [];
        mid = [];
        measuredUFTop = [];
        top = [];
        histoData = [];
    end
else
    measuredUF = [];
    measuredUFBottom = [];
    bottom = [];
    measuredUFMid = [];
    mid = [];
    measuredUFTop = [];
    top = [];
    histoData = [];
end



%% Determine the structure of the cost function:
% Determine the intervals of the top 10% and top 2% amplitude cost functions:
if isempty(mid)
    Ampcs = cumsum(sum(targetEvents2D))/sum(sum(targetEvents2D));
    bottom = find(Ampcs > .5, 1);
    mid = find(Ampcs > .9, 1);
    top = find(Ampcs > .98, 1);
end

% Determine the basis vector of the cost function:
if isempty(measuredUF)
    costBasis = [zeros(1,15) fliplr(0: 2*optimisationParameters.SAD :30*optimisationParameters.SAD)];
else
    costBasis = fliplr(0: 2*optimisationParameters.SAD :60*optimisationParameters.SAD);
end

% Determine the cost scaling and bound vectors:
if strcmpi(firstCost, 'TwoDs')
    boundVector = [optimisationParameters.SAD optimisationParameters.AmpSAD optimisationParameters.RTSAD optimisationParameters.MAD...
        optimisationParameters.AmpMAD optimisationParameters.RTMAD optimisationParameters.AmpBottomSAD .5*optimisationParameters.AmpBottomSAD...
        optimisationParameters.AmpBottomMAD optimisationParameters.AmpMidSAD .5*optimisationParameters.AmpMidSAD optimisationParameters.AmpMidMAD...
        optimisationParameters.AmpTopSAD .5*optimisationParameters.AmpTopSAD optimisationParameters.AmpTopMAD (SD(1)-SD(2))/2];
elseif strcmpi(firstCost, 'Amps')
    boundVector = [optimisationParameters.AmpSAD optimisationParameters.RTSAD optimisationParameters.SAD optimisationParameters.AmpMAD...
        optimisationParameters.RTMAD optimisationParameters.MAD optimisationParameters.AmpBottomSAD .5*optimisationParameters.AmpBottomSAD...
        optimisationParameters.AmpBottomMAD optimisationParameters.AmpMidSAD .5*optimisationParameters.AmpMidSAD optimisationParameters.AmpMidMAD...
        optimisationParameters.AmpTopSAD .5*optimisationParameters.AmpTopSAD optimisationParameters.AmpTopMAD (SD(1)-SD(2))/2];
elseif strcmpi(firstCost, 'RTs')
    boundVector = [optimisationParameters.RTSAD optimisationParameters.AmpSAD optimisationParameters.SAD optimisationParameters.RTMAD...
        optimisationParameters.AmpMAD optimisationParameters.MAD optimisationParameters.AmpBottomSAD .5*optimisationParameters.AmpBottomSAD...
        optimisationParameters.AmpBottomMAD optimisationParameters.AmpMidSAD .5*optimisationParameters.AmpMidSAD optimisationParameters.AmpMidMAD...
        optimisationParameters.AmpTopSAD .5*optimisationParameters.AmpTopSAD optimisationParameters.AmpTopMAD (SD(1)-SD(2))/2];
end
if isempty(measuredUF)
    boundVector = [zeros(1,15) boundVector];
    if strcmpi(firstCost, 'TwoDs')
        costScale = 1./(boundVector/boundVector(16));
    elseif strcmpi(firstCost, 'Amps') || strcmpi(firstCost, 'RTs')
        costScale = 1./(boundVector/boundVector(18));
    end
else
    boundVector = [boundVector(1:end-1) boundVector];
    if firstCost
        costScale = 1./(boundVector/boundVector(16));
    elseif strcmpi(firstCost, 'Amps') || strcmpi(firstCost, 'RTs')
        costScale = 1./(boundVector/boundVector(18));
    end
    boundVector(1:15) = zeros(1,15);
end
costScale(isinf(costScale)) = 0;
boundVector(end) = 0;

costFuncStruct = struct('bottom', bottom, 'mid', mid, 'top', top, 'costBasis', costBasis, 'costScale', costScale, 'boundVector', boundVector,...
    'firstCost', firstCost);



%% Resume an earlier optimisation:
if ~isfield(optimisationParameters.options, 'resumeOptimisation') || strcmpi(optimisationParameters.options.resumeOptimisation, '')
    if isfield(optimisationParameters.options, 'headless') && optimisationParameters.options.headless
        button = 'No';
    else
        button = questdlg('Do you want to resume an earlier optimisation?','Resume optimisation','Yes','No','No');
    end
else
    button = 'Yes';
end
if strcmpi(button, 'Yes')
    if ~isfield(optimisationParameters.options, 'resumeOptimisation') || strcmpi(optimisationParameters.options.resumeOptimisation, '')
        [popFilename, popPathname, filterIndex] = uigetfile({'*.mat','MAT files (*.mat)'},'Population file', 'last_population.mat');
    else
        [popPathname, popFilename] = fileparts(optimisationParameters.options.resumeOptimisation);
        filterIndex = true;
    end
    if filterIndex
        popFilename = fullfile(popPathname, popFilename);
        warning('off', 'MATLAB:load:variableNotFound');
        load(popFilename, 'distParameters', 'population', 'simulationParameters');
        warning('on', 'MATLAB:load:variableNotFound');
        if ~exist('population','var')
            initialPopulation = zeros(populationSize1,nVars);
            initialPopulation(1,:) = distParameters;
            base = distParameters - .25*(upBound-loBound);
            base = max([base; loBound]);
            ceiling = distParameters + .25*(upBound-loBound);
            ceiling = min([ceiling; upBound]);
            initialPopulation(2:end,:) = repmat(base,populationSize1-1,1) + repmat(ceiling-base,populationSize1-1,1)...
                .* rand(populationSize1-1,nVars);
        else
            if isfield(optimisationParameters.options, 'headless') && optimisationParameters.options.headless
                button = 'Yes';
            else
                button = questdlg('Do you want to reset the population diversity?','Reset Diversity','Yes','No','No');
            end
            if strcmpi(button, 'Yes')
                initialPopulation = zeros(populationSize1,nVars);
                initialPopulation(1,:) = distParameters;
                base = distParameters - .25*(upBound-loBound);
                base = max([base; loBound]);
                ceiling = distParameters + .25*(upBound-loBound);
                ceiling = min([ceiling; upBound]);
                initialPopulation(2:end,:) = repmat(base,populationSize1-1,1) + repmat(ceiling-base,populationSize1-1,1)...
                    .* rand(populationSize1-1,nVars);
            else
                initialPopulation = population; %#ok<*NODEF>
            end
        end
    end
else
    initialPopulation = repmat(loBound,populationSize1,1) + repmat(upBound-loBound,populationSize1,1) .* rand(populationSize1,nVars);
end










%% Set the stopping criteria for optimisation:
if ~isempty(optimisationParameters.SDlobound) || ~isempty(optimisationParameters.SDupbound)
    fitnessLimit = costFuncStruct.costBasis(end) + costFuncStruct.costScale(end)*costFuncStruct.costVector(end);
else
    fitnessLimit = costFuncStruct.costBasis(end-1) + costFuncStruct.costScale(end-1)*costFuncStruct.boundVector(end-1);
end



%% Define options and the problem space of GA:
filtering.state = filtT.state;
filtering.nTargetSweeps = filtT.nSweeps;
filtering.nNoiseSweeps = filtN.nSweeps;
filtering.targetExcludedTimes = filtT.excludedTimes;
filtering.noiseExcludedTimes = filtN.excludedTimes;
filtering = rmfield(filtering, 'excludedTimes');

if ~cluster
    cd(dataDir);
    save('initVar', 'dataDir', 'baseline', 'targetEvents1D', 'targetEvents1D_RT', 'targetEvents2D', 'targetDataLength', 'noiseProperties',...
        'noiseExcludedTimes', 'distributionType', 'detectionParameters', 'classificationParameters', 'optimisationParameters', 'simulationParameters',...
        'filtering', 'tSpectrum', 'tauRange', 'cliff', 'measuredUF', 'measuredUFBottom', 'measuredUFMid', 'measuredUFTop', 'histoData',...
        'costFuncStruct', 'filtN');
    cd ..
end

f = @(x0)fitnessMinis(x0, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength, noiseProperties, noiseExcludedTimes,...
    distributionType, detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, tauRange,...
    cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fullParallel, draw, parallelCores);

eliteProportion = 0.05;
if fullParallel
    if cluster
        if draw
            options = gaoptimset('PopInitRange', initRange, 'Generations', noGenerations, 'EliteCount', ceil(eliteProportion*populationSize1),...
                'InitialPopulation', initialPopulation, 'MutationFcn', {@mutationadaptfeasible},...
                'PlotFcns',{@gaplotbestf @gaplotscorediversity @gaplotdistance @gaplotselection}, 'populationSize', populationSize1,...
                'UseParallel', 'always', 'Vectorized', 'off', 'FitnessLimit', fitnessLimit, 'StallGenLimit', noGenerations); %#ok<*GAOPT>
        else
            options = gaoptimset('PopInitRange', initRange, 'Generations', noGenerations, 'EliteCount', ceil(eliteProportion*populationSize1),...
                'InitialPopulation', initialPopulation, 'MutationFcn', {@mutationadaptfeasible},...
                'populationSize', populationSize1, 'UseParallel', 'always', 'Vectorized', 'off', 'FitnessLimit', fitnessLimit,...
                'StallGenLimit', noGenerations); %#ok<*UNRCH>
        end
    else
        if draw
            options = gaoptimset('PopInitRange', initRange, 'Generations', noGenerations, 'EliteCount', ceil(eliteProportion*populationSize1),...
                'InitialPopulation', initialPopulation, 'MutationFcn', {@mutationadaptfeasible},...
                'PlotFcns',{@gaplotbestf @gaplotscorediversity @gaplotdistance @gaplotselection}, 'populationSize', populationSize1,...
                'UseParallel', 'always', 'Vectorized', 'off', 'FitnessLimit', fitnessLimit, 'StallGenLimit', noGenerations,...
                'OutputFcn', outputFcn);
        else
            options = gaoptimset('PopInitRange', initRange, 'Generations', noGenerations, 'EliteCount', ceil(eliteProportion*populationSize1),...
                'InitialPopulation', initialPopulation, 'MutationFcn', {@mutationadaptfeasible},...
                'populationSize', populationSize1, 'UseParallel', 'always', 'Vectorized', 'off', 'FitnessLimit', fitnessLimit,...
                'StallGenLimit', noGenerations, 'OutputFcn', outputFcn);
        end
    end
else
    if draw
        options = gaoptimset('PopInitRange', initRange, 'Generations', noGenerations, 'EliteCount', ceil(eliteProportion*populationSize1),...
            'InitialPopulation', initialPopulation, 'MutationFcn', {@mutationadaptfeasible},...
            'PlotFcns',{@gaplotbestf @gaplotscorediversity @gaplotdistance @gaplotselection}, 'populationSize', populationSize1,...
            'FitnessLimit', fitnessLimit, 'StallGenLimit', noGenerations, 'OutputFcn', outputFcn);
    else
        options = gaoptimset('PopInitRange', initRange, 'Generations', noGenerations, 'EliteCount', ceil(eliteProportion*populationSize1),...
            'InitialPopulation', initialPopulation, 'MutationFcn', {@mutationadaptfeasible}, 'populationSize', populationSize1,...
            'FitnessLimit', fitnessLimit, 'StallGenLimit', noGenerations, 'OutputFcn', outputFcn);
    end
end

problem = struct('fitnessfcn', f, 'nvars', nVars, 'lb', loBound, 'ub', upBound, 'options', options);



%% Start the optimisation:
disp('Genetic algorithm has started...');
if cluster
    c = parcluster(clusterProfile);
    parpool(c);
    cjob = createJob(c);
    createTask(cjob, @ga, 1, {problem});
    submit(cjob);
    wait(cjob);
    data = cjob.fetchOutputs;
    population = data.population;
else
    clear SD TwoDs TwoDsTemp V baseline boundVector bottom button classificationParameters cliff clusterProfile costBasis costFuncStruct costScale...
        dataFilename dataPathname detectionParameters diffAmps diffAmpsDev diffAmpsDevMid diffAmpsDevTop diffAmpsMid diffAmpsNeg diffAmpsNegDev...
        diffAmpsNegDevMid diffAmpsNegDevTop diffAmpsNegMid diffAmpsNegTop diffAmpsTop diffRTs diffRTsDev diffRTsNeg diffRTsNegDev diffTwoDs...
        diffTwoDsDev diffTwoDsNeg diffTwoDsNegDev distributionType crossoverFraction dt estimated f filtN filtT filterIndex filtering filtfs...
        histoData iPlane iVar initialPopulation iterCount loBound maxAmps maxAmpsDev maxAmpsDevMid maxAmpsDevTop maxAmpsMid maxAmpsNeg maxAmpsNegDev...
        maxAmpsNegDevMid maxAmpsNegDevTop maxAmpsNegMid maxAmpsNegTop maxAmpsTop maxRTs maxRTsDev maxRTsNeg maxRTsNegDev maxTwoDs maxTwoDsDev...
        maxTwoDsNeg maxTwoDsNegDev measuredUF measuredUFBottom measuredUFMid measuredUFTop mid minAmps minAmpsDev minAmpsDevMid minAmpsDevTop...
        minAmpsMid minAmpsNeg minAmpsNegDev minAmpsNegDevMid minAmpsNegDevTop minAmpsNegMid minAmpsNegTop minAmpsTop minRTs minRTsDev minRTsNeg...
        minRTsNegDev minTwoDs minTwoDsDev minTwoDsNeg minTwoDsNegDev nVars noGenerations noiseExcludedTimes noiseFilename noiseProperties optimData...
        optimisationParameters options outPath outputFcn populationSize1 tSpectrum sweepDuration targetDataLength targetEvents targetEvents1D...
        targetEvents1D_RT targetEvents2D targetExcludedTimes targetFilename targetProperties tauRange top upBound waveformT wd
    [~, fitness, ~, ~, population] = ga(problem);
end



%% Save the results of optimisation:
cd(dataDir);
if draw
    saveas(gcf, 'evolution.fig');
end
if ~(fullParallel || cluster) || ~draw
    proceed = false;
    fileOrder = 1;
    while ~proceed
        filename = findFile('name', fileOrder);
        if length(filename) > 4 && strcmpi(filename(end-3:end), '.mat')
            save(filename,'population','-append');
            proceed = true;
        else
            fileOrder = fileOrder + 1;
        end
    end
end
cd ..
fclose all;

if draw
    if fullParallel || cluster
        cd(dataDir);
        load('initVar'); %#ok<*LOAD>
        proceed = false;
        fileOrder = 1;
        while ~proceed
            filename = findFile('name', fileOrder);
            if length(filename) > 4 && strcmpi(filename(end-3:end), '.mat')
                load(filename);
                proceed = true;
            else
                fileOrder = fileOrder + 1;
            end
        end
        
        simulationParameters.tau_sy1 = tau_sy1;
        simulationParameters.tau_sy2 = 0.1*tau_sy1;
        simulationParameters.tau_m = tau_m;
        
        filtN.state = 'spectrum';
        if draw
            waveformS.estimate = 3;
            waveformS.tau_m = simulationParameters.tau_m;
            waveformS.riseTimeArrayExt = classificationParameters.riseTimeArrayExt;
        else
            waveformS.estimate = 0;
        end
        options.SD = [0 1 0 0];
        if fitness == fitnessTail
            V = noiseProperties.sweep + simV + simVTail;
        else
            V = noiseProperties.sweep + simV;
        end
        [detectForSD, ~, nSpectrum, waveform] = detectMinis(V, noiseExcludedTimes, detectionParameters, filtN, waveformS, parallelCores, options);
        SD = detectForSD(1,22);
        nSpectrum = nSpectrum(2,:);
        if draw
            saveas(waveform.F(1), strcat(filename(1:end-4), '_waveform.jpg'));
            saveas(waveform.F(1), strcat(filename(1:end-4), '_waveform.fig'));
            saveas(waveform.F(2), strcat(filename(1:end-4), '_90prctRTs.jpg'));
            saveas(waveform.F(2), strcat(filename(1:end-4), '_90prctRTs.fig'));
        end
        cd ..
        
        if exist('population','var')
            plotSaveFP(draw, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength, simulatedEvents2D,...
                simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail, noiseExcludedTimes, SD, distributionType,...
                detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, nSpectrum, freq,...
                tauRange, cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fitness, fitnessTail, evalVector,...
                evalVectorTail, bestEvalVector, bestEvalVectorTail, population);
        else
            plotSaveFP(draw, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength, simulatedEvents2D,...
                simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail, noiseExcludedTimes, SD, distributionType,...
                detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, nSpectrum, freq,...
                tauRange, cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fitness, fitnessTail, evalVector,...
                evalVectorTail, bestEvalVector, bestEvalVectorTail);
        end
    else
        cd(dataDir);
        filename = findFile('name', 1);
        dataProperties = loadABF(strcat(filename(1:end-4),'.abf'));
        filtN.state = 'off';
        waveformS.estimate = 3;
        waveformS.tau_m = simulationParameters.tau_m;
        waveformS.riseTimeArrayExt = classificationParameters.riseTimeArrayExt;
        [~, ~, ~, waveform] = detectMinis(dataProperties.sweep, noiseExcludedTimes, detectionParameters, filtN, waveformS, parallelCores);
        saveas(waveform.F(1), strcat(filename(1:end-4), '_waveform.jpg'));
        saveas(waveform.F(1), strcat(filename(1:end-4), '_waveform.fig'));
        saveas(waveform.F(2), strcat(filename(1:end-4), '_90prctRTs.jpg'));
        saveas(waveform.F(2), strcat(filename(1:end-4), '_90prctRTs.fig'));
        cd ..
    end
end
end