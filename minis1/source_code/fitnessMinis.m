function fitness = fitnessMinis(distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength,...
    noiseProperties, excludedTimes, distributionType, detectionParameters, classificationParameters, optimisationParameters, simulationParameters,...
    filtering, targetSpectrum, tauRange, cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fullParallel,...
    draw, parallelCores)

persistent bestFitness individual addStrTail f1 f2 f3 f4 f5 f6 f7 f8 f9 f10

if isempty(individual)
    individual = 1;
end
if ~fullParallel
    fprintf('\r\nGenerating individual no.: %g\n',individual);
end



%% Initialise input variables:
amplitudeArraySim = classificationParameters.amplitudeArray;
amplitudeArraySim(amplitudeArraySim < simulationParameters.loSimAmp) = [];

if tauRange
    simulationParameters.R_m = distParameters(end);
    simulationParameters.tau_sy1 = distParameters(end-1);
else
    simulationParameters.R_m = simulationParameters.tau_m;
    simulationParameters.tau_sy1 = distParameters(end);
end
simulationParameters.tau_sy2 = .1*simulationParameters.tau_sy1;



%% Simulate and classify minis:
smoothWindow = max([floor((amplitudeArraySim(1)/(classificationParameters.amplitudeArray(2) - classificationParameters.amplitudeArray(1)))/2) 3]);
if tauRange
    [simV, V, minis2D, shapes] = simulateMinis(baseline, simulationParameters, distributionType, distParameters(1:end-2), noiseProperties.dt,...
        noiseProperties.baseline, noiseProperties.sweep, excludedTimes, smoothWindow, parallelCores, classificationParameters.amplitudeArray,...
        amplitudeArraySim, classificationParameters.riseTimeArray);
else
    [simV, V, minis2D, shapes] = simulateMinis(baseline, simulationParameters, distributionType, distParameters(1:end-1), noiseProperties.dt,...
        noiseProperties.baseline, noiseProperties.sweep, excludedTimes, smoothWindow, parallelCores, classificationParameters.amplitudeArray,...
        amplitudeArraySim, classificationParameters.riseTimeArray);
end
minis1D = sum(minis2D,1);
minis1D_RT = sum(minis2D,2)';
noiseDataLength = length(V) - length(excludedTimes);



%% Detect and classify simulated events:
waveform.estimate = false;
filtN.state = 'spectrum';
filtN.nSweeps = filtering.nNoiseSweeps;
filtN.excludedTimes = filtering.noiseExcludedTimes;
detectionParametersSim = detectionParameters;

% Adapt the Maximum time to peak window to a particular rise time
% distribution:
minis1D_RTprc = cumsum(minis1D_RT)/sum(minis1D_RT);
minis1D_RT75prc = find(minis1D_RTprc > .75, 1);
if ~isempty(minis1D_RT75prc)
    minis1D_RT75prc = classificationParameters.riseTimeArray(minis1D_RT75prc)*(.75/minis1D_RTprc(minis1D_RT75prc));
    SWstart = 2*minis1D_RT75prc + detectionParameters.BLduration;
    detectionParametersSim.SWstart = max([round(SWstart) detectionParameters.SWstart]);
end

if costFuncStruct.costScale(end)
    options.SD = [0 1 0 0];
    [simulatedEvents, ~, noiseSpectrum] = detectMinis(V, excludedTimes, detectionParametersSim, filtN, waveform, parallelCores, options);
    SD = simulatedEvents(1,22);
else
    [simulatedEvents, ~, noiseSpectrum] = detectMinis(V, excludedTimes, detectionParametersSim, filtN, waveform, parallelCores);
    SD = 0;
end
freq = noiseSpectrum(4,:);
noiseSpectrum = noiseSpectrum(2,:);
[simulatedEvents1D, simulatedEvents1D_RT, simulatedEvents2D] = classifyMinis(simulatedEvents(:,4), simulatedEvents(:,12), classificationParameters);



%% Scale simulated minis-dependent distributions:
lengthRatio = targetDataLength/noiseDataLength;
minis1D = lengthRatio*minis1D;
minis2D = lengthRatio*minis2D;
minis1D_RT = lengthRatio*minis1D_RT;
simulatedEvents1D = lengthRatio*simulatedEvents1D;
simulatedEvents2D = lengthRatio*simulatedEvents2D;
simulatedEvents1D_RT = lengthRatio*simulatedEvents1D_RT;



%% Evaluate Fitness:
[fitness, evalVector, bestEvalVector, meanMinisRTs, medianMinisRTs, firstCost, tailMinis, constraintFitness, allCosts, cliffEval] = fitnessEval( ...
    targetEvents1D, simulatedEvents1D, targetEvents1D_RT, simulatedEvents1D_RT, targetEvents2D, simulatedEvents2D, SD, costFuncStruct, ...
    measuredUF, histoData, shapes, cliff, minis2D, optimisationParameters);



%% Re-evaluate fitness:
if (~isempty(costFuncStruct.costBasis) && firstCost > 9) || firstCost > 24 || (~isempty(addStrTail) && addStrTail)
    [VTail, simVTail, shapesTail, tailMinis, simulatedEvents1DTail, simEvents2DTail, simulatedEvents1D_RTTail, SDTail, freqTail,...
        noiseSpectrumTail] = resimulateMinis(tailMinis, shapes, lengthRatio, simulationParameters, detectionParameters, classificationParameters,...
        baseline, distributionType, noiseProperties, V, excludedTimes, smoothWindow, amplitudeArraySim, filtN, waveform, costFuncStruct, parallelCores);
    
    [fitnessTail, evalVectorTail, bestEvalVectorTail, meanMinisRTsTail, medianMinisRTsTail] = fitnessEval(targetEvents1D, simulatedEvents1DTail,...
        targetEvents1D_RT, simulatedEvents1D_RTTail, targetEvents2D, simEvents2DTail, SDTail, costFuncStruct, measuredUF, histoData, shapesTail,...
        cliff, minis2D, optimisationParameters);
    
    if (fitnessTail < fitness) || (~isempty(addStrTail) && addStrTail)
        addStrTail = true;
        fitness = fitnessTail;
        evalVector = evalVectorTail;
        bestEvalVector = bestEvalVectorTail;
        meanMinisRTs = meanMinisRTsTail;
        medianMinisRTs = medianMinisRTsTail;
        
        V = VTail;
        shapes = shapesTail;
        minis1D = minis1D + sum(tailMinis,1);
        minis2D = minis2D + tailMinis;
        minis1D_RT = minis1D_RT + sum(tailMinis,2)';
        simulatedEvents1D = simulatedEvents1DTail;
        simulatedEvents2D = simEvents2DTail;
        simulatedEvents1D_RT = simulatedEvents1D_RTTail;
        SD = SDTail;
        freq = freqTail;
        noiseSpectrum = noiseSpectrumTail;
    end
else
    fitnessTail = inf;
    simVTail = zeros(size(V));
    shapesTail = zeros(1,size(shapes,2));
    simulatedEvents1DTail = zeros(size(simulatedEvents1D));
    simEvents2DTail = zeros(size(simulatedEvents2D));
    simulatedEvents1D_RTTail = zeros(size(simulatedEvents1D_RT));
    evalVectorTail = zeros(size(evalVector));
    bestEvalVectorTail = zeros(size(bestEvalVector));
    tailMinis = zeros(size(minis2D));
end



%% Generating a name for the files:
cd(dataDir);
filename = num2str(floor(fitness));
if length(filename) < 6
    if length(filename) == 4
        addStr = '00';
    elseif length(filename) == 5
        addStr = '0';
    end
    filename = strcat(addStr, num2str(fitness));
else
    filename = num2str(fitness);
end
filename = strrep(filename, '.', '_');



%% Recording the characteristics of the best distributions:
if ~fullParallel && (individual == 1 || (individual > 1 && fitness < bestFitness))
    bestFitness = fitness;
    if ~SD
        SD = stdMinis(15, noiseProperties.dt, V, excludedTimes, detectionParameters);
    end
    meanMinisAmplitudes = mean(shapes(:,2));
    medianMinisAmplitudes = median(shapes(:,2));
    targetSweepDuration = noiseProperties.dt*targetDataLength/1000;         % seconds
    count = sum(sum(minis2D));
    minisFrequency = count/targetSweepDuration;                             % Hz
    
    if ~exist('log.txt','file')
        fileID = fopen('log.txt','w');
    else
        fileID = fopen('log.txt','a');
    end
    simFilename = strcat(filename,'.abf');
    simsimFilename = strcat(filename,'sim','.abf');
    
    
    
    % Print:
    if individual > 1
        fprintf(fileID, '\r\n');
    end
    fprintf(fileID, '\r\nIndividual number:                                        %g', individual);
    fprintf(fileID, '\r\nFitness:                                                  %g', fitness);
    fprintf(fileID, '\r\nFile name containing simulations superimposed on data:    %s', simFilename);
    fprintf(fileID, '\r\nFile name containing simulations only:                    %s', simsimFilename);
    fprintf(fileID, '\r\nFile name containing addStred minis tail:                    %s', strcat(simsimFilename, 'tail'));
    fprintf(fileID, '\r\n* * *');
    fprintf(fileID, '\r\nTarget-relative SAD:                                       %g', evalVector(16));
    fprintf(fileID, '\r\nTarget-relative Amplitude SAD:                             %g', evalVector(17));
    fprintf(fileID, '\r\nTarget-relative Rise time SAD:                             %g', evalVector(18));
    fprintf(fileID, '\r\nTarget-relative MAD:                                       %g', evalVector(19));
    fprintf(fileID, '\r\nTarget-relative Amplitude MAD:                             %g', evalVector(20));
    fprintf(fileID, '\r\nTarget-relative Rise time MAD:                             %g', evalVector(21));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 50%% SAD:                     %g', evalVector(22));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 50%% lower SAD:               %g', evalVector(23));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 50%% MAD:                     %g', evalVector(24));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 10%% SAD:                     %g', evalVector(25));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 10%% lower SAD:               %g', evalVector(26));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 10%% MAD:                     %g', evalVector(27));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 2%% SAD:                      %g', evalVector(28));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 2%% lower SAD:                %g', evalVector(29));
    fprintf(fileID, '\r\nTarget-relative Amplitude top 2%% MAD:                      %g', evalVector(30));
    fprintf(fileID, '\r\n* * *');
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        fprintf(fileID, '\r\nGroup-relative (outlying) SAD:                             %g', evalVector(1));
        fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude SAD:                   %g', evalVector(2));
        fprintf(fileID, '\r\nGroup-relative (outlying) Rise time SAD:                   %g', evalVector(3));
        fprintf(fileID, '\r\nGroup-relative (outlying) MAD:                             %g', evalVector(4));
        fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude MAD:                   %g', evalVector(5));
        fprintf(fileID, '\r\nGroup-relative (outlying) Rise time MAD:                   %g', evalVector(6));
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        fprintf(fileID, '\r\nGroup-relative (outlying) SAD:                             %g', evalVector(3));
        fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude SAD:                   %g', evalVector(1));
        fprintf(fileID, '\r\nGroup-relative (outlying) Rise time SAD:                   %g', evalVector(2));
        fprintf(fileID, '\r\nGroup-relative (outlying) MAD:                             %g', evalVector(6));
        fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude MAD:                   %g', evalVector(4));
        fprintf(fileID, '\r\nGroup-relative (outlying) Rise time MAD:                   %g', evalVector(5));
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        fprintf(fileID, '\r\nGroup-relative (outlying) SAD:                             %g', evalVector(3));
        fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude SAD:                   %g', evalVector(2));
        fprintf(fileID, '\r\nGroup-relative (outlying) Rise time SAD:                   %g', evalVector(1));
        fprintf(fileID, '\r\nGroup-relative (outlying) MAD:                             %g', evalVector(6));
        fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude MAD:                   %g', evalVector(5));
        fprintf(fileID, '\r\nGroup-relative (outlying) Rise time MAD:                   %g', evalVector(4));
    end
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 50%% SAD:           %g', evalVector(7));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 50%% lower SAD:     %g', evalVector(8));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 50%% MAD:           %g', evalVector(9));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 10%% SAD:           %g', evalVector(10));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 10%% lower SAD:     %g', evalVector(11));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 10%% MAD:           %g', evalVector(12));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 2%% SAD:            %g', evalVector(13));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 2%% lower SAD:      %g', evalVector(14));
    fprintf(fileID, '\r\nGroup-relative (outlying) Amplitude top 2%% MAD:            %g', evalVector(15));
    fprintf(fileID, '\r\n* * *');
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        fprintf(fileID, '\r\nBest possible target-relative SAD:                         %g', bestEvalVector(1));
        fprintf(fileID, '\r\nBest possible target-relative Amplitude SAD:               %g', bestEvalVector(2));
        fprintf(fileID, '\r\nBest possible target-relative Rise time SAD:               %g', bestEvalVector(3));
        fprintf(fileID, '\r\nBest possible target-relative MAD:                         %g', bestEvalVector(4));
        fprintf(fileID, '\r\nBest possible target-relative Amplitude MAD:               %g', bestEvalVector(5));
        fprintf(fileID, '\r\nBest possible target-relative Rise time MAD:               %g', bestEvalVector(6));
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        fprintf(fileID, '\r\nBest possible target-relative SAD:                         %g', bestEvalVector(3));
        fprintf(fileID, '\r\nBest possible target-relative Amplitude SAD:               %g', bestEvalVector(1));
        fprintf(fileID, '\r\nBest possible target-relative Rise time SAD:               %g', bestEvalVector(2));
        fprintf(fileID, '\r\nBest possible target-relative MAD:                         %g', bestEvalVector(6));
        fprintf(fileID, '\r\nBest possible target-relative Amplitude MAD:               %g', bestEvalVector(4));
        fprintf(fileID, '\r\nBest possible target-relative Rise time MAD:               %g', bestEvalVector(5));
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        fprintf(fileID, '\r\nBest possible target-relative SAD:                         %g', bestEvalVector(3));
        fprintf(fileID, '\r\nBest possible target-relative Amplitude SAD:               %g', bestEvalVector(2));
        fprintf(fileID, '\r\nBest possible target-relative Rise time SAD:               %g', bestEvalVector(1));
        fprintf(fileID, '\r\nBest possible target-relative MAD:                         %g', bestEvalVector(6));
        fprintf(fileID, '\r\nBest possible target-relative Amplitude MAD:               %g', bestEvalVector(5));
        fprintf(fileID, '\r\nBest possible target-relative Rise time MAD:               %g', bestEvalVector(4));
    end
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 50%% SAD:       %g', bestEvalVector(7));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 50%% lower SAD: %g', bestEvalVector(8));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 50%% MAD:       %g', bestEvalVector(9));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 10%% SAD:       %g', bestEvalVector(10));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 10%% lower SAD: %g', bestEvalVector(11));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 10%% MAD:       %g', bestEvalVector(12));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 2%% SAD:        %g', bestEvalVector(13));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 2%% lower SAD:  %g', bestEvalVector(14));
    fprintf(fileID, '\r\nBest possible target-relative Amplitude top 2%% MAD:        %g', bestEvalVector(15));
    fprintf(fileID, '\r\n* * *');
    fprintf(fileID, '\r\nNumber of simulated minis:                                 %g', count);
    fprintf(fileID, '\r\nMinis frequency:                                           %g', minisFrequency);
    fprintf(fileID, '\r\nMean (amplitude, rise time):                               %g %g', meanMinisAmplitudes, meanMinisRTs);
    fprintf(fileID, '\r\nMedian (amplitude, rise time):                             %g %g', medianMinisAmplitudes, medianMinisRTs);
    fprintf(fileID, '\r\nStandard deviation (smoothed 15ms-window):                 %g', SD);
    fprintf(fileID, '\r\nMinis amplitude generation threshold:                      %g', amplitudeArraySim(1));
    fprintf(fileID, '\r\nTau_m:                                                     %g', simulationParameters.R_m);
    fprintf(fileID, '\r\nTau_sy1:                                                   %g', simulationParameters.tau_sy1);
    fprintf(fileID, '\r\nTau_sy2:                                                   %g', simulationParameters.tau_sy2);
    fprintf(fileID, '\r\nRise time interval:                                        %s', detectionParameters.RTinterval);
    if tauRange
        fprintf(fileID, '\r\nRanging taus:                                              Yes');
    else
        fprintf(fileID, '\r\nRanging taus:                                              No');
    end
    if cliff
        fprintf(fileID, '\r\nCliff constraint:                                          Yes');
    else
        fprintf(fileID, '\r\nCliff constraint:                                          No');
    end
    fprintf(fileID, '\r\n* * *');
    fprintf(fileID, '\r\nMinis source distribution type:                            %s', distributionType);
    if fitnessTail == fitness
        fprintf(fileID, '\r\nAmplitude distribution tail addStred:                         Yes');
    else
        fprintf(fileID, '\r\nAmplitude distribution tail addStred:                         No');
    end
    
    if strcmpi(distributionType, 'Normal') || strcmpi(distributionType, 'Bimodal Normal') || strcmpi(distributionType, 'Trimodal Normal')...
            || strcmpi(distributionType, 'Quadrimodal Normal')
        fprintf(fileID,'\r\nmu1_1:                          %g', distParameters(1));
        fprintf(fileID,'\r\nsigma1_1:                       %g', distParameters(2));
        fprintf(fileID,'\r\nscale_1:                        %g', distParameters(3));
        fprintf(fileID,'\r\nmu2_1:                          %g', distParameters(4));
        fprintf(fileID,'\r\nsigma2_1:                       %g', distParameters(5));
        fprintf(fileID,'\r\nrho_1:                          %g', distParameters(6));
        if strcmpi(distributionType, 'Bimodal Normal') || strcmpi(distributionType, 'Trimodal Normal') || strcmpi(distributionType, 'Quadrimodal Normal')
            fprintf(fileID,'\r\nmu1_2:                          %g', distParameters(7));
            fprintf(fileID,'\r\nsigma1_2:                       %g', distParameters(8));
            fprintf(fileID,'\r\nscale_2:                        %g', distParameters(9));
            fprintf(fileID,'\r\nmu2_2:                          %g', distParameters(10));
            fprintf(fileID,'\r\nsigma2_2:                       %g', distParameters(11));
            fprintf(fileID,'\r\nrho_2:                          %g', distParameters(12));
            if strcmpi(distributionType, 'Trimodal Normal') || strcmpi(distributionType, 'Quadrimodal Normal')
                fprintf(fileID,'\r\nmu1_3:                          %g', distParameters(13));
                fprintf(fileID,'\r\nsigma1_3:                       %g', distParameters(14));
                fprintf(fileID,'\r\nscale_3:                        %g', distParameters(15));
                fprintf(fileID,'\r\nmu2_3:                          %g', distParameters(16));
                fprintf(fileID,'\r\nsigma2_3:                       %g', distParameters(17));
                fprintf(fileID,'\r\nrho_3:                          %g', distParameters(18));
                if strcmpi(distributionType, 'Quadrimodal Normal')
                    fprintf(fileID,'\r\nmu1_4:                          %g', distParameters(19));
                    fprintf(fileID,'\r\nsigma1_4:                       %g', distParameters(20));
                    fprintf(fileID,'\r\nscale_4:                        %g', distParameters(21));
                    fprintf(fileID,'\r\nmu2_4:                          %g', distParameters(22));
                    fprintf(fileID,'\r\nsigma2_4:                       %g', distParameters(23));
                    fprintf(fileID,'\r\nrho_4:                          %g', distParameters(24));
                end
            end
        end
    elseif strcmpi(distributionType, 'Log Normal') || strcmpi(distributionType, 'Bimodal Log Normal') || strcmpi(distributionType, 'Trimodal Log Normal')
        fprintf(fileID,'\r\nmu1_1:                          %g', distParameters(1));
        fprintf(fileID,'\r\nA_1:                            %g', distParameters(2));
        fprintf(fileID,'\r\nalpha_1:                        %g', distParameters(3));
        fprintf(fileID,'\r\nscale_1:                        %g', distParameters(4));
        fprintf(fileID,'\r\nmu2_1:                          %g', distParameters(5));
        fprintf(fileID,'\r\nB_1:                            %g', distParameters(6));
        fprintf(fileID,'\r\nbeta_1:                         %g', distParameters(7));
        fprintf(fileID,'\r\nrho_1:                          %g', distParameters(8));
        if strcmpi(distributionType, 'Bimodal Log Normal') || strcmpi(distributionType, 'Trimodal Log Normal')
            fprintf(fileID,'\r\nmu1_2:                          %g', distParameters(9));
            fprintf(fileID,'\r\nA_2:                            %g', distParameters(10));
            fprintf(fileID,'\r\nalpha_2:                        %g', distParameters(11));
            fprintf(fileID,'\r\nscale_2:                        %g', distParameters(12));
            fprintf(fileID,'\r\nmu2_2:                          %g', distParameters(13));
            fprintf(fileID,'\r\nB_2:                            %g', distParameters(14));
            fprintf(fileID,'\r\nbeta_2:                         %g', distParameters(15));
            fprintf(fileID,'\r\nrho_2:                          %g', distParameters(16));
            if strcmpi(distributionType, 'Trimodal Log Normal')
                fprintf(fileID,'\r\nmu1_3:                          %g', distParameters(17));
                fprintf(fileID,'\r\nA_3:                            %g', distParameters(18));
                fprintf(fileID,'\r\nalpha_3:                        %g', distParameters(19));
                fprintf(fileID,'\r\nscale_3:                        %g', distParameters(20));
                fprintf(fileID,'\r\nmu2_3:                          %g', distParameters(21));
                fprintf(fileID,'\r\nB_3:                            %g', distParameters(22));
                fprintf(fileID,'\r\nbeta_3:                         %g', distParameters(23));
                fprintf(fileID,'\r\nrho_3:                          %g', distParameters(24));
            end
        end
    elseif strcmpi(distributionType, 'Gaussian') || strcmpi(distributionType, 'Bimodal Gaussian') || strcmpi(distributionType, 'Trimodal Gaussian')...
            || strcmpi(distributionType, 'Quadrimodal Gaussian')
        fprintf(fileID,'\r\nmu1_1:                          %g', distParameters(1));
        fprintf(fileID,'\r\nsigma1_1:                       %g', distParameters(2));
        fprintf(fileID,'\r\nscale_1:                        %g', distParameters(3));
        fprintf(fileID,'\r\nmu2_1:                          %g', distParameters(4));
        fprintf(fileID,'\r\nsigma2_1:                       %g', distParameters(5));
        fprintf(fileID,'\r\ntheta_1:                        %g', distParameters(6));
        if strcmpi(distributionType, 'Bimodal Gaussian') || strcmpi(distributionType, 'Trimodal Gaussian') || strcmpi(distributionType, 'Quadrimodal Gaussian')
            fprintf(fileID,'\r\nmu1_2:                          %g', distParameters(7));
            fprintf(fileID,'\r\nsigma1_2:                       %g', distParameters(8));
            fprintf(fileID,'\r\nscale_2:                        %g', distParameters(9));
            fprintf(fileID,'\r\nmu2_2:                          %g', distParameters(10));
            fprintf(fileID,'\r\nsigma2_2:                       %g', distParameters(11));
            fprintf(fileID,'\r\ntheta_2:                        %g', distParameters(12));
            if strcmpi(distributionType, 'Trimodal Gaussian') || strcmpi(distributionType, 'Quadrimodal Gaussian')
                fprintf(fileID,'\r\nmu1_3:                          %g', distParameters(13));
                fprintf(fileID,'\r\nsigma1_3:                       %g', distParameters(14));
                fprintf(fileID,'\r\nscale_3:                        %g', distParameters(15));
                fprintf(fileID,'\r\nmu2_3:                          %g', distParameters(16));
                fprintf(fileID,'\r\nsigma2_3:                       %g', distParameters(17));
                fprintf(fileID,'\r\ntheta_3:                        %g', distParameters(18));
                if strcmpi(distributionType, 'Quadrimodal Gaussian')
                    fprintf(fileID,'\r\nmu1_4:                          %g', distParameters(19));
                    fprintf(fileID,'\r\nsigma1_4:                       %g', distParameters(20));
                    fprintf(fileID,'\r\nscale_4:                        %g', distParameters(21));
                    fprintf(fileID,'\r\nmu2_4:                          %g', distParameters(22));
                    fprintf(fileID,'\r\nsigma2_4:                       %g', distParameters(23));
                    fprintf(fileID,'\r\ntheta_4:                        %g', distParameters(24));
                end
            end
        end
    elseif strcmpi(distributionType, 'Skew-normal')
        fprintf(fileID,'\r\nmu1:                            %g', distParameters(1));
        fprintf(fileID,'\r\nsigma1:                         %g', distParameters(2));
        fprintf(fileID,'\r\nscale:                          %g', distParameters(3));
        fprintf(fileID,'\r\nmu2:                            %g', distParameters(4));
        fprintf(fileID,'\r\nsigma2:                         %g', distParameters(5));
        fprintf(fileID,'\r\nrho:                            %g', distParameters(6));
        fprintf(fileID,'\r\nskew1:                          %g', distParameters(7));
        fprintf(fileID,'\r\nskew2:                          %g', distParameters(8));
    elseif strcmpi(distributionType, 'Bimodal Skew-normal')
        fprintf(fileID,'\r\nmu1_1:                          %g', distParameters(1));
        fprintf(fileID,'\r\nsigma1_1:                       %g', distParameters(2));
        fprintf(fileID,'\r\nscale1:                         %g', distParameters(3));
        fprintf(fileID,'\r\nmu2_1:                          %g', distParameters(4));
        fprintf(fileID,'\r\nsigma2_1:                       %g', distParameters(5));
        fprintf(fileID,'\r\nrho1:                           %g', distParameters(6));
        fprintf(fileID,'\r\nmu1_2:                          %g', distParameters(7));
        fprintf(fileID,'\r\nsigma1_2:                       %g', distParameters(8));
        fprintf(fileID,'\r\nscale2:                         %g', distParameters(9));
        fprintf(fileID,'\r\nmu2_2:                          %g', distParameters(10));
        fprintf(fileID,'\r\nsigma2_2:                       %g', distParameters(11));
        fprintf(fileID,'\r\nrho2:                           %g', distParameters(12));
        fprintf(fileID,'\r\nskew1_1:                        %g', distParameters(13));
        fprintf(fileID,'\r\nskew2_1:                        %g', distParameters(14));
        fprintf(fileID,'\r\nskew1_2:                        %g', distParameters(15));
        fprintf(fileID,'\r\nskew2_2:                        %g', distParameters(16));
    elseif strcmpi(distributionType, 'Trimodal Skew-normal')
        fprintf(fileID,'\r\nmu1_1:                          %g', distParameters(1));
        fprintf(fileID,'\r\nsigma1_1:                       %g', distParameters(2));
        fprintf(fileID,'\r\nscale1:                         %g', distParameters(3));
        fprintf(fileID,'\r\nmu2_1:                          %g', distParameters(4));
        fprintf(fileID,'\r\nsigma2_1:                       %g', distParameters(5));
        fprintf(fileID,'\r\nrho1:                           %g', distParameters(6));
        fprintf(fileID,'\r\nmu1_2:                          %g', distParameters(7));
        fprintf(fileID,'\r\nsigma1_2:                       %g', distParameters(8));
        fprintf(fileID,'\r\nscale2:                         %g', distParameters(9));
        fprintf(fileID,'\r\nmu2_2:                          %g', distParameters(10));
        fprintf(fileID,'\r\nsigma2_2:                       %g', distParameters(11));
        fprintf(fileID,'\r\nrho2:                           %g', distParameters(12));
        fprintf(fileID,'\r\nmu1_3:                          %g', distParameters(13));
        fprintf(fileID,'\r\nsigma1_3:                       %g', distParameters(14));
        fprintf(fileID,'\r\nscale3:                         %g', distParameters(15));
        fprintf(fileID,'\r\nmu2_3:                          %g', distParameters(16));
        fprintf(fileID,'\r\nsigma2_3:                       %g', distParameters(17));
        fprintf(fileID,'\r\nrho3:                           %g', distParameters(18));
        fprintf(fileID,'\r\nskew1_1:                        %g', distParameters(19));
        fprintf(fileID,'\r\nskew2_1:                        %g', distParameters(20));
        fprintf(fileID,'\r\nskew1_2:                        %g', distParameters(21));
        fprintf(fileID,'\r\nskew2_2:                        %g', distParameters(22));
        fprintf(fileID,'\r\nskew1_3:                        %g', distParameters(23));
        fprintf(fileID,'\r\nskew2_3:                        %g', distParameters(24));
    end
    
    
    
    % Save:
    dataFilename = strcat(filename,'.mat');
    varsToSave.simulatedEvents = simulatedEvents;
    varsToSave.distParameters = distParameters;
    varsToSave.baseline = baseline;
    varsToSave.targetEvents1D = targetEvents1D;
    varsToSave.targetEvents1D_RT = targetEvents1D_RT;
    varsToSave.targetEvents2D = targetEvents2D;
    varsToSave.targetDataLength = targetDataLength;
    varsToSave.excludedTimes = excludedTimes;
    varsToSave.distributionType = distributionType;
    varsToSave.detectionParameters = detectionParameters;
    varsToSave.classificationParameters = classificationParameters;
    varsToSave.optimisationParameters = optimisationParameters;
    varsToSave.simulationParameters = simulationParameters;
    varsToSave.filtering = filtering;
    varsToSave.targetSpectrum = targetSpectrum;
    varsToSave.tauRange = tauRange;
    varsToSave.cliff = cliff;
    varsToSave.measuredUF = measuredUF;
    varsToSave.histoData = histoData;
    varsToSave.costFuncStruct = costFuncStruct;
    varsToSave.individual = individual;
    varsToSave.simFilename = simFilename;
    varsToSave.simsimFilename = simsimFilename;
    varsToSave.fitness = fitness;
    varsToSave.fitnessTail = fitnessTail;
    varsToSave.evalVector = evalVector;
    varsToSave.bestEvalVector = bestEvalVector;
    varsToSave.evalVectorTail = evalVectorTail;
    varsToSave.bestEvalVectorTail = bestEvalVectorTail;
    varsToSave.meanMinisAmplitudes = meanMinisAmplitudes;
    varsToSave.meanMinisRTs = meanMinisRTs;
    varsToSave.medianMinisAmplitudes = medianMinisAmplitudes;
    varsToSave.medianMinisRTs = medianMinisRTs;
    varsToSave.SD = SD;
    varsToSave.count = count;
    varsToSave.minisFrequency = minisFrequency;
    varsToSave.minis1D = minis1D;
    varsToSave.minis1D_RT = minis1D_RT;
    varsToSave.minis2D = minis2D;
    varsToSave.tailMinis = tailMinis;
    varsToSave.shapes = shapes;
    varsToSave.shapesTail = shapesTail;
    varsToSave.simulatedEvents1D = simulatedEvents1D;
    varsToSave.simulatedEvents1D_RT = simulatedEvents1D_RT;
    varsToSave.simulatedEvents2D = simulatedEvents2D;
    varsToSave.simulatedEvents1DTail = simulatedEvents1DTail;
    varsToSave.simulatedEvents1D_RTTail = simulatedEvents1D_RTTail;
    varsToSave.simEvents2DTail = simEvents2DTail;
    varsToSave.constraintFitness = constraintFitness;
    varsToSave.allCosts = allCosts;
    varsToSave.criticalCosts = [2 5 13:15 17 20 28:30];
    %varsToSave.criticalCosts = [2 5 13:15 17 20 30];
    varsToSave.cliffEval = cliffEval;
    if any(ismember(find(allCosts), varsToSave.criticalCosts))
        varsToSave.satisfactoryFit = false;
    else
        varsToSave.satisfactoryFit = true;
    end
    if ~constraintFitness(1) && cliffEval(1) <= 0.99 && cliffEval(2) <= 0.5
        varsToSave.satisfactoryConstraints = true;
    else
        varsToSave.satisfactoryConstraints = false;
    end
    save(dataFilename, '-struct', 'varsToSave');
    
    if noiseProperties.nchans_to_save == 1
        writeABF(V,simFilename,1000/noiseProperties.dt,{'mV'});
        writeABF(simV,simsimFilename,1000/noiseProperties.dt,{'mV'});
        writeABF(simVTail,[simsimFilename(1:end-4) 'tail.abf'],1000/noiseProperties.dt,{'mV'});
    elseif noiseProperties.nchans_to_save == 2
        writeABF([V; noiseProperties.current],simFilename,1000/noiseProperties.dt,{'mV';'pA'});
        writeABF([simV; noiseProperties.current],simsimFilename,1000/noiseProperties.dt,{'mV';'pA'});
        writeABF([simVTail; zeros(size(noiseProperties.current))],[simsimFilename(1:end-4) 'tail.abf'],1000/noiseProperties.dt,{'mV';'pA'});
    else
        disp('error in noiseProperties.nchans_to_save');
    end
    
    
    
    % Plot:
    if draw
        if tauRange
            [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10] = plotMinis(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, targetEvents2D,...
                simulatedEvents2D, minis2D, shapes, distributionType, distParameters(1:end-2), tailMinis, amplitudeArraySim, targetSpectrum,...
                noiseSpectrum, freq, detectionParameters.RTinterval, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, evalVector,...
                bestEvalVector, costFuncStruct, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10);
        else
            [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10] = plotMinis(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, targetEvents2D,...
                simulatedEvents2D, minis2D, shapes, distributionType, distParameters(1:end-1), tailMinis, amplitudeArraySim, targetSpectrum,...
                noiseSpectrum, freq, detectionParameters.RTinterval, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, evalVector,...
                bestEvalVector, costFuncStruct, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10);
        end
        saveas(f1, strcat(filename,'_amplitudes.jpg'));
        saveas(f1, strcat(filename,'_amplitudes.fig'));
        saveas(f2, strcat(filename,'_rise_times.jpg'));
        saveas(f2, strcat(filename,'_rise_times.fig'));
        saveas(f3, strcat(filename,'_2D.jpg'));
        saveas(f3, strcat(filename,'_2D.fig'));
        saveas(f4, strcat(filename,'_minis.jpg'));
        saveas(f4, strcat(filename,'_minis.fig'));
        saveas(f5, strcat(filename,'_source_distributions.jpg'));
        saveas(f5, strcat(filename,'_source_distributions.fig'));
        saveas(f6, strcat(filename,'_fft.jpg'));
        saveas(f6, strcat(filename,'_fft.fig'));
        if evalVector(1)
            saveas(f7, strcat(filename,'_mult_amplitudes.jpg'));
            saveas(f7, strcat(filename,'_mult_amplitudes.fig'));
            saveas(f8, strcat(filename,'_mult_rise_times.jpg'));
            saveas(f8, strcat(filename,'_mult_rise_times.fig'));
            saveas(f9, strcat(filename,'_single_comparison_summary.jpg'));
            saveas(f9, strcat(filename,'_single_comparison_summary.fig'));
            saveas(f10, strcat(filename,'_group_comparison_summary.jpg'));
            saveas(f10, strcat(filename,'_group_comparison_summary.fig'));
        end
    end
    
elseif (individual == 1 || (~isempty(bestFitness) && (individual > 1 && fitness < bestFitness))) && ~isempty(shapes) && ~isempty(simulatedEvents2D)
    clear AmpsMAD AmpsMADBest AmpsMADMulti AmpsBottomMAD AmpsBottomMADBest AmpsBottomMADMulti AmpsBottomSAD AmpsBottomSADBest AmpsBottomSADMulti...
        AmpsBottomSADlowMulti AmpsBottomSADlow AmpsBottomSADlowBest AmpsMidMAD AmpsMidMADBest AmpsMidMADMulti AmpsMidSAD AmpsMidSADBest AmpsMidSADMulti...
        AmpsMidSADlowMulti AmpsMidSADlow AmpsMidSADlowBest AmpsSAD AmpsSADBest AmpsSADMulti AmpsTopMAD AmpsTopMADBest AmpsTopMADMulti AmpsTopSAD...
        AmpsTopSADBest AmpsTopSADMulti AmpsTopSADlowMulti AmpsTopSADlow AmpsTopSADlowBest MAD MADBest MADMulti RTsMAD RTsMADBest RTsMADMulti RTsSAD...
        RTsSADBest RTsSADMulti SAD SADBest SADMulti SD SWstart amplitudeArraySim baseline classificationParameters cliff costFuncStruct dataDir...
        detectionParameters distributionType excludedTimes filename1 filename10 filename11 filename2 filename3 filename4 filename5 filename6 filename7...
        filename8 filename9 filtN filtering fullParallel histoData lengthRatio measuredUF measuredUFBottom measuredUFMid measuredUFTop minis1D minis1D_RT...
        minis1D_RT75prc minis1D_RTprc minis2D multiEvalVector noiseDataLength noiseProperties noiseSpectrum optimisationParameters parallelCores...
        simulatedEvents simulatedEvents1D simulatedEvents1D_RT smoothWindow targetDataLength targetEvents1D targetEvents1D_RT targetEvents2D...
        targetSpectrum tauRange waveform
    bestFitness = fitness;
    proceed = false;
    fileCount = 1;
    while ~proceed
        if exist(strcat(filename, '.mat'), 'file')
            fileCount = fileCount + 1;
            fileCountStr = num2str(fileCount);
            if length(fileCountStr) == 1 %#ok<ISCL>
                fileCountStr = strcat('_0', fileCountStr);
            else
                fileCountStr = strcat('_', fileCountStr);
            end
            if fileCount > 2
                filename = strcat(filename(1:end-3), fileCountStr);
            else
                filename = strcat(filename, fileCountStr);
            end
        else
            proceed = true;
        end
    end
    dataFilename = strcat(filename, '.mat');
    
    save(dataFilename,'-v7.3');
    matObj = matfile(dataFilename, 'writable', true);
    matObj.simV(1:size(simV,1), 1:size(simV,2)) = simV;
    matObj.simVTail(1:size(simVTail,1), 1:size(simVTail,2)) = simVTail;
    matObj.distParameters(1:1,1:length(distParameters)) = distParameters;
    matObj.shapes(1:size(shapes,1),1:size(shapes,2)) = shapes;
    matObj.shapesTail(1:size(shapesTail,1),1:size(shapesTail,2)) = shapesTail;
    matObj.simulatedEvents2D(1:size(simulatedEvents2D,1),1:size(simulatedEvents2D,2)) = simulatedEvents2D;
    matObj.simEvents2DTail(1:size(simEvents2DTail,1),1:size(simEvents2DTail,2)) = simEvents2DTail;
    matObj.tailMinis(1:size(tailMinis,1),1:size(tailMinis,2)) = tailMinis;
    matObj.tau_sy1(1:1,1:1) = simulationParameters.tau_sy1;
    matObj.tau_m(1:1,1:1) = simulationParameters.R_m;
    matObj.freq(1:1,1:length(freq)) = freq;
    matObj.fitness(1:1,1:1) = fitness;
    matObj.fitnessTail(1:1,1:1) = fitnessTail;
    matObj.evalVector(1:1,1:length(evalVector)) = evalVector;
    matObj.bestEvalVector(1:1,1:length(bestEvalVector)) = bestEvalVector;
    matObj.evalVectorTail(1:1,1:length(evalVectorTail)) = evalVectorTail;
    matObj.bestEvalVectorTail(1:1,1:length(bestEvalVectorTail)) = bestEvalVectorTail;
end
cd ..

individual = individual + 1;
fclose all;
end