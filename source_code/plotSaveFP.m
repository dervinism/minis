function [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10] = plotSaveFP(draw, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT,...
    targetEvents2D, targetDataLength, current, simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail,...
    excludedTimes, SD, distributionType, detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering,...
    targetSpectrum, noiseSpectrum, freq, tauRange, cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct,...
    fitness, fitnessTail, evalVector, evalVectorTail, bestEvalVector, bestEvalVectorTail, varargin)

if nargin == 43
    population = varargin{1};
else
    population = distParameters;
end
f1 = []; f2 = []; f3 = []; f4 = []; f5 = []; f6 = []; f7 = []; f8 = []; f9 = []; f10 = [];

[minis1D, minis1D_RT, minis2D] = classifyMinis(shapes(:,2), shapes(:,3), classificationParameters);



%% Classify simulated minis+noise events:
simulatedEvents1D = sum(current,1);
simulatedEvents2D = current;
simulatedEvents1D_RT = sum(current,2)';



%% Scale simulated minis-dependent distributions:
noiseDataLength = length(V) - length(excludedTimes);
lengthRatio = targetDataLength/noiseDataLength;
minis1D = lengthRatio*minis1D;
minis2D = lengthRatio*minis2D;
minis1D_RT = lengthRatio*minis1D_RT;
minis1D(1) = 0;
minis2D(1,:) = zeros(1,size(minis2D,2));
minis1D_RT(1) = 0;



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
meanMinisAmplitudes = mean(shapes(:,2));
medianMinisAmplitudes = median(shapes(:,2));
meanMinisRTs = mean(shapes(:,3));
medianMinisRTs = median(shapes(:,3));
targetSweepDuration = noiseProperties.dt*targetDataLength/1000;             % seconds
count = sum(sum(minis2D));
minisFrequency = count/targetSweepDuration;                                 % Hz

if ~exist('log.txt','file')
    fileID = fopen('log.txt','w');
else
    fileID = fopen('log.txt','a');
end
simFilename = strcat(filename,'.abf');
simsimFilename = strcat(filename,'sim','.abf');



% Print:
fprintf(fileID, '\r\nFitness:                                                   %g', fitness);
fprintf(fileID, '\r\nFile name containing simulations superimposed on data:     %s', simFilename);
fprintf(fileID, '\r\nFile name containing simulations only:                     %s', simsimFilename);
fprintf(fileID, '\r\nFile name containing addStred minis tail:                     %s', strcat(simsimFilename, 'tail'));
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
fprintf(fileID, '\r\nMinis amplitude generation threshold:                      %g', simulationParameters.loSimAmp);
fprintf(fileID, '\r\nTau_m:                                                     %g', simulationParameters.tau_m);
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
varsToSave.population = population;
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
varsToSave.simEvents2DTail = simEvents2DTail; %#ok<*STRNU>
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
    amplitudeArraySim = classificationParameters.amplitudeArray;
    amplitudeArraySim(amplitudeArraySim < simulationParameters.loSimAmp) = [];
    if tauRange
        [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10] = plotMinis(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, targetEvents2D,...
            simulatedEvents2D, minis2D, shapes, distributionType, distParameters(1:end-2), tailMinis, amplitudeArraySim, targetSpectrum,...
            noiseSpectrum, freq, detectionParameters.RTinterval, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, evalVector,...
            bestEvalVector, costFuncStruct, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10);
    else
        [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10] = plotMinis(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, targetEvents2D,...
            simulatedEvents2D, minis2D, shapes, distributionType, distParameters(1:end-2), tailMinis, amplitudeArraySim, targetSpectrum,...
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
cd ..