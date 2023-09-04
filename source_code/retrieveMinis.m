function retrieveMinis(addTail)

cd C:\Users\Martynas\Desktop\Phd\Programs\minis_detection_and_simulation\simulations_multiple_minis\minis_program\functions
addpath(pwd);
dataDir = 'data';
cd(dataDir);
load('initvar');
proceed = false;
fileOrder = 1;
while ~proceed
    filename = findFile('name', fileOrder);
    if length(filename) > 4 && strcmpi(filename(end-3:end), '.mat')
        load(filename) %#ok<*LOAD>
        proceed = true;
    else
        fileOrder = fileOrder + 1;
    end
end

if exist('tau_sy1','var') && exist('tau_m','var')
    simulationParameters.tau_sy1 = tau_sy1;
    simulationParameters.tau_sy2 = 0.1*tau_sy1;
    simulationParameters.tau_m = tau_m;
end
if ~exist('simV','var')
    dataPropertiesSim = loadABF(strcat(filename(1:end-4),'sim.abf'));
    simV = dataPropertiesSim.sweep;
end
if ~exist('simVTail','var')
    dataPropertiesSimTail = loadABF(strcat(filename(1:end-4),'simtail.abf'));
    simVTail = dataPropertiesSimTail.sweep;
end

filtN.state = 'spectrum';
waveformS.estimate = 3;
waveformS.tau_m = simulationParameters.tau_m;
waveformS.riseTimeArrayExt = classificationParameters.riseTimeArrayExt;
options.SD = [0 1 0 0];
if (fitness == fitnessTail) || addTail
    V = noiseProperties.sweep + simV + simVTail; %#ok<*NODEF>
else
    V = noiseProperties.sweep + simV;
end
[detectForSD, ~, nSpectrum, waveform] = detectMinis(V, noiseExcludedTimes, detectionParameters, filtN, waveformS, 8, options);
SD = detectForSD(1,22);
if ~exist('freq', 'var')
    freq = nSpectrum(4,:);
end
nSpectrum = nSpectrum(2,:);
saveas(waveform.F(1), strcat(filename(1:end-4), '_waveform.jpg'));
saveas(waveform.F(1), strcat(filename(1:end-4), '_waveform.fig'));
saveas(waveform.F(2), strcat(filename(1:end-4), '_90prctRTs.jpg'));
saveas(waveform.F(2), strcat(filename(1:end-4), '_90prctRTs.fig'));
cd ..

if exist('population','var')
    if addTail
        shapes = [shapes; shapesTail];
        plotSaveFP(true, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength,...
            simEvents2DTail, simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail, noiseExcludedTimes, SD, distributionType,...
            detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, nSpectrum, freq, tauRange,...
            cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fitnessTail, fitnessTail, evalVectorTail,...
            evalVectorTail, bestEvalVectorTail, bestEvalVectorTail, population);
    else
        plotSaveFP(true, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength,...
            simulatedEvents2D, simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail, noiseExcludedTimes, SD, distributionType,...
            detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, nSpectrum, freq, tauRange,...
            cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fitness, fitnessTail, evalVector, evalVectorTail,...
            bestEvalVector, bestEvalVectorTail, population);
    end
else
    if addTail
        shapes = [shapes; shapesTail];
        plotSaveFP(true, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength,...
            simEvents2DTail, simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail, noiseExcludedTimes, SD, distributionType,...
            detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, nSpectrum, freq, tauRange,...
            cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fitnessTail, fitnessTail, evalVectorTail,...
            evalVectorTail, bestEvalVectorTail, bestEvalVectorTail);
    else
        plotSaveFP(true, distParameters, dataDir, baseline, targetEvents1D, targetEvents1D_RT, targetEvents2D, targetDataLength,...
            simulatedEvents2D, simEvents2DTail, tailMinis, noiseProperties, shapes, shapesTail, V, simV, simVTail, noiseExcludedTimes, SD, distributionType,...
            detectionParameters, classificationParameters, optimisationParameters, simulationParameters, filtering, tSpectrum, nSpectrum, freq, tauRange,...
            cliff, measuredUF, measuredUFBottom, measuredUFMid, measuredUFTop, histoData, costFuncStruct, fitness, fitnessTail, evalVector, evalVectorTail,...
            bestEvalVector, bestEvalVectorTail);
    end
end

fclose all;
disp('Task completed');