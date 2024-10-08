function [outputData, outputData2, figures] = dataMinis(dataDirectory, excludedTimes, detectionParameters, classificationParameters, pulseAmplitude, APinfuse, APblock,...
    minisInfuse, initFile, parallelCores)

global ratio



%% Load abf data files from the data directory:
cd(dataDirectory);
files = dir('*.abf');
dataProperties = struct([]);
excludedTimesArray = struct([]);
dataLength = 0;
excludedLength = 0;
nSweeps = zeros(length(files),1);
fileDuration = zeros(length(files),2);
sweepDuration = zeros(length(files),1);
startPulse = 1000*excludedTimes.startPulse;
endPulse = 1000*excludedTimes.endPulse;
progress = 100/length(files): 100/length(files) :100;
h = waitbar(0,'Loading files - 0%');
waitText1 = 'Loading files - ';
waitText2 = '%';
for i = 1:length(files)
    dataProperties{i} = loadABF(files(i).name);
    dataLength = dataLength + length(dataProperties{i}.sweep);
    fileDuration(i,1) = length(dataProperties{i}.sweep);
    fileDuration(i,2) = dataLength;
    nSweeps(i) = dataProperties{i}.hd.lActualEpisodes;
    sweepDuration(i) = floor((length(dataProperties{i}.sweep)*dataProperties{i}.dt - dataProperties{i}.dt)/nSweeps(i));
    excludedTimesArray{i} = calcExcludedTimes(sweepDuration(i), nSweeps(i), startPulse, endPulse, [], [], dataProperties{i}.dt);
    excludedLength = excludedLength + length(excludedTimesArray{i});
    waitText = [waitText1 num2str(progress(i)) waitText2];
    waitbar(progress(i)/100, h, waitText);
end
close(h);



%% Detect minis:
fileDuration = [0 0; fileDuration];
detectionParameters.sampleInterval = dataProperties{1}.dt;
detectionParameters.smoothWindow = round(detectionParameters.smoothWindow/detectionParameters.sampleInterval);
detectionParameters.smoothWindowLite = 8;
filtering.state = 'off';
waveform.estimate = 0;
minis = {[]};
V = minis;
fileAxis = zeros(1, fileDuration(end,2));
h = waitbar(0,'Detecting minis-like events - 0%');
waitText1 = 'Detecting minis-like events - ';
for i = 1:length(files)
    [minis{i}, ~, ~, ~, ~, V{i}] = detectMinis(dataProperties{i}.sweep, excludedTimesArray{i}, detectionParameters, filtering, waveform, parallelCores);
    fileAxis(fileDuration(i,2)+1:fileDuration(i+1,2)) = initFile + (i-1) + (1:fileDuration(i+1,1))*(1/fileDuration(i+1,1)) - 1/fileDuration(i+1,1);
    waitText = [waitText1 num2str(progress(i)) waitText2];
    waitbar(progress(i)/100, h, waitText);
end
close(h);
% %            1               2          3                     4           5          6-7         8-9         10                11
% coreMinis = [finalisedPeaks, peakTimes, finalisedPeakIndices, amplitudes, baselines, tBaselines, nBaselines, elementRiseTimes, riseTimes,...
%     riseTimes1090, t10s, t50s, t90s, n10s, n50s, n90s, v10s, v50s, v90s];
% %   12             13    14    15    16    17    18    19    20    21



%% Estimate tau_m, capacitance, and the neuron's input resistance:
tau_m = zeros(2*length(files), 1);
tauPulse = tau_m;
tauPulseEff = tau_m;
capacitance = tau_m;
capacitanceEff = tau_m;
Rseries = tau_m;
nHalfF = zeros(1,2*length(files));
progress = 34/(2*length(files)): 34/(2*length(files)) :34;
h = waitbar(0,'Estimating parameters - 0%');
waitText1 = 'Estimating parameters - ';
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
for i = 1:2*length(files)
    file = ceil(i/2);
    minisCurrent = minis{file};
    dataLength = length(dataProperties{file}.sweep);
    sweepDuration = dataLength/dataProperties{file}.hd.lActualEpisodes;
    if rem(i,2)
        halfSweeps = ceil(dataProperties{file}.hd.lActualEpisodes/2);
        minisCurrent = minisCurrent(minisCurrent(:,3) <= round(fileDuration(file+1,1)/2), :);
        nHalfF(i) = fileDuration(file,2) + halfSweeps*sweepDuration;
        iData = 1:halfSweeps*sweepDuration;
    else
        halfSweeps = floor(dataProperties{file}.hd.lActualEpisodes/2);
        minisCurrent = minisCurrent(minisCurrent(:,3) > round(fileDuration(file+1,1)/2), :);
        nHalfF(i) = fileDuration(file+1,2);
        iData = dataLength - halfSweeps*sweepDuration + 1 : dataLength;
    end
    % Estimate effective tau_m based on PSPs and the average and median amplitude of top 10% minis:
    [~, ~, ~, ~, parameters] = averageWaveEffect(minisCurrent, dataProperties{file}.sweep, dataProperties{file}.dt, 10, excludedTimesArray{file}, 10,...
        detectionParameters.RTinterval, classificationParameters.riseTimeArrayExt, false);
    if ~isempty(parameters.tau_m)
        tau_m(i) = parameters.tau_m;
    else
        tau_m(i) = NaN;
    end
    % Estimate tau_m, total capacitance, and series resistance based on impulses:
    if ~isempty(startPulse) && length(startPulse) > 1
        [tauPulse(i), tauPulseEff(i), capacitance(i), capacitanceEff(i), Rseries(i)] = peel(1, detectionParameters.pulseDuration,...
            dataProperties{file}.sweep(iData), dataProperties{file}.current(iData), halfSweeps, dataProperties{file}.dt, false, startPulse(1),...
            endPulse(1), startPulse(2), endPulse(2), pulseAmplitude);
        [tauPulseNeg, tauPulseEffNeg, capacitanceNeg, capacitanceEffNeg] = peel(-1, detectionParameters.pulseDuration, dataProperties{file}.sweep(iData),...
            dataProperties{file}.current(iData), halfSweeps, dataProperties{file}.dt, false, [], [], startPulse(2), endPulse(2), pulseAmplitude);
        tauPulse(i) = mean([tauPulse(i) tauPulseNeg]);
        tauPulseEff(i) = mean([tauPulseEff(i) tauPulseEffNeg]);
        capacitance(i) = mean([capacitance(i) capacitanceNeg]);
        capacitanceEff(i) = mean([capacitanceEff(i) capacitanceEffNeg]);
    end
    waitText = [waitText1 num2str(progress(i)) waitText2];
    waitbar(progress(i)/100, h, waitText);
end
warning('on', 'MATLAB:plot:IgnoreImaginaryXYPart');



%% Estimate the mean and median amplitude of top 10% minis:
averageAmp20 = zeros(sum(nSweeps), 1);
medianAmp20 = averageAmp20;
averageAmp10 = averageAmp20;
medianAmp10 = averageAmp20;
averageAmp5 = averageAmp20;
medianAmp5 = averageAmp20;
averageAmp2 = averageAmp20;
medianAmp2 = averageAmp20;
averageAmp1 = averageAmp20;
medianAmp1 = averageAmp20;
nSweepEnd = averageAmp20;
sweepCount = 0;
progress = 34+33/sum(nSweeps): 33/sum(nSweeps) :34+33;
for iFile = 1:length(files)
    minisFile = minis{iFile};
    sweepDuration = floor(fileDuration(iFile+1,1)/nSweeps(iFile));
    for iSweep = 1:nSweeps(iFile)
        sweepCount = sweepCount + 1;
        minisSweep = minisFile(minisFile(:,3)>(iSweep-1)*sweepDuration+1 & minisFile(:,3)<iSweep*sweepDuration,:);
        % Find minis larger than the 80th percentile of minis amplitudes:
        prct80 = prctile(minisSweep(:,4),80);
        iPrct80 = minisSweep(:,4) >= prct80;
        averageAmp20(sweepCount) = mean(minisSweep(iPrct80,4));
        medianAmp20(sweepCount) = median(minisSweep(iPrct80,4));
        % Find minis larger than the 90th percentile of minis amplitudes:
        prct90 = prctile(minisSweep(:,4),90);
        iPrct90 = minisSweep(:,4) >= prct90;
        averageAmp10(sweepCount) = mean(minisSweep(iPrct90,4));
        medianAmp10(sweepCount) = median(minisSweep(iPrct90,4));
        % Find minis larger than the 95th percentile of minis amplitudes:
        prct95 = prctile(minisSweep(:,4),95);
        iPrct95 = minisSweep(:,4) >= prct95;
        averageAmp5(sweepCount) = mean(minisSweep(iPrct95,4));
        medianAmp5(sweepCount) = median(minisSweep(iPrct95,4));
        % Find minis larger than the 98th percentile of minis amplitudes:
        prct98 = prctile(minisSweep(:,4),98);
        iPrct98 = minisSweep(:,4) >= prct98;
        averageAmp2(sweepCount) = mean(minisSweep(iPrct98,4));
        medianAmp2(sweepCount) = median(minisSweep(iPrct98,4));
        % Find minis larger than the 99th percentile of minis amplitudes:
        prct99 = prctile(minisSweep(:,4),99);
        iPrct99 = minisSweep(:,4) >= prct99;
        averageAmp1(sweepCount) = mean(minisSweep(iPrct99,4));
        medianAmp1(sweepCount) = median(minisSweep(iPrct99,4));
        
        nSweepEnd(sweepCount) = fileDuration(iFile,2) + iSweep*sweepDuration;
        
        waitText = [waitText1 num2str(progress(sweepCount)) waitText2];
        waitbar(progress(sweepCount)/100, h, waitText);
    end
end



%% Estimate SD:
initSTD100 = struct([]);
initSTDsmooth100 = initSTD100;
initMean100 = initSTD100;
initMeanSmooth100 = initSTD100;
initSTD15 = initSTD100;
initSTDsmooth15 = initSTD100;
initMean15 = initSTD100;
initMeanSmooth15 = initSTD100;
inittSTD100 = initSTD100;
inittSTD15 = initSTD100;
dataLength100 = 0;
dataLengths100 = zeros(1,length(files)+1);
dataLength15 = 0;
dataLengths15 = dataLengths100;
progress = 34+33+17/length(files): 17/length(files) :34+33+17;
for iFile = 1:length(files)
    data = dataProperties{iFile}.sweep;
    dataSmooth = V{iFile};
    
    % Calculate excluded times:
    excludedT = excludedTimesArray{iFile};
    excludedIndices = round(excludedT/dataProperties{iFile}.dt) + 1;
    exclIndLogical = zeros(1, length(data));
    exclIndLogical(excludedIndices) = 1;
    circ = detectionParameters.SWstart + round(detectionParameters.BLduration/dataProperties{iFile}.dt);
    response = ones(1, 1 + 2*circ);
    excludedIndices = logical(conv(exclIndLogical, response, 'same'));
    excludedIndices(1 : (detectionParameters.SWstart + detectionParameters.refractoryPeriod)/dataProperties{iFile}.dt - 1) = true;
    
    [initSTD100{iFile}, initSTDsmooth100{iFile}, inittSTD100{iFile}, initMean100{iFile}, initMeanSmooth100{iFile}] = dataSD(100, dataProperties{iFile}.dt,...
        data, dataSmooth, excludedIndices, fileDuration(iFile,2));
    dataLength100 = dataLength100 + length(initSTD100{iFile});
    dataLengths100(iFile+1) = dataLength100;
    [initSTD15{iFile}, initSTDsmooth15{iFile}, inittSTD15{iFile}, initMean15{iFile}, initMeanSmooth15{iFile}] = dataSD(15, dataProperties{iFile}.dt, data,...
        dataSmooth, excludedIndices, fileDuration(iFile,2));
    dataLength15 = dataLength15 + length(initSTD15{iFile});
    dataLengths15(iFile+1) = dataLength15;
    
    waitText = [waitText1 num2str(progress(iFile)) waitText2];
    waitbar(progress(iFile)/100, h, waitText);
end

STD100 = zeros(1,dataLength100);
STDsmooth100 = STD100;
mean100 = STD100;
meanSmooth100 = STD100;
tSTD100 = STD100;
STD15 = zeros(1,dataLength15);
STDsmooth15 = STD15;
mean15 = STD15;
meanSmooth15 = STD15;
tSTD15 = STD15;
fileSTD100 = zeros(length(files),1);
fileSTD15 = zeros(length(files),1);
progress = 34+33+17+8/length(files): 8/length(files) :34+33+17+8;
for iFile = 1:length(files)
    STD100(dataLengths100(iFile)+1:dataLengths100(iFile+1)) = initSTD100{iFile};
    STDsmooth100(dataLengths100(iFile)+1:dataLengths100(iFile+1)) = initSTDsmooth100{iFile};
    mean100(dataLengths100(iFile)+1:dataLengths100(iFile+1)) = initMean100{iFile};
    meanSmooth100(dataLengths100(iFile)+1:dataLengths100(iFile+1)) = initMeanSmooth100{iFile};
    tSTD100(dataLengths100(iFile)+1:dataLengths100(iFile+1)) = inittSTD100{iFile};
    STD15(dataLengths15(iFile)+1:dataLengths15(iFile+1)) = initSTD15{iFile};
    STDsmooth15(dataLengths15(iFile)+1:dataLengths15(iFile+1)) = initSTDsmooth15{iFile};
    mean15(dataLengths15(iFile)+1:dataLengths15(iFile+1)) = initMean15{iFile};
    meanSmooth15(dataLengths15(iFile)+1:dataLengths15(iFile+1)) = initMeanSmooth15{iFile};
    tSTD15(dataLengths15(iFile)+1:dataLengths15(iFile+1)) = inittSTD15{iFile};
    fileSTD100(iFile) = mean(initSTD100{iFile});
    fileSTD15(iFile) = mean(initSTD15{iFile});
    
    waitText = [waitText1 num2str(progress(iFile)) waitText2];
    waitbar(progress(iFile)/100, h, waitText);
end

sweepSTD100 = zeros(sum(nSweeps), 1);
sweepSTDsmooth100 = sweepSTD100;
sweepMean100 = sweepSTD100;
sweepMeanSmooth100 = sweepSTD100;
sweepSTD15 = sweepSTD100;
sweepSTDsmooth15 = sweepSTD100;
sweepSTD15med = sweepSTD100;
sweepSTDsmooth15med = sweepSTD100;
sweepMean15 = sweepSTD15;
sweepMeanSmooth15 = sweepSTD15;
sweepCount = 0;
progress = 34+33+17+8+8/sum(nSweeps): 8/sum(nSweeps) :100;
for iFile = 1:length(files)
    fSTD100 = initSTD100{iFile};
    fSTDsmooth100 = initSTDsmooth100{iFile};
    fMean100 = initMean100{iFile};
    fMeansmooth100 = initMeanSmooth100{iFile};
    fSTD15 = initSTD15{iFile};
    fSTDsmooth15 = initSTDsmooth15{iFile};
    fMean15 = initMean15{iFile};
    fMeansmooth15 = initMeanSmooth15{iFile};
    sweepDuration = floor(fileDuration(iFile+1,1)/nSweeps(iFile));
    nSTD100 = round(inittSTD100{iFile}/dataProperties{iFile}.dt) + 1 - fileDuration(iFile,2);
    nSTD15 = round(inittSTD15{iFile}/dataProperties{iFile}.dt) + 1 - fileDuration(iFile,2);
    for iSweep = 1:nSweeps(iFile)
        sweepCount = sweepCount + 1;
        nSTD100sweep = find(nSTD100>(iSweep-1)*sweepDuration+1 & nSTD100<iSweep*sweepDuration);
        nSTD15sweep = find(nSTD15>(iSweep-1)*sweepDuration+1 & nSTD15<iSweep*sweepDuration);
        sweepSTD100(sweepCount) = mean(fSTD100(nSTD100sweep));
        sweepSTDsmooth100(sweepCount) = mean(fSTDsmooth100(nSTD100sweep));
        sweepMean100(sweepCount) = mean(fMean100(nSTD100sweep));
        sweepMeanSmooth100(sweepCount) = mean(fMeansmooth100(nSTD100sweep));
        sweepSTD15(sweepCount) = mean(fSTD15(nSTD15sweep));
        sweepSTDsmooth15(sweepCount) = mean(fSTDsmooth15(nSTD15sweep));
        sweepSTD15med(sweepCount) = median(fSTD15(nSTD15sweep));
        sweepSTDsmooth15med(sweepCount) = median(fSTDsmooth15(nSTD15sweep));
        sweepMean15(sweepCount) = mean(fMean15(nSTD15sweep));
        sweepMeanSmooth15(sweepCount) = mean(fMeansmooth15(nSTD15sweep));
        
        waitText = [waitText1 num2str(progress(iFile)) waitText2];
        waitbar(progress(iFile)/100, h, waitText);
    end
end
sweepSTDmean = [sweepSTD100 sweepSTDsmooth100 sweepSTD15 sweepSTDsmooth15 sweepSTD15med sweepSTDsmooth15med sweepMean100 sweepMeanSmooth100 sweepMean15...
    sweepMeanSmooth15];
close(h);










%% Plot the data:
% Plot the electrophysiological recording data:
tEnd = 0.001*fileDuration(end,2)*detectionParameters.sampleInterval - detectionParameters.sampleInterval;
tSweep = 0.001*(nSweepEnd*detectionParameters.sampleInterval - detectionParameters.sampleInterval);
fSweep = fileAxis(nSweepEnd);
tHalfF = 0.001*(nHalfF*detectionParameters.sampleInterval - detectionParameters.sampleInterval);
fHalfF = fileAxis(nHalfF);
swpSTD100 = fileAxis(round(tSTD100/detectionParameters.sampleInterval) + 1);
swpSTD15 = fileAxis(round(tSTD15/detectionParameters.sampleInterval) + 1);
ratio = (tEnd+detectionParameters.sampleInterval)/length(files);
clear fileAxis minis minisCurrent data dataSmooth excludedT excludedIndices exclIndLogical window rows vert horz inds arrayData arrayDataSmooth arrayExcl

button = questdlg('Save the electrophysiological recording figures?','Save Figures','Yes','No','Yes');
if strcmpi(button, 'Yes')
    warning('off', 'MATLAB:gui:latexsup:BadTeXString');
    screenSize = get(0, 'ScreenSize');
    screenSize = [screenSize(3) screenSize(4)];
    f1 = figure('Position', [17*screenSize(1)/2560 481*screenSize(2)/1024 2560*screenSize(1)/2560 440*screenSize(2)/1024]);
    set(f1, 'NumberTitle', 'off');
    set(f1, 'Name', 'Examine Data');
    
    uiwait(msgbox(...
        'Note that the graphs containing the electrophysiological data from each file will be saved using Matlab "fig" format in a designated folder "figures".',...
        'modal'));
    dataDir = [pwd filesep 'figures'];
    if ~exist(dataDir,'dir')
        mkdir(dataDir);
    end
    cd(dataDir);
    lastSweep = 0;
    progress = 100/length(files): 100/length(files) :100;
    h = waitbar(0,'Saving figures - 0%');
    waitText1 = 'Saving figures - ';
    for iFile = 1:length(files)
        vFile = V{iFile};
        tFile = (1:length(vFile))*dataProperties{iFile}.dt - dataProperties{iFile}.dt;
        figure(f1);
        p1 = plot(tFile,vFile, 'k', 'markersize',4);
        hold on
        p2 = plot(excludedTimesArray{iFile},vFile(round(excludedTimesArray{iFile}/dataProperties{iFile}.dt) + 1), '.y', 'markersize',5);
        ax = gca;
        set(ax,'Position' , [0.025 0.1 0.965 0.8]);
        set(ax,'TickLength', [0.002 0.01]);
        xlabel('Time (ms)');
        ylabel('mV or nA');
        
        nameStr = [files(iFile).name(1:end-4) '___no' num2str(iFile+initFile-1) '_sw' num2str(lastSweep+1) '-' num2str(lastSweep+nSweeps(iFile))];
        titleStr = ['Smoothed data from file ' nameStr];
        title(titleStr,'Interpreter','none');
        legend([p1 p2],'Data','Excluded times', 'Location','NorthEast');
        hold off
        
        saveas(f1, nameStr, 'fig');
        lastSweep = lastSweep + nSweeps(iFile);
        waitText = [waitText1 num2str(progress(iFile)) waitText2];
        figure(h);
        waitbar(progress(iFile)/100, h, waitText);
    end
    warning('on', 'MATLAB:gui:latexsup:BadTeXString');
    close(f1);
    close(h);
    cd ..
end

clear screenSize h dataDir lastSweep iFile excludedTimes vFile V fileDuration tFile dataProperties p1 p2 excludedTimesArray ax titleStr fileLength...
   fileDuration sweepDuration startPulse endPulse startGlitch endGlitch progress waitText1 waitText2 filtering waveform i parameters initSTD100...
   initSTDsmooth100 inittSTD100 dataLength100 dataLengths100 initSTD15 initSTDsmooth15 inittSTD15 dataLength15 dataLengths15 circ response f1 fig



% Save the workspace:
button = questdlg('Save the current workspace for later re-use?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.mat','MAT files (*.mat)'},'Save workspace as', dataDirectory);
    if filterIndex
        eventFilename = fullfile(eventPathname, eventFilename);
        save(eventFilename);
    end
end



% Plot the average and median amplitude value of the top 10% of detected events:
f2 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(tSweep,medianAmp10', 'LineStyle',':', 'LineWidth',1.5, 'Color','r');
ax1 = gca;
xlabel('Recoeding time (s)');
ylabel('mV or nA');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
l2 = line(fSweep,averageAmp10', 'LineStyle',':', 'LineWidth',1.5, 'Color','g', 'Parent',ax2);
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l3 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f2, 'NumberTitle', 'off');
set(f2, 'Name', 'Examine Amplitudes');
title('Amplitudes of the top 10% of detected events');
if nTimeMarks
    legend([l1 l2 l3(1)],'Median','Mean','Drug timing', 'Location','NorthEast');
else
    legend([l1 l2],'Median','Mean', 'Location','NorthEast');
end

% Create the zoom object for the figure f2:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f2:
figurePanHandle = pan(f2);
set(figurePanHandle,'ActionPostCallback',@minisPan);


% Plot tau_m:
f3 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(tHalfF,tau_m, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('Time constant (ms)');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l2 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l2(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f3, 'NumberTitle', 'off');
set(f3, 'Name', 'Examine PSP-based Effective Membrane Time Constant');
title('Effective passive membrane time constant of the top 10% of detected events');
if nTimeMarks
    legend([l1 l2(1)],'Effective membrane time constant','Drug timing', 'Location','NorthEast');
else
    legend(l1,'Effective membrane time constant', 'Location','NorthEast');
end

% Create the zoom object for the figure f3:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f3:
figurePanHandle = pan(f3);
set(figurePanHandle,'ActionPostCallback',@minisPan);



% Plot SD (100ms window):
f4 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(0.001*tSTD100, STD100, 'Color','b');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('mV or nA');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
l2 = line(swpSTD100,STDsmooth100, 'Color','g', 'Parent',ax2);
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l3 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f4, 'NumberTitle', 'off');
set(f4, 'Name', 'Examine Standard Deviation (100 ms Window)');
title('Standard deviation of the electrophysiological recording data (100 ms window)');
if nTimeMarks
    legend([l1 l2 l3(1)],'Raw data','Smoothed data','Drug timing', 'Location','NorthEast');
else
    legend([l1 l2],'Raw data','Smoothed data', 'Location','NorthEast');
end

% Create the zoom object for the figure f4:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f4:
figurePanHandle = pan(f4);
set(figurePanHandle,'ActionPostCallback',@minisPan);



% Plot SD (15ms window):
f5 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(0.001*tSTD15, STD15, 'Color','b');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('mV or nA');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
l2 = line(swpSTD15,STDsmooth15, 'Color','g', 'Parent',ax2);
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l3 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f5, 'NumberTitle', 'off');
set(f5, 'Name', 'Examine Standard Deviation (15 ms Window)');
title('Standard deviation of the electrophysiological recording data (15 ms window)');
if nTimeMarks
    legend([l1 l2 l3(1)],'Raw data','Smoothed data','Drug timing', 'Location','NorthEast');
else
    legend([l1 l2],'Raw data','Smoothed data', 'Location','NorthEast');
end

% Create the zoom object for the figure f5:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f5:
figurePanHandle = pan(f5);
set(figurePanHandle,'ActionPostCallback',@minisPan);



% Plot the baseline:
f6 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(0.001*tSTD100, mean100, 'Color','b');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('mV or nA');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
l2 = line(swpSTD100,meanSmooth100, 'Color','g', 'Parent',ax2);
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l3 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f6, 'NumberTitle', 'off');
set(f6, 'Name', 'Examine Baseline Fluctuations (100 ms Window)');
title('Baseline of the electrophysiological recording data (100 ms window)');
if nTimeMarks
    legend([l1 l2 l3(1)],'Raw data','Smoothed data','Drug timing', 'Location','NorthEast');
else
    legend([l1 l2],'Raw data','Smoothed data', 'Location','NorthEast');
end

% Create the zoom object for the figure f6:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f6:
figurePanHandle = pan(f6);
set(figurePanHandle,'ActionPostCallback',@minisPan);



% Plot pulse-based tau_m:
f7 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(tHalfF,tauPulseEff, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('Time constant (ms)');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
l2 = line(fHalfF,tauPulse, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g', 'Parent',ax2);
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l3 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f7, 'NumberTitle', 'off');
set(f7, 'Name', 'Examine Impulse-based Passive Membrane Time Constant');
title('Impulse-based passive membrane time constant');
if nTimeMarks
    legend([l2 l1 l3(1)],'Membrane time constant','Effective time constant','Drug timing', 'Location','NorthEast');
else
    legend([l2 l1],'Membrane time constant','Effective time constant', 'Location','NorthEast');
end

% Create the zoom object for the figure f7:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f7:
figurePanHandle = pan(f7);
set(figurePanHandle,'ActionPostCallback',@minisPan);



% Plot capacitance:
f8 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(tHalfF,capacitanceEff, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('Capacitance (pF)');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
l2 = line(fHalfF,capacitance, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g', 'Parent',ax2);
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l3 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l3(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f8, 'NumberTitle', 'off');
set(f8, 'Name', 'Examine Total Capacitance');
title('Total capacitance');
if nTimeMarks
    legend([l2 l1 l3(1)],'Capacitance','Effective capacitance','Drug timing', 'Location','NorthEast');
else
    legend([l2 l1],'Capacitance','Effective capacitance', 'Location','NorthEast');
end

% Create the zoom object for the figure f8:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f8:
figurePanHandle = pan(f8);
set(figurePanHandle,'ActionPostCallback',@minisPan);



% Plot neuron's input resistance:
f9 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
axes('Position', [.04 .07 .92 .82]);
l1 = line(tHalfF,Rseries, 'LineStyle',':', 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
ax1 = gca;
xlabel('Recording time (s)');
ylabel('Resistance (Megaohms)');

ax2 = axes('Position', get(ax1,'Position'), 'XAxisLocation','top', 'YAxisLocation','left', 'Color','none', 'XColor','k', 'YColor','k');
xlabel('Files');
axisvals1 = get(ax1,'ylim');
axisvals2 = get(ax2,'ylim');
yLimit = [min([axisvals1(1) axisvals2(1)]) max([axisvals1(2) axisvals2(2)])];
set(ax1,'xlim',[0 tEnd]);
set(ax2,'xlim',[initFile initFile+length(files)]);
set(ax1,'ylim',yLimit);
set(ax2,'ylim',yLimit);

xTimeMarks = [APinfuse APblock minisInfuse];
yTimeMarks1 = [yLimit(1) yLimit(1) yLimit(1)];
nTimeMarks = length(xTimeMarks);
yTimeMarks1(nTimeMarks+1 : end) = [];
yTimeMarks2 = [yLimit(2) yLimit(2) yLimit(2)];
yTimeMarks2(nTimeMarks+1 : end) = [];
line(xTimeMarks,yTimeMarks1, 'LineStyle',':', 'Marker','^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
line(xTimeMarks,yTimeMarks2, 'LineStyle',':', 'Marker','v', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'Parent',ax1);
l2 = zeros(1,nTimeMarks);
for iMark = 1:nTimeMarks
    l2(iMark) = line([xTimeMarks(iMark) xTimeMarks(iMark)],[yLimit(1) yLimit(2)], 'Color','r', 'Parent',ax1);
end

set(f9, 'NumberTitle', 'off');
set(f9, 'Name', 'Examine neuron''s input resistance');
title('Neuron''s input resistance (assuming electrode''s resistance is balanced)');
if nTimeMarks
    legend([l1 l2(1)],'R_N','Drug timing', 'Location','NorthEast');
else
    legend(l1,'R_N', 'Location','NorthEast');
end

% Create the zoom object for the figure f9:
figureZoomHandle = zoom;
set(figureZoomHandle,'ActionPostCallback',@minisZoom);

% Create the pan object for the figure f9:
figurePanHandle = pan(f9);
set(figurePanHandle,'ActionPostCallback',@minisPan);



figures = [f2 f3 f4 f5 f6 f7 f8 f9];
outputData = [fSweep' averageAmp20 medianAmp20 averageAmp10 medianAmp10 averageAmp5 medianAmp5 averageAmp2 medianAmp2 averageAmp1 medianAmp1 sweepSTDmean];
outputData2 = [fHalfF' tau_m tauPulse tauPulseEff capacitance capacitanceEff Rseries];



fclose all;
disp('Task completed');



% %% Choose sweeps to include:
% options.WindowStyle = 'Normal';
% answer = inputdlg('Enter sweep numbers', 'Choose Sweeps to Include', 1, [], options);
% if ~answer
% endPulseStr = strrep(endPulseStr, ' ', '');
% endPulseStrSplit = regexp(endPulseStr, ',*', 'split');
% endPulse = strarray2numarray(endPulseStrSplit');
% if sum(isnan(endPulse))
%     errmsg = 'Error: at least one of the entries is NaN';
%     msgbox(errmsg,'Error','Error')
% end
% end