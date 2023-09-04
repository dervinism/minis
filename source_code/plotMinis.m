function [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10] = plotMinis(amplitudeArray, riseTimeArray, target, current, minis, shapes, distributionType,...
    distributionParameters, tailMinis, amplitudeArraySim, targetSpectrum, noiseSpectrum, freq, RTtype, measuredUF, measuredUFBottom, measuredUFMid,...
    measuredUFTop, histoData, evalVector, bestEvalVector, costFuncStruct, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10)
% PLOMTMINIS is a helper function of the MANUALFITTINGMINIS and
% optimiseMinis functions used in MINIS GUI. For more information see the
% code.
%

draw = true;

%% Initialise figure windows:
if isempty(f1)
    f1 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
end
if isempty(f2)
    f2 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
end
if isempty(f3)
    f3 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
end
if isempty(f4)
    f4 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
end
if isempty(f5)
    f5 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
end
if isempty(f6)
    f6 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
end
ampLim2D = (amplitudeArray(2) - amplitudeArray(1))/2;
RTLim2D = (riseTimeArray(2) - riseTimeArray(1))/2;



%% Calculate the mean and the median:
shapes(:,shapes(:,2) > amplitudeArray(end) | shapes(:,3) > riseTimeArray(end)) = [];
meanMinisAmplitudes = mean(shapes(:,2));
meanMinisRT = mean(shapes(:,3));
medianMinisAmplitudes = median(shapes(:,2));
medianMinisRT = median(shapes(:,3));



%% Drawing uni-dimensional fits:
figure(f1);
targetAmpcs = cumsum(sum(target))/sum(sum(target));
targetAmp99prc = find(targetAmpcs > .99, 1);
currentAmpcs = cumsum(sum(current))/sum(sum(current));
currentAmp99prc = find(currentAmpcs > .99, 1);
iEndAmp = min([max([targetAmp99prc currentAmp99prc])+1 length(amplitudeArray)]);
plot([0 amplitudeArray(end)], [0 0], 'k-');
hold on
p1 = plot(amplitudeArray, sum(target,1), 'b.-');
p2 = plot(amplitudeArray, sum(current,1),'r.-');
p3 = plot(amplitudeArray, sum(target,1)-sum(current,1), 'g.-');
xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
set(f1, 'NumberTitle', 'off');
set(f1, 'Name', 'One-dimensional Amplitude Histograms');
title('Amplitude distribution fit');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
hold off

figure(f2);legend([p1 p2 p3],'Target','Simulated','Target-simulated','Location','NorthEast');
targetRTcs = cumsum(sum(target,2))/sum(sum(target));
targetRT99prc = find(targetRTcs > .99, 1);
currentRTcs = cumsum(sum(current,2))/sum(sum(current));
currentRT99prc = find(currentRTcs > .99, 1);
iEndRT = min([max([targetRT99prc currentRT99prc])+1 length(riseTimeArray)]);
plot([0 riseTimeArray(end)], [0 0], 'k-');
hold on
p1 = plot(riseTimeArray, sum(target,2),'b.-');
p2 = plot(riseTimeArray, sum(current,2),'r.-');
p3 = plot(riseTimeArray, sum(target,2)-sum(current,2), 'g.-');
xlim([riseTimeArray(1) riseTimeArray(find(sum(target,2),1,'last')+1)]);
set(f2, 'NumberTitle', 'off');
set(f2, 'Name', 'One-dimensional Rise time Histograms');
if strcmpi(RTtype, '10-90%')
    title('10-90% rise time distribution fit');
    xlabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    title('20-80% rise time distribution fit');
    xlabel('20-80% rise times (ms)');
end
ylabel('Number of Events');
legend([p1 p2 p3],'Target','Simulated','Target-simulated','Location','NorthEast');
hold off



%% Drawing two-dimensional fits:
cmin = min([min(min(target)) min(min(current))]);
cmax = max([max(max(target)) max(max(current))]);
clim = [cmin cmax];
figure(f3);
subplot(2,2,1,'replace');
hold on
imagesc(amplitudeArray, riseTimeArray, target, clim);
set(gca,'YDir','normal');
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
set(f3, 'NumberTitle', 'off');
set(f3, 'Name', '2-dimensional Histograms');
title('Target 2D event distribution');
xlabel('Amplitude(mV)');
if strcmpi(RTtype, '10-90%')
    ylabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    ylabel('20-80% rise times (ms)');
end
colorbar;
hold off

subplot(2,2,2,'replace');
hold on
imagesc(amplitudeArray, riseTimeArray, current, clim);
set(gca,'YDir','normal');
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
title('Current 2D event distribution');
xlabel('Amplitude(mV)');
if strcmpi(RTtype, '10-90%')
    ylabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    ylabel('20-80% rise times (ms)');
end
colorbar;
hold off

subplot(2,2,3,'replace');
hold on
imagesc(amplitudeArray, riseTimeArray, target-current);
set(gca,'YDir','normal');
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
title('Discrepancy of 2D event distributions');
xlabel('Amplitude(mV)');
if strcmpi(RTtype, '10-90%')
    ylabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    ylabel('20-80% rise times (ms)');
end
colorbar;
hold off

subplot(2,2,4,'replace');
hold on
surf(amplitudeArray, riseTimeArray, target-current);
shading interp;
set(gca,'YDir','normal');
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
title('Discrepancy of 2D event distributions');
xlabel('Amplitude(mV)');
if strcmpi(RTtype, '10-90%')
    ylabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    ylabel('20-80% rise times (ms)');
end
grid on
colorbar;
view(3);
hold off



%% Drawing minis distribution graphs:
figure(f4);
subplot(2,2,1,'replace');
hold on
minisAmpcs = cumsum(sum(minis))/sum(sum(minis));
minisAmp99prc = find(minisAmpcs > .999, 1);
iEndAmp = min([minisAmp99prc+1 length(amplitudeArray)]);
if isnan(iEndAmp)
    iEndAmp = length(amplitudeArray);
end
plot(amplitudeArray, sum(minis,1), 'b.-');
xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
set(f4, 'NumberTitle', 'off');
set(f4, 'Name', 'Simulated minis Histograms');
title('Minis amplitude distribution');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
s1Str(1) = {['mean: ' num2str(meanMinisAmplitudes)]};
s1Str(2) = {['median: ' num2str(medianMinisAmplitudes)]};
text(.7*amplitudeArray(iEndAmp),.9*max(sum(minis,1)),s1Str);
hold off

subplot(2,2,2,'replace');
hold on
minisRTcs = cumsum(sum(minis,2))/sum(sum(minis));
minisRT99prc = find(minisRTcs > .999, 1);
iEndRT = min([minisRT99prc+1 length(riseTimeArray)]);
if isnan(iEndRT)
    iEndRT = length(riseTimeArray);
end
plot(riseTimeArray, sum(minis, 2), 'b.-');
xlim([riseTimeArray(1) riseTimeArray(end)]);
title('Minis rise time distribution');
if strcmpi(RTtype, '10-90%')
    xlabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    xlabel('20-80% rise times (ms)');
end
ylabel('Number of Events');
s2Str(1) = {['mean: ' num2str(meanMinisRT)]};
s2Str(2) = {['median: ' num2str(medianMinisRT)]};
text(.7*riseTimeArray(iEndRT),.9*max(sum(minis,2)),s2Str);
hold off

subplot(2,2,3,'replace');
hold on
imagesc(amplitudeArray, riseTimeArray, minis);
set(gca,'YDir','normal');
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(end)+RTLim2D]);
title('minis 2D distribution');
xlabel('Amplitude(mV)');
if strcmpi(RTtype, '10-90%')
    ylabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    ylabel('20-80% rise times (ms)');
end
colorbar;
hold off

minis = round(minis);
nEvents = sum(sum(minis));
events = zeros(nEvents, 2);
bins = length(amplitudeArray)*length(riseTimeArray);
eventCounter = 1;
for iBin = 1:bins
    [row, col] = ind2sub([length(riseTimeArray) length(amplitudeArray)], iBin);
    events(eventCounter: eventCounter + minis(iBin)-1, :) = repmat([amplitudeArray(col) riseTimeArray(row)], minis(iBin), 1);
    eventCounter = eventCounter + minis(iBin);
end
subplot(2,2,4,'replace');
if ~isempty(events)
    hold on
    hist3([events(:,1) events(:,2)], {amplitudeArray riseTimeArray}, 'FaceAlpha', .65);
    set(gcf, 'renderer', 'opengl');
    xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
    ylim([riseTimeArray(1)-RTLim2D riseTimeArray(end)+RTLim2D]);
    title('2D distribution of simulated minis');
    xlabel('Amplitude(mV)');
    if strcmpi(RTtype, '10-90%')
        ylabel('10-90% rise times (ms)');
    elseif strcmpi(RTtype, '20-80%')
        ylabel('20-80% rise times (ms)');
    end
    grid on
    view(3);
    hold off
end



%% Drawing components of the minis source distribution:
if draw
    [eventMatrix, constMat] = minisDistribution(distributionType, distributionParameters, amplitudeArraySim, riseTimeArray(2:end)); %#ok<*UNRCH>
    eventMatrix = [zeros(1, size(eventMatrix,2)); round(eventMatrix)];
    eventMatrix = [zeros(size(eventMatrix,1), length(amplitudeArray)-length(amplitudeArraySim)) eventMatrix];
    constMat = cat(1, zeros(1, size(constMat,2), size(constMat,3)), round(constMat));
    constMat = cat(2, zeros(size(constMat,1), length(amplitudeArray)-length(amplitudeArraySim), size(constMat,3)), constMat);
    numStr = ['^s^t'; '^n^d'; '^r^d'; '^t^h'];
    
    figure(f5);
    set(f5, 'NumberTitle', 'off');
    set(f5, 'Name', 'Components of Minis Source Distribution');
    for iMat = 1:size(constMat,3)+2
        subplot(2,3,iMat,'replace');
        hold on
        if iMat < 5
            plotConst(constMat(:,:,iMat), amplitudeArray, riseTimeArray, ampLim2D, iEndAmp, RTLim2D, RTtype, strcat(num2str(iMat), numStr(iMat,:),...
                ' minis source distribution component (rounded)'));
        elseif iMat == 5
            plotConst(eventMatrix, amplitudeArray, riseTimeArray, ampLim2D, iEndAmp, RTLim2D, RTtype, 'Minis source distribution (non-smoothed)');
        else
            plotConst(round(tailMinis), amplitudeArray, riseTimeArray, ampLim2D, iEndAmp, RTLim2D, RTtype, 'Minis added tail (rounded)');
        end
        hold off
    end
    
    
    
    %% Drawing fft graphs:
    smoothWindow = 40;
    if length(noiseSpectrum) < length(targetSpectrum)
        targetSpectrum = targetSpectrum(1:length(noiseSpectrum));
    elseif length(targetSpectrum) < length(noiseSpectrum)
        noiseSpectrum = noiseSpectrum(1:length(targetSpectrum));
        freq = freq(1:length(targetSpectrum));
    end
    figure(f6);
    subplot(2,2,1,'replace');
    loglog(freq, noiseSpectrum);
    hold on
    loglog(freq, targetSpectrum,'r');
    filestring = 'Power spectra of data';
    set(f6, 'NumberTitle', 'off');
    set(f6, 'Name', 'Fourier Transform (time to frequency domain)');
    title(filestring);
    xlabel('Frequency (Hz)');
    ylabel('Power (mV^2/Hz)');
    legend('Simulated','Target','Location','NorthEast');
    hold off
    
    subplot(2,2,2,'replace');
    nSpecSmooth = gsmooth(noiseSpectrum, smoothWindow, 3);
    mnSpecSmooth = gsmooth(targetSpectrum, smoothWindow, 3);
    loglog(freq, abs(nSpecSmooth));
    hold on
    loglog(freq, abs(mnSpecSmooth),'r');
    filestring = 'Power spectra of data (smoothed)';
    title(filestring);
    xlabel('Frequency (Hz)');
    ylabel('Power (mV^2/Hz)');
    legend('Simulated','Target','Location','NorthEast');
    hold off
    
    subplot(2,2,3,'replace');
    loglog(freq, abs(targetSpectrum - noiseSpectrum), 'g');
    hold on
    filestring = 'Difference of power spectra';
    title(filestring);
    xlabel('Frequency (Hz)');
    ylabel('Power (mV^2/Hz)');
    legend('Difference modulus','Location','NorthEast');
    hold off
    
    subplot(2,2,4,'replace');
    diffSmooth = gsmooth(targetSpectrum - noiseSpectrum, smoothWindow, 3);
    loglog(freq, abs(diffSmooth), 'g');
    hold on
    filestring = 'Difference of power spectra (smoothed)';
    title(filestring);
    xlabel('Frequency (Hz)');
    ylabel('Power (mV^2/Hz)');
    legend('Difference modulus','Location','NorthEast');
    hold off
end



%% Drawing additional visualisations:
if evalVector(1)
    % Amplitude histogram line-ups:
    if isempty(f7)
        f7 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    end
    figure(f7);
    Amps = [histoData.Amps; sum(current,1)];
    Amp99prc = zeros(1,size(Amps,1));
    pAmp = Amp99prc;
    for iData = 1:size(Amps,1)
        Ampcs = cumsum(Amps(iData,:))/sum(Amps(iData,:));
        Amp99prc(iData) = find(Ampcs > .99, 1);
    end
    iEndAmp = min([max(Amp99prc)+1 length(amplitudeArray)]);
    plot([0 amplitudeArray(end)], [0 0], 'k-');
    hold on
    colour = ([zeros(size(Amps,1)-1,2) ones(size(Amps,1)-1,1)]).*repmat((1/(size(Amps,1)-1): 1/(size(Amps,1)-1) :1)',1,3);
    colour = [colour; [1 0 0]];
    for iData = 1:size(Amps,1)
        pAmp(iData) = plot(amplitudeArray, Amps(iData,:), '.-', 'color', colour(iData,:));
    end
    xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
    set(f7, 'NumberTitle', 'off');
    set(f7, 'Name', 'One-dimensional Amplitude Histograms (Real + Simulated)');
    title('Real and simulated amplitude distributions');
    xlabel('Amplitude(mV)');
    ylabel('Number of Events');
    legend([pAmp(floor(iData/2)) pAmp(iData)],'Real data','Simulated','Location','NorthEast');
    hold off
    
    
    
    % Rise time histogram line-ups:
    if isempty(f8)
        f8 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    end
    figure(f8);
    RTs = [histoData.RTs; sum(current,2)'];
    RT99prc = zeros(1,size(RTs,1));
    pRT= RT99prc;
    for iData = 1:size(RTs,1)
        RTcs = cumsum(RTs(iData,:))/sum(RTs(iData,:));
        RT99prc(iData) = find(RTcs > .99, 1);
    end
    iEndRT = min([max(RT99prc)+1 length(riseTimeArray)]);
    plot([0 riseTimeArray(end)], [0 0], 'k-');
    hold on
    colour = ([zeros(size(Amps,1)-1,2) ones(size(Amps,1)-1,1)]).*repmat((1/(size(Amps,1)-1): 1/(size(Amps,1)-1) :1)',1,3);
    colour = [colour; [1 0 0]];
    for iData = 1:size(RTs,1)
        pRT(iData) = plot(riseTimeArray, RTs(iData,:), '.-', 'color', colour(iData,:));
    end
    xlim([riseTimeArray(1) riseTimeArray(iEndRT)]);
    set(f8, 'NumberTitle', 'off');
    set(f8, 'Name', 'One-dimensional Rise Time Histograms (Real + Simulated)');
    title('Real and simulated amplitude distributions');
    xlabel('Amplitude(mV)');
    ylabel('Number of Events');
    legend([pRT(floor(iData/2)) pRT(iData)],'Real data','Simulated','Location','NorthEast');
    hold off
    
    
    
    
    
    % Single comparison optimisation summary plot:
    if isempty(f9)
        f9 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    end
    figure(f9);
    set(f9, 'NumberTitle', 'off');
    set(f9, 'Name', 'Single comparison optimisation summary');
    sp = subplot(2,8,1,'replace');
    set(sp,'position',[0.04 0.55 0.083 0.4]);
    hold on
    p1 = plot([0.5 1.5], [0.5 0.5], 'y', 'MarkerSize', 10);
    p2 = plot([0.5 1.5], [0.5 0.5], 'MarkerSize', 10);
    p3 = plot([0.5 1.5], [0.5 0.5], 'r', 'MarkerSize', 10);
    p4 = plot(1.1, 0.5, 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    p5 = plot(1.1, 0.5, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    xlim([0 2]);
    ylim([0 2]);
    title('Legend');
    ylabel('Event count');
    legend([p5 p4 p3 p2 p1],'Current','Min','Target','50cent','Range','Location','SouthWest');
    XLim = xlim;
    xpos = XLim(1) + (XLim(2)-XLim(1))*.5;
    YLim = ylim;
    ypos = YLim(2) - (YLim(2)-YLim(1))*.14;
    text(xpos, ypos, {'Single'}, 'color', 'k', 'HorizontalAlignment','center', 'Fontweight','bold');
    ypos = YLim(2) - (YLim(2)-YLim(1))*.21;
    text(xpos, ypos, {'Comparison'}, 'color', 'k', 'HorizontalAlignment','center', 'Fontweight','bold');
    ypos = YLim(2) - (YLim(2)-YLim(1))*.28;
    text(xpos, ypos, {'Optimisation'}, 'color', 'k', 'HorizontalAlignment','center', 'Fontweight','bold');
    ypos = YLim(2) - (YLim(2)-YLim(1))*.35;
    text(xpos, ypos, {'Summary'}, 'color', 'k', 'HorizontalAlignment','center', 'Fontweight','bold');
    hold off
    
    
    sp = subplot(2,8,2,'replace');
    set(sp,'position',[0.155 0.55 0.083 0.4]);
    hold on
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        plot([0.5 1.5], [measuredUF(21) measuredUF(21)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(33) measuredUF(33)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(9) measuredUF(9)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        plot([0.5 1.5], [measuredUF(13) measuredUF(13)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(25) measuredUF(25)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(1) measuredUF(1)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        plot([0.5 1.5], [measuredUF(17) measuredUF(17)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(29) measuredUF(29)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(5) measuredUF(5)], 'MarkerSize', 10);
    end
    plot([0.5 1.5], [costFuncStruct.boundVector(16) costFuncStruct.boundVector(16)], 'r', 'MarkerSize', 10);
    if bestEvalVector(1) > costFuncStruct.boundVector(16)
        plot(1, bestEvalVector(1), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(1), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(16) > costFuncStruct.boundVector(16)
        plot(1, evalVector(16), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(16), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('Combined SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Amplitude SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Rise time SAD');
    end
    hold off
    
    
    sp = subplot(2,8,3,'replace');
    set(sp,'position',[0.27 0.55 0.083 0.4]);
    hold on
    if strcmpi(costFuncStruct.firstCost, 'TwoDs') || strcmpi(costFuncStruct.firstCost, 'RTs')
        plot([0.5 1.5], [measuredUF(13) measuredUF(13)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(25) measuredUF(25)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(1) measuredUF(1)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        plot([0.5 1.5], [measuredUF(17) measuredUF(17)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(29) measuredUF(29)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(5) measuredUF(5)], 'MarkerSize', 10);
    end
    plot([0.5 1.5], [costFuncStruct.boundVector(17) costFuncStruct.boundVector(17)], 'r', 'MarkerSize', 10);
    if bestEvalVector(2) > costFuncStruct.boundVector(17)
        plot(1, bestEvalVector(2), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(2), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(17) > costFuncStruct.boundVector(17)
        plot(1, evalVector(17), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(17), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Amplitude SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Rise time SAD');
    end
    hold off
    
    
    sp = subplot(2,8,4,'replace');
    set(sp,'position',[0.385 0.55 0.083 0.4]);
    hold on
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        plot([0.5 1.5], [measuredUF(17) measuredUF(17)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(29) measuredUF(29)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(5) measuredUF(5)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'Amps') || strcmpi(costFuncStruct.firstCost, 'RTs')
        plot([0.5 1.5], [measuredUF(21) measuredUF(21)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(33) measuredUF(33)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(9) measuredUF(9)], 'MarkerSize', 10);
    end
    plot([0.5 1.5], [costFuncStruct.boundVector(18) costFuncStruct.boundVector(18)], 'r', 'MarkerSize', 10);
    if bestEvalVector(3) > costFuncStruct.boundVector(18)
        plot(1, bestEvalVector(3), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(3), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(18) > costFuncStruct.boundVector(18)
        plot(1, evalVector(18), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(18), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('Rise time SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Combined SAD');
    end
    hold off
    
    
    sp = subplot(2,8,5,'replace');
    set(sp,'position',[0.5 0.55 0.083 0.4]);
    hold on
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        plot([0.5 1.5], [measuredUF(22) measuredUF(22)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(34) measuredUF(34)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(10) measuredUF(10)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        plot([0.5 1.5], [measuredUF(14) measuredUF(14)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(26) measuredUF(26)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(2) measuredUF(2)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        plot([0.5 1.5], [measuredUF(18) measuredUF(18)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(30) measuredUF(30)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(6) measuredUF(6)], 'MarkerSize', 10);
    end
    plot([0.5 1.5], [costFuncStruct.boundVector(19) costFuncStruct.boundVector(19)], 'r', 'MarkerSize', 10);
    if bestEvalVector(4) > costFuncStruct.boundVector(19)
        plot(1, bestEvalVector(4), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(4), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(19) > costFuncStruct.boundVector(19)
        plot(1, evalVector(19), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(19), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('Combined MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Amplitude MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Rise time MAD');
    end
    hold off
    
    
    sp = subplot(2,8,6,'replace');
    set(sp,'position',[0.615 0.55 0.083 0.4]);
    hold on
    if strcmpi(costFuncStruct.firstCost, 'TwoDs') || strcmpi(costFuncStruct.firstCost, 'RTs')
        plot([0.5 1.5], [measuredUF(14) measuredUF(14)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(26) measuredUF(26)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(2) measuredUF(2)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        plot([0.5 1.5], [measuredUF(18) measuredUF(18)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(30) measuredUF(30)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(6) measuredUF(6)], 'MarkerSize', 10);
    end
    plot([0.5 1.5], [costFuncStruct.boundVector(20) costFuncStruct.boundVector(20)], 'r', 'MarkerSize', 10);
    if bestEvalVector(5) > costFuncStruct.boundVector(20)
        plot(1, bestEvalVector(5), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(5), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(20) > costFuncStruct.boundVector(20)
        plot(1, evalVector(20), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(20), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Amplitude MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Rise time MAD');
    end
    hold off
    
    
    sp = subplot(2,8,7,'replace');
    set(sp,'position',[0.73 0.55 0.083 0.4]);
    hold on
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        plot([0.5 1.5], [measuredUF(18) measuredUF(18)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(30) measuredUF(30)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(6) measuredUF(6)], 'MarkerSize', 10);
    elseif strcmpi(costFuncStruct.firstCost, 'Amps') || strcmpi(costFuncStruct.firstCost, 'RTs')
        plot([0.5 1.5], [measuredUF(22) measuredUF(22)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(34) measuredUF(34)], 'y', 'MarkerSize', 10);
        plot([0.5 1.5], [measuredUF(10) measuredUF(10)], 'MarkerSize', 10);
    end
    plot([0.5 1.5], [costFuncStruct.boundVector(21) costFuncStruct.boundVector(21)], 'r', 'MarkerSize', 10);
    if bestEvalVector(6) > costFuncStruct.boundVector(21)
        plot(1, bestEvalVector(6), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(6), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(21) > costFuncStruct.boundVector(21)
        plot(1, evalVector(21), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(21), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('Rise time MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Combined MAD');
    end
    hold off
    
    
    sp = subplot(2,8,8,'replace');
    set(sp,'position',[0.845 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFBottom(5) measuredUFBottom(5)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFBottom(9) measuredUFBottom(9)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFBottom(1) measuredUFBottom(1)], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(22) costFuncStruct.boundVector(22)], 'r', 'MarkerSize', 10);
    if bestEvalVector(7) > costFuncStruct.boundVector(22)
        plot(1, bestEvalVector(7), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(7), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(22) > costFuncStruct.boundVector(22)
        plot(1, evalVector(22), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(22), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top50% SAD');
    hold off
    
    
    sp = subplot(2,8,9,'replace');
    set(sp,'position',[0.04 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFBottom(1)/2-abs(measuredUFBottom(1)-measuredUFBottom(5))/2 measuredUFBottom(1)/2-abs(measuredUFBottom(1)-measuredUFBottom(5))/2], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFBottom(1)/2+abs(measuredUFBottom(1)-measuredUFBottom(9))/2 measuredUFBottom(1)/2+abs(measuredUFBottom(1)-measuredUFBottom(9))/2], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFBottom(1)/2 measuredUFBottom(1)/2], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(23) costFuncStruct.boundVector(23)], 'r', 'MarkerSize', 10);
    if bestEvalVector(8) > costFuncStruct.boundVector(23)
        plot(1, bestEvalVector(8), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(8), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(23) > costFuncStruct.boundVector(23)
        plot(1, evalVector(23), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(23), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top50% low SAD');
    ylabel('Event count');
    hold off
    
    
    sp = subplot(2,8,10,'replace');
    set(sp,'position',[0.155 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFBottom(6) measuredUFBottom(6)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFBottom(10) measuredUFBottom(10)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFBottom(2) measuredUFBottom(2)], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(24) costFuncStruct.boundVector(24)], 'r', 'MarkerSize', 10);
    if bestEvalVector(9) > costFuncStruct.boundVector(24)
        plot(1, bestEvalVector(9), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(9), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(24) > costFuncStruct.boundVector(24)
        plot(1, evalVector(24), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(24), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top50% MAD');
    hold off
    
    
    sp = subplot(2,8,11,'replace');
    set(sp,'position',[0.27 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFMid(5) measuredUFMid(5)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFMid(9) measuredUFMid(9)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFMid(1) measuredUFMid(1)], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(25) costFuncStruct.boundVector(25)], 'r', 'MarkerSize', 10);
    if bestEvalVector(10) > costFuncStruct.boundVector(25)
        plot(1, bestEvalVector(10), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(10), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(25) > costFuncStruct.boundVector(25)
        plot(1, evalVector(25), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(25), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top10% SAD');
    hold off
    
    
    sp = subplot(2,8,12,'replace');
    set(sp,'position',[0.385 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFMid(1)/2-abs(measuredUFMid(1)-measuredUFMid(5))/2 measuredUFMid(1)/2-abs(measuredUFMid(1)-measuredUFMid(5))/2], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFMid(1)/2+abs(measuredUFMid(1)-measuredUFMid(9))/2 measuredUFMid(1)/2+abs(measuredUFMid(1)-measuredUFMid(9))/2], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFMid(1)/2 measuredUFMid(1)/2], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(26) costFuncStruct.boundVector(26)], 'r', 'MarkerSize', 10);
    if bestEvalVector(11) > costFuncStruct.boundVector(26)
        plot(1, bestEvalVector(11), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(11), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(26) > costFuncStruct.boundVector(26)
        plot(1, evalVector(26), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(26), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top10% low SAD');
    hold off
    
    
    sp = subplot(2,8,13,'replace');
    set(sp,'position',[0.5 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFMid(6) measuredUFMid(6)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFMid(10) measuredUFMid(10)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFMid(2) measuredUFMid(2)], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(27) costFuncStruct.boundVector(27)], 'r', 'MarkerSize', 10);
    if bestEvalVector(12) > costFuncStruct.boundVector(27)
        plot(1, bestEvalVector(12), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(12), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(27) > costFuncStruct.boundVector(27)
        plot(1, evalVector(27), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(27), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top10% MAD');
    hold off
    
    
    sp = subplot(2,8,14,'replace');
    set(sp,'position',[0.615 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFTop(5) measuredUFTop(5)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFTop(9) measuredUFTop(9)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFTop(1) measuredUFTop(1)], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(28) costFuncStruct.boundVector(28)], 'r', 'MarkerSize', 10);
    if bestEvalVector(13) > costFuncStruct.boundVector(28)
        plot(1, bestEvalVector(13), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(13), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(28) > costFuncStruct.boundVector(28)
        plot(1, evalVector(28), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(28), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top2% SAD');
    hold off
    
    
    sp = subplot(2,8,15,'replace');
    set(sp,'position',[0.73 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFTop(1)/2-abs(measuredUFTop(1)-measuredUFTop(5))/2 measuredUFTop(1)/2-abs(measuredUFTop(1)-measuredUFTop(5))/2], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFTop(1)/2+abs(measuredUFTop(1)-measuredUFTop(9))/2 measuredUFTop(1)/2+abs(measuredUFTop(1)-measuredUFTop(9))/2], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFTop(1)/2 measuredUFTop(1)/2], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(29) costFuncStruct.boundVector(29)], 'r', 'MarkerSize', 10);
    if bestEvalVector(14) > costFuncStruct.boundVector(29)
        plot(1, bestEvalVector(14), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(14), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(29) > costFuncStruct.boundVector(29)
        plot(1, evalVector(29), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(29), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top2% low SAD');
    hold off
    
    
    sp = subplot(2,8,16,'replace');
    set(sp,'position',[0.845 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [measuredUFTop(6) measuredUFTop(6)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFTop(10) measuredUFTop(10)], 'y', 'MarkerSize', 10);
    plot([0.5 1.5], [measuredUFTop(2) measuredUFTop(2)], 'MarkerSize', 10);
    plot([0.5 1.5], [costFuncStruct.boundVector(30) costFuncStruct.boundVector(30)], 'r', 'MarkerSize', 10);
    if bestEvalVector(15) > costFuncStruct.boundVector(30)
        plot(1, bestEvalVector(15), 'v', 'MarkerEdgeColor',[1 0.5 0.2]);
    else
        plot(1, bestEvalVector(15), 'v', 'MarkerEdgeColor','g');
    end
    if evalVector(30) > costFuncStruct.boundVector(30)
        plot(1, evalVector(30), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(30), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top2% MAD');
    hold off
    
    
    
    
    
    % Group comparison optimisation summary plot:
    if isempty(f10)
        f10 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    end
    figure(f10);
    set(f10, 'NumberTitle', 'off');
    set(f10, 'Name', 'Group comparison optimisation summary');
    sp = subplot(2,8,1,'replace');
    set(sp,'position',[0.04 0.55 0.083 0.4]);
    hold on
    p1 = plot([0.5 1.5], [0.2 0.2], 'r', 'MarkerSize', 10);
    p2 = plot(1, 0.2, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    xlim([0 2]);
    ylim([0 2]);
    title('Legend');
    ylabel('Event count');
    legend([p2 p1],'Current','Target', 'Location','SouthWest');
    XLim = xlim;
    xpos = XLim(1) + (XLim(2)-XLim(1))*.5;
    YLim = ylim;
    ypos = YLim(2) - (YLim(2)-YLim(1))*.14;
    text(xpos, ypos, {'Group'}, 'color','k', 'HorizontalAlignment','center', 'Fontweight','bold');
    ypos = YLim(2) - (YLim(2)-YLim(1))*.21;
    text(xpos, ypos, {'Comparison'}, 'color','k', 'HorizontalAlignment','center', 'Fontweight','bold');
    ypos = YLim(2) - (YLim(2)-YLim(1))*.28;
    text(xpos, ypos, {'Optimisation'}, 'color','k', 'HorizontalAlignment','center', 'Fontweight','bold');
    ypos = YLim(2) - (YLim(2)-YLim(1))*.35;
    text(xpos, ypos, {'Summary'}, 'color','k', 'HorizontalAlignment','center', 'Fontweight','bold');
    hold off
    
    
    sp = subplot(2,8,2,'replace');
    set(sp,'position',[0.155 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(1) > 0
        plot(1, evalVector(1), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Amplitude SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Rise time SAD');
    end
    hold off
    
    
    sp = subplot(2,8,3,'replace');
    set(sp,'position',[0.27 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(2) > 0
        plot(1, evalVector(2), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(2), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Amplitude SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Rise time SAD');
    end
    hold off
    
    
    sp = subplot(2,8,4,'replace');
    set(sp,'position',[0.385 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(3) > 0
        plot(1, evalVector(3), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(3), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('Rise time SAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('SAD');
    end
    hold off
    
    
    sp = subplot(2,8,5,'replace');
    set(sp,'position',[0.5 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(4) > 0
        plot(1, evalVector(4), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(4), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Amplitude MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Rise time MAD');
    end
    hold off
    
    
    sp = subplot(2,8,6,'replace');
    set(sp,'position',[0.615 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(5) > 0
        plot(1, evalVector(5), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(5), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('Amplitude MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        title('Rise time MAD');
    end
    hold off
    
    
    sp = subplot(2,8,7,'replace');
    set(sp,'position',[0.73 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(6) > 0
        plot(1, evalVector(6), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(6), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        title('Rise time MAD');
    elseif strcmpi(costFuncStruct.firstCost, 'Amps') || strcmpi(costFuncStruct.firstCost, 'RTs')
        title('MAD');
    end
    hold off
    
    
    sp = subplot(2,8,8,'replace');
    set(sp,'position',[0.845 0.55 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(7) > 0
        plot(1, evalVector(7), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(7), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top50% SAD');
    hold off
    
    
    sp = subplot(2,8,9,'replace');
    set(sp,'position',[0.04 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(8) > 0
        plot(1, evalVector(8), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(8), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top50% low SAD');
    ylabel('Event count');
    hold off
    
    
    sp = subplot(2,8,10,'replace');
    set(sp,'position',[0.155 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(9) > 0
        plot(1, evalVector(9), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(9), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top50% MAD');
    hold off
    
    
    sp = subplot(2,8,11,'replace');
    set(sp,'position',[0.27 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(10) > 0
        plot(1, evalVector(10), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(10), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top10% SAD');
    hold off
    
    
    sp = subplot(2,8,12,'replace');
    set(sp,'position',[0.385 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(11) > 0
        plot(1, evalVector(11), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(11), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top10% low SAD');
    hold off
    
    
    sp = subplot(2,8,13,'replace');
    set(sp,'position',[0.5 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 13);
    if evalVector(12) > 0
        plot(1, evalVector(12), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(12), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top10% MAD');
    hold off
    
    
    sp = subplot(2,8,14,'replace');
    set(sp,'position',[0.615 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(13) > 0
        plot(1, evalVector(13), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(13), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top2% SAD');
    hold off
    
    
    sp = subplot(2,8,15,'replace');
    set(sp,'position',[0.73 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(14) > 0
        plot(1, evalVector(14), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(14), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top2% low SAD');
    hold off
    
    
    sp = subplot(2,8,16,'replace');
    set(sp,'position',[0.845 0.04 0.083 0.4]);
    hold on
    plot([0.5 1.5], [0 0], 'r', 'MarkerSize', 10);
    if evalVector(15) > 0
        plot(1, evalVector(15), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    else
        plot(1, evalVector(15), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    end
    xlim([0 2]);
    title('Amplitude top2% MAD');
    hold off
else
    f7 = [];
    f8 = [];
    f9 = [];
    f10 = [];
end





function plotConst(constMat, amplitudeArray, riseTimeArray, ampLim2D, iEndAmp, RTLim2D, RTtype, titleStr)

if sum(sum(constMat)) < 0
    binColour = 'r';
    constMat = abs(constMat);
end
nEvents = sum(sum(constMat));
events = zeros(nEvents, 2);
bins = length(amplitudeArray)*length(riseTimeArray);
eventCounter = 1;
for iBin = 1:bins
    [row, col] = ind2sub([length(riseTimeArray) length(amplitudeArray)], iBin);
    events(eventCounter: eventCounter + constMat(iBin)-1, :) = repmat([amplitudeArray(col) riseTimeArray(row)], constMat(iBin), 1);
    eventCounter = eventCounter + constMat(iBin);
end
title(titleStr);
if ~isempty(events)
    if exist('binColour','var')
        hist3([events(:,1) events(:,2)],{amplitudeArray riseTimeArray}, 'FaceAlpha',.5, 'FaceColor',binColour);
    else
        hist3([events(:,1) events(:,2)],{amplitudeArray riseTimeArray}, 'FaceAlpha',.65);
    end
    set(gcf, 'renderer', 'opengl');
    xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
    ylim([riseTimeArray(1)-RTLim2D riseTimeArray(end)+RTLim2D]);
    xlabel('Amplitude(mV)');
    if strcmpi(RTtype, '10-90%')
        ylabel('10-90% rise times (ms)');
    elseif strcmpi(RTtype, '20-80%')
        ylabel('20-80% rise times (ms)');
    end
    grid on
    view(3);
end