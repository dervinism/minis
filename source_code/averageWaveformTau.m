function [waveform, parameters] = averageWaveformTau(minis, V, dt, SWduration, excludedT, tBLduration, classificationParameters, peelDuration, draw)
% AVERAGEWAVEFORM is a helper subfunction of detect minis. It estimates the
% average waveform of the minis detected using the detectMinis program.
%
%   [WAVEFORM PARAMETERS] = AVERAGEWAVEFORM(minis, V, dt, SWduration,...
%       excludedT, tBLduration, peelDuration, cutoff, draw)
%   Estimates the average waveform. MINIS is the matrix produced by the
%   minis detection algorithm of the detectMinis program. V is the
%   electrophysiological recording trace in milivolts (mV) or nanoamperes
%   (nA). DT is the sampling interval in miliseconds (ms). SWDURATION is
%   the beginning of the search window that corresponds to
%   searchParameters.SWstart as used by the minis detection algorithm, ms.
%   EXCLUDEDT is the vector of excluded times. TBLDURATION is the duration
%   of the baseline period, ms. PEELDURATION is a scalar variable
%   representing the duration of the period used for estimating the
%   membrane time constant, tau, in ms (or, in other words, for fitting the
%   regression line to the logarithm of the minis decay phase). DRAW is a
%   string variable controling the dislpay of the average waveform and
%   regression line fitting fingures. Setting 'on' turns on the figure
%   displays, whereas setting to 'off' turns them off.
%   WAVEFORM is a vector containg the average waveform, mV or nA.
%   PARAMETERS is a structure variable containing fields with average
%   waveform parameters 'peak' (mV or nA), 'BL' (baseline, mV or nA), 'Amp'
%   (amplitude, mV or nA), 'risetime' (ms), 'tau' (dendritic membrane time
%   constant, ms).
%
%   [WAVEFORM PARAMETERS F] = AVERAGEWAVEFORM(...)
%   In addition outputs the handles of the average waveform and regression
%   line fitting fingures. For this option to work, the user has to set the
%   input DRAW variable to 'on'.
%


dAmp = (classificationParameters.amplitudeArray(2) - classificationParameters.amplitudeArray(1));
amplitudeArray = classificationParameters.amplitudeArray(1): dAmp :classificationParameters.amplitudeArray(end);
tau = zeros(1,length(amplitudeArray));
for iAmpBin = 1:length(amplitudeArray)
    if iAmpBin ~= 1
        loBound = amplitudeArray(iAmpBin) - dAmp/2;
        upBound = amplitudeArray(iAmpBin) + dAmp/2;
    else
        loBound = 0;
        upBound = amplitudeArray(iAmpBin) + dAmp/2;
    end
    minisBin = minis(minis(:,4)<upBound & minis(:,4)>=loBound,:);
    if ~isempty(minisBin)
        standardWaveformIndices = 0:12*round(SWduration/dt)-1;
        standardWaveformIndices = repmat(standardWaveformIndices, length(minisBin(:,17)), 1);
        waveformStartIndices = minisBin(:,17) - 3*SWduration/dt;
        
        % Eliminate minis that contain excluded times:
        midWave = ceil(size(standardWaveformIndices,2)/2);
        waveformMidIndices = waveformStartIndices + midWave;
        rowsToExclude = zeros(length(waveformMidIndices), 1);
        excludedInd = [1 round(excludedT/dt)];
        for iRow = 1:length(rowsToExclude)
            distances = abs(excludedInd - waveformMidIndices(iRow));
            if min(distances) <= midWave
                rowsToExclude(iRow) = 1;
            end
        end
        waveformStartIndices(logical(rowsToExclude)) = [];
        
        % Estimate the average waveform:
        waveformArray = repmat(waveformStartIndices, 1, size(standardWaveformIndices, 2));
        waveformArray = waveformArray + standardWaveformIndices(1:length(waveformStartIndices),:);
        % waveformArray = V(waveformArray) - repmat(minis(logical(~rowsToExclude),20),1,size(waveformArray, 2)); % subtract the 'baseline'
        waveformArray = V(waveformArray);
        if size(waveformArray,1) > 1
            waveform = mean(waveformArray);
        else
            waveform = waveformArray;
        end
        waveTimeline = 0:dt:length(waveform)*dt-dt;
        
        % Estimate the peak, the baseline, and the amplitude:
        [parameters.peak, iPeak] = max(waveform(1:midWave));
        [~, iTrough] = min(waveform(1:iPeak));
        iBLend = iTrough + round(.5*tBLduration/dt);
        iBLstart = 1;
        parameters.BL = mean(waveform(iBLstart:iBLend));
        [parameters.risetime, ~, i10] = calcRT(waveTimeline(1:iPeak), waveform(1:iPeak), parameters.BL);
        % Shift the baseline closer to the 10% rise time mark if necessary:
        if i10 - iBLend > ceil(.75*tBLduration/dt)
            iSWstart = i10 - round(tBLduration/dt);
            [~, iTrough] = min(waveform(iSWstart:iPeak));
            iTrough = iSWstart - 1 + iTrough;
            iBLend = min([iTrough + round(.5*tBLduration/dt) iPeak]);
            parameters.BL = mean(waveform(iBLstart:iBLend));
            [parameters.risetime] = calcRT(waveTimeline(1:iPeak), waveform(1:iPeak), parameters.BL);
        end
        parameters.Amp = parameters.peak - parameters.BL;
        
        % Calculate the rise-time:
        [parameters.risetime, t10, i10, t50, i50, t90, i90] = calcRT(waveTimeline(1:iPeak), waveform(1:iPeak), parameters.BL);
        
        % Visualise the waveform:
        if strcmpi(draw,'on')
            f1 = figure('position', [50 50 1200 600], 'Name', 'Averaged Mini-like Events', 'NumberTitle', 'off');
            plot(waveTimeline, waveform, 'k', 'markersize', 4);
            hold on
            plot(waveTimeline(iPeak), waveform(iPeak), '^', 'markersize', 10, 'markeredgecolor', 'r');
            plot(t10, waveform(i10), 'r.', 'markersize', 10);
            plot(t50, waveform(i50), 'r.', 'markersize', 10);
            plot(t90, waveform(i90), 'r.', 'markersize', 10);
            plot(waveTimeline(iBLstart), parameters.BL, 'c.', 'markersize', 10);
            plot(waveTimeline(iBLend), parameters.BL, 'b.', 'markersize', 10);
            title('Smoothed Data');
            xlabel('Time (ms)');
            ylabel('Voltage (mV)');
            hold off
        end
        
        
        % Estimate the membrane time constant (tau):
        subtrWaveform = waveform(iPeak:end) - parameters.BL;
        iStartPeel = 1;
        iEndPeel = round(peelDuration/dt);
        lnPeel = real(log(subtrWaveform));
        lnTimeline = (0:1:length(lnPeel)-1)*dt;
        if exist('draw','var') && strcmpi(draw,'on')
            f2 = figure('position', [50 50 1200 600], 'Name', 'Natural Logarithm of the Decay Period', 'NumberTitle', 'off');
            plot(lnTimeline, lnPeel, 'k.', 'markersize', 6);
            hold on
            title('Smoothed Data');
            xlabel('Time (ms)');
            ylabel('Voltage (mV)');
            varargout = {[f1 f2]}; %#ok<*NASGU>
        end
        
        independent = [lnTimeline(iStartPeel:iEndPeel)' ones(iEndPeel-iStartPeel+1, 1)];
        coefficients = regress(lnPeel(iStartPeel:iEndPeel)', independent);
        regressionLine = coefficients(1)*lnTimeline(iStartPeel:iEndPeel) + coefficients(2);
        if exist('draw','var') && strcmpi(draw,'on')
            plot(lnTimeline(iStartPeel:iEndPeel), regressionLine, 'r');
            hold off
        end
        tau(iAmpBin) = -1/coefficients(1);
    else
        tau(iAmpBin) = 0;
    end
    clear waveformArray
end
if strcmpi(draw,'on')
    figure('position', [50 50 1200 600], 'Name', 'Range of Taus', 'NumberTitle', 'off');
    plot(amplitudeArray, tau, 'k', 'markersize', 4);
    hold on
end

% Find the relationship between the amplitude size and tau:
iStart = find(amplitudeArray >= .07, 1);
iEnd = find(amplitudeArray <= .24, 1, 'last');
independent = [amplitudeArray(iStart:iEnd)' ones(iEnd-iStart+1, 1)];
coefficients = regress(tau(iStart:iEnd)', independent);
regressionLine = coefficients(1)*amplitudeArray(iStart:iEnd) + coefficients(2);
if strcmpi(draw,'on')
    plot(amplitudeArray(iStart:iEnd), regressionLine, 'r');
    hold off
end

tauEst = coefficients(2) + coefficients(1)*amplitudeArray;
parameters.tau = mean(tauEst);
unrealTau = find(tauEst<=5,1);
if ~isempty(unrealTau)
    tauEst(unrealTau:end) = tauEst(unrealTau-1)*ones(1,length(tauEst)-unrealTau+1);
end
parameters.tauEst = tauEst;
%parameters.tauEst = 6:20/length(tauEst):26;