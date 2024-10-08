function [tau, effTau, capacitance, capacitanceEff, R_s] = peel(sign, pulseDuration, V, I, nSweeps, dt, draw, startTime1, endTime1,...
    startTime2, endTime2, varargin)
% PEEL estimates the dendritic membrane time constant tau given the
% dendritic membrane potential recording data.
%
%   [TAU CAPACITANCE R_S] = PEEL(SIGN, PEELDURATION)
%   estimates the dendritic membrane time constant tau , total dendritic
%   capacitance, and series resistance given three estimation parameters.
%   SIGN is a scalar value indicating whether to assess a positive or
%   negative pulse for estimating tau. 1 indicates positive and -1
%   indicates negative. PEELDURATION is the duration of the peeling period
%   in ms. V is an actual recording data
%   vector in millivolts (mV) or nanoamperes (nA). NSWEEPS is the number of
%   data sweeps comprising the entire data vetor V. DT is the data sampling
%   interval in ms. DATAPROPERTIES is a structure variable that is output
%   by the LOADABF function. DRAW is a string variable that can be set to
%   'on' for displaying figures or 'off' otherwise. The default state is
%   'off'. STARTTIME1 is a scalar that sets the start time of the recording
%   region containing the first pulse. The pulse should be a long one used
%   for estimating the pseudo series resitance. If you do not wish to
%   estimate the series resistance, set this and the following input
%   variables to be empty. ENDTIME1 is a scalar that sets the end time of
%   the recording region containing the first pulse. The pulse should be a
%   long one used for estimating the pseudo series resitance. If you do not
%   wish to estimate the series resistance, set this and the previous input
%   variables to be empty. STARTTIME2 is a scalar that sets the start time
%   of the recording region containing the second pulse. The pulse should
%   be a short one used for estimating the passive membrane time constant
%   and capacitance. ENDTIME2 is a scalar that sets the end time of the
%   recording region containing the second pulse. The pulse should be a
%   short one used for estimating the passive membrane time constant and
%   capacitance.
%   The output variable TAU is a scalar representing the passive dendritic
%   membrane time constant in milliseconds (ms). EFFTAU is a scalar
%   representing the effective passive dendritic membrane time constant.
%   CAPACITANCE is a scalar variable representing the total membrane
%   capacitance of the recorded neuron in picofarads (pF). R_S is a scalar
%   corresponding to the series resistance in Ohms. This output variable is
%   set to zero if the STARTTIME1 is empty.
%


% persistent figNum
%% Allign the sweeps:
sweepDuration = floor(length(V)/nSweeps);
stackedSweeps = sign*reshape(V(1:sweepDuration*nSweeps), sweepDuration, nSweeps)';
if isempty(endTime2)
    endTime2 = size(stackedSweeps,2)*dt - dt;
end
stackedCurrents = sign*reshape(I(1:sweepDuration*nSweeps), sweepDuration, nSweeps)';
stackedCurrents2 = stackedCurrents(:, round(startTime2/dt)+1 : round(endTime2/dt)+1);

midPoint = ceil(size(stackedCurrents2,2)/3);
if sign > 0
    forwardShift = stackedCurrents2(:,1:midPoint) - circshift(stackedCurrents2(:,1:midPoint), [0 1]);
elseif sign < 0
    forwardShift = stackedCurrents2(:,midPoint+1:end) - circshift(stackedCurrents2(:,midPoint+1:end), [0 1]);
end
[~, currentStart] = max(forwardShift,[],2);
iStart = currentStart - min(currentStart) + 1;
iEnd = repmat(size(stackedSweeps,2),size(stackedSweeps,1),1) - max(iStart) + iStart;

shiftedSweeps = zeros(size(stackedSweeps,1),size(stackedSweeps,2)-max(iStart)+1);
shiftedCurrents = shiftedSweeps;
for iTrace = 1:size(stackedSweeps,1)
    shiftedSweeps(iTrace,:) = stackedSweeps(iTrace, iStart(iTrace):iEnd(iTrace));
    shiftedCurrents(iTrace,:) = stackedCurrents(iTrace, iStart(iTrace):iEnd(iTrace));
end



%% Average the sweeps for the second pulse:
shiftedSweeps2 = shiftedSweeps(:, round(startTime2/dt)+1 : round(endTime2/dt)+1);
shiftedCurrents2 = shiftedCurrents(:, round(startTime2/dt)+1 : round(endTime2/dt)+1);
averageSweep = mean(shiftedSweeps2,1);
averageCurrent = mean(shiftedCurrents2,1);



%% Plot the data:
t = dt*(1:length(averageSweep)) - dt;
if draw
    figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83], 'Name', 'Averaged Raw Data', 'NumberTitle', 'off');
    plot(t, averageSweep, 'k', 'markersize', 4);
    hold on
    title('Raw Data');
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
end



%% Estimate the injected current and injected pulse duration:
if sum(averageCurrent)
    if sign > 0
        forwardShift = averageCurrent(:,1:midPoint) - circshift(averageCurrent(:,1:midPoint), [0 1]);
        backwardShift = averageCurrent(:,1:midPoint) - circshift(averageCurrent(:,1:midPoint), [0 -1]);
    elseif sign < 0
        forwardShift = averageCurrent(:,midPoint+1:end) - circshift(averageCurrent(:,midPoint+1:end), [0 1]);
        backwardShift = averageCurrent(:,midPoint+1:end) - circshift(averageCurrent(:,midPoint+1:end), [0 -1]);
    end
    [~, currentStart] = max(forwardShift);
    currentStart = currentStart - 1;
    [~, currentEnd] = max(backwardShift);
    blStart = currentStart - (currentEnd - (currentStart-1));
    bl = blStart : currentStart-1;
    if nargin == 13
        injCurr = varargin{1};
    else
        if sign > 0
            injCurr = mean(averageCurrent(currentStart+1:currentEnd)) - mean(averageCurrent(bl));
        elseif sign < 0
            injCurr = mean(averageCurrent(midPoint+currentStart+1:midPoint+currentEnd)) - mean(averageCurrent(bl));
        end
    end
else
    injCurr = 0;
end



%% Estimate tau_m:
% Find and mark the first peak:
[peak, iPeak] = max(averageSweep);
if exist('draw','var') && strcmpi(draw,'on')
    plot(t(iPeak), peak, '^', 'markersize', 10, 'markeredgecolor', 'g');
end

% Locate the baseline:
decayTime = 50;
endBL = iPeak*dt-pulseDuration;
iEndBL = round(endBL/dt);
startBL = endBL - 12;
iStartBL = round(startBL/dt);

% Subtract and plot the baseline:
BL = mean(averageSweep(iStartBL:iEndBL));
subtractedTrace = averageSweep(iPeak: iPeak + round(decayTime/dt)) - BL;
subtractedTimeline = t(iPeak: iPeak + round(decayTime/dt));
if draw
    plot(t(iStartBL:iEndBL), BL, 'r', 'markersize', 4);
    hold off
end

% Locate the peeling period:
lnPeel = log(subtractedTrace);
diffLnPeel = diff(lnPeel);
smoothDiff = smooth(diffLnPeel, round(1/dt));
indexArray = 1:round(1/dt);
indexArray = repmat(indexArray, length(diffLnPeel)-round(1/dt)+1, 1);
indexArray = indexArray + repmat((0:length(diffLnPeel)-round(1/dt))', 1, size(indexArray,2));
stdDiff = std(diffLnPeel(indexArray), 1, 2)';
iStartPeel = max([round(3/dt) find(smoothDiff>=1.05*mean(diffLnPeel(round(2/dt)-1:round(25/dt)-1)),1) + round(.5/dt) + 1]);
if isempty(iStartPeel)
    iStartPeel = round(4/dt);
end
minStdDiff = min(stdDiff(iStartPeel - round(.5/dt) : end));
iEndPeel = iStartPeel + find(stdDiff(iStartPeel+1:end)>4*minStdDiff,1) + round(.5/dt);

% Plot the peeling period:
if draw
    figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83], 'Name', 'Averaged Peeling Period', 'NumberTitle', 'off');
    plot(subtractedTimeline, subtractedTrace, 'k', 'markersize', 4);
    hold on
    title('Raw Data');
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    hold off
end

% Plot the natural logarithm of the peeling period:
lnTimeline = 0: dt :(length(lnPeel)-1)*dt;
if draw
    figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83], 'Name', 'Natural Logarithm of the Peeling Period', 'NumberTitle', 'off');
    plot(lnTimeline, lnPeel, 'k.', 'markersize', 6);
    hold on
    title('Raw Data');
    xlabel('Time (ms)');
    ylabel('ln(voltage (mV) - baseline voltage (mV))');
end

% Fit and plot the regression line:
independent = [lnTimeline(iStartPeel:iEndPeel)' ones(iEndPeel-iStartPeel+1, 1)];
coefficients = regress(lnPeel(iStartPeel:iEndPeel)', independent);
regressionLine = coefficients(1)*t(iStartPeel:iEndPeel) + coefficients(2);
if draw
    plot(lnTimeline(iStartPeel:iEndPeel), regressionLine, 'r');
    hold off
    %     if isempty(figNum)
    %         figNum = 1;
    %     else
    %         figNum = figNum + 1;
    %     end
    %     saveas(gcf, num2str(figNum), 'fig');
    %     close(gcf);
end

% Estimate tau_m:
tau = -1/coefficients(1);



%% Estimate effective tau_m
waveformNorm = subtractedTrace(iStartPeel:end);
waveformNorm = gsmooth(waveformNorm, round(2/dt), 3);
iStart = 1;
Vstart = waveformNorm(iStart);
Vend = Vstart*1/exp(1);
iEnd = find(waveformNorm <= Vend, 1);
if isempty(iEnd)
    iEnd = length(waveformNorm);
end
effTau = (iEnd - iStart)*dt;                                                % Effective tau_m (1/e^th life)



%% Estimate capacitance:
% Ref: Major et al., 1993 Solutions for Transients in Arbitrarily Branching
% Cables: I. Voltage Recording with a Somatic Shunt, eq36:
if injCurr
    amplitude = exp(coefficients(2));
    A0_unit_impulse = amplitude/(injCurr*tau*(1 - exp(-pulseDuration/tau)));    % mV/(nA*ms) = mV/pC
    capacitance = 1/(0.001*A0_unit_impulse);                                    % pF = us/MOhm
    
    A0_unit_impulseEff = amplitude/(injCurr*tau*(1 - exp(-pulseDuration/effTau)));
    capacitanceEff = 1/(0.001*A0_unit_impulseEff);
else
    capacitance = [];
    capacitanceEff = [];
end



%% Estimate bridge balance error:
if ~isempty(startTime1)
    if isempty(endTime1)
        endTime1 = startTime2;
    end
    shiftedSweeps1 = shiftedSweeps(:, round(startTime1/dt)+1 : round(endTime1/dt)+1);
    shiftedCurrents1 = shiftedCurrents(:, round(startTime1/dt)+1 : round(endTime1/dt)+1);
    currArray = max(shiftedCurrents1,[],2);
    [~, iCurr] = max(currArray);
    largestTrace = shiftedCurrents1(iCurr,:);
    forwardShift = largestTrace - circshift(largestTrace, [0 1]);
    backwardShift = largestTrace - circshift(largestTrace, [0 -1]);
    [~, currentStart] = max(forwardShift);
    currentStart = max([1 currentStart - 1]);
    [~, currentEnd] = max(backwardShift);
    currentEnd = min([currentEnd, length(backwardShift)]);
    
    % Estimate voltage amplitudes:
    ampStart = currentStart + 2;
    ampEnd = ampStart + 2;
    amplitudes = shiftedSweeps1(:,ampEnd) - shiftedSweeps1(:,ampStart);
    
    % Estimate current amplitudes:
    currPulses = mean(shiftedCurrents1(:,currentStart:currentEnd),2);
    blStart = currentStart - (currentEnd - (currentStart-1));
    bl = blStart : currentStart-1;
    if blStart <= 0
        blStart = currentEnd + 1;
        bl = blStart : currentEnd + (currentEnd - (currentStart-1));
        if bl(end) > size(shiftedCurrents1,2)
            blPre = currentStart - 1;
            blPost = size(shiftedCurrents1,2) - currentEnd;
            if blPre >= blPost
                bl = 1 : currentStart-1;
            else
                bl = currentEnd + 1 : size(shiftedCurrents1,2);
            end
        end
    end
    currBaselines = mean(shiftedCurrents1(:,bl),2);
    currAmplitudes = currPulses - currBaselines;
    
    % Estimate pseudo series resistance (R_series):
    R_s = mean(amplitudes./currAmplitudes);                                 % Megaohms
else
    R_s = [];
end