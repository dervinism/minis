function [simV, vNoiseMinis, eventMatrix, shapes] = simulateMinis(baseline, parameters, distributionType, distributionParameters, dt, basis, sweep,...
    excludedTimes, smoothWindow, parallelCores, amplitudeArray, amplitudeArraySim, riseTimeArray)






%% Simulate a sample of minis:
tDuration = 20; % duration of a simulated mini, tau
[V, ~, sitesLengths] = sampleMinis(parameters, tDuration, dt, parallelCores, riseTimeArray(2:end));
clear parameters tDuration





%% Obtaining minis distribution:
if size(distributionParameters,1) == 1
    eventMatrix = minisDistribution(distributionType, distributionParameters, amplitudeArraySim, riseTimeArray(2:end));
    eventMatrix = [zeros(1, size(eventMatrix,2)); round(eventMatrix)];
    if sum(sum(eventMatrix)) == 0
        eventMatrix(1,1) = 1;
    end
    
    if strcmpi(baseline, 'Subtracted')
        if size(basis,2) - size(eventMatrix,2) > 0
            eventMatrix = [zeros(size(eventMatrix,1), size(basis,2) - size(eventMatrix,2)) eventMatrix];
            basis(:,1:size(basis,2) - size(eventMatrix,2)) = zeros(size(eventMatrix,1),...
                size(basis,2) - size(eventMatrix,2));
            eventMatrix = basis + eventMatrix;
            eventMatrix = round(smoothMinis(eventMatrix, smoothWindow, 3));
        else
            eventMatrix = eventMatrix + basis;
        end
    else
        if basis - size(eventMatrix,2) > 0
            eventMatrix = [zeros(size(eventMatrix,1), basis - size(eventMatrix,2)) eventMatrix];
            eventMatrix = round(smoothMinis(eventMatrix, smoothWindow, 3));
        end
    end
elseif size(distributionParameters,1) > 1
    eventMatrix = distributionParameters;
end
clear distributionParameters baseline parallelCores distributionType amplitudeArraySim





%% Simulating minis on top of the noise:
t = dt*(1:1:length(sweep)) - dt;
[simV, ~, shapes] = positionMinis(V, t, excludedTimes, eventMatrix, amplitudeArray, riseTimeArray, sitesLengths);
vNoiseMinis = sweep + simV;
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%                            Local functions                              %
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, riseTimeArray, varargout] = sampleMinis(parameters, tDuration, dt, parallelCores, varargin)
% The local function sampleMinis simulates a sample of minis based on a
% user supplied uni-dimensional amplitude distribution or a two-dimensional
% amplitude and 10-90% rise time distribution.

% Obtaining certain dendritic parameters that are not contained in the input:
tau_m = parameters.tau_m; % passive membrane time constant, ms
R_m = tau_m*parameters.C_m*1000;
lambda = sqrt((R_m*parameters.d)/(4*parameters.R_i)); % length constant of cylindrical core conductor, cm

% Setting up single mini simulation parameters:
simPoints = round(tDuration*tau_m/dt); % number of data points composing the simulated mini
t = dt:dt:dt*simPoints; % range of the time variable for simulating a single mini, ms

% Simulating samples of minis:
riseTimeArray = varargin{1};
V = zeros(length(riseTimeArray),length(t)+1);
inputSites = zeros(length(riseTimeArray),1);
dendriticLengths = zeros(length(riseTimeArray),1);
L = parameters.L;
t10 = zeros(length(riseTimeArray),1);
t90 = zeros(length(riseTimeArray),1);
if parallelCores == 1
    for iRT = 1:length(riseTimeArray)
        [V(iRT,:), t10(iRT), t90(iRT), inputSites(iRT), dendriticLengths(iRT)] = singleMini(parameters, tau_m, lambda, riseTimeArray(iRT), t, L);
    end
else
    parfor iRT = 1:length(riseTimeArray)
        [V(iRT,:), t10(iRT), t90(iRT), inputSites(iRT), dendriticLengths(iRT)] = singleMini(parameters, tau_m, lambda, riseTimeArray(iRT), t, L); %#ok<PFOUS>
    end
end
V(11:end,1:3) = zeros(length(riseTimeArray)-10,3);
V = [zeros(1,size(V,2)); V];
varargout = {[[0; inputSites] [0; dendriticLengths]]};

% Drawing:
draw = 0;
if draw == 1
    f1 = figure('position', [5 300 600 300]);
    figure(f1);
    for iRT = 1:length(riseTimeArray)
        plot([0 t], V(iRT,:),'b-');
        hold on
        plot(t10(iRT), .1,'b.');
        plot(t90(iRT), .9,'b.');
    end
    title('Minis shapes');
    xlabel('time, ms');
    ylabel('membrane potential, mV');
    hold off
end
end



function [V, t10, t90, inputSite, dendriticLength] = singleMini(parameters, tau_m, lambda, iriseTimeArray, t, L)
L = max([L parameters.L]);
dX = L/2;
X = dX;
tau_sy1 = parameters.tau_sy1;
tau_sy2 = parameters.tau_sy2;
riseTime = 0;
bisectX = 1;
stepping = 0;
%f1 = figure('position', [5 300 600 300]);
while bisectX == 1
    while riseTime < iriseTimeArray
        [V,riseTime,t10,t90] = finiteCableLumped(parameters.nSeries, L, lambda, parameters.d, parameters.C_m, t, tau_m, X, parameters.Q, tau_sy1, tau_sy2);
        %         figure(f1);
        %         plot([0 t], V,'g-');
        %         hold on
        %         plot(t10, .1,'g.');
        %         plot(t90, .9,'g.');
        if ~stepping && tau_sy1 == parameters.tau_sy1 && riseTime > iriseTimeArray && X ~= 0
            X = 0;
            riseTime = 0;
        elseif tau_sy1 == parameters.tau_sy1 && riseTime > iriseTimeArray && X == 0 && L <= 3.5
            L = 1.15*L;
            riseTime = 0;
        elseif riseTime > iriseTimeArray && X == 0 && L > 3.5
            tau_sy1 = .8*tau_sy1;
            tau_sy2 = .8*tau_sy2;
            riseTime = 0;
        elseif X == L && riseTime < iriseTimeArray
            L = 1.15*L;
            dX = L/2;
            X = dX;
            riseTime = 0;
        elseif riseTime < iriseTimeArray
            X = X + dX;
            stepping = 1;
        end
    end
    
    riseTime = round(riseTime*10^15)/10^15;
    if riseTime == iriseTimeArray
        bisectX = 0;
        inputSite = X;
        dendriticLength = L;
        X = 0;
    elseif exist('previousX','var') && abs(X - previousX) <= 0.01
        X = max([0 X - dX/2]);
        [V, ~, t10, t90] = finiteCableLumped(parameters.nSeries, L, lambda, parameters.d, parameters.C_m, t, tau_m, X, parameters.Q, tau_sy1, tau_sy2);
        %         figure(f1);
        %         plot([0 t], V,'b-');
        %         hold on
        %         plot(t10, .1,'b.');
        %         plot(t90, .9,'b.');
        bisectX = 0;
        inputSite = max([0 X - dX/2]);
        dendriticLength = L;
        X = 0;
    else
        riseTime = 0;
    end
    X = max([X - dX 0]);
    dX = dX/2;
    previousX = X;
end
end



function [simV, eventMatrix, shapes] = positionMinis(V, t, excludedT, nEvents, amplitudeArray, riseTimeArray, varargin)
simV = zeros(1,length(t));
dt = t(2) - t(1);
adjt = t;
adjt(round(excludedT/dt)+1) = [];
sitesLengths = varargin{1};
nMinis = sum(sum(nEvents)); % the number of minis that will be added over the total recording sweep
shapes = zeros(nMinis,9); % for the purpose of keeping a record of which particular waveforms were chosen
countMinis = 0;
for iAmplitude = 1:length(amplitudeArray)
    nMinis_amp = sum(nEvents(:,iAmplitude)); % the number of minis of particular amplitude that will be added over the total recording sweep
    cumnMinis = cumsum(nEvents(:,iAmplitude));
    for iMini = 1:nMinis_amp
        onsetMinis = round(length(adjt)*rand); % uniform probability of start time along the file
        if onsetMinis == 0
            onsetMinis = length(adjt);
        end
        onsetMinisT = adjt(onsetMinis);
        onsetMinis = round(onsetMinisT/dt);
        countMinis = countMinis + 1;
        iMinisComparisonArray = iMini*ones(length(riseTimeArray),1);
        compareBins = le(iMinisComparisonArray, cumnMinis);
        rtBin = length(riseTimeArray) - sum(compareBins) + 1;
        try
            riseTime = riseTimeArray(rtBin);
        catch err
            disp(err)
        end
        durationMinis = length(V(rtBin,:));
        offsetMinis = onsetMinis + durationMinis - 1;
        if offsetMinis <= length(t) % checking whether the overlayed mini is within the length of recording sweep
            simV(onsetMinis:offsetMinis) = simV(onsetMinis:offsetMinis) + amplitudeArray(iAmplitude)*V(rtBin,:);
        else
            durationMinis1 = length(t) - onsetMinis + 1;
            simV(onsetMinis:end) = simV(onsetMinis:end) + amplitudeArray(iAmplitude)*V(rtBin,1:durationMinis1);
            durationMinis2 = durationMinis - durationMinis1;
            simV(1:durationMinis2) = simV(1:durationMinis2) + amplitudeArray(iAmplitude)*V(rtBin,durationMinis1+1:end);
        end
        [~, iPeak] = max(V(rtBin,:));
        peakMinis = onsetMinis + iPeak - 1;
        if peakMinis > numel(t)
          peakMinis = peakMinis - numel(t);
        end
        peakMinisT = t(peakMinis);
        shapes(countMinis,:) = [countMinis amplitudeArray(iAmplitude) riseTime sitesLengths(rtBin,:)...
            onsetMinis onsetMinisT peakMinis peakMinisT];
    end
end
eventMatrix = nEvents;

draw = 0;
if draw == 1
    f3 = figure('position', [5 300 600 300]);
    figure(f3);
    if dimension == 1
        %pcolor(amplitudeArray, sortedRiseTimeArray, eventMatrix);
        imagesc(amplitudeArray, sortedRiseTimeArray, eventMatrix);
    elseif dimension == 1.5 || dimension == 2
        %pcolor(amplitudeArray, riseTimeArray, eventMatrix);
        imagesc(amplitudeArray, riseTimeArray, eventMatrix);
    end
    set(gca,'YDir','normal');
    title('Two-dimensional distribution of simulated minis');
    xlabel('Amplitude(mV)');
    ylabel('10-90% rise times (ms)');
    colorbar;
end
end