function output = simulateDetectEvaluate(filename, excludedTimesInit, detectionParameters, simulationParameters, optimisationParameters, classificationParameters, filtering, parallelCores, varargin)
% simulateDetectEvaluate(filename, excludedTimes, detectionParameters, simulationParameters, optimisationParameters, classificationParameters, filtering, parallelCores, outputDirName)
%
% Function simulates minis drawn from a specified minis distribution,
% superimposes these simulated events on a noise file (if one is
% specified), carries out detection on the simulated data, and evaluates
% the detection performance based on signal detection theory measures.



rng('shuffle');

if nargin > 8
  outputDirName = varargin{1};
end



%% Parameters
draw = false;
nAddedEvents = round(mean(optimisationParameters.options.bounds(:,3))); %[500 1000 2000 4000 8000 16000 32000 64000 128000];
% if contains(filename,'p103a') || contains(filename,'p106b')
%     %nAddedEvents = nAddedEvents*2.5;
%     nAddedEvents = nAddedEvents*4;
% else
%     %nAddedEvents = nAddedEvents*1.25;
%     nAddedEvents = nAddedEvents*2;
% end
scaleNoise = true;
noiseScaleFactor = 1; %[1 1.2 1.4 1.8 2.6 4.2];
amplitudeScaleFactor = 1; %[1/6 1/3 1/2 2/3 5/6 1 7/6 8/6 9/6 10/6];
nRuns = 1; %4;
hitWindow = 10; % ms
noiseWindow = 20; % ms


for iPopulation = 1:numel(nAddedEvents)
    optimisationParameters.options.bounds(1,3) = nAddedEvents(iPopulation);
    optimisationParameters.options.bounds(2,3) = nAddedEvents(iPopulation);
    for iFactor = 1:numel(noiseScaleFactor)
        if iPopulation > 1 && iFactor > 1
            continue
        end
        
        
        %% Process the file:
        if ~exist('noiseProperties', 'var')
            if isempty(filename)
                noiseProperties.filename = [];
                noiseProperties.hd = [];
                noiseProperties.nchans_to_save = 1;
                noiseProperties.sweep = -75.*ones(1,4000000);
                noiseProperties.dt = 0.0500; % smpling rate of 20 kHz
                noiseProperties.tOffset = 0;
                noiseProperties.current = [];
                filtN.state = filtering.state;
                excludedTimesInit = [];
            else
                % Load the file:
                if ischar(filename)
                    noiseProperties = loadABF(filename);
                else
                    noiseProperties = filename;
                    noiseProperties.hd.lActualEpisodes = noiseProperties.lActualEpisodes;
                end
                dt = noiseProperties.dt;
                
                % Determine filtering mode:
                filtN.state = filtering.state;
                filtN.nSweeps = noiseProperties.hd.lActualEpisodes;
                filtN.excludedTimes = excludedTimesInit;
                if isfield(filtering, 'filtfs')
                  filtN.filtfs = filtering.filtfs;
                end
                
                % Determine excluded times for a noise file:
                sweepDuration = (length(noiseProperties.sweep)*dt - dt)/noiseProperties.hd.lActualEpisodes;
                excludedTimes = calcExcludedTimes(sweepDuration, noiseProperties.hd.lActualEpisodes, 1000*excludedTimesInit.startPulse,...
                    1000*excludedTimesInit.endPulse, 1000*excludedTimesInit.startGlitch, 1000*excludedTimesInit.endGlitch, dt);
            end
            
            if strcmpi(filtN.state, 'on')
                [noiseProperties.sweep, ~, f2, filtering.filtfs] = filterMinis(noiseProperties.sweep, noiseProperties.dt, filtN, true);
                close(f2);
            end
            initSweep = noiseProperties.sweep;
        end
        
        if scaleNoise
            noiseProperties.sweep = (initSweep - mean(initSweep)) * noiseScaleFactor(iFactor);
        end
        
        
        
        %% Initialise input variables:
        distributionType = optimisationParameters.distType;
        distParameters = (optimisationParameters.options.bounds(1,:)...
            + optimisationParameters.options.bounds(2,:))./2;
        
        noiseProperties.baseline = length(classificationParameters.amplitudeArray);
        
        amplitudeArraySim = classificationParameters.amplitudeArray;
        amplitudeArraySim(amplitudeArraySim < simulationParameters.loSimAmp) = [];
        
        simulationParameters.R_m = simulationParameters.tau_m;
        simulationParameters.tau_sy1 = 2;
        simulationParameters.tau_sy2 = .1;
        
        
        
        %% Initialise output variables
        sensitivity = cell(nRuns,1);
        specificity = cell(nRuns,1);
        FPR = cell(nRuns,1);
        dPrime = cell(nRuns,1);
        performance = cell(nRuns,1);
        simulatedEventInfo = cell(nRuns,1);
        
        
        
        %% Pick up output directory
        if ~exist('outputDirName', 'var')
            outputDirName = uigetdir(pwd,...
                'Choose output directory for saving simulated voltage traces and detection algorithm performance data');
            if ~exist('outputDirName', 'file')
                mkdir(outputDirName);
            end
        end
        %if ~exist('outputDirNameRaw', 'var')
        %    outputDirNameRaw = [outputDirName '_raw'];
        %    if ~exist('outputDirNameRaw', 'file')
        %        mkdir(outputDirNameRaw);
        %   end
        %end
        
        
        
        %% Generate output filename
        currentTime = datestr(now);
        outputFilename = ['algorithm_performance_data__amp' num2str(optimisationParameters.options.bounds(1,1))...
            '_sAmp' num2str(optimisationParameters.options.bounds(1,2)) '_n' num2str(optimisationParameters.options.bounds(1,3))...
            '_RT' num2str(optimisationParameters.options.bounds(1,4)) '_sRT' num2str(optimisationParameters.options.bounds(1,5))...
            '_rho' num2str(optimisationParameters.options.bounds(1,6)) '_thr' num2str(detectionParameters.Amplobound)...
            '_noiseScaleFactor' num2str(noiseScaleFactor(iFactor)), '_smoothWindow' num2str(detectionParameters.smoothWindow),...
            '__' currentTime];
        outputFilename = strrep(outputFilename,'-','_');
        outputFilename = strrep(outputFilename,' ','_');
        outputFilename = strrep(outputFilename,':','_');
        outputFilename = strrep(outputFilename,'.','p');
        
        
        
        %% Detect noise events
        detectionParametersSim = detectionParameters;
        detectionParametersSim.sampleInterval = noiseProperties.dt;
        detectionParametersSim.smoothWindow = round(detectionParametersSim.smoothWindow/detectionParametersSim.sampleInterval);
        waveform.estimate = false;
        filtN.state = 'spectrum';
        options.summaryPlot = false;
        options.edit = false;
        [~, noiseVfilt, ~, ~, ~, noiseV] = detectMinis(noiseProperties.sweep, excludedTimes, detectionParametersSim, filtN, waveform, parallelCores, options); %#ok<*ASGLU>
        %writeABF(single([noiseV; zeros(size(noiseProperties.current))]), [outputDirName filesep 'noiseFiltSmooth.abf'], 1000/dt, {'mV';'pA'});
        %writeABF(single([noiseVfilt; zeros(size(noiseProperties.current))]), [outputDirName filesep 'noiseFilt.abf'], 1000/dt, {'mV';'pA'});
        % excludedInds = zeros(size(noiseProperties.sweep));
        % excludedInds(round(excludedTimes(excludedTimes > 0)./noiseProperties.dt)) = 1;
        % noiseSD = std(noiseV(~logical(excludedInds)));
        minPeakWidth = 0.5/noiseProperties.dt;
        minPeakAmp = 0.01;
        % [~, falseI] = findpeaks(noiseV, 'MinPeakWidth',minPeakWidth, 'MinPeakProminence',minPeakAmp);
        filtNoiseV = movmean(noiseV,noiseWindow/noiseProperties.dt);
        [~, falseI] = findpeaks(filtNoiseV, 'MinPeakWidth',minPeakWidth, 'MinPeakProminence',minPeakAmp);
        [~, adjustedFalseI] = max(noiseV(falseI - noiseWindow/2 + 1 : falseI + noiseWindow/2));
        falseI = falseI - noiseWindow/2 + adjustedFalseI;
        falseT = falseI.*dt; % False (noise) events
        falseI(logical(ismember(round(falseT./dt),round((excludedTimes+dt)./dt)))) = [];
        falseT(logical(ismember(round(falseT./dt),round((excludedTimes+dt)./dt)))) = [];
        
        
        
        %parfor iRun = 1:nRuns
        for iRun = 1:nRuns
            
            %% Simulate minis:
            smoothWindow = max([floor((amplitudeArraySim(1)/(classificationParameters.amplitudeArray(2) - classificationParameters.amplitudeArray(1)))/2) 3]); %#ok<*PFBNS>
            for iAmp = 1:numel(amplitudeScaleFactor)
                if iPopulation == 1 && iFactor == 1
                    distParametersScaled = distParameters;
                    distParametersScaled(1) = distParametersScaled(1)*amplitudeScaleFactor(iAmp);
                    if iAmp == 1
                        [simV, V, minis2D, shapes] = simulateMinis('Zero', simulationParameters, distributionType, distParametersScaled, noiseProperties.dt,...
                            noiseProperties.baseline, noiseProperties.sweep, excludedTimes, smoothWindow, parallelCores, classificationParameters.amplitudeArray,...
                            amplitudeArraySim, classificationParameters.riseTimeArray); % simV: simulated events only; V: noise + simulated events
                        minis1D = sum(minis2D,1); %#ok<*NASGU>
                        minis1D_RT = sum(minis2D,2)';
                    else
                        [simVAmp, V, minis2D, shapesAmp] = simulateMinis('Zero', simulationParameters, distributionType, distParametersScaled, noiseProperties.dt,...
                            noiseProperties.baseline, V, excludedTimes, smoothWindow, parallelCores, classificationParameters.amplitudeArray,...
                            amplitudeArraySim, classificationParameters.riseTimeArray); % simV: simulated events only; V: noise + simulated events
                        simV = simV + simVAmp;
                        shapes = [shapes; shapesAmp];
                        minis1D = minis1D + sum(minis2D,1); %#ok<*NASGU>
                        minis1D_RT = minis1D_RT + sum(minis2D,2)';
                    end
                else
                    [simV, V, minis2D, shapes] = simulateMinis('Zero', simulationParameters, distributionType, distParametersScaled, noiseProperties.dt,...
                        noiseProperties.baseline, noiseProperties.sweep, excludedTimes, smoothWindow, parallelCores, classificationParameters.amplitudeArray,...
                        amplitudeArraySim, classificationParameters.riseTimeArray); % simV: simulated events only; V: noise + simulated events
                    minis1D = sum(minis2D,1); %#ok<*NASGU>
                    minis1D_RT = sum(minis2D,2)';
                end
            end
            simulatedEventInfo{iRun} = shapes;
            
            
            
            %% Detect simulated events:
            % Adapt the Maximum time to peak window to a particular rise time
            % distribution:
            detectionParametersSim2 = detectionParametersSim;
            minis1D_RTprc = cumsum(minis1D_RT)/sum(minis1D_RT);
            minis1D_RT75prc = find(minis1D_RTprc > .75, 1);
            if ~isempty(minis1D_RT75prc)
                minis1D_RT75prc = classificationParameters.riseTimeArray(minis1D_RT75prc)*(.75/minis1D_RTprc(minis1D_RT75prc));
                SWstart = 2*minis1D_RT75prc + detectionParameters.BLduration;
                detectionParametersSim2.SWstart = max([round(SWstart) detectionParameters.SWstart]);
            end
            
            [simulatedEvents, filterV, ~, ~, ~, smoothV] = detectMinis(V, excludedTimes, detectionParametersSim2, filtN, waveform, parallelCores, options);
            
            
            
            %% Save simulated trace
            if iRun < 10
                simFilename = strcat(outputFilename,'_000',num2str(iRun),'.abf');
                simFilenameRaw = strcat(outputFilename,'_000',num2str(iRun),'_raw.abf');
            elseif iRun < 100
                simFilename = strcat(outputFilename,'_00',num2str(iRun),'.abf');
                simFilenameRaw = strcat(outputFilename,'_00',num2str(iRun),'_raw.abf');
            elseif iRun < 1000
                simFilename = strcat(outputFilename,'_0',num2str(iRun),'.abf');
                simFilenameRaw = strcat(outputFilename,'_0',num2str(iRun),'_raw.abf');
            end
            if ~isempty(outputDirName)
                writeABF(single([smoothV; zeros(size(noiseProperties.current))]), [outputDirName filesep simFilename], 1000/dt, {'mV';'pA'});
                %writeABF(single([simV; zeros(size(noiseProperties.current))]), [outputDirName filesep simFilename(1:end-4) 'simOnly.abf'], 1000/dt, {'mV';'pA'});
                %writeABF(single([filterV; zeros(size(noiseProperties.current))]), [outputDirNameRaw filesep simFilename], 1000/dt, {'mV';'pA'});
                writeABF(single([filterV; zeros(size(noiseProperties.current))]), [outputDirName filesep simFilenameRaw], 1000/dt, {'mV';'pA'});

                %writeBinary(int16(smoothV), [outputDirName filesep simFilename(1:end-4) '.dat']);
                %writeBinary(int16(smoothI), [outputDirName filesep simFilename(1:end-4) '_current.dat']);
            end
            
            
            
            %% Evaluate detection performance
            % 1. Take true events
            % 2. Pick their closest detected events.
            % 3. Include true events not exceeding -5 and 5 ms relative to their closest detected events as true positives.
            % 4. Otherwise class them as false negatives.
            % 5. Class detected events that have not been marked as true positives as false alarms
            % 6. Include false events not exceeding -5 and 5 ms relative to their closest detected events as correct rejections.
            % 7. Calculate sensitivity, specificity, and other measures of detection performance.
            
            % Associate detected events with true (simulated minis/ ground true events) and false (peaks in the noise data) events
            crWindow = hitWindow; % ms
            [~, sortInd] = sort(shapes(:,9));
            sortedShapes = shapes(sortInd,:);
            trueI = sortedShapes(:,8); % True events
            trueT = sortedShapes(:,9);
            positivesI = simulatedEvents(:,3); % Positive identifications (includes both hits and false alarms)
            positivesT = simulatedEvents(:,2);
            positivesAssociated2true = zeros(size(positivesI)); % Positive identifications associated to the nearest true events
            positivesAssociated2false = zeros(size(positivesI)); % Positive identifications associated to the nearest false events
            for iPositive = 1:numel(positivesI)
                trueDist = abs(trueT - positivesT(iPositive));
                [~, nearestTrueI] = min(trueDist);
                positivesAssociated2true(iPositive) = trueI(nearestTrueI);
                falseDist = abs(falseT - positivesT(iPositive));
                [~, nearestFalseI] = min(falseDist);
                positivesAssociated2false(iPositive) = falseI(nearestFalseI);
            end
            
            % Calculate the distance to the nearest neighbour
            distanceToTheRight = abs([trueT(2:end); inf] - trueT);
            distanceToTheLeft = abs([inf; trueT(1:end-1)] - trueT);
            distance2neighbour{iRun} = min([distanceToTheRight; distanceToTheLeft],[],1); %#ok<*AGROW>
            
            % Locate hits and misses
            truePositives = zeros(1,numel(simV));
            falseNegatives = zeros(1,numel(simV));
            for iMini = 1:numel(trueI)
                detectedPositivesT = positivesT(trueI(iMini) == positivesAssociated2true);
                detectedPositivesI = positivesI(trueI(iMini) == positivesAssociated2true);
                if ~isempty(detectedPositivesT)
                    detectedPositivesDist = abs(detectedPositivesT - trueT(iMini));
                    [minDist, minDistI] = min(detectedPositivesDist);
                    if minDist <= hitWindow/2
                        truePositives(detectedPositivesI(minDistI)) = 1;
                    else
                        falseNegatives(trueI(iMini)) = 1;
                    end
                else
                    falseNegatives(trueI(iMini)) = 1;
                end
            end
            
            % Locate false alarms
            falsePositives = zeros(1,numel(simV));
            falsePositives(positivesI) = 1;
            falsePositives(logical(truePositives)) = 0;
            
            % Locate correct rejections
            positivesAssociated2false(logical(ismember(positivesI, find(truePositives)))) = [];
            trueNegatives = zeros(1,numel(simV));
            trueNegatives(falseI) = 1;
            for iMini = 1:numel(falseI)
              iPositivesAssociated2false = positivesAssociated2false(falseI(iMini) == positivesAssociated2false);
              tPositivesAssociated2false = falseT(ismember(falseI, iPositivesAssociated2false));
              if ~isempty(tPositivesAssociated2false)
                detectedPositivesDist = abs(tPositivesAssociated2false - positivesT);
                minDist = min(detectedPositivesDist);
                if minDist <= crWindow/2
                  trueNegatives(falseI(iMini)) = 0;
                end
              end
            end
            trueNegatives(logical(falseNegatives) | logical(truePositives) | logical(falsePositives)) = 0; % just as an insurance
            
            sensitivity{iRun} = sum(truePositives)/(sum(truePositives) + sum(falseNegatives)); %#ok<*PFOUS> % True positive rate
            specificity{iRun} = sum(trueNegatives)/(sum(trueNegatives) + sum(falsePositives)); % Correct rejection rate
            FPR{iRun} = sum(falsePositives)/(sum(trueNegatives) + sum(falsePositives)); % False positive rate
            if sensitivity{iRun} == 1
                sensitivityApprox = 1-(1e-6);
            else
                sensitivityApprox = sensitivity{iRun};
            end
            if FPR{iRun} == 0
                FPRapprox = 1e-6;
            else
                FPRapprox = FPR{iRun};
            end
            dPrime{iRun} = dprime_simple(sensitivityApprox, FPRapprox);
            
            allTrue = zeros(1,numel(simV));
            allTrue(trueI) = 1;
            allTrue2 = zeros(1,numel(simV));
            allTrue2(logical(truePositives) | logical(falseNegatives)) = 1;
            allPositive = zeros(1,numel(simV));
            allPositive(positivesI) = 1;
            performance{iRun} = sparse([allTrue; allTrue2; allPositive; truePositives; falseNegatives; falsePositives; trueNegatives]);
            % simulated event positions; hits (detected positions) + misses; hits (detected positions) + false alarms;
            % hits (detected positions); misses; false alarms; correct rejections
            
            if draw
                time = (1:numel(smoothV)).*dt;
                figure; plot(time, smoothV); hold on;
                plot(time(falseI), smoothV(falseI), '.r', 'MarkerSize',10);
                plot(time(trueI), smoothV(trueI), '.g', 'MarkerSize',10);
                plot(time(logical(truePositives)), smoothV(logical(truePositives)), 'og', 'MarkerSize',10);
                plot(time(logical(falsePositives)), smoothV(logical(falsePositives)), 'oy', 'MarkerSize',10);
                plot(time(logical(trueNegatives)), smoothV(logical(trueNegatives)), 'or', 'MarkerSize',10);
                plot(time(logical(falseNegatives)), smoothV(logical(falseNegatives)), 'ob', 'MarkerSize',10);
                legend('V_m','False events','True events','True positive','False positive','Correct rejection','False rejection');
            end
            
            
            
            %% Display detected and true event histograms:
            if draw
                [simulatedEvents1D, simulatedEvents1D_RT, simulatedEvents2D] = classifyMinis(simulatedEvents(:,4),...
                    simulatedEvents(:,12), classificationParameters); %#ok<*UNRCH>
                amplitudeArray = classificationParameters.amplitudeArray;
                ampLim2D = (amplitudeArray(2) - amplitudeArray(1))/2;
                riseTimeArray = classificationParameters.riseTimeArray;
                RTLim2D = (riseTimeArray(2) - riseTimeArray(1))/2;
                
                % Drawing the uni-dimensional amplitude histogram:
                f1 = figure('position', [50 50 1200 600]);
                figure(f1);
                subplot(2,1,1,'replace');
                hold on
                targetAmpcs = cumsum(simulatedEvents1D)/sum(simulatedEvents1D);
                targetAmp99prc = find(targetAmpcs > .99, 1);
                iEndAmp = min([targetAmp99prc+1 length(amplitudeArray)]);
                xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
                plot([0 amplitudeArray(end)], [0 0], 'k-');
                p1 = plot(amplitudeArray, simulatedEvents1D, 'b.-');
                p2 = plot(amplitudeArray, minis1D, 'g.-');
                set(f1, 'NumberTitle', 'off');
                set(f1, 'Name', 'One-dimensional histograms');
                title('Amplitude distribution');
                xlabel('Amplitude(mV)');
                ylabel('Number of Events');
                legend([p1 p2], {'Detected','Simulated'});
                hold off
                
                % Drawing the uni-dimensional rise time histogram:
                subplot(2,1,2,'replace');
                hold on
                targetRTcs = cumsum(simulatedEvents1D_RT)/sum(simulatedEvents1D_RT);
                targetRT99prc = find(targetRTcs > .99, 1);
                iEndRT = min([targetRT99prc+1 length(riseTimeArray)]);
                xlim([riseTimeArray(1) riseTimeArray(iEndRT)]);
                plot([0 riseTimeArray(end)], [0 0], 'k-');
                p1 = plot(riseTimeArray, simulatedEvents1D_RT,'b.-');
                p2 = plot(riseTimeArray, minis1D_RT, 'g.-');
                if strcmpi(detectionParameters.RTinterval, '10-90%')
                    title('10-90% rise time distribution');
                    xlabel('10-90% rise times (ms)');
                elseif strcmpi(detectionParameters.RTinterval, '20-80%')
                    title('20-80% rise time distribution');
                    xlabel('20-80% rise times (ms)');
                end
                ylabel('Number of Events');
                legend([p1 p2], {'Detected','Simulated'});
                hold off
                
                % Drawing the two-dimensional histogram:
                f2 = figure('position', [50 50 1200 600]);
                figure(f2);
                hold on
                xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
                ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
                imagesc(amplitudeArray, riseTimeArray, simulatedEvents2D);
                set(gca,'YDir','normal');
                set(f2, 'NumberTitle', 'off');
                set(f2, 'Name', 'Two-dimensional histogram');
                title('2D detected event distribution');
                xlabel('Amplitude(mV)');
                if strcmpi(detectionParameters.RTinterval, '10-90%')
                    ylabel('10-90% rise times (ms)');
                elseif strcmpi(detectionParameters.RTinterval, '20-80%')
                    ylabel('20-80% rise times (ms)');
                end
                colorbar;
                hold off
            end
            
            
            
            %% Display detected and true event cumulative histograms:
            if draw
                simulatedEvents1D_cs = cumsum(simulatedEvents1D)/sum(simulatedEvents1D);
                simulatedEvents1D_RT_cs = cumsum(simulatedEvents1D_RT)/sum(simulatedEvents1D_RT);
                minis1D_cs = cumsum(minis1D)/sum(minis1D);
                minis1D_RT_cs = cumsum(minis1D_RT)/sum(minis1D_RT);
                
                % Drawing the uni-dimensional cumulative amplitude histogram:
                f1 = figure('position', [50 50 1200 600]);
                figure(f1);
                subplot(2,1,1,'replace');
                hold on
                xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
                plot([0 amplitudeArray(end)], [0 0], 'k-');
                p1 = plot(amplitudeArray, simulatedEvents1D_cs, 'b.-');
                p2 = plot(amplitudeArray, minis1D_cs, 'g.-');
                set(f1, 'NumberTitle', 'off');
                set(f1, 'Name', 'One-dimensional histograms');
                title('Cumulative amplitude distribution');
                xlabel('Amplitude(mV)');
                ylabel('Number of Events');
                legend([p1 p2], {'Detected','Simulated'});
                hold off
                
                % Drawing the uni-dimensional cumulative rise time histogram:
                subplot(2,1,2,'replace');
                hold on
                xlim([riseTimeArray(1) riseTimeArray(iEndRT)]);
                plot([0 riseTimeArray(end)], [0 0], 'k-');
                p1 = plot(riseTimeArray, simulatedEvents1D_RT_cs,'b.-');
                p2 = plot(riseTimeArray, minis1D_RT_cs, 'g.-');
                if strcmpi(detectionParameters.RTinterval, '10-90%')
                    title('Cumulative 10-90% rise time distribution');
                    xlabel('10-90% rise times (ms)');
                elseif strcmpi(detectionParameters.RTinterval, '20-80%')
                    title('Cumulative 20-80% rise time distribution');
                    xlabel('20-80% rise times (ms)');
                end
                ylabel('Number of Events');
                legend([p1 p2], {'Detected','Simulated'});
                hold off
            end
            
            
            
            %% Assign output variables
            simData{iRun} = single([smoothV; zeros(size(noiseProperties.current))]);
            simDataRaw{iRun} = single([filterV; zeros(size(noiseProperties.current))]);
        end
        
        
        
        %% Save algorithm performance data
        save([outputDirName filesep outputFilename '.mat'], 'sensitivity','specificity','FPR','dPrime','performance','falseI','falseT',...
            'noiseScaleFactor','iFactor','dt','filename','excludedTimes','detectionParameters','simulationParameters',...
            'optimisationParameters','classificationParameters','filtering','distance2neighbour','simulatedEventInfo', '-v7.3');
        
          
          
        %% Assign output variable
        if iPopulation == 1 && iFactor == 1 && nRuns == 1
            output.simData = simData{1};
            output.simDataRaw = simDataRaw{1};
            detectionPerformance.sensitivity = sensitivity{1};
            detectionPerformance.specificity = specificity{1};
            detectionPerformance.FPR = FPR{1};
            detectionPerformance.dPrime = dPrime{1};
            detectionPerformance.performance = performance{1};
            detectionPerformance.falseI = falseI;
            detectionPerformance.falseT = falseT;
            output.detectionParameters = detectionParameters;
            output.simulationParameters = simulationParameters;
            output.optimisationParameters = optimisationParameters;
            output.classificationParameters = classificationParameters;
            output.detectionPerformance = detectionPerformance;
            if isfield(filtering, 'filtfs')
                output.filtfs = filtering.filtfs;
            else
                output.filtfs = {};
            end
            output.simulatedEventInfo = simulatedEventInfo{1};
        else
            output = [];
        end
    end
end