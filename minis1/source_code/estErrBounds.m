function [Amps, AmpsNeg, RTs, RTsNeg, twoDs, twoDsNeg, fileNames, fileSweeps, SD, fileAmps, fileRTs, fileAmpsNeg, fileRTsNeg, sLength, sweepcount, ld] = estErrBounds(...
    excludedTimes, detectionParameters, classificationParameters, ld, parallelCores)

%% Load the files and divide them into sweeps:
filtering = noiseFilterDlg;
optionsP.SD = [1 1 1 0];
filt.state = 'off';
waveform.estimate = false;

cd(ld);
files = dir('*.abf');
startPulse = 1000*excludedTimes.startPulse;
endPulse = 1000*excludedTimes.endPulse;
%
fileNames = struct([]);
fileSweeps = zeros(length(files),1);
SD = zeros(length(files),3);
fileAmps = zeros(length(files),length(classificationParameters.amplitudeArray));
fileRTs = zeros(length(files),length(classificationParameters.riseTimeArray));
fileAmpsNeg = fileAmps;
fileRTsNeg = fileRTs;
%
Amps = [];
AmpsNeg = [];
RTs = [];
RTsNeg = [];
twoDs = struct([]);
twoDsNeg = struct([]);
sLength = [];
sweepcount = 0;
%
progress = 100/length(files): 100/length(files) :100;
h = waitbar(0,'Loading files and detecting minis-like events - 0%');
waitText1 = 'Loading files and detecting minis-like events - ';
waitText2 = '%';
for iFile = 1:length(files)
    filenameSplit = regexp(files(iFile).name, '\\*', 'split');
    filenameShort = char(filenameSplit(end));
    fileNames{iFile} = filenameShort;
    dataProperties = loadABF(files(iFile).name);
    detectionParameters2 = detectionParameters;
    detectionParameters2.sampleInterval = dataProperties.dt;
    detectionParameters2.smoothWindow = round(detectionParameters2.smoothWindow/detectionParameters2.sampleInterval);
    detectionParameters2.smoothWindowLite = 8;
    nSweepDuration = round(length(dataProperties.sweep)/dataProperties.hd.lActualEpisodes);
    sweepDuration = nSweepDuration*dataProperties.dt - dataProperties.dt;
    excludedTimes2 = calcExcludedTimes(sweepDuration, dataProperties.hd.lActualEpisodes, startPulse, endPulse, [], [], dataProperties.dt);
    excludedInd = round(excludedTimes2/dataProperties.dt)+1;
    
    % Filter the noise:
    filtering.nSweeps = dataProperties.hd.lActualEpisodes;
    filtering.excludedTimes = excludedTimes;
    if strcmpi(filtering.state,'on')
        if iFile == 1
            [dataProperties.sweep, ~, f0, filtfs] = filterMinis(dataProperties.sweep, dataProperties.dt, filtering, true);
            close(f0);
        else
            dataProperties.sweep = filterMinis(dataProperties.sweep, dataProperties.dt, filtering, false, [], filtfs);
        end
    end
    
    % Detect mini-like events and summarise the detection:
    detectedEvents = detectMinis(dataProperties.sweep, excludedTimes2, detectionParameters2, filt, waveform, parallelCores, optionsP);
    detectedEventsNeg = detectMinis(-dataProperties.sweep, excludedTimes2, detectionParameters2, filt, waveform, parallelCores);
    detectedEventsNeg(:,[1,4,5,19:21]) = -detectedEventsNeg(:,[1,4,5,19:21]);
    SD(iFile,:) = detectedEvents(1,22:24);
    fileSweeps(iFile) = dataProperties.hd.lActualEpisodes;
    % Divide into sweeps:
    for iSweep = 1:dataProperties.hd.lActualEpisodes
        sweepcount = sweepcount + 1;
        
        % Positive events:
        sweepEvents = detectedEvents(detectedEvents(:,2) > sweepDuration*(iSweep-1) & detectedEvents(:,2) <= sweepDuration*(iSweep),:);
        [sweepEvents1D, sweepEvents1D_RT, sweepEventsTwoD] = classifyMinis(sweepEvents(:,4), sweepEvents(:,12), classificationParameters);
        Amps = [Amps; sweepEvents1D]; %#ok<AGROW>
        RTs = [RTs; sweepEvents1D_RT]; %#ok<AGROW>
        twoDs{sweepcount} = sweepEventsTwoD;
        
        % Negative events:
        sweepEvents = detectedEventsNeg(detectedEventsNeg(:,2) > sweepDuration*(iSweep-1) & detectedEventsNeg(:,2) <= sweepDuration*(iSweep),:);
        [sweepEvents1DNeg, sweepEvents1D_RTNeg, sweepEventsTwoDNeg] = classifyMinisNeg(sweepEvents(:,4), sweepEvents(:,12), classificationParameters);
        AmpsNeg = [AmpsNeg; sweepEvents1DNeg]; %#ok<AGROW>
        RTsNeg = [RTsNeg; sweepEvents1D_RTNeg]; %#ok<AGROW>
        twoDsNeg{sweepcount} = sweepEventsTwoDNeg;
        
        if iSweep > 1
            sweepTime = (iSweep-1)*sweepDuration+dataProperties.dt: dataProperties.dt :iSweep*sweepDuration;
        else
            sweepTime = (iSweep-1)*sweepDuration: dataProperties.dt :iSweep*sweepDuration;
        end
        excludedIndSweep = excludedInd(excludedInd >= sweepTime(1)/dataProperties.dt & excludedInd <= sweepTime(end)/dataProperties.dt);
        sLength = [sLength; (length(sweepTime) - length(excludedIndSweep))*dataProperties.dt]; %#ok<AGROW>
    end
    fileAmps(iFile,:) = sum(Amps(end-dataProperties.hd.lActualEpisodes+1:end,:),1);
    fileRTs(iFile,:) = sum(RTs(end-dataProperties.hd.lActualEpisodes+1:end,:),1);
    fileAmpsNeg(iFile,:) = sum(AmpsNeg(end-dataProperties.hd.lActualEpisodes+1:end,:),1);
    fileRTsNeg(iFile,:) = sum(RTsNeg(end-dataProperties.hd.lActualEpisodes+1:end,:),1);
    
    waitText = [waitText1 num2str(progress(iFile)) waitText2];
    waitbar(progress(iFile)/100, h, waitText);
end
close(h);

if exist('files', 'var') && isempty(files)
    Amps = [];
    AmpsNeg = [];
    RTs = [];
    RTsNeg = [];
    twoDs = [];
    twoDsNeg = [];
    fileNames = [];
    fileSweeps = [];
    SD = [];
    fileAmps = [];
    fileRTs = [];
    fileAmpsNeg = [];
    fileRTsNeg = [];
    sLength = [];
    sweepcount = [];
    msgbox('Error: The provided path directory does not contain abf files. Please provide another path directory or copy abf files to the existing one.',...
        'Error','Error');
end