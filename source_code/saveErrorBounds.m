function ld = saveErrorBounds(AmpsMean, RTsMean, TwoDsMean, AmpsMedian, RTsMedian, TwoDsMedian, AmpsMin, RTsMin, TwoDsMin, AmpsMax, RTsMax, TwoDsMax,...
    AmpsPrct, RTsPrct, TwoDsPrct, boundAmps, boundRTs, boundTwoDs, boundAmpsLinear, boundRTsLinear, boundTwoDsLinear, AmpsMeanBottom, AmpsMedianBottom,...
    AmpsMinBottom, AmpsMaxBottom, AmpsPrctBottom, boundAmpsBottom, boundAmpsLinearBottom, AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid,...
    boundAmpsMid, boundAmpsLinearMid, AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, boundAmpsTop, boundAmpsLinearTop, sLength,...
    slFactor, fileNames, fileSweeps, flFactor, pow, powWgh, estCount, UFest, AmpsLSE, RTsLSE, TwoDsLSE, AmpsLinearLSE, RTsLinearLSE, TwoDsLinearLSE,...
    AmpsLSEBottom, AmpsLinearLSEBottom, AmpsLSEMid, AmpsLinearLSEMid, AmpsLSETop, AmpsLinearLSETop, SD, optimData, dataF, cntF, F, G, H, A, B, C, FBottom,...
    ABottom, FMid, AMid, FTop, ATop, detectionParameters, classificationParameters, ld, wd, expand)

dataDir = 'error_bounds';

button = questdlg('Save the error bound estimates in a text file?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    cd(wd);
    if ~exist(dataDir,'dir')
        mkdir(dataDir);
    end
    cd(dataDir);
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.txt','Text files (*.txt)'},'Save Error Bound Estimates as', ld);
    if filterIndex
        iEndShort = min([size(AmpsMean,2) 90]);
        iEnd = min([size(boundAmps,2) 90]);
        AmpsMean = AmpsMean(:,1:iEndShort);
        RTsMean = RTsMean(:,1:iEndShort);
        TwoDsMean = TwoDsMean(:,1:iEndShort);
        AmpsMedian = AmpsMedian(:,1:iEndShort);
        RTsMedian = RTsMedian(:,1:iEndShort);
        TwoDsMedian = TwoDsMedian(:,1:iEndShort);
        AmpsMin = AmpsMin(:,1:iEndShort);
        RTsMin = RTsMin(:,1:iEndShort);
        TwoDsMin = TwoDsMin(:,1:iEndShort);
        AmpsMax = AmpsMax(:,1:iEndShort);
        RTsMax = RTsMax(:,1:iEndShort);
        TwoDsMax = TwoDsMax(:,1:iEndShort);
        AmpsPrct = AmpsPrct(:,1:iEndShort);
        RTsPrct = RTsPrct(:,1:iEndShort);
        TwoDsPrct = TwoDsPrct(:,1:iEndShort);
        AmpsMeanBottom  = AmpsMeanBottom(:,1:iEndShort);
        AmpsMedianBottom = AmpsMedianBottom(:,1:iEndShort);
        AmpsMinBottom = AmpsMinBottom(:,1:iEndShort);
        AmpsMaxBottom = AmpsMaxBottom(:,1:iEndShort);
        AmpsPrctBottom = AmpsPrctBottom(:,1:iEndShort);
        AmpsMeanMid  = AmpsMeanMid(:,1:iEndShort);
        AmpsMedianMid = AmpsMedianMid(:,1:iEndShort);
        AmpsMinMid = AmpsMinMid(:,1:iEndShort);
        AmpsMaxMid = AmpsMaxMid(:,1:iEndShort);
        AmpsPrctMid = AmpsPrctMid(:,1:iEndShort);
        AmpsMeanTop = AmpsMeanTop(:,1:iEndShort);
        AmpsMedianTop = AmpsMedianTop(:,1:iEndShort);
        AmpsMinTop = AmpsMinTop(:,1:iEndShort);
        AmpsMaxTop = AmpsMaxTop(:,1:iEndShort);
        AmpsPrctTop = AmpsPrctTop(:,1:iEndShort);
        boundAmps = boundAmps(:,1:iEnd);
        boundRTs = boundRTs(:,1:iEnd);
        boundTwoDs = boundTwoDs(:,1:iEnd);
        boundAmpsBottom = boundAmpsBottom(:,1:iEnd);
        boundAmpsMid = boundAmpsMid(:,1:iEnd);
        boundAmpsTop = boundAmpsTop(:,1:iEnd);
        boundAmpsLinear = boundAmpsLinear(:,1:iEnd);
        boundRTsLinear = boundRTsLinear(:,1:iEnd);
        boundTwoDsLinear = boundTwoDsLinear(:,1:iEnd);
        boundAmpsLinearBottom = boundAmpsLinearBottom(:,1:iEnd);
        boundAmpsLinearMid = boundAmpsLinearMid(:,1:iEnd);
        boundAmpsLinearTop = boundAmpsLinearTop(:,1:iEnd);
        ld = eventPathname;
        headerSize = repmat('%8.8g\t', 1, size(boundAmps,2));
        headerSize = strcat([headerSize(1:end-1) 'n']);
        headerSizeShort = repmat('%8.8g\t', 1, size(AmpsMean,2));
        headerSizeShort = strcat([headerSizeShort(1:end-1) 'n']);
        headerSize2 = repmat('%8.8g\t', 1, length(fileSweeps));
        headerSize2 = strcat([headerSize2(1:end-1) 'n']);
        headerSize3 = repmat('%8.8g\t', 1, size(SD,1));
        headerSize3 = strcat([headerSize3(1:end-1) 'n']);
        fid = fopen(fullfile(eventPathname,eventFilename),'wt+');
        
        fprintf(fid,'\r\nMaximum time to peak (ms):        %8.8g', detectionParameters.SWstart);
        fprintf(fid,'\r\nBaseline duration (ms):           %8.8g', detectionParameters.BLduration);
        fprintf(fid,'\r\nPeak integration period (ms):     %8.8g', detectionParameters.refractoryPeriod);
        fprintf(fid,'\r\nAmplitude lower bound (mV or nA): %8.8g', detectionParameters.Amplobound);
        fprintf(fid,'\r\nAmplitude upper bound (mV or nA): %8.8g', detectionParameters.Ampupbound');
        fprintf(fid,'\r\nGaussian smoothing window (ms):   %8.8g', detectionParameters.smoothWindow);
        fprintf(fid,'\r\nRise time interval:               %8.8s', detectionParameters.RTinterval);
        fprintf(fid,'\r\nRise time bin size:               %8.8g', classificationParameters.riseTimeArray(2));
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Mini EPSP sum of absolute deviations (SAD)');
        fprintf(fid,'\r\nLine 3:  Largest mini EPSP overall deviation (maximum absolute deviation - MAD)');
        fprintf(fid,'\r\nLine 4:  Mini IPSP SAD');
        fprintf(fid,'\r\nLine 5:  Mini IPSP MAD');
        fprintf(fid,'\r\n');
        
        if expand
            fprintf(fid,'\r\n6-score combined SAD 50th centile (corresponds to single score %2.2gth-EPSPs and %2.2gth-IPSPs centiles):\n',...
                [optimData.prct optimData.prctNeg]);
        else
            fprintf(fid,'\r\n3-score combined SAD 50th centile (corresponds to single score %2.2gth-EPSPs and %2.2gth-IPSPs centiles):\n',...
                [optimData.prct optimData.prctNeg]);
        end
        
        fprintf(fid,'\r\n\nAmplitudes:\n');
        fprintf(fid, headerSizeShort, AmpsPrct');
        
        fprintf(fid,'\r\n\nRise times:\n');
        fprintf(fid, headerSizeShort, RTsPrct');
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times:\n');
        fprintf(fid, headerSizeShort, TwoDsPrct');
        
        fprintf(fid,'\r\n\nData medians:\n');
        
        fprintf(fid,'\r\n\nAmplitudes:\n');
        fprintf(fid, headerSizeShort, AmpsMedian');
        
        fprintf(fid,'\r\n\nRise times:\n');
        fprintf(fid, headerSizeShort, RTsMedian');
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times:\n');
        fprintf(fid, headerSizeShort, TwoDsMedian');
        
        fprintf(fid,'\r\n\nData means:\n');
        
        fprintf(fid,'\r\n\nAmplitudes:\n');
        fprintf(fid, headerSizeShort, AmpsMean');
        
        fprintf(fid,'\r\n\nRise times:\n');
        fprintf(fid, headerSizeShort, RTsMean');
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times:\n');
        fprintf(fid, headerSizeShort, TwoDsMean');
        
        fprintf(fid,'\r\n\nData minima:\n');
        
        fprintf(fid,'\r\n\nAmplitudes:\n');
        fprintf(fid, headerSizeShort, AmpsMin');
        
        fprintf(fid,'\r\n\nRise times:\n');
        fprintf(fid, headerSizeShort, RTsMin');
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times:\n');
        fprintf(fid, headerSizeShort, TwoDsMin');
        
        fprintf(fid,'\r\n\nData maxima:\n');
        
        fprintf(fid,'\r\n\nAmplitudes:\n');
        fprintf(fid, headerSizeShort, AmpsMax');
        
        fprintf(fid,'\r\n\nRise times:\n');
        fprintf(fid, headerSizeShort, RTsMax');
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times:\n');
        fprintf(fid, headerSizeShort, TwoDsMax');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Top 50%% mini EPSP amplitude SAD');
        fprintf(fid,'\r\nLine 3:  Top 50%% largest mini EPSP amplitude MAD');
        fprintf(fid,'\r\nLine 4:  Top 50%% mini IPSP amplitude SAD');
        fprintf(fid,'\r\nLine 5:  Top 50%% mini IPSP amplitude MAD');
        fprintf(fid,'\r\n');
        
        if expand
            fprintf(fid,'\r\n6-score combined SAD 50th centile:\n');
        else
            fprintf(fid,'\r\n3-score combined SAD 50th centile:\n');
        end
        
        fprintf(fid, headerSizeShort, AmpsPrctBottom');
        
        fprintf(fid,'\r\n\nData medians:\n');
        
        fprintf(fid, headerSizeShort, AmpsMedianBottom');
        
        fprintf(fid,'\r\n\nData means:\n');
        
        fprintf(fid, headerSizeShort, AmpsMeanBottom');
        
        fprintf(fid,'\r\n\nData minima:\n');
        
        fprintf(fid, headerSizeShort, AmpsMinBottom');
        
        fprintf(fid,'\r\n\nData maxima:\n');
        
        fprintf(fid, headerSizeShort, AmpsMaxBottom');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Top 10%% mini EPSP amplitude SAD');
        fprintf(fid,'\r\nLine 3:  Top 10%% largest mini EPSP amplitude MAD');
        fprintf(fid,'\r\nLine 4:  Top 10%% mini IPSP amplitude SAD');
        fprintf(fid,'\r\nLine 5:  Top 10%% mini IPSP amplitude MAD');
        fprintf(fid,'\r\n');
        
        if expand
            fprintf(fid,'\r\n6-score combined SAD 50th centile:\n');
        else
            fprintf(fid,'\r\n3-score combined SAD 50th centile:\n');
        end
        
        fprintf(fid, headerSizeShort, AmpsPrctMid');
        
        fprintf(fid,'\r\n\nData medians:\n');
        
        fprintf(fid, headerSizeShort, AmpsMedianMid');
        
        fprintf(fid,'\r\n\nData means:\n');
        
        fprintf(fid, headerSizeShort, AmpsMeanMid');
        
        fprintf(fid,'\r\n\nData minima:\n');
        
        fprintf(fid, headerSizeShort, AmpsMinMid');
        
        fprintf(fid,'\r\n\nData maxima:\n');
        
        fprintf(fid, headerSizeShort, AmpsMaxMid');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  top 2%% mini EPSP amplitude SAD');
        fprintf(fid,'\r\nLine 3:  top 2%% largest mini EPSP amplitude MAD');
        fprintf(fid,'\r\nLine 4:  top 2%% mini IPSP amplitude SAD');
        fprintf(fid,'\r\nLine 5:  top 2%% mini IPSP amplitude MAD');
        fprintf(fid,'\r\n');
        
        if expand
            fprintf(fid,'\r\n6-score combined SAD 50th centile:\n');
        else
            fprintf(fid,'\r\n3-score combined SAD 50th centile:\n');
        end
        
        fprintf(fid, headerSizeShort, AmpsPrctTop');
        
        fprintf(fid,'\r\n\nData medians:\n');
        
        fprintf(fid, headerSizeShort, AmpsMedianTop');
        
        fprintf(fid,'\r\n\nData means:\n');
        
        fprintf(fid, headerSizeShort, AmpsMeanTop');
        
        fprintf(fid,'\r\n\nData minima:\n');
        
        fprintf(fid, headerSizeShort, AmpsMinTop');
        
        fprintf(fid,'\r\n\nData maxima:\n');
        
        fprintf(fid, headerSizeShort, AmpsMaxTop');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Mini EPSP SAD bound (model)');
        fprintf(fid,'\r\nLine 3:  Mini EPSP MAD (model)');
        fprintf(fid,'\r\nLine 4:  Mini IPSP SAD bound (model)');
        fprintf(fid,'\r\nLine 5:  Mini IPSP MAD (model)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  Mini EPSP SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 8:  Mini EPSP MAD (confidence-weighted model)');
        fprintf(fid,'\r\nLine 9:  Mini IPSP SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 10: Mini IPSP MAD (confidence-weighted model)');
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\nAmplitudes (model):\n');
        fprintf(fid, headerSize, boundAmps(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmps([1 6:9],:)');
        LSEheader = ['%2.2g\t' '%8.8g\n'];
        lines = [2:5, 7:10];
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLSE]);
        
        fprintf(fid,'\r\n\nRise times (model):\n');
        fprintf(fid, headerSize, boundRTs(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundRTs([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; RTsLSE]);
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times (model):\n');
        fprintf(fid, headerSize, boundTwoDs(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundTwoDs([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; TwoDsLSE]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Top 50%% mini EPSP amplitude SAD bound (model)');
        fprintf(fid,'\r\nLine 3:  Top 50%% mini EPSP amplitude MAD (model)');
        fprintf(fid,'\r\nLine 4:  Top 50%% mini IPSP amplitude SAD bound (model)');
        fprintf(fid,'\r\nLine 5:  Top 50%% mini IPSP amplitude MAD (model)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  Top 50%% mini EPSP amplitude SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 8:  Top 50%% mini EPSP amplitude MAD (confidence-weighted model)');
        fprintf(fid,'\r\nLine 9:  Top 50%% mini IPSP amplitude SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 10: Top 50%% mini IPSP amplitude MAD (confidence-weighted model)');
        fprintf(fid,'\r\n\n');
        
        fprintf(fid, headerSize, boundAmpsBottom(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsBottom([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLSEBottom]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Top 10%% mini EPSP amplitude SAD bound (model)');
        fprintf(fid,'\r\nLine 3:  Top 10%% mini EPSP amplitude MAD (model)');
        fprintf(fid,'\r\nLine 4:  Top 10%% mini IPSP amplitude SAD bound (model)');
        fprintf(fid,'\r\nLine 5:  Top 10%% mini IPSP amplitude MAD (model)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  Top 10%% mini EPSP amplitude SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 8:  Top 10%% mini EPSP amplitude MAD (confidence-weighted model)');
        fprintf(fid,'\r\nLine 9:  Top 10%% mini IPSP amplitude SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 10: Top 10%% mini IPSP amplitude MAD (confidence-weighted model)');
        fprintf(fid,'\r\n\n');
        
        fprintf(fid, headerSize, boundAmpsMid(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsMid([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLSEMid]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  top 2%% mini EPSP amplitude SAD bound (model)');
        fprintf(fid,'\r\nLine 3:  top 2%% mini EPSP amplitude MAD (model)');
        fprintf(fid,'\r\nLine 4:  top 2%% mini IPSP amplitude SAD bound (model)');
        fprintf(fid,'\r\nLine 5:  top 2%% mini IPSP amplitude MAD (model)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  top 2%% mini EPSP amplitude SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 8:  top 2%% mini EPSP amplitude MAD (confidence-weighted model)');
        fprintf(fid,'\r\nLine 9:  top 2%% mini IPSP amplitude SAD bound (confidence-weighted model)');
        fprintf(fid,'\r\nLine 10: top 2%% mini IPSP amplitude MAD (confidence-weighted model)');
        fprintf(fid,'\r\n\n');
        
        fprintf(fid, headerSize, boundAmpsTop(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsTop([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLSETop]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Mini EPSP SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 3:  Mini EPSP MAD (linear fit)');
        fprintf(fid,'\r\nLine 4:  Mini IPSP SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 5:  Mini IPSP MAD (linear fit)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  Mini EPSP SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 8:  Mini EPSP MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 9:  Mini IPSP SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 10: Mini IPSP MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\nAmplitudes (linear fit):\n');
        fprintf(fid, headerSize, boundAmpsLinear(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsLinear([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLinearLSE]);
        
        fprintf(fid,'\r\n\nRise times (linear fit):\n');
        fprintf(fid, headerSize, boundRTsLinear(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundRTsLinear([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; RTsLinearLSE]);
        
        fprintf(fid,'\r\n\nCombined amplitudes and rise times (linear fit):\n');
        fprintf(fid, headerSize, boundTwoDsLinear(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundTwoDsLinear([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; TwoDsLinearLSE]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Top 50%% mini EPSP amplitude SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 3:  Top 50%% mini EPSP amplitude MAD (linear fit)');
        fprintf(fid,'\r\nLine 4:  Top 50%% mini IPSP amplitude SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 5:  Top 50%% mini IPSP amplitude MAD (linear fit)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  Top 50%% mini EPSP amplitude SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 8:  Top 50%% mini EPSP amplitude MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 9:  Top 50%% mini IPSP amplitude SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 10: Top 50%% mini IPSP amplitude MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\n\n');
        
        fprintf(fid, headerSize, boundAmpsLinearBottom(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsLinearBottom([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLinearLSEBottom]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  Top 10%% mini EPSP amplitude SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 3:  Top 10%% mini EPSP amplitude MAD (linear fit)');
        fprintf(fid,'\r\nLine 4:  Top 10%% mini IPSP amplitude SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 5:  Top 10%% mini IPSP amplitude MAD (linear fit)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  Top 10%% mini EPSP amplitude SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 8:  Top 10%% mini EPSP amplitude MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 9:  Top 10%% mini IPSP amplitude SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 10: Top 10%% mini IPSP amplitude MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\n\n');
        
        fprintf(fid, headerSize, boundAmpsLinearMid(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsLinearMid([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLinearLSEMid]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nLine 1:  File size in sweeps');
        fprintf(fid,'\r\nLine 2:  top 2%% mini EPSP amplitude SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 3:  top 2%% mini EPSP amplitude MAD (linear fit)');
        fprintf(fid,'\r\nLine 4:  top 2%% mini IPSP amplitude SAD bound (linear fit)');
        fprintf(fid,'\r\nLine 5:  top 2%% mini IPSP amplitude MAD (linear fit)');
        fprintf(fid,'\r\n');
        fprintf(fid,'\r\nLine 6:  File size in sweeps');
        fprintf(fid,'\r\nLine 7:  top 2%% mini EPSP amplitude SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 8:  top 2%% mini EPSP amplitude MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 9:  top 2%% mini IPSP amplitude SAD bound (confidence-weighted linear fit)');
        fprintf(fid,'\r\nLine 10: top 2%% mini IPSP amplitude MAD (confidence-weighted linear fit)');
        fprintf(fid,'\r\n\n');
        
        fprintf(fid, headerSize, boundAmpsLinearTop(1:5,:)');
        fprintf(fid,'\r\n');
        fprintf(fid, headerSize, boundAmpsLinearTop([1 6:9],:)');
        fprintf(fid,'\r\n\nFitting least square error:\n');
        fprintf(fid, LSEheader, [lines; AmpsLinearLSETop]);
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nMini EPSP amplitude SAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundAmps(2,1) pow]);
        fprintf(fid,'\r\nMini EPSP amplitude MAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundAmps(3,1) pow]);
        fprintf(fid,'\r\nMini IPSP amplitude SAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundAmps(4,1) pow]);
        fprintf(fid,'\r\nMini IPSP amplitude MAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundAmps(5,1) pow]);
        fprintf(fid,'\r\nMini EPSP amplitude SAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundAmps(6,1) powWgh]);
        fprintf(fid,'\r\nMini EPSP amplitude MAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundAmps(7,1) powWgh]);
        fprintf(fid,'\r\nMini IPSP amplitude SAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundAmps(8,1) powWgh]);
        fprintf(fid,'\r\nMini IPSP amplitude MAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundAmps(9,1) powWgh]);
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\nMini EPSP rise time SAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundRTs(2,1) pow]);
        fprintf(fid,'\r\nMini EPSP rise time MAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundRTs(3,1) pow]);
        fprintf(fid,'\r\nMini IPSP rise time SAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundRTs(4,1) pow]);
        fprintf(fid,'\r\nMini IPSP rise time MAD bound equation (model):                                            %16.14g*((n)^%1.14g)',[boundRTs(5,1) pow]);
        fprintf(fid,'\r\nMini EPSP rise time SAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundRTs(6,1) powWgh]);
        fprintf(fid,'\r\nMini EPSP rise time MAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundRTs(7,1) powWgh]);
        fprintf(fid,'\r\nMini IPSP rise time SAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundRTs(8,1) powWgh]);
        fprintf(fid,'\r\nMini IPSP rise time MAD bound equation (confidence-weighted model):                        %16.14g*((n)^%1.14g)',[boundRTs(9,1) powWgh]);
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\nMini EPSP combined amplitude and rise time SAD bound equation (model):                     %16.14g*((n)^%1.14g)',[boundTwoDs(2,1) pow]);
        fprintf(fid,'\r\nMini EPSP combined amplitude and rise time MAD bound equation (model):                     %16.14g*((n)^%1.14g)',[boundTwoDs(3,1) pow]);
        fprintf(fid,'\r\nMini IPSP combined amplitude and rise time SAD bound equation (model):                     %16.14g*((n)^%1.14g)',[boundTwoDs(4,1) pow]);
        fprintf(fid,'\r\nMini IPSP combined amplitude and rise time MAD bound equation (model):                     %16.14g*((n)^%1.14g)',[boundTwoDs(5,1) pow]);
        fprintf(fid,'\r\nMini EPSP combined amplitude and rise time SAD bound equation (confidence-weighted model): %16.14g*((n)^%1.14g)',[boundTwoDs(6,1) powWgh]);
        fprintf(fid,'\r\nMini EPSP combined amplitude and rise time MAD bound equation (confidence-weighted model): %16.14g*((n)^%1.14g)',[boundTwoDs(7,1) powWgh]);
        fprintf(fid,'\r\nMini IPSP combined amplitude and rise time SAD bound equation (confidence-weighted model): %16.14g*((n)^%1.14g)',[boundTwoDs(8,1) powWgh]);
        fprintf(fid,'\r\nMini IPSP combined amplitude and rise time MAD bound equation (confidence-weighted model): %16.14g*((n)^%1.14g)',[boundTwoDs(9,1) powWgh]);
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\nWhere n is the length of the file in sweeps.\n');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsBottom(2,1) pow]);
        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsBottom(3,1) pow]);
        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsBottom(4,1) pow]);
        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsBottom(5,1) pow]);
        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsBottom(6,1) powWgh]);
        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsBottom(7,1) powWgh]);
        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsBottom(8,1) powWgh]);
        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsBottom(9,1) powWgh]);
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsMid(2,1) pow]);
        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsMid(3,1) pow]);
        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsMid(4,1) pow]);
        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsMid(5,1) pow]);
        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsMid(6,1) powWgh]);
        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsMid(7,1) powWgh]);
        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsMid(8,1) powWgh]);
        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsMid(9,1) powWgh]);
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsTop(2,1) pow]);
        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsTop(3,1) pow]);
        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsTop(4,1) pow]);
        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound equation (model):                                   %16.14g*((n)^%1.14g)',[boundAmpsTop(5,1) pow]);
        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsTop(6,1) powWgh]);
        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsTop(7,1) powWgh]);
        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsTop(8,1) powWgh]);
        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound equation (confidence-weighted model):               %16.14g*((n)^%1.14g)',[boundAmpsTop(9,1) powWgh]);
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\nWhere n is the length of the file in sweeps.\n');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nSweep length (ms) and length compensation factor:\n');
        fprintf(fid,'%16.14g\t %16.14g\n',[sLength slFactor]');
        
        fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
        
        fprintf(fid,'\r\nFile names:');
        for nf = 1:size(fileSweeps,1)
            fprintf(fid,'\r\n%s',fileNames{nf});
        end
        fprintf(fid,'\r\n');
        
        fprintf(fid,'\r\nFile sweep count:\n');
        fprintf(fid, headerSize2, fileSweeps);
        
        fprintf(fid,'\r\nCompensation factors:\n');
        fprintf(fid, headerSize2, flFactor);
        
        fprintf(fid,'\r\nTotal standard deviation of the data files:\n');
        fprintf(fid, headerSize3, SD(:,1));
        fprintf(fid,'\r\nStandard deviation of the data files (15ms-window average):\n');
        fprintf(fid, headerSize3, SD(:,2));
        fprintf(fid,'\r\nStandard deviation of the data files (30ms-window average):\n');
        fprintf(fid, headerSize3, SD(:,3));
        
        if estCount
            for est = 1:estCount
                
                fprintf(fid,'\r\n\n****************************************************************************************************************************************************************\n');
                
                fprintf(fid,'\r\nFile sizes:');
                fprintf(fid,'%16.14g\t %16.14g\n', UFest(est).estimate.UF);
                
                if UFest(est).estimate.estimated
                    fprintf(fid,'\r\nSAD and MAD bound estimates and predicted values for scaled files (when noise and target files differ in size)\n');
                    fprintf(fid,'\r\nEstimated values based on data:');
                    if expand
                        fprintf(fid,'\r\nMini EPSP amplitude SAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(1));
                        fprintf(fid,'\r\nMini EPSP amplitude MAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(2));
                        fprintf(fid,'\r\nMini IPSP amplitude SAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(3));
                        fprintf(fid,'\r\nMini IPSP amplitude MAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(4));
                        fprintf(fid,'\r\nMini EPSP rise time SAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(5));
                        fprintf(fid,'\r\nMini EPSP rise time MAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(6));
                        fprintf(fid,'\r\nMini IPSP rise time SAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(7));
                        fprintf(fid,'\r\nMini IPSP rise time MAD bound (6-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(8));
                        fprintf(fid,'\r\nMini EPSP combined SAD bound (6-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(9));
                        fprintf(fid,'\r\nMini EPSP combined MAD bound (6-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(10));
                        fprintf(fid,'\r\nMini IPSP combined SAD bound (6-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(11));
                        fprintf(fid,'\r\nMini IPSP combined MAD bound (6-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(12));
                    else
                        fprintf(fid,'\r\nMini EPSP amplitude SAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(1));
                        fprintf(fid,'\r\nMini EPSP amplitude MAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(2));
                        fprintf(fid,'\r\nMini IPSP amplitude SAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(3));
                        fprintf(fid,'\r\nMini IPSP amplitude MAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(4));
                        fprintf(fid,'\r\nMini EPSP rise time SAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(5));
                        fprintf(fid,'\r\nMini EPSP rise time MAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(6));
                        fprintf(fid,'\r\nMini IPSP rise time SAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(7));
                        fprintf(fid,'\r\nMini IPSP rise time MAD bound (3-score combined SAD 50th percentile):         %16.14g', UFest(est).estimate.measuredUF(8));
                        fprintf(fid,'\r\nMini EPSP combined SAD bound (3-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(9));
                        fprintf(fid,'\r\nMini EPSP combined MAD bound (3-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(10));
                        fprintf(fid,'\r\nMini IPSP combined SAD bound (3-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(11));
                        fprintf(fid,'\r\nMini IPSP combined MAD bound (3-score combined SAD 50th percentile):          %16.14g', UFest(est).estimate.measuredUF(12));
                    end
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nMini EPSP amplitude SAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(13));
                    fprintf(fid,'\r\nMini EPSP amplitude MAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(14));
                    fprintf(fid,'\r\nMini IPSP amplitude SAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(15));
                    fprintf(fid,'\r\nMini IPSP amplitude MAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(16));
                    fprintf(fid,'\r\nMini EPSP rise time SAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(17));
                    fprintf(fid,'\r\nMini EPSP rise time MAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(18));
                    fprintf(fid,'\r\nMini IPSP rise time SAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(19));
                    fprintf(fid,'\r\nMini IPSP rise time MAD bound (min):                                          %16.14g', UFest(est).estimate.measuredUF(20));
                    fprintf(fid,'\r\nMini EPSP combined SAD bound (min):                                           %16.14g', UFest(est).estimate.measuredUF(21));
                    fprintf(fid,'\r\nMini EPSP combined MAD bound (min):                                           %16.14g', UFest(est).estimate.measuredUF(22));
                    fprintf(fid,'\r\nMini IPSP combined SAD bound (min):                                           %16.14g', UFest(est).estimate.measuredUF(21));
                    fprintf(fid,'\r\nMini IPSP combined MAD bound (min):                                           %16.14g', UFest(est).estimate.measuredUF(22));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nMini EPSP amplitude SAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(23));
                    fprintf(fid,'\r\nMini EPSP amplitude MAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(24));
                    fprintf(fid,'\r\nMini IPSP amplitude SAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(25));
                    fprintf(fid,'\r\nMini IPSP amplitude MAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(26));
                    fprintf(fid,'\r\nMini EPSP rise time SAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(27));
                    fprintf(fid,'\r\nMini EPSP rise time MAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(28));
                    fprintf(fid,'\r\nMini IPSP rise time SAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(29));
                    fprintf(fid,'\r\nMini IPSP rise time MAD bound (max):                                          %16.14g', UFest(est).estimate.measuredUF(30));
                    fprintf(fid,'\r\nMini EPSP combined SAD bound (max):                                           %16.14g', UFest(est).estimate.measuredUF(31));
                    fprintf(fid,'\r\nMini EPSP combined MAD bound (max):                                           %16.14g', UFest(est).estimate.measuredUF(32));
                    fprintf(fid,'\r\nMini IPSP combined SAD bound (max):                                           %16.14g', UFest(est).estimate.measuredUF(33));
                    fprintf(fid,'\r\nMini IPSP combined MAD bound (max):                                           %16.14g', UFest(est).estimate.measuredUF(34));
                    fprintf(fid,'\r\n');
                    
                    if expand
                        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(1));
                        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(2));
                        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(3));
                        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(4));
                    else
                        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(1));
                        fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(2));
                        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(3));
                        fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFBottom(4));
                    end
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFBottom(5));
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFBottom(6));
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFBottom(7));
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFBottom(8));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFBottom(9));
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFBottom(10));
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFBottom(11));
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFBottom(12));
                    fprintf(fid,'\r\n');
                    
                    if expand
                        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(1));
                        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(2));
                        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(3));
                        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound (6-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(4));
                    else
                        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(1));
                        fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(2));
                        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(3));
                        fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound (3-score combined SAD 50th percentile): %16.14g', UFest(est).estimate.measuredUFMid(4));
                    end
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFMid(5));
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFMid(6));
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFMid(7));
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound (min):                                  %16.14g', UFest(est).estimate.measuredUFMid(8));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFMid(9));
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFMid(10));
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFMid(11));
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound (max):                                  %16.14g', UFest(est).estimate.measuredUFMid(12));
                    fprintf(fid,'\r\n');
                    
                    if expand
                        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound (6-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(1));
                        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound (6-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(2));
                        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound (6-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(3));
                        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound (6-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(4));
                    else
                        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound (3-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(1));
                        fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound (3-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(2));
                        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound (3-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(3));
                        fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound (3-score combined SAD 50th percentile):  %16.14g', UFest(est).estimate.measuredUFTop(4));
                    end
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound (min):                                   %16.14g', UFest(est).estimate.measuredUFTop(5));
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound (min):                                   %16.14g', UFest(est).estimate.measuredUFTop(6));
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound (min):                                   %16.14g', UFest(est).estimate.measuredUFTop(7));
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound (min):                                   %16.14g', UFest(est).estimate.measuredUFTop(8));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound (max):                                   %16.14g', UFest(est).estimate.measuredUFTop(9));
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound (max):                                   %16.14g', UFest(est).estimate.measuredUFTop(10));
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound (max):                                   %16.14g', UFest(est).estimate.measuredUFTop(11));
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound (max):                                   %16.14g', UFest(est).estimate.measuredUFTop(12));
                    fprintf(fid,'\r\n');
                end
                if ~isempty(UFest(est).estimate.predictedUF)
                    fprintf(fid,'\r\nPredicted values based on theory:\n');
                    fprintf(fid,'\r\nMini EPSP amplitude SAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(1));
                    fprintf(fid,'\r\nMini EPSP amplitude MAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(2));
                    fprintf(fid,'\r\nMini IPSP amplitude SAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(3));
                    fprintf(fid,'\r\nMini IPSP amplitude MAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(4));
                    fprintf(fid,'\r\nMini EPSP rise time SAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(5));
                    fprintf(fid,'\r\nMini EPSP rise time MAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(6));
                    fprintf(fid,'\r\nMini IPSP rise time SAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(7));
                    fprintf(fid,'\r\nMini IPSP rise time MAD bound:                                                %16.14g', UFest(est).estimate.predictedUF(8));
                    fprintf(fid,'\r\nMini EPSP combined SAD bound:                                                 %16.14g', UFest(est).estimate.predictedUF(9));
                    fprintf(fid,'\r\nMini EPSP combined MAD bound:                                                 %16.14g', UFest(est).estimate.predictedUF(10));
                    fprintf(fid,'\r\nMini IPSP combined SAD bound:                                                 %16.14g', UFest(est).estimate.predictedUF(11));
                    fprintf(fid,'\r\nMini IPSP combined MAD bound:                                                 %16.14g', UFest(est).estimate.predictedUF(12));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nMini EPSP amplitude SAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(1));
                    fprintf(fid,'\r\nMini EPSP amplitude MAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(2));
                    fprintf(fid,'\r\nMini IPSP amplitude SAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(3));
                    fprintf(fid,'\r\nMini IPSP amplitude MAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(4));
                    fprintf(fid,'\r\nMini EPSP rise time SAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(5));
                    fprintf(fid,'\r\nMini EPSP rise time MAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(6));
                    fprintf(fid,'\r\nMini IPSP rise time SAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(7));
                    fprintf(fid,'\r\nMini IPSP rise time MAD bound (confidence-weighted):                          %16.14g', UFest(est).estimate.predictedUFWgh(8));
                    fprintf(fid,'\r\nMini EPSP combined SAD bound (confidence-weighted):                           %16.14g', UFest(est).estimate.predictedUFWgh(9));
                    fprintf(fid,'\r\nMini EPSP combined MAD bound (confidence-weighted):                           %16.14g', UFest(est).estimate.predictedUFWgh(10));
                    fprintf(fid,'\r\nMini IPSP combined SAD bound (confidence-weighted):                           %16.14g', UFest(est).estimate.predictedUFWgh(11));
                    fprintf(fid,'\r\nMini IPSP combined MAD bound (confidence-weighted):                           %16.14g', UFest(est).estimate.predictedUFWgh(12));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound:                                        %16.14g', UFest(est).estimate.predictedUFBottom(1));
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound:                                        %16.14g', UFest(est).estimate.predictedUFBottom(2));
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound:                                        %16.14g', UFest(est).estimate.predictedUFBottom(3));
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound:                                        %16.14g', UFest(est).estimate.predictedUFBottom(4));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound:                                        %16.14g', UFest(est).estimate.predictedUFMid(1));
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound:                                        %16.14g', UFest(est).estimate.predictedUFMid(2));
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound:                                        %16.14g', UFest(est).estimate.predictedUFMid(3));
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound:                                        %16.14g', UFest(est).estimate.predictedUFMid(4));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound:                                         %16.14g', UFest(est).estimate.predictedUFTop(1));
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound:                                         %16.14g', UFest(est).estimate.predictedUFTop(2));
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound:                                         %16.14g', UFest(est).estimate.predictedUFTop(3));
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound:                                         %16.14g', UFest(est).estimate.predictedUFTop(4));
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTheoretical model:\n');
                    scaleFactor = UFest(est).estimate.UF(1)/UFest(est).estimate.UF(2);
                    fprintf(fid,'\r\nTheoretical model equations:');
                    fprintf(fid,'\r\nMini EPSP amplitude SAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(2,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini EPSP amplitude MAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(3,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini IPSP amplitude SAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(4,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini IPSP amplitude MAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(5,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini EPSP rise time SAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(2,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini EPSP rise time MAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(3,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini IPSP rise time SAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(4,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini IPSP rise time MAD bound:                                                1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(5,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini EPSP combined SAD bound:                                                 1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(2,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini EPSP combined MAD bound:                                                 1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(3,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini IPSP combined SAD bound:                                                 1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(4,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nMini IPSP combined MAD bound:                                                 1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(5,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nMini EPSP amplitude SAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(6,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini EPSP amplitude MAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(7,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini IPSP amplitude SAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(8,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini IPSP amplitude MAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmps(9,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini EPSP rise time SAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(6,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini EPSP rise time MAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(7,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini IPSP rise time SAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(8,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini IPSP rise time MAD bound (confidence-weighted):                          1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundRTs(9,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini EPSP combined SAD bound (confidence-weighted):                           1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(6,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini EPSP combined MAD bound (confidence-weighted):                           1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(7,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini IPSP combined SAD bound (confidence-weighted):                           1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(8,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\nMini IPSP combined MAD bound (confidence-weighted):                           1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundTwoDs(9,1) UFest(est).estimate.UF(1) powWgh scaleFactor UFest(est).estimate.UF(2) powWgh]);
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude SAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsBottom(2,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nTop 50%% mini EPSP amplitude MAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsBottom(3,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude SAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsBottom(4,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nTop 50%% mini IPSP amplitude MAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsBottom(5,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude SAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsMid(2,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nTop 10%% mini EPSP amplitude MAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsMid(3,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude SAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsMid(4,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\nTop 10%% mini IPSP amplitude MAD bound:                                        1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsMid(5,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\n');
                    
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude SAD bound:                                         1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsTop(2,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\ntop 2%% mini EPSP amplitude MAD bound:                                         1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsTop(3,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude SAD bound:                                         1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsTop(4,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\ntop 2%% mini IPSP amplitude MAD bound:                                         1/2*%5.14g*((%1.14g)^%1.14g + %2.14g*(%1.14g)^%1.14g)',[boundAmpsTop(5,1) UFest(est).estimate.UF(1) pow scaleFactor UFest(est).estimate.UF(2) pow]);
                    fprintf(fid,'\r\n');
                end
                
                fprintf(fid,'\r\n****************************************************************************************************************************************************************\n');
                
                fprintf(fid,'\nPrint-outs for the Excel sheet:\n');
                
                headerSizeExcel = repmat('%16.8g\t', 1, 7);
                headerSizeExcel = strcat([headerSizeExcel(1:end-1) 'n']);
                
                
                % EPSPs:
                fprintf(fid,'\nEPSPs:\n');
                fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n', '', 'Data 50th centile', 'Data min', 'Data max',...
                    'Model', 'Weighted model', 'Linear regular', 'Linear weighted');
                
                % EPSPs.Combined SAD:
                fprintf(fid, '%16s\t', 'Combined SAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(9), UFest(est).estimate.measuredUF(21),...
                    UFest(est).estimate.measuredUF(33), UFest(est).estimate.predictedUF(9), UFest(est).estimate.predictedUFWgh(9),...
                    boundTwoDsLinear(2,UFest(est).estimate.UF(1)), boundTwoDsLinear(6,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, TwoDsLSE(1), TwoDsLSE(5), TwoDsLinearLSE(1), TwoDsLinearLSE(5)]);
                
                % EPSPs.Amplitude SAD:
                fprintf(fid, '%16s\t', 'Amplitude SAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(1), UFest(est).estimate.measuredUF(13),...
                    UFest(est).estimate.measuredUF(25), UFest(est).estimate.predictedUF(1), UFest(est).estimate.predictedUFWgh(1),...
                    boundAmpsLinear(2,UFest(est).estimate.UF(1)), boundAmpsLinear(6,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSE(1), AmpsLSE(5), AmpsLinearLSE(1), AmpsLinearLSE(5)]);
                
                % EPSPs.Amplitude top 50% SAD:
                fprintf(fid, '%16s\t', 'Amplitude top 50% SAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFBottom(1), UFest(est).estimate.measuredUFBottom(5),...
                    UFest(est).estimate.measuredUFBottom(9), UFest(est).estimate.predictedUFBottom(1), UFest(est).estimate.predictedUFWghBottom(1),...
                    boundAmpsLinearBottom(2,UFest(est).estimate.UF(1)), boundAmpsLinearBottom(6,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSEBottom(1), AmpsLSEBottom(5), AmpsLinearLSEBottom(1), AmpsLinearLSEBottom(5)]);
                
                % EPSPs.Amplitude top 10% SAD:
                fprintf(fid, '%16s\t', 'Amplitude top 10% SAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFMid(1), UFest(est).estimate.measuredUFMid(5),...
                    UFest(est).estimate.measuredUFMid(9), UFest(est).estimate.predictedUFMid(1), UFest(est).estimate.predictedUFWghMid(1),...
                    boundAmpsLinearMid(2,UFest(est).estimate.UF(1)), boundAmpsLinearMid(6,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSEMid(1), AmpsLSEMid(5), AmpsLinearLSEMid(1), AmpsLinearLSEMid(5)]);
                
                % EPSPs.Amplitude top 2% SAD:
                fprintf(fid, '%16s\t', 'Amplitude top 2% SAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFTop(1), UFest(est).estimate.measuredUFTop(5),...
                    UFest(est).estimate.measuredUFTop(9), UFest(est).estimate.predictedUFTop(1), UFest(est).estimate.predictedUFWghTop(1),...
                    boundAmpsLinearTop(2,UFest(est).estimate.UF(1)), boundAmpsLinearTop(6,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSETop(1), AmpsLSETop(5), AmpsLinearLSETop(1), AmpsLinearLSETop(5)]);
                
                % EPSPs.Rise time SAD:
                fprintf(fid, '%16s\t', 'Rise time SAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(5), UFest(est).estimate.measuredUF(17),...
                    UFest(est).estimate.measuredUF(29), UFest(est).estimate.predictedUF(5), UFest(est).estimate.predictedUFWgh(5),...
                    boundRTsLinear(2,UFest(est).estimate.UF(1)), boundRTsLinear(6,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, RTsLSE(1), RTsLSE(5), RTsLinearLSE(1), RTsLinearLSE(5)]);
                
                % EPSPs.Combined MAD:
                fprintf(fid, '%16s\t', 'Combined MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(10), UFest(est).estimate.measuredUF(22),...
                    UFest(est).estimate.measuredUF(34), UFest(est).estimate.predictedUF(10), UFest(est).estimate.predictedUFWgh(10),...
                    boundTwoDsLinear(3,UFest(est).estimate.UF(1)), boundTwoDsLinear(7,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, TwoDsLSE(2), TwoDsLSE(6), TwoDsLinearLSE(2), TwoDsLinearLSE(6)]);
                
                % EPSPs.Amplitude MAD:
                fprintf(fid, '%16s\t', 'Amplitude MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(2), UFest(est).estimate.measuredUF(14),...
                    UFest(est).estimate.measuredUF(26), UFest(est).estimate.predictedUF(2), UFest(est).estimate.predictedUFWgh(2),...
                    boundAmpsLinear(3,UFest(est).estimate.UF(1)), boundAmpsLinear(7,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSE(2), AmpsLSE(6), AmpsLinearLSE(2), AmpsLinearLSE(6)]);
                
                % EPSPs.Amplitude top 50% MAD:
                fprintf(fid, '%16s\t', 'Amplitude top 50% MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFBottom(2), UFest(est).estimate.measuredUFBottom(6),...
                    UFest(est).estimate.measuredUFBottom(10), UFest(est).estimate.predictedUFBottom(2), UFest(est).estimate.predictedUFWghBottom(2),...
                    boundAmpsLinearBottom(3,UFest(est).estimate.UF(1)), boundAmpsLinearBottom(7,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSEBottom(2), AmpsLSEBottom(6), AmpsLinearLSEBottom(2), AmpsLinearLSEBottom(6)]);
                
                % EPSPs.Amplitude top 10% MAD:
                fprintf(fid, '%16s\t', 'Amplitude top 10% MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFMid(2), UFest(est).estimate.measuredUFMid(6),...
                    UFest(est).estimate.measuredUFMid(10), UFest(est).estimate.predictedUFMid(2), UFest(est).estimate.predictedUFWghMid(2),...
                    boundAmpsLinearMid(3,UFest(est).estimate.UF(1)), boundAmpsLinearMid(7,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSEMid(2), AmpsLSEMid(6), AmpsLinearLSEMid(2), AmpsLinearLSEMid(6)]);
                
                % EPSPs.Amplitude top 2% MAD:
                fprintf(fid, '%16s\t', 'Amplitude top 2% MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFTop(2), UFest(est).estimate.measuredUFTop(6),...
                    UFest(est).estimate.measuredUFTop(10), UFest(est).estimate.predictedUFTop(2), UFest(est).estimate.predictedUFWghTop(2),...
                    boundAmpsLinearTop(3,UFest(est).estimate.UF(1)), boundAmpsLinearTop(7,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSETop(2), AmpsLSETop(6), AmpsLinearLSETop(2), AmpsLinearLSETop(6)]);
                
                % EPSPs.Rise time MAD:
                fprintf(fid, '%16s\t', 'Rise time MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(6), UFest(est).estimate.measuredUF(18),...
                    UFest(est).estimate.measuredUF(30), UFest(est).estimate.predictedUF(6), UFest(est).estimate.predictedUFWgh(6),...
                    boundRTsLinear(3,UFest(est).estimate.UF(1)), boundRTsLinear(7,UFest(est).estimate.UF(1))]);
                fprintf(fid, '%16s\t', '>> fitting LSE');
                fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, RTsLSE(2), RTsLSE(6), RTsLinearLSE(2), RTsLinearLSE(6)]);
                
                
                % IPSPs:
                %                             fprintf(fid,'\nIPSPs:\n');
                %                             fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n', '', 'Data 50th centile', 'Data min', 'Data max',...
                %                                 'Model', 'Weighted model', 'Linear regular', 'Linear weighted');
                %
                %                             % IPSPs.Combined SAD:
                %                             fprintf(fid, '%16s\t', 'SAD');
                %                             fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(11), UFest(est).estimate.measuredUF(23),...
                %                                 UFest(est).estimate.measuredUF(35), UFest(est).estimate.predictedUF(11), UFest(est).estimate.predictedUFWgh(11),...
                %                                 boundTwoDsLinear(4,UFest(est).estimate.UF(1)), boundTwoDsLinear(8,UFest(est).estimate.UF(1))]);
                %                             fprintf(fid, '%16s\t', '>> fitting LSE');
                %                             fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, TwoDsLSE(3), TwoDsLSE(7), TwoDsLinearLSE(3), TwoDsLinearLSE(7)]);
                %
                %                             % IPSPs.Amplitude SAD:
                %                             fprintf(fid, '%16s\t', 'SAD');
                %                             fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(3), UFest(est).estimate.measuredUF(15),...
                %                                 UFest(est).estimate.measuredUF(27), UFest(est).estimate.predictedUF(3), UFest(est).estimate.predictedUFWgh(3),...
                %                                 boundAmpsLinear(4,UFest(est).estimate.UF(1)), boundAmpsLinear(8,UFest(est).estimate.UF(1))]);
                %                             fprintf(fid, '%16s\t', '>> fitting LSE');
                %                             fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSE(3), AmpsLSE(7), AmpsLinearLSE(3), AmpsLinearLSE(7)]);
                %
                %                             % IPSPs.Rise time SAD:
                %                             fprintf(fid, '%16s\t', 'SAD');
                %                             fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(7), UFest(est).estimate.measuredUF(19),...
                %                                 UFest(est).estimate.measuredUF(31), UFest(est).estimate.predictedUF(7), UFest(est).estimate.predictedUFWgh(7),...
                %                                 boundRTsLinear(4,UFest(est).estimate.UF(1)), boundRTsLinear(8,UFest(est).estimate.UF(1))]);
                %                             fprintf(fid, '%16s\t', '>> fitting LSE');
                %                             fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, RTsLSE(3), RTsLSE(7), RTsLinearLSE(3), RTsLinearLSE(7)]);
                %
                %                             % IPSPs.Combined MAD:
                %                             fprintf(fid, '%16s\t', 'SAD');
                %                             fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(12), UFest(est).estimate.measuredUF(24),...
                %                                 UFest(est).estimate.measuredUF(36), UFest(est).estimate.predictedUF(12), UFest(est).estimate.predictedUFWgh(12),...
                %                                 boundTwoDsLinear(5,UFest(est).estimate.UF(1)), boundTwoDsLinear(9,UFest(est).estimate.UF(1))]);
                %                             fprintf(fid, '%16s\t', '>> fitting LSE');
                %                             fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, TwoDsLSE(4), TwoDsLSE(8), TwoDsLinearLSE(4), TwoDsLinearLSE(8)]);
                %
                %                             % IPSPs.Amplitude MAD:
                %                             fprintf(fid, '%16s\t', 'SAD');
                %                             fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(4), UFest(est).estimate.measuredUF(16),...
                %                                 UFest(est).estimate.measuredUF(28), UFest(est).estimate.predictedUF(4), UFest(est).estimate.predictedUFWgh(4),...
                %                                 boundAmpsLinear(5,UFest(est).estimate.UF(1)), boundAmpsLinear(9,UFest(est).estimate.UF(1))]);
                %                             fprintf(fid, '%16s\t', '>> fitting LSE');
                %                             fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, AmpsLSE(4), AmpsLSE(8), AmpsLinearLSE(4), AmpsLinearLSE(8)]);
                %
                %                             % IPSPs.Rise time MAD:
                %                             fprintf(fid, '%16s\t', 'SAD');
                %                             fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(8), UFest(est).estimate.measuredUF(20),...
                %                                 UFest(est).estimate.measuredUF(32), UFest(est).estimate.predictedUF(8), UFest(est).estimate.predictedUFWgh(8),...
                %                                 boundRTsLinear(5,UFest(est).estimate.UF(1)), boundRTsLinear(9,UFest(est).estimate.UF(1))]);
                %                             fprintf(fid, '%16s\t', '>> fitting LSE');
                %                             fprintf(fid, headerSizeExcel, [NaN, NaN, NaN, RTsLSE(4), RTsLSE(8), RTsLinearLSE(4), RTsLinearLSE(8)]);


                fprintf(fid,'\nPrint-outs for settings entries:\n');
                
                headerSizeExcel = repmat('%16.8g\t', 1, 2);
                headerSizeExcel = strcat([headerSizeExcel(1:end-1) 'n']);
                
                
                % EPSPs:
                fprintf(fid,'\nEPSPs:\n');
                fprintf(fid, '%16s\t %16s\t %16s\n', '', 'SAD', 'MAD');
                
                % EPSPs.Combined SAD & MAD:
                fprintf(fid, '%16s\t', 'Combined SAD & MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(9), UFest(est).estimate.measuredUF(10)]);
                
                % EPSPs.Amplitude SAD & MAD:
                fprintf(fid, '%16s\t', 'Amplitude SAD & MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(1), UFest(est).estimate.measuredUF(2)]);
                
                % EPSPs.Rise time SAD & MAD:
                fprintf(fid, '%16s\t', 'Rise time SAD & MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUF(5), UFest(est).estimate.measuredUF(6)]);

                % EPSPs.Amplitude top 50% SAD & MAD:
                fprintf(fid, '%16s\t', 'Amplitude top 50% SAD & MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFBottom(1), UFest(est).estimate.measuredUFBottom(2)]);
                
                % EPSPs.Amplitude top 10% SAD & MAD:
                fprintf(fid, '%16s\t', 'Amplitude top 10% SAD & MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFMid(1), UFest(est).estimate.measuredUFMid(2)]);
                
                % EPSPs.Amplitude top 2% SAD & MAD:
                fprintf(fid, '%16s\t', 'Amplitude top 2% SAD & MAD');
                fprintf(fid, headerSizeExcel, [UFest(est).estimate.measuredUFTop(1), UFest(est).estimate.measuredUFTop(2)]);
            end
        end
    end
    cd ..
end

button = questdlg('Save the SAD and MAD bound and histogram plots?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    uiwait(msgbox('Note that the figure files will be saved using Matlab "fig" format in a designated folder "error_bounds".','modal'));
    cd(wd);
    if ~exist(dataDir,'dir')
        mkdir(dataDir);
    end
    cd(dataDir);
    % Histograms:
    saveas(dataF(1), 'Mini_EPSP_amplitude_histograms', 'fig');
    saveas(dataF(2), 'Mini_EPSP_rise_time_histograms', 'fig');
    saveas(dataF(3), 'Mini_IPSP_amplitude_histograms', 'fig');
    saveas(dataF(4), 'Mini_IPSP_rise_time_histograms', 'fig');
    % Centiles:
    if expand
        saveas(cntF(1), 'EPSP 6-score combined SAD 50th centile', 'fig');
        saveas(cntF(2), 'IPSP 6-score combined SAD 50th centile', 'fig');
    else
        saveas(cntF(1), 'EPSP 3-score combined SAD 50th centile', 'fig');
        saveas(cntF(2), 'IPSP 3-score combined SAD 50th centile', 'fig');
    end
    if ~isempty(F)
        % Amplitude error fits:
        for iF = 1:length(F)
            filename = strcat(['amplitude_non-linear_fits_' num2str(iF)]);
            saveas(F(iF), filename, 'fig');
        end
        % Top 50% amplitude error fits:
        for iFBottom = 1:length(FBottom)
            filename = strcat(['top_50%_amplitude_non-linear_fits_' num2str(iFBottom)]);
            saveas(FBottom(iFBottom), filename, 'fig');
        end
        % Top 10% amplitude error fits:
        for iFMid = 1:length(FMid)
            filename = strcat(['top_10%_amplitude_non-linear_fits_' num2str(iFMid)]);
            saveas(FMid(iFMid), filename, 'fig');
        end
        % top 2% amplitude error fits:
        for iFTop = 1:length(FTop)
            filename = strcat(['top_2%_amplitude_non-linear_fits_' num2str(iFTop)]);
            saveas(FTop(iFTop), filename, 'fig');
        end
        % Rise time error fits:
        for iG = 1:length(G)
            filename = strcat(['rise_time_non-linear_fits_' num2str(iG)]);
            saveas(G(iG), filename, 'fig');
        end
        % Combined amplitude and rise time error fits:
        for iH = 1:length(H)
            filename = strcat(['combined_non-linear_fits_' num2str(iH)]);
            saveas(H(iH), filename, 'fig');
        end
        
        % Amplitude linear error fits:
        for iA = 1:length(A)
            filename = strcat(['amplitude_linear_fits' num2str(iA)]);
            saveas(A(iA), filename, 'fig');
        end
        % Top 10% amplitude linear error fits:
        for iABottom = 1:length(ABottom)
            filename = strcat(['top_50%_amplitude_linear_fits' num2str(iABottom)]);
            saveas(ABottom(iABottom), filename, 'fig');
        end
        % Top 10% amplitude linear error fits:
        for iAMid = 1:length(AMid)
            filename = strcat(['top_10%_amplitude_linear_fits' num2str(iAMid)]);
            saveas(AMid(iAMid), filename, 'fig');
        end
        % top 2% amplitude linear error fits:
        for iATop = 1:length(ATop)
            filename = strcat(['top_2%_amplitude_linear_fits' num2str(iATop)]);
            saveas(ATop(iATop), filename, 'fig');
        end
        % Rise time linear error fits through the origin:
        for iB = 5:8
            filename = strcat(['rise_time_linear_fits' num2str(iB-4)]);
            saveas(B(iB), filename, 'fig');
        end
        % Combined amplitude and rise time linear error fits through the origin:
        for iC = 5:8
            filename = strcat(['combined_linear_fits' num2str(iC-4)]);
            saveas(C(iC), filename, 'fig');
        end
    end
    cd ..
end

button = questdlg('Close the figure windows?','Close Figure','Yes','No','Yes');
if strcmpi(button, 'Yes')
    if isempty(F)
        close([dataF cntF]);
    else
        close([dataF cntF F G H A B C FBottom ABottom FMid AMid FTop ATop]);
    end
end

button = questdlg('Save the error bound data for optimisation?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    cd(wd);
    [eventFilename, eventPathname, filterIndex] = uiputfile({'*.mat','MAT files (*.mat)'},'Save error bound data for optimisation as', ld);
    if filterIndex
        ld = eventPathname;
        eventFilename = fullfile(eventPathname, eventFilename);
        save(eventFilename,'optimData');
    end
end