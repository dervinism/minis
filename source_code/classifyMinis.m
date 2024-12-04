function [ampHisto, RTHisto, fullHisto] = classifyMinis(amplitudes, riseTimes, classificationParameters)

fullHisto = hist2d(amplitudes, riseTimes, classificationParameters.amplitudeArrayExt(2:end), classificationParameters.riseTimeArrayExt(2:end));
fullHisto = fullHisto(1:end-1,1:end-1);
fullHisto = [zeros(1,size(fullHisto,2)); fullHisto];
fullHisto = [zeros(size(fullHisto,1),1) fullHisto];
ampHisto = sum(fullHisto,1);
RTHisto = sum(fullHisto,2)';
end