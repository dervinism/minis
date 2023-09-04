function [ampHisto, RTHisto, fullHisto] = classifyMinisNeg(amplitudes, riseTimes, classificationParameters)

fullHisto = fliplr(hist2d(amplitudes, riseTimes, -fliplr(classificationParameters.amplitudeArrayExt), classificationParameters.riseTimeArrayExt));
fullHisto = fullHisto(1:end-1,1:end-1);
ampHisto = sum(fullHisto,1);
RTHisto = sum(fullHisto,2)';
end