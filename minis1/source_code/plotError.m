function [f1, f2] = plotError(detectedEvents1D, detectedEvents1D_RT, classificationParameters, varargin)

f1 = figure('position', [50 50 1200 600]);
f2 = figure('position', [50 50 1200 600]);

figure(f1);
hold on
plot(classificationParameters.amplitudeArray, detectedEvents1D, '.-');
set(f1, 'NumberTitle', 'off');
set(f1, 'Name', 'One-dimensional histograms');
title('Amplitude distribution');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
legend('toggle','NorthEast');
hold off

figure(f2);
hold on
plot(classificationParameters.riseTimeArray, detectedEvents1D_RT,'.-');
title('10-90% rise time distribution');
xlabel('10-90% rise times (ms)');
ylabel('Number of Events');
legend('toggle','NorthEast');
hold off