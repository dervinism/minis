function [f1, f2, f3, f4] = plotCompare(amplitudeArray, riseTimeArray, target, current, RTtype, mnSpectra, nSpectra, f, varargin)
% PLOTCOMPARE is a helper function of the COMPAREMINIS function used in
% MINIS GUI. For more information see the code.
%



if nargin == 11
    f1 = varargin{1};
    f2 = varargin{2};
    f3 = varargin{3};
	f4 = varargin{4};
else
    f1 = figure('position', [50 50 1200 600]);
    f2 = figure('position', [50 50 1200 600]);
    f3 = figure('position', [50 50 1200 600]);
	f4 = figure('position', [50 50 1200 600]);
end
ampLim2D = (amplitudeArray(2) - amplitudeArray(1))/2;
RTLim2D = (riseTimeArray(2) - riseTimeArray(1))/2;



%% Drawing uni-dimensional histos:
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
p3 = plot(amplitudeArray, sum(target-current,1),'g.-');
xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
set(f1, 'NumberTitle', 'off');
set(f1, 'Name', 'One-dimensional Amplitude Histograms');
title('Amplitude distributions');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
legend([p1 p2 p3],'Target','Noise','Target-noise','Location','NorthEast');
hold off

figure(f2)
targetRTcs = cumsum(sum(target,2))/sum(sum(target));
targetRT99prc = find(targetRTcs > .99, 1);
currentRTcs = cumsum(sum(current,2))/sum(sum(current));
currentRT99prc = find(currentRTcs > .99, 1);
iEndRT = min([max([targetRT99prc currentRT99prc])+1 length(riseTimeArray)]);
plot([0 riseTimeArray(end)], [0 0], 'k-');
hold on
p1 = plot(riseTimeArray, sum(target,2),'b.-');
p2 = plot(riseTimeArray, sum(current,2),'r.-');
p3 = plot(riseTimeArray, sum(target-current,2),'g.-');
xlim([riseTimeArray(1) riseTimeArray(iEndRT)]);
set(f2, 'NumberTitle', 'off');
set(f2, 'Name', 'One-dimensional Rise Time Histograms');
if strcmpi(RTtype, '10-90%')
    title('10-90% rise time distributions');
    xlabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    title('20-80% rise time distributions');
    xlabel('20-80% rise times (ms)');
end
ylabel('Number of Events');
legend([p1 p2 p3],'Target','Noise','Target-noise','Location','NorthEast');
hold off



%% Drawing two-dimensional histos:
figure(f3);
subplot(2,2,1,'replace');
hold on
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
imagesc(amplitudeArray, riseTimeArray, target);
set(gca,'YDir','normal');
set(f3, 'NumberTitle', 'off');
set(f3, 'Name', 'Two-dimensional Histograms');
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
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
imagesc(amplitudeArray, riseTimeArray, current);
set(gca,'YDir','normal');
title('Noise 2D event distribution');
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
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
imagesc(amplitudeArray, riseTimeArray, target-current);
set(gca,'YDir','normal');
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
xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
surf(amplitudeArray, riseTimeArray, target-current);
shading interp;
set(gca,'YDir','normal');
title('Discrepancy of 2D event distributions: 3D view');
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



%% Drawing power and phase spectra:
smoothWindow = 40;
if length(nSpectra(1,:)) < length(mnSpectra(1,:))
    mnSpectra = mnSpectra(:,1:length(nSpectra(1,:)));
elseif length(mnSpectra(1,:)) < length(nSpectra(1,:))
    nSpectra = nSpectra(:,1:length(mnSpectra(1,:)));
    f = f(1:length(mnSpectra(1,:)));
end
nSpectrum = nSpectra(1,:);
mnSpectrum = mnSpectra(1,:);
% nSpectrumPh = nSpectra(2,:);
% mnSpectrumPh = mnSpectra(2,:);

% Drawing power spectra:
figure(f4);
subplot(2,2,1,'replace');
loglog(f, nSpectrum);
hold on
loglog(f, mnSpectrum,'r');
filestring = 'Power spectra of data';
set(f4, 'NumberTitle', 'off');
set(f4, 'Name', 'Fourier Transform (time to frequency domain)');
title(filestring);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');
legend('Noise','Target','Location','SouthWest');
hold off

subplot(2,2,2,'replace');
nSpecSmooth = gsmooth(nSpectrum, smoothWindow, 3);
mnSpecSmooth = gsmooth(mnSpectrum, smoothWindow, 3);
loglog(f, abs(nSpecSmooth));
hold on
loglog(f, abs(mnSpecSmooth),'r');
filestring = 'Power spectra of data (smoothed)';
title(filestring);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');
legend('Noise','Target','Location','SouthWest');
hold off

subplot(2,2,3,'replace');
loglog(f, abs(mnSpectrum - nSpectrum), 'g');
hold on
filestring = 'Difference of power spectra';
title(filestring);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');
legend('Difference modulus','Location','SouthWest');
hold off

subplot(2,2,4,'replace');
diffSmooth = gsmooth(mnSpectrum - nSpectrum, smoothWindow, 3);
loglog(f, abs(diffSmooth), 'g');
hold on
filestring = 'Difference of power spectra (smoothed)';
title(filestring);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');
legend('Difference modulus','Location','SouthWest');
hold off

% Drawing phase spectra:
% figure(f5);
% subplot(2,2,1,'replace');
% loglog(f, nSpectrumPh);
% hold on
% loglog(f, mnSpectrumPh,'r');
% filestring = 'Phase spectra of data';
% set(f5, 'NumberTitle', 'off');
% set(f5, 'Name', 'Fourier Transform (time to frequency domain)');
% title(filestring);
% xlabel('Frequency (Hz)');
% ylabel('Angle (rad)');
% legend('Simulated','Target','Location','SouthWest');
% hold off
% 
% subplot(2,2,2,'replace');
% nSpecSmooth = gsmooth(nSpectrumPh, smoothWindow, 3);
% mnSpecSmooth = gsmooth(mnSpectrumPh, smoothWindow, 3);
% loglog(f, abs(nSpecSmooth));
% hold on
% loglog(f, abs(mnSpecSmooth),'r');
% filestring = 'Phase spectra of data (smoothed)';
% title(filestring);
% xlabel('Frequency (Hz)');
% ylabel('Angle (rad)');
% legend('Simulated','Target','Location','SouthWest');
% hold off
% 
% subplot(2,2,3,'replace');
% loglog(f, abs(mnSpectrumPh - nSpectrumPh), 'g');
% hold on
% filestring = 'Difference of phase spectra';
% title(filestring);
% xlabel('Frequency (Hz)');
% ylabel('Angle (rad)');
% legend('Difference modulus','Location','SouthWest');
% hold off
% 
% subplot(2,2,4,'replace');
% diffSmooth = gsmooth(mnSpectrumPh - nSpectrumPh, smoothWindow, 3);
% loglog(f, abs(diffSmooth), 'g');
% hold on
% filestring = 'Difference of phase spectra (smoothed)';
% title(filestring);
% xlabel('Frequency (Hz)');
% ylabel('Angle (rad)');
% legend('Difference modulus','Location','SouthWest');
% hold off