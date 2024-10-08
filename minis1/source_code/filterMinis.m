function [V, spectrum, f1, filtfsOut] = filterMinis(V, dt, filterParameters, summaryPlot, varargin)

if nargin > 4
    filename = varargin{1};
    if nargin == 6
        filtfs = varargin{2};
    end
end

if ~exist('filtfs','var') && isfield(filterParameters,'filtfs')
    filtfs = filterParameters.filtfs;
end

% FFT:
notchPlot = false;
if ~isempty(filterParameters.excludedTimes.endPulse)
    endPulse = filterParameters.excludedTimes.endPulse(end);
else
    endPulse = 0;
end
if summaryPlot
    if exist('filename','var') && ~isempty(filename)
        [spectrum, filtSpectrum, ~, f1] = fftMinis(V, dt, filterParameters.nSweeps, endPulse, summaryPlot, notchPlot, filename);
    else
        [spectrum, filtSpectrum, ~, f1] = fftMinis(V, dt, filterParameters.nSweeps, endPulse, summaryPlot, notchPlot);
    end
else
    [spectrum, filtSpectrum] = fftMinis(V, dt, filterParameters.nSweeps, endPulse, summaryPlot, notchPlot);
    f1 = [];
end

% Obtain band-stop frequencies:
if ~strcmpi(filterParameters.state, 'spectrum')
    if ~notchPlot
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        if ~exist('filtfs','var') || isempty(filtfs)
            filtfs = inputdlg('Choose stop band frequencies separated by commas:','FFT',1,{'50, 150'},options);
        end
        filtfsOut = filtfs;
        if ~isempty(filtfs)
            filtfs = cell2mat(filtfs);
            filtfs = strrep(filtfs, ' ', '');
            filtfs = regexp(filtfs, ',*', 'split');
            filtSpectrum(3,1:length(filtfs)) = strarray2numarray(filtfs');
        end
    end
    
    sweepDuration = floor(length(V)/filterParameters.nSweeps);
    excludedTimes = [filterParameters.excludedTimes.startPulse; filterParameters.excludedTimes.endPulse];
    if ~isempty(excludedTimes)
        exclIndsArray = round(1000*excludedTimes/dt) + 1;
    end
    
    % Butterworth filter:
    Rp = .01;                                                               % Passband ripple, dB
    Rs = 3;                                                                 % Stopband attenuation, dB
    for filtF = 1:length(filtfs)
        Wp = [(filtSpectrum(3,filtF)-3)/spectrum(4,end)...                  % Low and High passbands, normalised frequency
            (filtSpectrum(3,filtF)+3)/spectrum(4,end)];
        Ws = [(filtSpectrum(3,filtF)-.5)/spectrum(4,end)...                 % Stopband, normalised frequency
            (filtSpectrum(3,filtF)+.5)/spectrum(4,end)];
        [n, Wn] = buttord(Wp, Ws, Rp, Rs);                                  % n is a filter order
        if n > 2
            n = 2;
        end
        [b, a] = butter(n, Wn, 'stop');
        
        for iSweep = 1:filterParameters.nSweeps
            if ~isempty(excludedTimes)
                for iPart = 1:size(exclIndsArray,2)+1
                    if iPart == 1
                        V((iSweep-1)*sweepDuration+1 : exclIndsArray(1,1)-1) = filtfilt(double(b), double(a),...
                            double(V((iSweep-1)*sweepDuration+1 : exclIndsArray(1,1)-1)));
                    elseif iPart == size(exclIndsArray,2)+1
                        V((iSweep-1)*sweepDuration+exclIndsArray(2,end)+1 : end) = filtfilt(double(b), double(a),...
                            double(V((iSweep-1)*sweepDuration+exclIndsArray(2,end)+1 : end)));
                    else
                        V((iSweep-1)*sweepDuration+exclIndsArray(2,iPart-1)+1 : (iSweep-1)*sweepDuration+exclIndsArray(1,iPart)) = filtfilt(double(b),...
                            double(a), double(V((iSweep-1)*sweepDuration+exclIndsArray(2,iPart-1)+1 : (iSweep-1)*sweepDuration+exclIndsArray(1,iPart))));
                    end
                end
            else
                V((iSweep-1)*sweepDuration+1 : iSweep*sweepDuration) = filtfilt(double(b), double(a),...
                    double(V((iSweep-1)*sweepDuration+1 : iSweep*sweepDuration)));
            end
        end
    end
    
    if summaryPlot
        if exist('f1', 'var') && ~isempty(f1)
            close(f1);
        end
        if exist('filename','var')
            [spectrum, ~, ~, f1] = fftMinis(V, dt, filterParameters.nSweeps, endPulse, summaryPlot, notchPlot, filename);
        else
            [spectrum, ~, ~, f1] = fftMinis(V, dt, filterParameters.nSweeps, endPulse, summaryPlot, notchPlot);
        end
    else
        spectrum = fftMinis(V, dt, filterParameters.nSweeps, endPulse, summaryPlot, notchPlot);
    end
end

% Comment out this part if you don't want to use the high-pass filter:
%     Wp = 50/spectrum(3,end);
%     Ws = 40/spectrum(3,end);
%     [n, Wn] = buttord(Wp, Ws, Rp, Rs);
%     [b, a] = butter(n, Wn, 'high');
%     V = filtfilt(b, a, V);
end



function [fSpectrum, ROIspectrum, samplingF, varargout] = fftMinis(V, dt, nSweeps, endPulse, spectrumPlot, notchPlot, varargin)
% FFTMINIS is a helper subfunction of detectMinis. It performs the Fourier
% time-to-frequency transform of the electrophysiological recording and
% estimates the four strongest frequency components in the range of
% 40-360Hz.
%
%   [FSPECTRUM, ROISPECTRUM, SAMPLINGF] = FFTMINIS(V, dt, spectrumPlot)
%   performs Fourier transform. V is the recording trace after excluded
%   times analysis, mV or nA. DT is the sampling interval in miliseconds.
%   SPECTRUMPLOT is a string variable that controls the display of the
%   frequency vs. amplitude graph. Can be set to either 'on' for dislpaying
%   or 'off' otherwise. The default value is 'off'. NOTCHPLOT is a string
%   variable that controls the display of the stop band frequency marks.
%   Can be set to either 'on' for dislpaying or 'off' otherwise. The
%   default value is 'off'.FSPECTRUM is a matrix composed of four row
%   vectors. The first vector contains amplitudes. The second vector
%   contains a power spectrum which should correspond to decibels (dB). The
%   third vector contains phase information. The fourth vector contains
%   frequencies. ROISPECTRUM is a matrix with three row vectors that are
%   maximum four elements long. The row vectors have the same meaning as
%   for the fSpectrum variable. The columns correspond to the largest
%   amplitude frequencies in the range of 40-360Hz. SAMPLINGF is the
%   sampling frequency in Hertz (Hz).
%
%   [FSPECTRUM, ROISPECTRUM, SAMPLINGF] = FFTMINIS(..., filename)
%   You can suply the data file name as a string. The file name would apear
%   in the title of the frequency spectrum graph.
%
%   [FSPECTRUM, ROISPECTRUM, SAMPLINGF, F1] = FFTMINIS(...)
%   In addition outputs the handle of the frequency spectrum (frequency vs.
%   amplitude) graph.
%


% Initialise variables:
if nargin == 7
    filename = varargin{1};
    if isstruct(filename) && isfield(filename, 'filename')
        filename = filename.filename;
    end
end
V = V - mean(V);                                                            % Subtract the DC component (mean)
samplingF = 1000/dt;
endPulse = endPulse*1000/dt;

% FFT:
% The DC component of V is Amps(1), and Amps((length(V)+1)/2))> is the Nyquist frequency component of x. If nfft is odd, however,
% the Nyquist frequency component is not evaluated, and the number of unique points is (length(V)+1)/2 . This can be generalized
% for both cases to ceil((length(V)+1)/2).
sweepDuration = 2*floor(round(length(V)/nSweeps)/2);
fftDuration = 2*floor((sweepDuration - endPulse - 1)/2);
f = samplingF/2*linspace(0,1,fftDuration/2+1);
Amps = zeros(nSweeps,fftDuration/2+1);
power = Amps;
phase = Amps;
for iSweep = 1:nSweeps
    iAmps = fftshift(fft(V(iSweep*sweepDuration-fftDuration+1:iSweep*sweepDuration)))/fftDuration;
    Amps(iSweep,:) = fliplr(2*abs(iAmps(1:fftDuration/2+1)));               % Multiply by 2 to compensate for the energy loss due to truncation
    power(iSweep,:) = (Amps(iSweep,:).^2*length(f)^2)./(sum(f.^2)*f);
    power(iSweep,1) = 0;
    phase(iSweep,:) = angle(iAmps(1:fftDuration/2+1));
end
clear V
Amps = mean(Amps,1);
power = mean(power,1);
phase = mean(phase,1);
fSpectrum = [Amps; power; phase; f];                                        % Spectrum contains an amplitude, a decibel, a phase, and a frequency vectors

% Plot single-sided amplitude spectrum:
if spectrumPlot
    h = get(0,'CurrentFigure');
    if isempty(h)
        f1 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    else
        figurecount = length(findobj('Type','figure'));
        f1 = figure(figurecount + 1);
        set(f1, 'Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83]);
    end
    loglog(f, power);
    hold on
    if exist('filename','var')
        filestring = sprintf('Single-Sided Amplitude Spectrum of %s', filename);
    else
        filestring = 'Single-Sided Amplitude Spectrum of data';
    end
    set(f1, 'NumberTitle', 'off');
    set(f1, 'Name', 'Fourier Transform (time to frequency domain)');
    title(filestring);
    xlabel('Frequency (Hz)');
    ylabel('Power (mV^2/Hz)');
end

% Estimate large amplitude noise frequencies (within 40-360Hz band):
cutoff = zeros(1,2);
cutoff(1) = find(f>40,1);
cutoff(2) = find(f>360,1)-1;
ampsAdj = Amps(cutoff(1):cutoff(2));
ampsAdj = ampsAdj - 1E-3*(fliplr(1:length(ampsAdj))/length(ampsAdj));
ampsAdj = ampsAdj - .25E-3;
ampsAdj(ampsAdj<0) = 0;
fAdj = f(cutoff(1):cutoff(2));
frequencies = zeros(1,4);
magnitudes = zeros(1,4);
for iBand = 1:4
    [~, ind] = max(ampsAdj);
    if ind
        frequencies(iBand) = fAdj(ind);
        magnitudes(iBand) = Amps(cutoff(1) - 1 + ind);
        ampsAdj(fAdj>frequencies(iBand)-1 & fAdj<frequencies(iBand)+1) = 0;
    end
end
frequencies(frequencies <= 0) = [];
magnitudes(magnitudes <= 0) = [];
ROIspectrum = [magnitudes; 20*log10(magnitudes); frequencies];

% Mark the estimated frequencies:
if spectrumPlot
    if notchPlot
        figure(f1);
        for iFreq = 1:length(frequencies)
            axes1 = get(gcf,'CurrentAxes');
            set(axes1,'XLim',[f(cutoff(1)) f(cutoff(2))]);
            plot(frequencies(iFreq), magnitudes(iFreq), 'o', 'markersize', 10, 'markeredgecolor', 'r');
        end
        hold off
    end
    varargout = {f1};
end
end