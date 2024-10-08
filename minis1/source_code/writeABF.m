function err = writeABF(channelData, filename, sampleRateHz, units)
% err = writeABF(channelData, filename, sampleRateHz, units)
%
% Function takes row matrix data and saves it as an Axon binary file (ABF)
% version 1. Rows of the matrix must correspond to individual recording
% channel data.
% Input: channelData - a row matrix (each row is a channel). A row or a
%                      column vector is treated as single channel data.
%        filename - a string with a file name.
%        sampleRateHz - a scalar with data sampling frequency in Hz.
%        units - a cell array of unit strings. Default is {'mV'}.
% Output: err - a system-dependent error message if writeABF fails to
%               create a writable file. Otherwise, err is an empty
%               character vector.


%% Parse function input
assert(isnumeric(channelData));
if size(channelData,2) == 1
  channelData = channelData';
end

assert(ischar(filename));

assert(isscalar(sampleRateHz));

if nargin < 4
  units = {mV'};
elseif ~iscell(units)
  assert(ischar(units));
  units = {units};
else
  for i = 1:numel(units)
    assert(ischar(units{i}));
  end
end


%% Open a binary file for writing data into
[fid, err] = fopen(filename,'wb','ieee-le');
if (fid < 0)
  error(err);
end


%% Define constants for ABF1 files
blockSize = 512;
headerBlocks = 4;
maxNumCh = 16;


%% Determine dimensions of data
channelCount = size(channelData, 1);
channelPointCount = size(channelData, 2);
dataPointCount = channelPointCount*channelCount;


%% Calculate how large our file must be and create a byte array file of that size
bytesPerPoint = 2;
dataBlocks = dataPointCount * bytesPerPoint / blockSize + 1;
arraySize = (dataBlocks + headerBlocks) * blockSize;
fseek(fid,0,  'bof'); fwrite(fid,zeros(1,arraySize),'uint8');


%% Populate only the useful header data values
fseek(fid,0,  'bof'); fwrite(fid,'ABF ','uchar'); % fFileSignature
fseek(fid,4,  'bof'); fwrite(fid,1.3,'float'); % fFileVersionNumber
fseek(fid,8,  'bof'); fwrite(fid,5,'short'); % nOperationMode (5 is episodic)
fseek(fid,10, 'bof'); fwrite(fid,dataPointCount,'long'); % lActualAcqLength
fseek(fid,16, 'bof'); fwrite(fid,1,'long'); % lActualEpisodes
fseek(fid,40, 'bof'); fwrite(fid,headerBlocks,'long'); % lDataSectionPtr
fseek(fid,100,'bof'); fwrite(fid,0,'short'); % nDataFormat is 1 for float32
fseek(fid,122,'bof'); fwrite(fid,(1e6/sampleRateHz)/channelCount,'float'); % fADCSampleInterval
fseek(fid,138,'bof'); fwrite(fid,dataPointCount,'long'); % lNumSamplesPerEpisode


%% Populate header data values relating to channels
fseek(fid,120,'bof'); fwrite(fid,channelCount,'short'); % nADCNumChannels
for ch = 1:maxNumCh
	if ch <= channelCount
    fseek(fid, 378 + 2*(ch-1), 'bof'); fwrite(fid,ch-1,'short'); % nADCPtoLChannelMap
    fseek(fid, 410 + 2*(ch-1), 'bof'); fwrite(fid,ch-1,'short'); % nADCSamplingSeq
  else
    fseek(fid, 378 + 2*(ch-1), 'bof'); fwrite(fid,-1,'short');
    fseek(fid, 410 + 2*(ch-1), 'bof'); fwrite(fid,-1,'short');
  end
end


%% Determine the peak data deviation from zero
maxVal = zeros(1,channelCount);
for ch = 1:channelCount
  maxVal(ch) = max(abs(channelData(ch,:)));
  if maxVal(ch) == 0
    maxVal(ch) = 1e-9;
  end
end


%% ADC adjustments

% These ADC adjustments are used for integer conversion. It's a good idea
% to populate these with non-zero values even when using float32 notation
% to avoid divide-by-zero errors when loading ABFs.

fSignalGain = 1; % always 1
for ch = 1:maxNumCh
  if ch <= channelCount
    fADCProgrammableGain(ch) = (10/10^(floor(log10(maxVal(ch)))));
  else
    fADCProgrammableGain(ch) = 1;
  end
end
lADCResolution = 2^15; % 16-bit signed = +/- 32768


%% Set the scaling factor to be the biggest allowable to accommodate the data
fADCRange = 10;
valueScale = zeros(1,channelCount);
for ch = 1:channelCount
  fInstrumentScaleFactor = 1;
  for i = 1:10
    fInstrumentScaleFactor = fInstrumentScaleFactor/10;
    valueScaleI = fADCProgrammableGain(ch) * lADCResolution / fADCRange * fInstrumentScaleFactor;
    maxDeviationFromZero = 32767 / valueScaleI;
    if maxDeviationFromZero >= maxVal(ch)
      valueScale(ch) = valueScaleI;
      break
    end
  end
end


%% Prepare units as a space-padded 8-byte string
unitString = {};
for ch = 1:channelCount
  unitStringCh = units{ch};
    while numel(unitStringCh) < 8
      unitStringCh = [unitStringCh ' '];
    end
    unitString{ch} = unitStringCh; %#ok<*AGROW>
end


%% Store the scale data in the header
fseek(fid,252,'bof'); fwrite(fid,lADCResolution,'long');
fseek(fid,244,'bof'); fwrite(fid,fADCRange,'float');
for ch = 1:maxNumCh
  fseek(fid, 922+(ch-1)*4,  'bof'); fwrite(fid,fInstrumentScaleFactor,'float');
  fseek(fid, 1050+(ch-1)*4, 'bof'); fwrite(fid,fSignalGain,'float');
  fseek(fid, 730+(ch-1)*4,  'bof'); fwrite(fid,fADCProgrammableGain(ch),'float');
  fseek(fid, 602+(ch-1)*8,  'bof'); fwrite(fid,unitString{min([ch channelCount])},'uchar');
end


%% Interleave signals
channelData = reshape(channelData, [1 dataPointCount]);
    

%% Scale signals
valueScale = repmat(valueScale, [1 dataPointCount/numel(valueScale)]);
channelData = channelData.*valueScale;
    

%% Fill the rest of data with interleaved and scaled signals
dataByteOffset = blockSize * headerBlocks;
fseek(fid,dataByteOffset,'bof'); fwrite(fid,channelData,'short'); % signals

fclose(fid);