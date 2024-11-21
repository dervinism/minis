function noiseProperties = loadABF(varargin)
% The local function loadABF loads an Axon Binary Format (*.abf) file with
% electyrophysiological recording data and extracts its properties.

useNativeLoadFunction = false;

if nargin
  filename = varargin{1};
  if nargin == 2
    pathname = varargin{2};
  end
else
  [filename, pathname] = uigetfile({'*.abf','Axon ABF files (*.abf)'},'Choose an abf file', '*.abf');
  disp(['User selected file ' fullfile(pathname, filename) '  loading...' ]);
end


% Extract properties of a noise-alone file:
if strcmpi(filename(end-2:end),'abf')
  if exist('pathname', 'var')
    fileFullName = fullfile(pathname, filename);
  else
    fileFullName = filename;
  end
  if useNativeLoadFunction
    fid = fopen(fileFullName,'rb','ieee-le'); %#ok<*UNRCH> 
    fseek(fid,4,'bof');
    fileVersionNumber = fread(fid,1,'float'); %file format version number

    hd = get_abf_header_info(fileFullName);
    hd.nADCchans = length(hd.fADCProgrammableGain);
    desired_chan = hd.nADCSamplingSeq(1)+1;
    hd.sampled_chan_mask = 0*hd.sampled_chan_mask;
    hd.sampled_chan_mask(desired_chan) = 1;
    if length(hd.nADCSamplingSeq) > 1
      secondary_chan = hd.nADCSamplingSeq(2);
      hd.sampled_chan_mask(secondary_chan+1) = 1;
    end
    %     data_units = hd.sADCUnits(desired_chan+1,:);
    %     for ichar = length(data_units):-1:2  % strip off trailing spaces from units name
    %         if data_units(ichar) == ' '
    %             data_units = data_units(1:end-1);
    %         end
    %     end
    nchans_to_read = length(hd.sampled_chan_mask(hd.sampled_chan_mask>0));  % still would allow us to read >1 channels in future

    %if round(fileVersionNumber, 3) == 1.84
    [abf_data,abf_hdr] = loadabfnew(fileFullName,hd);
    %elseif round(fileVersionNumber, 3) < 2
    %    [abf_data, si] = abfload(fileFullName);
    %else
    %  error('ABF version 2 is not supported');
    %end
    if nchans_to_read > 1 && length(hd.nADCSamplingSeq) > 1
      currentChan = 2;
    else
      currentChan = desired_chan;
    end
  else
    [data,si,hd] = abf2load(fileFullName);

    %desired_chan = find(contains(hd.recChNames,'Vmem'));
    %currentChan = find(contains(hd.recChNames,'Iinj') || contains(hd.recChNames,'Ipatch'));
    desired_chan = find(startsWith(hd.recChNames,'V'));
    currentChan = find(startsWith(hd.recChNames,'I'));
    if isempty(desired_chan)
      warning('Cannot locate voltage recording channel. Expecting voltage channel name to start with letter V');
      if isscalar(currentChan)
        desired_chan = currentChan;
      elseif numel(currentChan) == 2
        desired_chan = currentChan(2);
        currentChan = currentChan(1);
      elseif isempty(currentChan)
        desired_chan = 1;
        currentChan = 2;
      end
    end
    if isempty(currentChan)
      warning('Cannot locate current recording channel. Expecting current channel name to start with letter I');
      if isscalar(desired_chan)
        currentChan = desired_chan;
      elseif numel(desired_chan) == 2
        currentChan = desired_chan(2);
        desired_chan = desired_chan(1);
      end
    end

    if isempty(currentChan) && isscalar(desired_chan)
      currentChan = desired_chan;
    end

    hd.sampled_chan_mask = zeros(size(hd.nADCSamplingSeq(:)'));
    hd.sampled_chan_mask(desired_chan) = 1;
    if length(hd.nADCSamplingSeq) > 1
      secondary_chan = hd.nADCSamplingSeq(2);
      hd.sampled_chan_mask(secondary_chan+1) = 1;
    end
    nchans_to_read = length(hd.sampled_chan_mask(hd.sampled_chan_mask>0));

    dataDims = size(data);
    if numel(size(data)) == 3
      abf_data = [reshape(data(:,desired_chan,:), [dataDims(1)*dataDims(3) 1]) ...
        reshape(data(:,currentChan,:), [dataDims(1)*dataDims(3) 1])];
    else
      abf_data = [reshape(data(:,desired_chan), [dataDims(1) 1]) ...
        reshape(data(:,currentChan), [dataDims(1) 1])];
    end

    abf_hdr = [NaN, 0, si]; 
  end
else
  disp('not an abf file!');
  return
end

noiseProperties.filename = filename;
if exist('pathname', 'var')
  noiseProperties.pathname = pathname;
end
noiseProperties.hd = hd;
noiseProperties.nchans_to_save = nchans_to_read;
noiseProperties.sweep = abf_data(:,desired_chan)';
% may want to store as integers with a scale factor to save memory?
noiseProperties.dt = abf_hdr(3)*0.001;
noiseProperties.tOffset = abf_hdr(2);
noiseProperties.current = abf_data(:,currentChan)';
end



function hd = get_abf_header_info(fname)
% This function reads an abf file.
% 6.5.2000  ssamp should be given as the subsample factor assigned to each ADC channel.
% so subsampling by 100 of eye data on channels 4 and 5 would be indicated
% ssamp=[0 0 0 0 100 100]
% note that all channels are subsampled by the same amount, the max subsample element.

[fid,err] = fopen(fname,'rb','ieee-le');
if (fid < 0)
  error(err);
end

%parameters
hd.nADCchannels=16;  %total number of acquision channels possible on DigiData 2000;
hd.sampled_chan_mask=zeros(1,hd.nADCchannels);
% limit ssamp to subsampling all channels by the same amount.

% make header structure by reading in dummy data
hd.comments = char(fread(fid,100,'char')');
hdr = fread(fid,1,'short');          % n_channels
hdr = [hdr; fread(fid,1,'short')];   % n_timebase
hdr = [hdr; fread(fid,1,'long')];    % n_samp_interval
hdr = [hdr; fread(fid,1,'long')];    % n_samples
hdr = [hdr; fread(fid,8,'short')];   % gain[8]
hdr = [hdr; fread(fid,1,'long')];    %#ok<*NASGU> % time

fseek(fid,0,'bof');
hd.lFileSignature=fread(fid,1,'long');

fseek(fid,4,'bof');
hd.fFileVersionNumber=fread(fid,1,'float'); %file format version number

fseek(fid,8,'bof');
hd.nOperationMode=fread(fid,1,'short'); %3=gap-free  5=episodic-stim

fseek(fid,10,'bof');
hd.lActualAcqLength=fread(fid,1,'long'); %total number of samples (lActualAcqLength)

fseek(fid,14,'bof');
hd.nNumPointsIgnored=fread(fid,1,'short'); %

fseek(fid,16,'bof');
hd.lActualEpisodes=fread(fid,1,'long'); %total number of sweeps

fseek(fid,32,'bof');
hd.fHeaderVersionNumber=fread(fid,1,'float'); %total number of sweeps incl. any average

fseek(fid,36,'bof');
hd.nFileType=fread(fid,1,'short'); %total number of sweeps incl. any average

fseek(fid,38,'bof');
hd.nMSBinFormat=fread(fid,1,'short'); %total number of sweeps incl. any average

fseek(fid,40,'bof'); % 40 bytes into the file
hd.lDataSectionPtr=fread(fid,1,'long'); %Block number for start of Data section (file starts at block zero)
hd.dataoffset=hd.lDataSectionPtr; % offset in blocks of 512 bytes is 40 into file

% now set to correct values (comments left undefined; gain and time never used);

fseek(fid,100,'bof');
hd.nDataFormat=fread(fid,1,'short');

fseek(fid,118,'bof');
hd.channel_count_acquired=fread(fid,1,'short');

fseek(fid,120,'bof');
hd.nADCNumChannels=fread(fid,1,'short');
%number of channels acquired (nADCNumChannels in Axon)
hd.nchan=hd.nADCNumChannels;

hd.nsamp=round(hd.lActualAcqLength/hd.nchan); %divide total data length in abf file by nchan to get per channel length
%By convention, set timebase to 0 to indicate that this is an ABF file
fseek(fid,122,'bof'); %read in fADCSampleInterval
hd.fADCSampleInterval=fread(fid,1,'float');  % us between samples
hd.sampint_per_chan_us=hd.fADCSampleInterval*hd.nchan; %convert to sample interval per channel, in microseconds

fseek(fid,138,'bof'); %
hd.lNumSamplesPerEpisode=fread(fid,1,'long');

fseek(fid,146,'bof'); %
hd.lEpisodesPerRun=fread(fid,1,'long');

fseek(fid,162,'bof'); %
hd.nFirstEpisodeInRun=fread(fid,1,'short');

fseek(fid,178,'bof'); %
hd.fEpisodeStartToStart=fread(fid,1,'float');

fseek(fid,182,'bof'); %
hd.fRunStartToStart=fread(fid,1,'float');

fseek(fid,186,'bof'); %
hd.fTrialStartToStart=fread(fid,1,'float');

fseek(fid,190,'bof'); %
hd.lAverageCount=fread(fid,1,'long');

fseek(fid,194,'bof'); %
hd.lClockChange=fread(fid,1,'long');

fseek(fid,206,'bof'); %
hd.nDataDisplayMode=fread(fid,1,'short');

fseek(fid,208,'bof');
hd.lDisplayAverageUpdate=fread(fid,1,'long'); %number of ADC counts corresponding to fADCRange

fseek(fid,212,'bof');
hd.nChannelStatsStrategy=fread(fid,1,'short'); %number of ADC counts corresponding to fADCRange

fseek(fid,218,'bof');
hd.lSamplesPerTrace=fread(fid,1,'long'); %number of ADC counts corresponding to fADCRange

fseek(fid,222,'bof');
hd.lStartDisplayNum=fread(fid,1,'long'); %number of ADC counts corresponding to fADCRange

fseek(fid,226,'bof');
hd.lFinishDisplayNum=fread(fid,1,'long'); %number of ADC counts corresponding to fADCRange

% Read range and resolution
fseek(fid,244,'bof');
hd.fADCRange=fread(fid,1,'float'); %positive full-scale input in volts

fseek(fid,252,'bof');
hd.lADCResolution=fread(fid,1,'long'); %number of ADC counts corresponding to fADCRange

fseek(fid,260,'bof');
hd.nExperimentType=fread(fid,1,'short'); % 0 = voltage clamp

% Determine sequence in which data was stored and establish a ssampget vector

fseek(fid,378,'bof');
hd.nADCPtoLChannelMap =  fread(fid,hd.nADCchannels,'short')'; %

fseek(fid,410,'bof');
for j=1:hd.nchan % nADCchannels
  hd.nADCSamplingSeq(j)=fread(fid,1,'short'); %list the order in which ADC channels were sampled
  hd.sampled_chan_mask(1 + hd.nADCSamplingSeq(j) )=1;
  % something like 0 2  6  ie the ADC channel numbers of all the channels actually sampled, in order
end

fseek(fid,442,'bof');
hd.sADCChannelName=[];
for j=1: hd.nADCchannels
  chan_name='';
  % could do all this a lot more efficiently e.g. abfload.m way
  chan_name =  char(fread(fid,10,'uchar'))'; %#ok<*FREAD> %
  hd.sADCChannelName= [hd.sADCChannelName; chan_name ];
end

fseek(fid,602,'bof');
hd.sADCUnits=[];
for j=1: hd.nADCchannels
  chan_unit='';
  chan_unit =  char(fread(fid,8,'uchar'))'; %
  hd.sADCUnits= [hd.sADCUnits; chan_unit ]; %could do all this a lot more efficiently e.g. abfload.m way
end

fseek(fid,730,'bof');
hd.fADCProgrammableGain =  fread(fid,hd.nADCchannels,'float')'; %

fseek(fid,794,'bof');
hd.fADCDisplayAmplification =  fread(fid,hd.nADCchannels,'float')'; %  !!!!!

fseek(fid,7858,'bof');
hd.fADCDisplayOffset =  fread(fid,hd.nADCchannels,'float')'; % !!!!!!

fseek(fid,922,'bof');
hd.fInstrumentScaleFactor =  fread(fid,hd.nADCchannels,'float')';
fseek(fid,986,'bof');
hd.fInstrumentOffset =  fread(fid,hd.nADCchannels,'float')'; %

fseek(fid,1050,'bof');
hd.fSignalGain =  fread(fid,hd.nADCchannels,'float')'; %

fseek(fid,1114,'bof');
hd.fSignalOffset =  fread(fid,hd.nADCchannels,'float')'; %

fseek(fid,1378,'bof');
hd.fDACScalefactor =  fread(fid,1,'float')';%

fseek(fid,1394,'bof');
hd.fDACHoldingLevel =  fread(fid,1,'float')';%

fseek(fid,2034,'bof');
hd.lHeaderSize =  fread(fid,1,'long')';%

fseek(fid,4576,'bof');
hd.fTelegraphAdditGain =  fread(fid,hd.nADCchannels,'float')'; %

hd.scalef=(hd.fADCRange./(hd.lADCResolution.*hd.fInstrumentScaleFactor.*hd.fADCProgrammableGain.*hd.fTelegraphAdditGain) );

fseek(fid,0,'bof'); % go back to start
hd.whole_header_as_bytes=uint8(fread(fid,hd.lHeaderSize,'uchar')); % read all the bytes in the header en masse
% so now have a byte-wise copy of the entire header, to use as a
% 'background' when writing abf files, if need be
fclose(fid);
end



function [x,hdr,comments] = loadabfnew(fname,hd)
% This function reads in a .abf file.
% 6.5.2000  ssamp should be given as the subsample factor assigned to each ADC channel.
% so subsampling by 100 of eye data on channels 4 and 5 would be indicated
% ssamp=[0 0 0 0 100 100]
% note that all channels are subsample by the same amount, the max subsample element.
% Modifications by EA, 6.6.2000
ssamp=hd.sampled_chan_mask;
[fid,err] = fopen(fname,'rb','ieee-le');
if (fid < 0)
  error(err);
end

%parameters
nADCchannels=16;  %total number of acquision channels possible on DigiData 2000;

% limit ssamp to subsampling all channels by the same amount.
ssamp=max(ssamp)*(ssamp>0);

% make header structure by reading in dummy data
comments = setstr(fread(fid,100,'char')'); %#ok<*DSTSTR>
hdr = fread(fid,1,'short');          % n_channels
hdr = [hdr; fread(fid,1,'short')];   % n_timebase
hdr = [hdr; fread(fid,1,'long')];    % n_samp_interval
hdr = [hdr; fread(fid,1,'long')];    % n_samples
hdr = [hdr; fread(fid,8,'short')];   % gain[8]
hdr = [hdr; fread(fid,1,'long')];    % time

% now set to correct values (comments left undefined; gain and time never used);
fseek(fid,120,'bof');
hdr(1)=fread(fid,1,'short');
nchan = hdr(1); %number of channels acquired (nADCNumChannels in Axon)
fseek(fid,10,'bof');
hdr(4)=fread(fid,1,'long'); %total number of samples (lActualAcqLength)
hdr(4)=hdr(4)/hdr(1); %divide total data length in abf file by nchan to get per channel length
nsamp = hdr(4);
%By convention, set timebase to 0 to indicate that this is an ABF file and does not
%necessarily conform to the timebase = 5 or 10 msec rule.
hdr(2)=0;
fseek(fid,122,'bof'); %read in fADCSampleInterval i.e. clock time between sample points in cosnecutive (interleaved) channels
hdr(3)=fread(fid,1,'float');
hdr(3)=hdr(3)*hdr(1); %convert to sample interval per channel, in microseconds

% Determine sequence in which data was stored and establish a ssampget vector
fseek(fid,410,'bof');
for j=1:nADCchannels
  nADCSamplingSeq(j)=fread(fid,1,'short'); %#ok<*AGROW> %list the order in which ADC channels were sampled
  % returns -1 after all sampled channels listed
end
jx=1;ChannelID=[];
for j=1:length(ssamp)
  if ssamp(j)>0
    ChannelID=find(nADCSamplingSeq==(j-1));  %j-1 since ssamp(1) corresponds to ADC channel zero
    if ~isempty(ChannelID)
      PositionStored(jx)=ChannelID; % identifies the location of the desired channel in the data storage section.
      jx=jx+1;
    else
      error(['ADC channel ' num2str(j) 'was not sampled. Restart.'])
    end
    ChannelID=[];
  end
end

% Set data matrix dimensions
ss = max(ssamp);
ns = ceil(nsamp/ss);
nc = sum(ssamp > 0); % number of channels being read
x = zeros(ns,nc,'single');

% Read range and resolution
fseek(fid,244,'bof');
fADCRange=fread(fid,1,'float'); %positive full-scale input in volts
fseek(fid,252,'bof');
lADCResolution=fread(fid,1,'long'); %number of ADC counts corresponding to fADCRange
% as of 6/5/2000, all files are 10.24 Volts fADCRange and 32768 (or 2^15) lADCResolution.
% aka 16 bit integer for +/- 10 Volts.

% Read in datablock offset for 512 byte blocks
fseek(fid,40,'bof');
dataoffset=fread(fid,1,'long'); %Block number for start of Data section (file starts at block zero)
%position at datablock start
fseek(fid,(dataoffset*512),'bof'); %blocks 512 in length.
% as of 6/5/2000, header version 1.5, header was 4 blocks (2048) in length, and dataoffset
% was 4

% Read data an scale factor information
% index i keeps track of the ADC channel number
% index ix keeps track of the ADC channels with data
% PositionStored takes into account the sampling sequence
ix = 1;

for i=1:length(ssamp) %walk through ADC channels 0:(length(ssamp)-1)
  if ssamp(i) > 0
    fseek(fid,(hd.dataoffset*512)+(2*(PositionStored(ix)-1)),'bof');
    x(:,ix) = fread(fid,ns,'short',(nchan*(ssamp(i)-1)*2)+((nchan-1)*2));
    fseek(fid,922+(4*(i-1)),'bof');
    fInstrumentScaleFactor=fread(fid,1,'float');
    fseek(fid,986+(4*(i-1)),'bof');
    fInstrumentOffset=fread(fid,1,'float');
    fseek(fid,1114+(4*(i-1)),'bof');
    fSignalOffset=fread(fid,1,'float');
    fseek(fid,730+(4*(i-1)),'bof');
    fADCProgrammableGain=fread(fid,1,'float');
    fseek(fid,4576+(4*(i-1)),'bof');
    fTelegraphAdditGain =  fread(fid,1,'float');
    if ~fTelegraphAdditGain
      scalef=(fADCRange/(lADCResolution*fInstrumentScaleFactor*fADCProgrammableGain) );
    else
      scalef=(fADCRange/(lADCResolution*fInstrumentScaleFactor*fADCProgrammableGain*fTelegraphAdditGain) );
    end
    x(:,ix) = x(:,ix)*scalef + fInstrumentOffset + fSignalOffset;
    ix = ix+1;
  end
end
fclose(fid);
end