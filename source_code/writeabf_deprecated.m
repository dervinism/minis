function err = writeabf(fname,hd,x)
x=x';
ssamp=hd.sampled_chan_mask;

[fid,err] = fopen(fname,'wb','ieee-le');
if (fid < 0)
    error(err);
end

%write out the whole original header as a template
count=fwrite(fid,hd.whole_header_as_bytes,'uchar'); % hdr = fread(fid,1,'short');          % n_channel
%parameters

% limit ssamp to subsampling all channels by the same amount.

n_chans_to_write=size(x,2);
n_pts_per_chan=size(x,1); %#ok<*NASGU>
hd.nchan=n_chans_to_write;

jx=1;ChannelID=[];
for j=1: length(ssamp) %n_chans_to_write
    if ssamp(j)>0
        ChannelID=find(hd.nADCSamplingSeq==(j-1));  %j-1 since ssamp(1) corresponds to ADC channel zero
        if ~isempty(ChannelID)
            PositionStored(jx)=ChannelID; %#ok<*AGROW> % identifies the location of the desired channel in the data storage section.
            jx=jx+1;
        else
            error(['ADC channel ' num2str(j) 'was not sampled. Restart.']);
        end
        ChannelID=[];
    end
end

ix=1;
for i=1:length(ssamp) %walk through ADC channels 0:(length(ssamp)-1)
    if ssamp(i) > 0
        nskip=((hd.nADCNumChannels-1)*2); % it skips *before* it writes, not after (i.e opposite of fread)
        %start_byte=(hd.dataoffset*512)+(2*(PositionStored(ix)-1)) - nskip;
        start_byte=6196+(2*(PositionStored(ix)-1)) - nskip;
        fseek(fid,start_byte,'bof'); %data is 2 byte short and multiplexed.
        % as in, interleaved
        % furthermore, skipped ADC channels are not stored as data.
        % debug this better: only works if v in 0th chan
        %xx(:,ix)=( round((x(:,ix) - hd.fSignalOffset(ix) - hd.fInstrumentOffset(ix) )/hd.scalef(ix) ) ); %???? GM 27/5/11 !!!!
        xx(:,ix)= (x(:,ix) - hd.fSignalOffset(ix) - hd.fInstrumentOffset(ix))/hd.scalef(ix)/0.01;
        intxx(:,ix)=int16(xx(:,ix));
        % 'int16'  is usuallly same as 'short' = 2byte =16bit integer +-302768
        count = fwrite(fid,intxx(:,ix),'short',nskip); % last arg is bytes to skip *before* each
        %count = fwrite(fid,intxx(:,ix),'int16',nskip);
        % int16 need to convert it first to short with all right gains and offsets
        ix = ix+1;
    end
end
fclose(fid);
end