function tau = doublePeel(pulseDuration, filename, draw, startTime, endTime)
% DOUBLEPEEL estimates the dendritic membrane time constant tau given the
% dendritic membrane potential recording data.
%
%   TAU = DOUBLEPEEL(PEELAFTERPEAK, PEELDURATION, FILENAME, DRAW,...
%       STARTTIME, ENDTIME)
%   estimates the dendritic membrane time constant tau given two estimation
%   parameters, the name of the data file, and the data time interval for
%   estimation. PEELAFTERPEEK is the gap between the peak the start of the
%   peeling period in milliseconds (ms). PEELDURATION is the duration of
%   the peeling period in ms. FILENAME is a string variable containing the
%   full name of the data file inculding the path. DRAW is a string
%   variable that can be set to 'on' or 'off' in order to display figures
%   or not, respectively. STARTTIME is your chosen start time of the
%   averaged sweep trace. ENDTIME is your chosen end time of the averaged
%   sweep trace.
%
%   The function can also be run in a standalone mode. In this case the
%   default parameter values are PEELAFTERPEAK = 5, PEELDURATION = 10,
%   DRAW = 'off'.
%


% Initialise the input variables:
if ~nargin
    draw = true;
    
    
    % Load an abf file:
    if ~exist('filename','var')
        [filename, filepath, filterIndex] = uigetfile({'*.abf','Axon ABF files (*.abf)'},'Choose an abf file', '*.abf');
        if filterIndex
            disp(['User selected file' fullfile(filepath, filename) '  loading...']);
            dataProperties = loadABF(filename, filepath);
        else
            tau = [];
            return
        end
    end
else
    if ischar(filename)
        dataProperties = loadABF(filename);
    else
        dataProperties = filename;
        dataProperties.hd.lActualEpisodes = dataProperties.lActualEpisodes;
    end
end

% Estimate tau:
tau1 = peel(1, pulseDuration, dataProperties.sweep, dataProperties.current, dataProperties.hd.lActualEpisodes, dataProperties.dt, draw, [], [], startTime,...
    endTime);
tau2 = peel(-1, pulseDuration, dataProperties.sweep, dataProperties.current, dataProperties.hd.lActualEpisodes, dataProperties.dt, draw, [], [], startTime,...
    endTime);
tau = (tau1 + tau2)/2;