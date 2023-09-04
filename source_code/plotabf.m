function [fV, fI] = plotabf(varargin)
% [fV, fI] = plotabf(filename)
%
% Function displays membrane potential data stored in an Axon Binary File.
% Input: filename - file name (optional). If no filename is supplied the
%                   function prompts you to navigate to the file.
% Output: fV - membrane potential figure handle.
%         fI - injected current figure handle.

% Load the data
if nargin
    dataProperties = loadABF(varargin{1});
else
    dataProperties = loadABF();
end

V = dataProperties.sweep;
I = dataProperties.current;
assert(numel(V) == numel(I));
t = dataProperties.dt:dataProperties.dt:dataProperties.dt*numel(V);

% Display the membrane potential data
options.nameString = [dataProperties.filename ': Membrane potential'];
options.titleString = options.nameString;
options.dataType = 'Membrane potential';
options.dataUnits = 'mV';
fV = plotData(t, V, options);

% Display injected current data
options.nameString = [dataProperties.filename ': Injected current'];
options.titleString = options.nameString;
options.dataType = 'Injected current';
options.dataUnits = 'nA';
fI = plotData(t, I, options);