function [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t, V, vBL, varargin)
% CALCRT calculates 10-90% or 20-80% rise time and associated indices given
% the time and data vectors and the baseline value.
%
%   T1090 = CALCRT(RTtype, t, data, baseline)
%   or T2080 = CALCRT(RTtype, t, data, baseline)
%   calculates 10-90% or 20-80% rise time, respectively. RTTYPE is a string
%   representing the rise time interval of choice: '10-90%' for a 10-90%
%   rise time and '20-80%' for a 20-80% rise time. T is the time vector in
%   milliseconds. DATA is the data vector that is either in millivolts or
%   in nanoamperes. BASELINE is the baseline value in millivolts or
%   nanoamperes, respectively.
%
%   T1090 = CALCRT(..., interactMode, pauseTime)
%   or T2080 = CALCRT(..., interactMode, pauseTime)
%   in addition displays the calculated measures in a figure. INTERACTMODE
%   should be set to true whereas PAUSETIME should be given in seconds.
%
%   [T1090, T10, I10, T50, I50, T90, I90] = CALCRT(...)
%   or [T2080, T20, I20, T50, I50, T80, I80] = CALCRT(...)
%   in addition calculated the associated times and indices of 10%, 50%,
%   and 90% or 20%, 50%, and 80% amplitude change, respectively.
%



if nargin == 6
    interactMode = varargin{1};
    pauseTime = varargin{2};
end


if strcmpi(RTtype, '10-90%')
    weights = [.9 .5 .1];
elseif strcmpi(RTtype, '20-80%')
    weights = [.8 .5 .2];
end


dataLength = length(V);
V = V - vBL;
if min(V) > 0
    V = V - min(V);
end

% Find the rise times backwards:
V90 = V - weights(1)*V(end);
V90original = V90;
V90(V90 > 0) = 0;
i90lastNegative = find(V90, 1, 'last');
[~, i90rev] = min(abs(fliplr(V90original(i90lastNegative:end))));
i90 = dataLength - i90rev + 1;
t90 = t(i90);

V50 = V - weights(2)*V(end);
V50original = V50;
V50(V50 > 0) = 0;
i50lastNegative = find(V50, 1, 'last');
[~, i50rev] = min(abs(fliplr(V50original(i50lastNegative:end))));
i50 = dataLength - i50rev + 1;
t50 = t(i50);

V10 = V - weights(3)*V(end);
V10original = V10;
V10(V10 > 0) = 0;
i10lastNegative = find(V10, 1, 'last');
[~, i10rev] = min(abs(fliplr(V10original(i10lastNegative:end))));
i10 = dataLength - i10rev + 1;
t10 = t(i10);
t1090 = t90 - t10;

% Test the function performance:
if exist('interactMode','var') && interactMode
    h = get(0,'CurrentFigure');
    if isempty(h)
        f1 = figure(1);
    else
        f1 = figure(h + 1);
    end
    
    plot(t, V, 'k');
    hold on
    plot(t90, V(i90), 'r.', 'markersize', 10);
    plot(t50, V(i50), 'r.', 'markersize', 10);
    plot(t10, V(i10), 'r.', 'markersize', 10);
    hold off
    pause(pauseTime*2);
    close(f1);
end