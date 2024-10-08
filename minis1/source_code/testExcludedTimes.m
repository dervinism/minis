function [startTimesOut, endTimesOut] = testExcludedTimes(startTimesStr, endTimesStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd)
% TESTEXCLUDEDTIMES tests whether the start and end excluded time vectors
% (1) contain only numerical values, (2) are of equal lengths, (3) contain
% end time values that are larger than start time values, and (4) converts
% the time vectors into a numerical format.
%
%   [STARTTIMESOUT ENDTIMESOUT] = TESTEXCLUDEDTIMES(startTimesStr,...
%       endTimesStr, errmsgNaN, errmsgVectorLengths, errmsgStartEnd)
%   applies the aforementioned tests and returns converted vector variables.
%   errmsgNaN is a string variable containing an error message displayed
%   when excluded time values are not numerical. errmsgVectorLengths is an
%   error message displayed when time vectors are not of equal lengths.
%   errmsgStartEnd is an error message displayed when at least one of the
%   end time values is smaller than its corresponding start value.
%

% Numerical values?
startTimesStr = strrep(startTimesStr, ' ', '');
startTimesStr = regexp(startTimesStr, ',*', 'split');
startTimes = strarray2numarray(startTimesStr');                             % Convert to numbers
if ~isempty(startTimesStr) && sum(~strcmpi(startTimesStr, '...')) && sum(isnan(startTimes))
    msgbox(errmsgNaN,'Error','Error');
    return
end

% Numerical values?
endTimesStr = strrep(endTimesStr, ' ', '');
endTimesStr = regexp(endTimesStr, ',*', 'split');
endTimes = strarray2numarray(endTimesStr');                                 % Convert to numbers
if ~isempty(endTimesStr) && sum(~strcmpi(endTimesStr, '...')) && sum(isnan(endTimes))
    msgbox(errmsgNaN,'Error','Error');
    return
end

% Equal lengths?
if length(startTimes) ~= length(endTimes)
    msgbox(errmsgVectorLengths,'Error','Error');
    return
end

% Start value larger than end value?
diff = endTimes - startTimes;
diff(diff > 0) = 0;
negativeValues = sum(abs(diff));
if negativeValues
    msgbox(errmsgStartEnd,'Error','Error');
    return
end

startTimesOut = startTimes;
endTimesOut = endTimes;