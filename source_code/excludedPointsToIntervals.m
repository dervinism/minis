function intervalsMs = excludedPointsToIntervals(excludedPointsMs, dt)
% EXCLUDEDPOINTSTOINTERVALS collapses a vector of individual excluded
% sample time points (ms, one entry per excluded sample -- the internal
% representation detectMinis uses, e.g. its excludedT variable) into an
% Nx2 matrix of contiguous [start end] intervals (ms), for use with
% INCIDENCERATE. DT is the sampling interval (ms).

points = sort(excludedPointsMs(:));
if isempty(points)
    intervalsMs = zeros(0,2);
    return
end

% A gap of more than ~1.5 sample intervals marks a break between
% separately-excluded regions:
breaks = find(diff(points) > 1.5*dt);
starts = [points(1); points(breaks+1)];
ends = [points(breaks); points(end)] + dt;
intervalsMs = [starts ends];
end
