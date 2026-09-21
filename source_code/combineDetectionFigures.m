function combinedFig = combineDetectionFigures(sourceFigs, panelTitle)
% COMBINEDETECTIONFIGURES merges the axes of several existing figures into
% a single figure of tiled panels, then closes the source figures.
%
%   COMBINEDFIG = COMBINEDETECTIONFIGURES(sourceFigs, panelTitle) accepts
%   SOURCEFIGS, an array of figure handles (invalid/already-closed handles
%   are silently skipped). Each source figure may itself contain more than
%   one axes (e.g. a 2-subplot histogram figure) -- every axes found
%   becomes its own panel, in the order the source figures were given.
%   Legends and colorbars attached to a copied axes are preserved.
%   PANELTITLE is used as both the combined figure's window title and its
%   overall tiledlayout title.
%
%   Returns an empty graphics array if no valid source figures/axes were
%   found (nothing is created in that case).

validFigs = sourceFigs(ishghandle(sourceFigs));

allAxes = gobjects(0);
for iFig = 1:numel(validFigs)
    % findobj enumerates a figure's axes newest-first; flip back to
    % creation order (matching left-to-right/top-to-bottom subplot order):
    axesInFig = findobj(validFigs(iFig), 'Type', 'axes');
    axesInFig = flipud(axesInFig(:));
    allAxes = [allAxes; axesInFig]; %#ok<AGROW>
end

if isempty(allAxes)
    combinedFig = gobjects(0);
    close(validFigs);
    return
end

nPanels = numel(allAxes);

combinedFig = figure('Units', 'Normalized', 'Position', [0.02 0.05 0.96 0.85],...
    'Name', panelTitle, 'NumberTitle', 'off');
% 'flow' (rather than a fixed nRows x nCols grid) arranges exactly as many
% tiles as there are axes, so a panel count that doesn't factor evenly
% (e.g. 5) never leaves an empty trailing tile:
tl = tiledlayout(combinedFig, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, panelTitle, 'Interpreter', 'none');

for iAx = 1:nPanels
    ax = allAxes(iAx);
    % A legend/colorbar is not a child of its axes in the graphics
    % hierarchy, so it must be copied together with the axes in the same
    % copyobj call for the association to survive:
    decorations = gobjects(0);
    leg = get(ax, 'Legend');
    if ~isempty(leg) && isgraphics(leg)
        decorations(end+1) = leg; %#ok<AGROW>
    end
    cb = get(ax, 'Colorbar');
    if ~isempty(cb) && isgraphics(cb)
        decorations(end+1) = cb; %#ok<AGROW>
    end
    copied = copyobj([ax decorations], tl);
    % copyobj does NOT auto-assign a tile in 'flow' mode the way nexttile
    % does -- every copied axes silently defaults to Layout.Tile = 1 (and
    % therefore identical position) unless set explicitly, which stacks
    % every panel on top of each other:
    copied(1).Layout.Tile = iAx;
    % Source titles/labels often contain the recording's file name (e.g.
    % the frequency-spectrum panel), which the default 'tex' interpreter
    % would otherwise render with underscores as subscripts:
    copied(1).Title.Interpreter = 'none';
end

close(validFigs);
end
