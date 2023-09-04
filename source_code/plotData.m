function f1 = plotData(t, V, options)
% Draws membrane potential data  and returns a handles to the figure. This
% subfunction is designated for drawing raw and smooth data in detectMinis
% and estimateNoiseMinis functions.
%
%   f1 = plotData(t, V, options) draws a data graph and returns a handle to
%   it. t is a time vector, V is a vector containing data, and options is a
%   structure variable following fields:
%   'nameString' - the name displayed on the figure window header;
%   'titleString' - the figure title.
%   'dataType' - the data type (i.e., 'Membrane potential' or 'Current
%       clamp');
%   'dataUnits' - data measurement units (i.e., 'mV' or 'nA');
%
%   If V is an array of two vectors, instead a figure of two types of
%   data superimposed on each other is drawn and a handle is returned.
%

screenSize = get(0, 'ScreenSize');
screenSize = [screenSize(3) screenSize(4)];

h = get(0,'CurrentFigure');
if isempty(h)
    f1 = figure(1);
else
    figurecount = length(findobj('Type','figure'));
    f1 = figure(figurecount + 1);
end

set(f1, 'position', [17*screenSize(1)/2560 481*screenSize(2)/1024 2560*screenSize(1)/2560 440*screenSize(2)/1024]);
set(f1, 'NumberTitle', 'off');
set(f1, 'Name', options.nameString);
set(gca,'Position' , [0.025 0.1 0.965 0.815]);
set(gca,'TickLength', [0.002 0.01]);
plot(t, V(1,:), 'k', 'markersize', 4);
hold on;
if size(V,1) == 2
    plot(t, V(2,:), 'r', 'markersize', 4);
end
xlabel('Time (ms)');
ylabel(sprintf('%s  (%s)', options.dataType, options.dataUnits));
title(options.titleString, 'Interpreter','none');
hold off;
end