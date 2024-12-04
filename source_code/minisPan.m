function minisPan(obj, ~)
% obj         handle to the figure that has been clicked on

global ratio initFile

ax = findobj(obj, 'type','axes');
if length(ax) == 3
    ax(1) = [];
end
axes1 = ax(1);
axes2 = ax(2);
axisvals1 = [get(axes1,'xlim'), get(axes1,'ylim')];

% Focus on bottom axes
set(axes2,'xlim', [ratio*(axisvals1(1)-initFile) ratio*(axisvals1(2)-initFile)]);
set(axes2,'ylim', axisvals1(3:4));