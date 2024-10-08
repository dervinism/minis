function state = customPlotGA(~,state,flag)

persistent f1 f2 f3

s = gca;
delete(s);
% s = get(gcf,'children');
% positions = get(s,'Position');
% for pos = 1:length(positions)
%     if isequal(round(10*positions{pos}), [1 1 3 2]);
%         delete(s(pos));
%     elseif isequal(round(10*positions{pos}), [1 7 3 2]);
%         set(s(pos), 'Units','normalized','Position',[.1 .64 .35 .27]);
%     elseif isequal(round(10*positions{pos}), [6 7 3 2]);
%         set(s(pos), 'Units','normalized','Position',[.585 .64 .35 .27]);
%     elseif isequal(round(10*positions{pos}), [1 4 3 2]);
%         set(s(pos), 'Units','normalized','Position',[.1 .18 .35 .27]);
%     elseif isequal(round(10*positions{pos}), [6 4 3 2]);
%         set(s(pos), 'Units','normalized','Position',[.585 .18 .35 .27]);
%     end
% end
f0 = gcf;
dataDir = 'data';
cd(dataDir);
filename = findFile('name', 1);
if strcmpi(filename(end-2:end),'mat')
    if ~strcmpi(filename(end-7:end-4),'full')
        load(filename); %#ok<*LOAD>
        saveas(f0, strcat(num2str(fitness),'_full_evolution.jpg'));
        load('initVar');
        cd ..
        simulationParameters.tau_sy1 = tau_sy1;
        simulationParameters.tau_sy2 = 0.1*tau_sy1;
        if ~strcmpi(flag,'init')
            [f1, f2, f3] = plotSaveFP(dataDir, distributionParameters, V, shapes, classificationParameters, noiseProperties, target, targetDataLength,...
                current, detectionParameters, simulationParameters, costFunction, distributionType, dimension, tauRange, tauSyEst, baseline, filtering,...
                f1, f2, f3);
        else
            [f1, f2, f3] = plotSaveFP(dataDir, distributionParameters, V, shapes, classificationParameters, noiseProperties, target, targetDataLength,...
                current, detectionParameters, simulationParameters, costFunction, distributionType, dimension, tauRange, tauSyEst, baseline, filtering);
        end
    else
        cd ..
    end
else
    cd ..
end