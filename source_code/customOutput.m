function [state, options, optchanged] = customOutput(options, state, ~, path)

persistent bestFitness lastImprovement zeroDiv

gNumber = state.Generation ;                                                % Get the current generation number
fitness = min(state.Score);                                                 % Best score in the current generation
population = state.Population;                                              % The current population

if gNumber == 0 || bestFitness > fitness
    bestFitness = fitness;
    lastImprovement = gNumber;
end

fprintf('\r\nCurrent generation number: %g\r\n', gNumber);
fprintf('\r\nLast improvement in generation: %g\r\n', lastImprovement);
fprintf('\r\nBest score so far: %g\r\n', bestFitness);
fprintf('\r\n');

currDir = pwd;
cd(path)
if ~exist('genNumber.txt','file')
    fileID = fopen('genNumber.txt','w');
else
    fileID = fopen('genNumber.txt','a');
end
fprintf(fileID, '\r\nCurrent generation number: %g\r\n', gNumber);
fprintf(fileID, '\r\nLast improvement in generation: %g\r\n', lastImprovement);
fprintf(fileID, '\r\nBest score so far: %g\r\n', bestFitness);
fprintf(fileID, '\r\n');
cd(currDir);
optchanged = false;


%% Stop ga if the population diversity is zero:
if sum(sum(population)) == sum(size(population,1)*population(1,:))
    if isempty(zeroDiv)
        zeroDiv = 1;
        fprintf('\r\nZero population diversity');
        fprintf(fileID, '\r\nZero population diversity');
    else
        state.StopFlag = 1;
        fprintf('\r\nOptimisation terminated: Zero population diversity');
        fprintf(fileID, '\r\nOptimisation terminated: Zero population diversity');
    end
end


%% Stop ga if there is no improvement in the fitness value for 100 generations in a row:
if gNumber - lastImprovement >= 100
    state.StopFlag = 1;
    fprintf('\r\nOptimisation terminated: No improvement for 100 generations');
    fprintf(fileID, '\r\nOptimisation terminated: No improvement for 100 generations');
end