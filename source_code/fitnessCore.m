function [fitness, constraints, firstCost, allCosts, cliff] = fitnessCore( ...
  evalVector, meanMinisRTs, medianMinisRTs, cliff, minis2D, costFuncStruct, optimisationParameters)

costFuncType = 'parallel';

% Complying with the skewness constraint:
iCostBasis = ceil(0.03*numel(costFuncStruct.costBasis));
firstCost = 0;
allCosts = zeros(size(evalVector));
constraints = [0 0 0];
skewness = meanMinisRTs/medianMinisRTs;
if skewness <= 1
    fitness = costFuncStruct.costBasis(iCostBasis) + (1/skewness)*costFuncStruct.costBasis(end-1);
else
    fitness = 0;
end
constraints(1) = fitness;

% Complying with the cliff constraint:
iCostBasis = ceil(0.03*numel(costFuncStruct.costBasis));
if cliff
    cliff = detectCliff(minis2D);
    if any(cliff)
        if cliff(1)
            constraints(2) = costFuncStruct.costBasis(iCostBasis) + cliff(1)*costFuncStruct.costBasis(end-1);
            fitness = fitness + constraints(2);
        end
        if cliff(2)
            constraints(3) = costFuncStruct.costBasis(iCostBasis) + cliff(2)*costFuncStruct.costBasis(end-1);
            fitness = fitness + constraints(3);
        end
    end
else
    cliff = [0 0];
end

if strcmpi(costFuncType, 'serial') && fitness
    fitness = fitness + costFuncStruct.costBasis(1) + costFuncStruct.costScale(1)*evalVector(1);
    return
end

% Evaluating fitness:
costVector = evalVector - costFuncStruct.boundVector;
allCosts = costVector > 0;
firstCost = find(allCosts, 1);
if ~isempty(firstCost) && firstCost == length(evalVector)
    if ~isempty(optimisationParameters.SDlobound) && evalVector(end) < optimisationParameters.SDlobound
        fitness = fitness + costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*(optimisationParameters.SDlobound - evalVector(end));
    elseif ~isempty(optimisationParameters.SDlobound) && evalVector(end) > optimisationParameters.SDupbound
        fitness = fitness + costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*(evalVector(end) - optimisationParameters.SDupbound);
    else
        fitness = fitness + costFuncStruct.costBasis(firstCost) + 0;
    end
elseif ~isempty(firstCost)
    if strcmpi(costFuncType, 'serial')
        fitness = fitness + costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*evalVector(firstCost);
    elseif strcmpi(costFuncType, 'parallel')
        if allCosts(end)
            if ~isempty(optimisationParameters.SDlobound) && evalVector(end) < optimisationParameters.SDlobound
                fitness = fitness + costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*(optimisationParameters.SDlobound - evalVector(end));
            elseif ~isempty(optimisationParameters.SDlobound) && evalVector(end) > optimisationParameters.SDupbound
                fitness = fitness + costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*(evalVector(end) - optimisationParameters.SDupbound);
            else
                fitness = fitness + costFuncStruct.costBasis(firstCost) + 0;
            end
            allCosts(end) = 0;
            fitness = fitness + sum(costFuncStruct.costBasis(allCosts) + costFuncStruct.costScale(allCosts).*evalVector(allCosts));
            allCosts(end) = 1;
        else
            fitness = fitness + sum(costFuncStruct.costBasis(allCosts) + costFuncStruct.costScale(allCosts).*evalVector(allCosts));
        end
    end
else
    firstCost = find(evalVector, 1, 'last');
    fitness = fitness + costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*evalVector(firstCost);
end