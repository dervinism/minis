function [fitness, firstCost] = fitnessCore(evalVector, meanMinisRTs, medianMinisRTs, cliff, minis2D, costFuncStruct, optimisationParameters)

% Complying with the skewness constraint:
skewness = meanMinisRTs/medianMinisRTs;
if skewness <= 1
    fitness = costFuncStruct.costBasis(1) + 8*costFuncStruct.costBasis(end-1);
    firstCost = 0;
    return
end

% Complying with the cliff constraint:
if cliff
    cliff = detectCliff(minis2D);
    if any(cliff)
      if cliff(1)
        fitness = costFuncStruct.costBasis(1) + 2*costFuncStruct.costBasis(end-1) + cliff(1)*costFuncStruct.costBasis(end-1);
      end
      if cliff(2)
        fitness = costFuncStruct.costBasis(1) + 2*costFuncStruct.costBasis(end-1) + cliff(2)*costFuncStruct.costBasis(end-1);
      end
      firstCost = 0;
      return
    end
end

% Evaluating fitness:
costVector = evalVector - costFuncStruct.boundVector;
firstCost = find(costVector > 0, 1);
if ~isempty(firstCost) && firstCost == length(evalVector)
    if ~isempty(optimisationParameters.SDlobound) && evalVector(end) < optimisationParameters.SDlobound
        fitness = costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*(optimisationParameters.SDlobound - evalVector(end));
    elseif ~isempty(optimisationParameters.SDlobound) && evalVector(end) > optimisationParameters.SDupbound
        fitness = costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*(evalVector(end) - optimisationParameters.SDupbound);
    else
        fitness = costFuncStruct.costBasis(firstCost) + 0;
    end
elseif ~isempty(firstCost)
    fitness = costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*evalVector(firstCost);
else
    firstCost = find(evalVector, 1, 'last');
    fitness = costFuncStruct.costBasis(firstCost) + costFuncStruct.costScale(firstCost)*evalVector(firstCost);
end