function [fitness, constraints, firstCost, allCosts, cliff] = fitnessCore( ...
  evalVector, meanMinisRTs, medianMinisRTs, stdMinisRTs, cliff, minis2D, ...
  shapes, costFuncStruct, optimisationParameters)

costFuncType = 'parallel'; % 'serial' or 'parallel'

% Complying with the skewness constraint:
iCostBasis = ceil(0.03*numel(costFuncStruct.costBasis));
firstCost = 0; %#ok<*NASGU>
allCosts = zeros(size(evalVector));
constraints = [0 0 0 0 0 0];
skewness = (meanMinisRTs - medianMinisRTs)/stdMinisRTs;
if skewness <= 0
    fitness = costFuncStruct.costBasis(iCostBasis) + (0.5/1e-9)*costFuncStruct.costBasis(end-1);
elseif skewness < 0.25
    fitness = costFuncStruct.costBasis(iCostBasis) + (0.5/skewness)*costFuncStruct.costBasis(end-1);
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

% Complying with the unimodality constraint
Thr = 0.45;
amps = shapes(:,2);
[~, bcAmp] = bimodalitycoeff(amps);
if bcAmp > Thr
  constraints(4) = costFuncStruct.costBasis(iCostBasis) + (bcAmp-Thr)*costFuncStruct.costBasis(end-1);
  fitness = fitness + constraints(4);
end

RTs = shapes(:,3);
[~, bcRTs] = bimodalitycoeff(RTs);
if bcRTs > Thr
  constraints(5) = costFuncStruct.costBasis(iCostBasis) + (bcRTs-Thr)*costFuncStruct.costBasis(end-1);
  fitness = fitness + constraints(5);
end

% Complying with rise time mean contraint
Thr = 5;
if mean(RTs) > 5
  constraints(6) = costFuncStruct.costBasis(iCostBasis) + (mean(RTs)-Thr)*costFuncStruct.costBasis(end-1);
  fitness = fitness + constraints(6);
end

%figure; plot(0:0.01:10, sum(minis2D));
%figure; plot(0:0.5:11, sum(minis2D,2));

if strcmpi(costFuncType, 'serial') && fitness
    fitness = fitness + costFuncStruct.costBasis(1) + costFuncStruct.costScale(1)*evalVector(1); %#ok<*UNRCH>
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