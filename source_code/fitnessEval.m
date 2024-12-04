function [fitness, evalVector, bestEvalVector, meanMinisRTs, medianMinisRTs, firstCost, tailMinis, constraintFitness, allCosts, cliff] = fitnessEval( ...
    targetEvents1D, simulatedEvents1D, targetEvents1D_RT, simulatedEvents1D_RT, targetEvents2D, simulatedEvents2D, SD, costFuncStruct, measuredUF, ...
    histoData, shapes, cliff, minis2D, optimisationParameters)

[SAD, AmpsSAD, RTsSAD, MAD, AmpsMAD, RTsMAD, AmpsBottomSAD, AmpsBottomSADlow, AmpsBottomMAD, AmpsMidSAD, AmpsMidSADlow, AmpsMidMAD, AmpsTopSAD, AmpsTopSADlow,...
    AmpsTopMAD, tailMinis] = fitnessSingle(targetEvents1D, simulatedEvents1D, targetEvents1D_RT, simulatedEvents1D_RT, targetEvents2D, simulatedEvents2D,...
    costFuncStruct.bottom, costFuncStruct.mid, costFuncStruct.top);
if strcmpi(costFuncStruct.firstCost, 'TwoDs')
    evalVector = [SAD AmpsSAD RTsSAD MAD AmpsMAD RTsMAD AmpsBottomSAD AmpsBottomSADlow AmpsBottomMAD AmpsMidSAD AmpsMidSADlow AmpsMidMAD AmpsTopSAD...
        AmpsTopSADlow AmpsTopMAD];
elseif strcmpi(costFuncStruct.firstCost, 'Amps')
    evalVector = [AmpsSAD RTsSAD SAD AmpsMAD RTsMAD MAD AmpsBottomSAD AmpsBottomSADlow AmpsBottomMAD AmpsMidSAD AmpsMidSADlow AmpsMidMAD AmpsTopSAD...
        AmpsTopSADlow AmpsTopMAD];
elseif strcmpi(costFuncStruct.firstCost, 'RTs')
    evalVector = [RTsSAD AmpsSAD SAD RTsMAD AmpsMAD MAD AmpsBottomSAD AmpsBottomSADlow AmpsBottomMAD AmpsMidSAD AmpsMidSADlow AmpsMidMAD AmpsTopSAD...
        AmpsTopSADlow AmpsTopMAD];
end
if ~isempty(measuredUF)
    [SADMulti, AmpsSADMulti, RTsSADMulti, MADMulti, AmpsMADMulti, RTsMADMulti, AmpsBottomSADMulti, AmpsBottomSADlowMulti, AmpsBottomMADMulti, AmpsMidSADMulti,...
        AmpsMidSADlowMulti, AmpsMidMADMulti, AmpsTopSADMulti, AmpsTopSADlowMulti, AmpsTopMADMulti, SADBest, AmpsSADBest, RTsSADBest, MADBest, AmpsMADBest,...
        RTsMADBest, AmpsBottomSADBest, AmpsBottomSADlowBest, AmpsBottomMADBest, AmpsMidSADBest, AmpsMidSADlowBest, AmpsMidMADBest, AmpsTopSADBest,...
        AmpsTopSADlowBest, AmpsTopMADBest] = fitnessMulti(simulatedEvents1D, simulatedEvents1D_RT, simulatedEvents2D, histoData, costFuncStruct.bottom,...
        costFuncStruct.mid, costFuncStruct.top);
    if strcmpi(costFuncStruct.firstCost, 'TwoDs')
        multiEvalVector = [SADMulti AmpsSADMulti RTsSADMulti MADMulti AmpsMADMulti RTsMADMulti AmpsBottomSADMulti AmpsBottomSADlowMulti AmpsBottomMADMulti...
            AmpsMidSADMulti AmpsMidSADlowMulti AmpsMidMADMulti AmpsTopSADMulti AmpsTopSADlowMulti AmpsTopMADMulti];
        bestEvalVector = [SADBest AmpsSADBest RTsSADBest MADBest AmpsMADBest RTsMADBest AmpsBottomSADBest AmpsBottomSADlowBest AmpsBottomMADBest...
            AmpsMidSADBest AmpsMidSADlowBest AmpsMidMADBest AmpsTopSADBest AmpsTopSADlowBest AmpsTopMADBest];
    elseif strcmpi(costFuncStruct.firstCost, 'Amps')
        multiEvalVector = [AmpsSADMulti RTsSADMulti SADMulti AmpsMADMulti RTsMADMulti MADMulti AmpsBottomSADMulti AmpsBottomSADlowMulti AmpsBottomMADMulti...
            AmpsMidSADMulti AmpsMidSADlowMulti AmpsMidMADMulti AmpsTopSADMulti AmpsTopSADlowMulti AmpsTopMADMulti];
        bestEvalVector = [AmpsSADBest RTsSADBest SADBest AmpsMADBest RTsMADBest MADBest AmpsBottomSADBest AmpsBottomSADlowBest AmpsBottomMADBest...
            AmpsMidSADBest AmpsMidSADlowBest AmpsMidMADBest AmpsTopSADBest AmpsTopSADlowBest AmpsTopMADBest];
    elseif strcmpi(costFuncStruct.firstCost, 'RTs')
        multiEvalVector = [RTsSADMulti AmpsSADMulti SADMulti RTsMADMulti AmpsMADMulti MADMulti AmpsBottomSADMulti AmpsBottomSADlowMulti AmpsBottomMADMulti...
            AmpsMidSADMulti AmpsMidSADlowMulti AmpsMidMADMulti AmpsTopSADMulti AmpsTopSADlowMulti AmpsTopMADMulti];
        bestEvalVector = [RTsSADBest AmpsSADBest SADBest RTsMADBest AmpsMADBest MADBest AmpsBottomSADBest AmpsBottomSADlowBest AmpsBottomMADBest...
            AmpsMidSADBest AmpsMidSADlowBest AmpsMidMADBest AmpsTopSADBest AmpsTopSADlowBest AmpsTopMADBest];
    end
else
    multiEvalVector = zeros(1,15);
    bestEvalVector = zeros(1,15);
end
evalVector = [multiEvalVector evalVector SD]; % multiEvalVector is simulation comparison to all target files
                                              % evalVector is simulation comprison to the single chosen target file

meanMinisRTs = mean(shapes(:,3));
medianMinisRTs = median(shapes(:,3));
stdMinisRTs = std(shapes(:,3));

[fitness, constraintFitness, firstCost, allCosts, cliff] = fitnessCore(evalVector, meanMinisRTs, medianMinisRTs, stdMinisRTs, cliff, minis2D, shapes, costFuncStruct, optimisationParameters);
criticalCosts = [2 5 13:15 17 20 28:30];
%criticalCosts = [2 5 13:15 17 20 30];
if any(ismember(find(allCosts), criticalCosts))
    criticalCosts = criticalCosts(ismember(criticalCosts, find(allCosts)));
else
    criticalCosts = [];
end
disp(['First cost: ' num2str(firstCost) '  Fitness: ' num2str(fitness) '  Cliffs: ' num2str(cliff) ...
  '  Constraints: ' num2str(constraintFitness) '  Costs: ' num2str(find(allCosts)) '  Critical costs: ' num2str(criticalCosts)]);