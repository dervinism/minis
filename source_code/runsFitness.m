function [fitAmp, fitRT, pAmp, pRT] = runsFitness(difference1D, lastBin, difference1D_RT, lastBinRT, fit)

diff1Drunt = difference1D(1:lastBin);
diff1Drunt(diff1Drunt == 0) = [];
diff1D_RTrunt = difference1D_RT(1:lastBinRT);
diff1D_RTrunt(diff1D_RTrunt == 0) = [];
[~, pAmp] = runstest(diff1Drunt,0,'tail','right');
[~, pRT] = runstest(diff1D_RTrunt,0,'tail','right');
if fit
    % Runs test fitness:
    if pAmp <= .1
        fitAmp = 1;
    else
        fitAmp = 1+pAmp;
    end
    if pRT <= .1
        fitRT = 1;
    else
        fitRT = 1+pRT;
    end
else
    fitAmp = [];
    fitRT = [];
end
end