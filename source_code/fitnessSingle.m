function [SAD, AmpSAD, RTSAD, MAD, AmpMAD, RTMAD, AmpBottomSAD, AmpBottomSADlow, AmpBottomMAD, AmpMidSAD, AmpMidSADlow, AmpMidMAD, AmpTopSAD, AmpTopSADlow, AmpTopMAD,...
    tailMinis] = fitnessSingle(targetEvents1D, simulatedEvents1D, targetEvents1D_RT, simulatedEvents1D_RT, targetEvents2D, simulatedEvents2D, bottom, mid, top)

deviations1Dinit = targetEvents1D - simulatedEvents1D;
deviations1D = abs(deviations1Dinit);
deviations1Dlow = deviations1Dinit;
deviations1Dlow(deviations1Dlow < 0) = 0;
deviations1D_RT = abs(targetEvents1D_RT - simulatedEvents1D_RT);
difference2D = targetEvents2D - simulatedEvents2D;
deviations2D = abs(difference2D);

SAD = sum(sum(deviations2D));
AmpSAD = sum(deviations1D);
RTSAD = sum(deviations1D_RT);
MAD = max(max(deviations2D));
AmpMAD = max(deviations1D);
RTMAD = max(deviations1D_RT);

AmpBottomSAD = sum(deviations1D(bottom:end));
AmpBottomSADlow = sum(deviations1Dlow(bottom:end));
AmpBottomMAD = max(deviations1D(bottom:end));

AmpMidSAD = sum(deviations1D(mid:end));
AmpMidSADlow = sum(deviations1Dlow(mid:end));
AmpMidMAD = max(deviations1D(mid:end));
tailMinis = difference2D;
tailMinis(:,1:mid-1) = zeros(size(tailMinis,1),mid-1);
tailMinis(tailMinis < 0) = 0;
tailMinis = tailMinis*(sum(AmpMidSADlow)/sum(sum(tailMinis)));

AmpTopSAD = sum(deviations1D(top:end));
AmpTopSADlow = sum(deviations1Dlow(top:end));
AmpTopMAD = max(deviations1D(top:end));