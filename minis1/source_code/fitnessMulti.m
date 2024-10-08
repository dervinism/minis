function [TwoDsSADMulti, AmpsSADMulti, RTsSADMulti, TwoDsMADMulti, AmpsMADMulti, RTsMADMulti, AmpsBottomSADMulti, AmpsBottomSADlowMulti, AmpsBottomMADMulti,...
    AmpsMidSADMulti, AmpsMidSADlowMulti, AmpsMidMADMulti, AmpsTopSADMulti, AmpsTopSADlowMulti, AmpsTopMADMulti, TwoDsSAD, AmpsSAD, RTsSAD, TwoDsMAD, AmpsMAD, RTsMAD,...
    AmpsBottomSAD, AmpsBottomSADlow, AmpsBottomMAD, AmpsMidSAD, AmpsMidSADlow, AmpsMidMAD, AmpsTopSAD, AmpsTopSADlow, AmpsTopMAD] = fitnessMulti(...
    simulatedEvents1D, simulatedEvents1D_RT, simulatedEvents2D, histoData, bottom, mid, top)



%% Calculate the multiple-error scores:
dataPoints = size(histoData.Amps,1);
AmpsSADinit = histoData.Amps - repmat(simulatedEvents1D, dataPoints, 1);
AmpsSAD = abs(AmpsSADinit);
RTsSAD = abs(histoData.RTs - repmat(simulatedEvents1D_RT, dataPoints, 1));
TwoDsSAD = abs(histoData.TwoDs - repmat(simulatedEvents2D, [1 1 dataPoints]));

AmpsBottomMAD = min(max(AmpsSAD(:,bottom:end),[],2));
AmpsBottomSAD = min(sum(AmpsSAD(:,bottom:end),2));
AmpsSADinit(AmpsSADinit < 0) = 0;
AmpsBottomSADlow = min(sum(AmpsSADinit(:,bottom:end),2));

AmpsMidMAD = min(max(AmpsSAD(:,mid:end),[],2));
AmpsMidSAD = min(sum(AmpsSAD(:,mid:end),2));
AmpsMidSADlow = min(sum(AmpsSADinit(:,mid:end),2));

AmpsTopMAD = min(max(AmpsSAD(:,top:end),[],2));
AmpsTopSAD = min(sum(AmpsSAD(:,top:end),2));
AmpsTopSADlow = min(sum(AmpsSADinit(:,top:end),2));

AmpsMAD = min(max(AmpsSAD,[],2));
RTsMAD = min(max(RTsSAD,[],2));
TwoDsMAD = min(max(max(TwoDsSAD,[],1),[],2));
AmpsSAD = min(sum(AmpsSAD,2));
RTsSAD = min(sum(RTsSAD,2));
TwoDsSAD = min(sum(sum(TwoDsSAD,1),2));



%% Calculate the 'sticking-out' area:
AmpsData = [histoData.Amps; simulatedEvents1D];
RTsData = [histoData.RTs; simulatedEvents1D_RT];
TwoDsData = cat(3, histoData.TwoDs, simulatedEvents2D);
AmpsSADMulti = zeros(dataPoints+1,1);
RTsSADMulti = AmpsSADMulti;
TwoDsSADMulti = AmpsSADMulti;
AmpsMADMulti = AmpsSADMulti;
RTsMADMulti = AmpsSADMulti;
TwoDsMADMulti = AmpsSADMulti;
AmpsBottomSADMulti = AmpsSADMulti;
AmpsBottomSADlowMulti = AmpsSADMulti;
AmpsBottomMADMulti = AmpsSADMulti;
AmpsMidSADMulti = AmpsSADMulti;
AmpsMidSADlowMulti = AmpsSADMulti;
AmpsMidMADMulti = AmpsSADMulti;
AmpsTopSADMulti = AmpsSADMulti;
AmpsTopSADlowMulti = AmpsSADMulti;
AmpsTopMADMulti = AmpsSADMulti;
for iData = 1:dataPoints+1
    if iData == 1
        AmpsSADMulti1 = AmpsData(iData,:) - max(AmpsData(iData+1:end,:));
        AmpsSADMulti1(AmpsSADMulti1 < 0) = 0;
        AmpsSADMulti2 = min(AmpsData(iData+1:end,:)) - AmpsData(iData,:);
        AmpsSADMulti2(AmpsSADMulti2 < 0) = 0;
        AmpsSADMulti(iData) = sum(AmpsSADMulti1) + sum(AmpsSADMulti2);
        AmpsMADMulti(iData) = max([max(AmpsSADMulti1) max(AmpsSADMulti2)]);
        AmpsBottomSADMulti(iData) = sum(AmpsSADMulti1(bottom:end)) + sum(AmpsSADMulti2(bottom:end));
        AmpsBottomSADlowMulti(iData) = sum(AmpsSADMulti2(bottom:end));
        AmpsBottomMADMulti(iData) = max([max(AmpsSADMulti1(bottom:end)) max(AmpsSADMulti2(bottom:end))]);
        AmpsMidSADMulti(iData) = sum(AmpsSADMulti1(mid:end)) + sum(AmpsSADMulti2(mid:end));
        AmpsMidSADlowMulti(iData) = sum(AmpsSADMulti2(mid:end));
        AmpsMidMADMulti(iData) = max([max(AmpsSADMulti1(mid:end)) max(AmpsSADMulti2(mid:end))]);
        AmpsTopSADMulti(iData) = sum(AmpsSADMulti1(top:end)) + sum(AmpsSADMulti2(top:end));
        AmpsTopSADlowMulti(iData) = sum(AmpsSADMulti2(top:end));
        AmpsTopMADMulti(iData) = max([max(AmpsSADMulti1(top:end)) max(AmpsSADMulti2(top:end))]);
        
        RTsSADMulti1 = RTsData(iData,:) - max(RTsData(iData+1:end,:));
        RTsSADMulti1(RTsSADMulti1 < 0) = 0;
        RTsSADMulti2 = min(RTsData(iData+1:end,:)) - RTsData(iData,:);
        RTsSADMulti2(RTsSADMulti2 < 0) = 0;
        RTsSADMulti(iData) = sum(RTsSADMulti1) + sum(RTsSADMulti2);
        RTsMADMulti(iData) = max([max(RTsSADMulti1) max(RTsSADMulti2)]);
        
        TwoDsSADMulti1 = TwoDsData(:,:,iData) - max(TwoDsData(:,:,iData+1:end),[],3);
        TwoDsSADMulti1(TwoDsSADMulti1 < 0) = 0;
        TwoDsSADMulti2 = min(TwoDsData(:,:,iData+1:end),[],3) - TwoDsData(:,:,iData);
        TwoDsSADMulti2(TwoDsSADMulti2 < 0) = 0;
        TwoDsSADMulti(iData) = sum(sum(TwoDsSADMulti1)) + sum(sum(TwoDsSADMulti2));
        TwoDsMADMulti(iData) = max([max(max(TwoDsSADMulti1)) max(max(TwoDsSADMulti2))]);
    elseif iData == dataPoints
        AmpsSADMulti1 = AmpsData(iData,:) - max(AmpsData(1:end-1,:));
        AmpsSADMulti1(AmpsSADMulti1 < 0) = 0;
        AmpsSADMulti2 = min(AmpsData(1:end-1,:)) - AmpsData(iData,:);
        AmpsSADMulti2(AmpsSADMulti2 < 0) = 0;
        AmpsSADMulti(iData) = sum(AmpsSADMulti1) + sum(AmpsSADMulti2);
        AmpsMADMulti(iData) = max([max(AmpsSADMulti1) max(AmpsSADMulti2)]);
        AmpsBottomSADMulti(iData) = sum(AmpsSADMulti1(bottom:end)) + sum(AmpsSADMulti2(bottom:end));
        AmpsBottomSADlowMulti(iData) = sum(AmpsSADMulti2(bottom:end));
        AmpsBottomMADMulti(iData) = max([max(AmpsSADMulti1(bottom:end)) max(AmpsSADMulti2(bottom:end))]);
        AmpsMidSADMulti(iData) = sum(AmpsSADMulti1(mid:end)) + sum(AmpsSADMulti2(mid:end));
        AmpsMidSADlowMulti(iData) = sum(AmpsSADMulti2(mid:end));
        AmpsMidMADMulti(iData) = max([max(AmpsSADMulti1(mid:end)) max(AmpsSADMulti2(mid:end))]);
        AmpsTopSADMulti(iData) = sum(AmpsSADMulti1(top:end)) + sum(AmpsSADMulti2(top:end));
        AmpsTopSADlowMulti(iData) = sum(AmpsSADMulti2(top:end));
        AmpsTopMADMulti(iData) = max([max(AmpsSADMulti1(top:end)) max(AmpsSADMulti2(top:end))]);
        
        RTsSADMulti1 = RTsData(iData,:) - max(RTsData(1:end-1,:));
        RTsSADMulti1(RTsSADMulti1 < 0) = 0;
        RTsSADMulti2 = min(RTsData(1:end-1:end,:)) - RTsData(iData,:);
        RTsSADMulti2(RTsSADMulti2 < 0) = 0;
        RTsSADMulti(iData) = sum(RTsSADMulti1) + sum(RTsSADMulti2);
        RTsMADMulti(iData) = max([max(RTsSADMulti1) max(RTsSADMulti2)]);
        
        TwoDsSADMulti1 = TwoDsData(:,:,iData) - max(TwoDsData(:,:,1:end-1),[],3);
        TwoDsSADMulti1(TwoDsSADMulti1 < 0) = 0;
        TwoDsSADMulti2 = min(TwoDsData(:,:,1:end-1),[],3) - TwoDsData(:,:,iData);
        TwoDsSADMulti2(TwoDsSADMulti2 < 0) = 0;
        TwoDsSADMulti(iData) = sum(sum(TwoDsSADMulti1)) + sum(sum(TwoDsSADMulti2));
        TwoDsMADMulti(iData) = max([max(max(TwoDsSADMulti1)) max(max(TwoDsSADMulti2))]);
    else
        AmpsSADMulti1 = AmpsData(iData,:) - max(AmpsData([1:iData-1 iData+1:end],:));
        AmpsSADMulti1(AmpsSADMulti1 < 0) = 0;
        AmpsSADMulti2 = min(AmpsData([1:iData-1 iData+1:end],:)) - AmpsData(iData,:);
        AmpsSADMulti2(AmpsSADMulti2 < 0) = 0;
        AmpsSADMulti(iData) = sum(AmpsSADMulti1) + sum(AmpsSADMulti2);
        AmpsMADMulti(iData) = max([max(AmpsSADMulti1) max(AmpsSADMulti2)]);
        AmpsBottomSADMulti(iData) = sum(AmpsSADMulti1(bottom:end)) + sum(AmpsSADMulti2(bottom:end));
        AmpsBottomSADlowMulti(iData) = sum(AmpsSADMulti2(bottom:end));
        AmpsBottomMADMulti(iData) = max([max(AmpsSADMulti1(bottom:end)) max(AmpsSADMulti2(bottom:end))]);
        AmpsMidSADMulti(iData) = sum(AmpsSADMulti1(mid:end)) + sum(AmpsSADMulti2(mid:end));
        AmpsMidSADlowMulti(iData) = sum(AmpsSADMulti2(mid:end));
        AmpsMidMADMulti(iData) = max([max(AmpsSADMulti1(mid:end)) max(AmpsSADMulti2(mid:end))]);
        AmpsTopSADMulti(iData) = sum(AmpsSADMulti1(top:end)) + sum(AmpsSADMulti2(top:end));
        AmpsTopSADlowMulti(iData) = sum(AmpsSADMulti2(top:end));
        AmpsTopMADMulti(iData) = max([max(AmpsSADMulti1(top:end)) max(AmpsSADMulti2(top:end))]);
        
        RTsSADMulti1 = RTsData(iData,:) - max(RTsData([1:iData-1 iData+1:end],:));
        RTsSADMulti1(RTsSADMulti1 < 0) = 0;
        RTsSADMulti2 = min(RTsData([1:iData-1 iData+1:end],:)) - RTsData(iData,:);
        RTsSADMulti2(RTsSADMulti2 < 0) = 0;
        RTsSADMulti(iData) = sum(RTsSADMulti1) + sum(RTsSADMulti2);
        RTsMADMulti(iData) = max([max(RTsSADMulti1) max(RTsSADMulti2)]);
        
        TwoDsSADMulti1 = TwoDsData(:,:,iData) - max(TwoDsData(:,:,[1:iData-1 iData+1:end]),[],3);
        TwoDsSADMulti1(TwoDsSADMulti1 < 0) = 0;
        TwoDsSADMulti2 = min(TwoDsData(:,:,[1:iData-1 iData+1:end]),[],3) - TwoDsData(:,:,iData);
        TwoDsSADMulti2(TwoDsSADMulti2 < 0) = 0;
        TwoDsSADMulti(iData) = sum(sum(TwoDsSADMulti1)) + sum(sum(TwoDsSADMulti2));
        TwoDsMADMulti(iData) = max([max(max(TwoDsSADMulti1)) max(max(TwoDsSADMulti2))]);
    end
end
AmpsSADMulti = AmpsSADMulti(end) - max(AmpsSADMulti(1:end-1));
RTsSADMulti = RTsSADMulti(end) - max(RTsSADMulti(1:end-1));
TwoDsSADMulti = TwoDsSADMulti(end) - max(TwoDsSADMulti(1:end-1));
AmpsMADMulti = AmpsMADMulti(end) - max(AmpsMADMulti(1:end-1));
RTsMADMulti = RTsMADMulti(end) - max(RTsMADMulti(1:end-1));
TwoDsMADMulti = TwoDsMADMulti(end) - max(TwoDsMADMulti(1:end-1));
AmpsBottomSADMulti = AmpsBottomSADMulti(end) - max(AmpsBottomSADMulti(1:end-1));
AmpsBottomSADlowMulti = AmpsBottomSADlowMulti(end) - max(AmpsBottomSADlowMulti(1:end-1));
AmpsBottomMADMulti = AmpsBottomMADMulti(end) - max(AmpsBottomMADMulti(1:end-1));
AmpsMidSADMulti = AmpsMidSADMulti(end) - max(AmpsMidSADMulti(1:end-1));
AmpsMidSADlowMulti = AmpsMidSADlowMulti(end) - max(AmpsMidSADlowMulti(1:end-1));
AmpsMidMADMulti = AmpsMidMADMulti(end) - max(AmpsMidMADMulti(1:end-1));
AmpsTopSADMulti = AmpsTopSADMulti(end) - max(AmpsTopSADMulti(1:end-1));
AmpsTopSADlowMulti = AmpsTopSADlowMulti(end) - max(AmpsTopSADlowMulti(1:end-1));
AmpsTopMADMulti = AmpsTopMADMulti(end) - max(AmpsTopMADMulti(1:end-1));