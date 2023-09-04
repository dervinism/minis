function [AmpsSADMulti, RTsSADMulti, TwoDsSADMulti, AmpsMADMulti, RTsMADMulti, TwoDsMADMulti, AmpsMidSADMulti, AmpsMidMADMulti, AmpsTopSADMulti,...
    AmpsTopMADMulti] = multiErrorBounds(histoData, Mid, Top)



%% Calculate the 'sticking-out' area:
AmpsData = histoData.Amps;
RTsData = histoData.RTs;
TwoDsData = histoData.TwoDs;
dataPoints = size(AmpsData,1);
AmpsSADMulti = zeros(dataPoints+1,1);
RTsSADMulti = AmpsSADMulti;
TwoDsSADMulti = AmpsSADMulti;
AmpsMADMulti = AmpsSADMulti;
RTsMADMulti = AmpsSADMulti;
TwoDsMADMulti = AmpsSADMulti;
AmpsMidSADMulti = AmpsSADMulti;
AmpsMidMADMulti = AmpsSADMulti;
AmpsTopSADMulti = AmpsSADMulti;
AmpsTopMADMulti = AmpsSADMulti;
for iData = 1:dataPoints
    if iData == 1
        AmpsSADMulti1 = AmpsData(iData,:) - max(AmpsData(iData+1:end,:));
        AmpsSADMulti1(AmpsSADMulti1 < 0) = 0;
        AmpsSADMulti2 = min(AmpsData(iData+1:end,:)) - AmpsData(iData,:);
        AmpsSADMulti2(AmpsSADMulti2 < 0) = 0;
        AmpsSADMulti(iData) = sum(AmpsSADMulti1) + sum(AmpsSADMulti2);
        AmpsMADMulti(iData) = max(max(AmpsSADMulti1) + max(AmpsSADMulti2));
        AmpsMidSADMulti(iData) = sum(AmpsSADMulti1(Mid:end)) + sum(AmpsSADMulti2(Mid:end));
        AmpsMidMADMulti(iData) = max(max(AmpsSADMulti1(Mid:end)) + max(AmpsSADMulti2(Mid:end)));
        AmpsTopSADMulti(iData) = sum(AmpsSADMulti1(Top:end)) + sum(AmpsSADMulti2(Top:end));
        AmpsTopMADMulti(iData) = max(max(AmpsSADMulti1(Top:end)) + max(AmpsSADMulti2(Top:end)));
        
        RTsSADMulti1 = RTsData(iData,:) - max(RTsData(iData+1:end,:));
        RTsSADMulti1(RTsSADMulti1 < 0) = 0;
        RTsSADMulti2 = min(RTsData(iData+1:end,:)) - RTsData(iData,:);
        RTsSADMulti2(RTsSADMulti2 < 0) = 0;
        RTsSADMulti(iData) = sum(RTsSADMulti1) + sum(RTsSADMulti2);
        RTsMADMulti(iData) = max(max(RTsSADMulti1) + max(RTsSADMulti2));
        
        TwoDsSADMulti1 = TwoDsData(:,:,iData) - max(TwoDsData(:,:,iData+1:end),[],3);
        TwoDsSADMulti1(TwoDsSADMulti1 < 0) = 0;
        TwoDsSADMulti2 = min(TwoDsData(:,:,iData+1:end),[],3) - TwoDsData(:,:,iData);
        TwoDsSADMulti2(TwoDsSADMulti2 < 0) = 0;
        TwoDsSADMulti(iData) = sum(sum(TwoDsSADMulti1)) + sum(sum(TwoDsSADMulti2));
        TwoDsMADMulti(iData) = max(max(max(TwoDsSADMulti1)) + max(max(TwoDsSADMulti2)));
    elseif iData == dataPoints
        AmpsSADMulti1 = AmpsData(iData,:) - max(AmpsData(1:end-1,:));
        AmpsSADMulti1(AmpsSADMulti1 < 0) = 0;
        AmpsSADMulti2 = min(AmpsData(1:end-1,:)) - AmpsData(iData,:);
        AmpsSADMulti2(AmpsSADMulti2 < 0) = 0;
        AmpsSADMulti(iData) = sum(AmpsSADMulti1) + sum(AmpsSADMulti2);
        AmpsMADMulti(iData) = max(max(AmpsSADMulti1) + max(AmpsSADMulti2));
        AmpsMidSADMulti(iData) = sum(AmpsSADMulti1(Mid:end)) + sum(AmpsSADMulti2(Mid:end));
        AmpsMidMADMulti(iData) = max(max(AmpsSADMulti1(Mid:end)) + max(AmpsSADMulti2(Mid:end)));
        AmpsTopSADMulti(iData) = sum(AmpsSADMulti1(Top:end)) + sum(AmpsSADMulti2(Top:end));
        AmpsTopMADMulti(iData) = max(max(AmpsSADMulti1(Top:end)) + max(AmpsSADMulti2(Top:end)));
        
        RTsSADMulti1 = RTsData(iData,:) - max(RTsData(1:end-1,:));
        RTsSADMulti1(RTsSADMulti1 < 0) = 0;
        RTsSADMulti2 = min(RTsData(1:end-1:end,:)) - RTsData(iData,:);
        RTsSADMulti2(RTsSADMulti2 < 0) = 0;
        RTsSADMulti(iData) = sum(RTsSADMulti1) + sum(RTsSADMulti2);
        RTsMADMulti(iData) = max(max(RTsSADMulti1) + max(RTsSADMulti2));
        
        TwoDsSADMulti1 = TwoDsData(:,:,iData) - max(TwoDsData(:,:,1:end-1),[],3);
        TwoDsSADMulti1(TwoDsSADMulti1 < 0) = 0;
        TwoDsSADMulti2 = min(TwoDsData(:,:,1:end-1),[],3) - TwoDsData(:,:,iData);
        TwoDsSADMulti2(TwoDsSADMulti2 < 0) = 0;
        TwoDsSADMulti(iData) = sum(sum(TwoDsSADMulti1)) + sum(sum(TwoDsSADMulti2));
        TwoDsMADMulti(iData) = max(max(max(TwoDsSADMulti1)) + max(max(TwoDsSADMulti2)));
    else
        AmpsSADMulti1 = AmpsData(iData,:) - max(AmpsData([1:iData-1 iData+1:end],:));
        AmpsSADMulti1(AmpsSADMulti1 < 0) = 0;
        AmpsSADMulti2 = min(AmpsData([1:iData-1 iData+1:end],:)) - AmpsData(iData,:);
        AmpsSADMulti2(AmpsSADMulti2 < 0) = 0;
        AmpsSADMulti(iData) = sum(AmpsSADMulti1) + sum(AmpsSADMulti2);
        AmpsMADMulti(iData) = max(max(AmpsSADMulti1) + max(AmpsSADMulti2));
        AmpsMidSADMulti(iData) = sum(AmpsSADMulti1(Mid:end)) + sum(AmpsSADMulti2(Mid:end));
        AmpsMidMADMulti(iData) = max(max(AmpsSADMulti1(Mid:end)) + max(AmpsSADMulti2(Mid:end)));
        AmpsTopSADMulti(iData) = sum(AmpsSADMulti1(Top:end)) + sum(AmpsSADMulti2(Top:end));
        AmpsTopMADMulti(iData) = max(max(AmpsSADMulti1(Top:end)) + max(AmpsSADMulti2(Top:end)));
        
        RTsSADMulti1 = RTsData(iData,:) - max(RTsData([1:iData-1 iData+1:end],:));
        RTsSADMulti1(RTsSADMulti1 < 0) = 0;
        RTsSADMulti2 = min(RTsData([1:iData-1 iData+1:end],:)) - RTsData(iData,:);
        RTsSADMulti2(RTsSADMulti2 < 0) = 0;
        RTsSADMulti(iData) = sum(RTsSADMulti1) + sum(RTsSADMulti2);
        RTsMADMulti(iData) = max(max(RTsSADMulti1) + max(RTsSADMulti2));
        
        TwoDsSADMulti1 = TwoDsData(:,:,iData) - max(TwoDsData(:,:,[1:iData-1 iData+1:end]),[],3);
        TwoDsSADMulti1(TwoDsSADMulti1 < 0) = 0;
        TwoDsSADMulti2 = min(TwoDsData(:,:,[1:iData-1 iData+1:end]),[],3) - TwoDsData(:,:,iData);
        TwoDsSADMulti2(TwoDsSADMulti2 < 0) = 0;
        TwoDsSADMulti(iData) = sum(sum(TwoDsSADMulti1)) + sum(sum(TwoDsSADMulti2));
        TwoDsMADMulti(iData) = max(max(max(TwoDsSADMulti1)) + max(max(TwoDsSADMulti2)));
    end
end
AmpsSADMulti = max(AmpsSADMulti);
RTsSADMulti = max(RTsSADMulti);
TwoDsSADMulti = max(TwoDsSADMulti);
AmpsMADMulti = max(AmpsMADMulti);
RTsMADMulti = max(RTsMADMulti);
TwoDsMADMulti = max(TwoDsMADMulti);
AmpsMidSADMulti = max(AmpsMidSADMulti);
AmpsMidMADMulti = max(AmpsMidMADMulti);
AmpsTopSADMulti = max(AmpsTopSADMulti);
AmpsTopMADMulti = max(AmpsTopMADMulti);