function cliff = detectCliff(minis2D)

cliff = [0 0];
cliff2D = false;
cliffAmps = true;

% Detect a cliff in the minis 2-D distribution
centiles = 99:-1:50;
cost = centiles/100;
if cliff2D
  minisVector = reshape(minis2D,1,size(minis2D,1)*size(minis2D,2)); %#ok<*UNRCH>
  minisVector(minisVector == 0) = [];
  zerosAdd1 = zeros(size(minis2D,1),1);
  zerosAdd2 = zeros(1,size(minis2D,2)+2);
  onesAdd2 = ones(1,size(minis2D,2)+2);
  minis2DExt = [zerosAdd1 minis2D zerosAdd1];
  minis2DExt = [onesAdd2; minis2DExt; zerosAdd2];
  for iCent = 1:numel(centiles)
    cent = prctile(minisVector,centiles(iCent));
    [I, J] = find(minis2DExt >= cent);
    if ~isempty(I) && numel(I) > 0
      for iMed = 1:numel(I)
        cliffMat = minis2DExt(max([1 I(iMed)-1]):I(iMed)+1, max([1 J(iMed)-1]):J(iMed)+1);
        iCliff = find(cliffMat == 0, 2);
        if ~isempty(iCliff) && numel(iCliff) == 2
          cliff(1) = cost(iCent); % Extra cost
          break
        end
      end
      if cliff(1)
        break
      end
    end
  end
end

% Detect a cliff in the minis 1-D amplitude distribution
if cliffAmps
  minis1D_Amp = sum(minis2D);
  minisVector = minis1D_Amp;
  minisVector(minisVector == 0) = [];
  minis1D_Amp = [0 minis1D_Amp 0];
  for iCent = 1:numel(centiles)
    cent = prctile(minisVector,centiles(iCent));
    J = find(minis1D_Amp >= cent);
    for iMed = 1:numel(J)
      cliffMat = minis1D_Amp(J(iMed)-1:J(iMed)+1);
      iCliff = find(cliffMat == 0, 1);
      if ~isempty(iCliff) && numel(iCliff) > 0
        cliff(2) = cost(iCent); % Extra cost
        return
      end
    end
  end
end
end