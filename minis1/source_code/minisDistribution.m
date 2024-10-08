function [nEvents, constMat] = minisDistribution(distributionType, distributionParameters, amplitudeArray, varargin)

ampLim2D = (amplitudeArray(2) - amplitudeArray(1))/2;
amplitudeArrayExt = -amplitudeArray(end): amplitudeArray(2)-amplitudeArray(1) :amplitudeArray(end);
riseTimeArray = varargin{1};
riseTimeArrayExt = -riseTimeArray(end): riseTimeArray(2)-riseTimeArray(1) :riseTimeArray(end);
RTLim2D = (riseTimeArray(2) - riseTimeArray(1))/2;
if strcmpi(distributionType, 'Normal')
    mu1 = distributionParameters(1);
    sigma1 = distributionParameters(2);
    scale = distributionParameters(3);
    mu2 = distributionParameters(4);
    sigma2 = distributionParameters(5);
    rho = distributionParameters(6);
    nEvents = normalBiv(mu1,sigma1,scale,mu2,sigma2,rho,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    constMat = nEvents;
elseif strcmpi(distributionType, 'Bimodal Normal')
    % First peak:
    mu1_1 = distributionParameters(1);
    sigma1_1 = distributionParameters(2);
    scale_1 = distributionParameters(3);
    mu2_1 = distributionParameters(4);
    sigma2_1 = distributionParameters(5);
    rho_1 = distributionParameters(6);
    nEvents_1 = normalBiv(mu1_1,sigma1_1,scale_1,mu2_1,sigma2_1,rho_1,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_1 < 0
        nEvents_1 = -nEvents_1;
    end
    % Second peak:
    mu1_2 = distributionParameters(7);
    sigma1_2 = distributionParameters(8);
    scale_2 = distributionParameters(9);
    mu2_2 = distributionParameters(10);
    sigma2_2 = distributionParameters(11);
    rho_2 = distributionParameters(12);
    nEvents_2 = normalBiv(mu1_2,sigma1_2,scale_2,mu2_2,sigma2_2,rho_2,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_2 < 0
        nEvents_2 = -nEvents_2;
    end
    nEvents = nEvents_1 + nEvents_2;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2);
elseif strcmpi(distributionType, 'Trimodal Normal')
    % First peak:
    mu1_1 = distributionParameters(1);
    sigma1_1 = distributionParameters(2);
    scale_1 = distributionParameters(3);
    mu2_1 = distributionParameters(4);
    sigma2_1 = distributionParameters(5);
    rho_1 = distributionParameters(6);
    nEvents_1 = normalBiv(mu1_1,sigma1_1,scale_1,mu2_1,sigma2_1,rho_1,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_1 < 0
        nEvents_1 = -nEvents_1;
    end
    % Second peak:
    mu1_2 = distributionParameters(7);
    sigma1_2 = distributionParameters(8);
    scale_2 = distributionParameters(9);
    mu2_2 = distributionParameters(10);
    sigma2_2 = distributionParameters(11);
    rho_2 = distributionParameters(12);
    nEvents_2 = normalBiv(mu1_2,sigma1_2,scale_2,mu2_2,sigma2_2,rho_2,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_2 < 0
        nEvents_2 = -nEvents_2;
    end
    % Third peak:
    mu1_3 = distributionParameters(13);
    sigma1_3 = distributionParameters(14);
    scale_3 = distributionParameters(15);
    mu2_3 = distributionParameters(16);
    sigma2_3 = distributionParameters(17);
    rho_3 = distributionParameters(18);
    nEvents_3 = normalBiv(mu1_3,sigma1_3,scale_3,mu2_3,sigma2_3,rho_3,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_3 < 0
        nEvents_3 = -nEvents_3;
    end
    nEvents = nEvents_1 + nEvents_2 + nEvents_3;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2, nEvents_3);
elseif strcmpi(distributionType, 'Quadrimodal Normal')
    % First peak:
    mu1_1 = distributionParameters(1);
    sigma1_1 = distributionParameters(2);
    scale_1 = distributionParameters(3);
    mu2_1 = distributionParameters(4);
    sigma2_1 = distributionParameters(5);
    rho_1 = distributionParameters(6);
    nEvents_1 = normalBiv(mu1_1,sigma1_1,scale_1,mu2_1,sigma2_1,rho_1,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_1 < 0
        nEvents_1 = -nEvents_1;
    end
    % Second peak:
    mu1_2 = distributionParameters(7);
    sigma1_2 = distributionParameters(8);
    scale_2 = distributionParameters(9);
    mu2_2 = distributionParameters(10);
    sigma2_2 = distributionParameters(11);
    rho_2 = distributionParameters(12);
    nEvents_2 = normalBiv(mu1_2,sigma1_2,scale_2,mu2_2,sigma2_2,rho_2,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_2 < 0
        nEvents_2 = -nEvents_2;
    end
    % Third peak:
    mu1_3 = distributionParameters(13);
    sigma1_3 = distributionParameters(14);
    scale_3 = distributionParameters(15);
    mu2_3 = distributionParameters(16);
    sigma2_3 = distributionParameters(17);
    rho_3 = distributionParameters(18);
    nEvents_3 = normalBiv(mu1_3,sigma1_3,scale_3,mu2_3,sigma2_3,rho_3,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_3 < 0
        nEvents_3 = -nEvents_3;
    end
    % Fourth peak:
    mu1_4 = distributionParameters(19);
    sigma1_4 = distributionParameters(20);
    scale_4 = distributionParameters(21);
    mu2_4 = distributionParameters(22);
    sigma2_4 = distributionParameters(23);
    rho_4 = distributionParameters(24);
    nEvents_4 = normalBiv(mu1_4,sigma1_4,scale_4,mu2_4,sigma2_4,rho_4,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D);
    if scale_4 < 0
        nEvents_4 = -nEvents_4;
    end
    nEvents = nEvents_1 + nEvents_2 + nEvents_3 + + nEvents_4;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2, nEvents_3, nEvents_4);
elseif strcmpi(distributionType, 'Log Normal')
    mu1 = distributionParameters(1);
    A = distributionParameters(2);
    alpha = distributionParameters(3);
    scale = distributionParameters(4);
    mu2 = distributionParameters(5);
    B = distributionParameters(6);
    beta = distributionParameters(7);
    rho = distributionParameters(8);
    amplitudeArrayExt = 0: amplitudeArray(2)-amplitudeArray(1) :amplitudeArray(end);
    riseTimeArrayExt = 0: riseTimeArray(2)-riseTimeArray(1) :riseTimeArray(end);
    nEvents = logNormalBiv(mu1,A,alpha,mu2,B,beta,rho,scale,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    constMat = nEvents;
elseif strcmpi(distributionType, 'Bimodal Log Normal')
    % First Peak:
    mu1_1 = distributionParameters(1);
    A_1 = distributionParameters(2);
    alpha_1 = distributionParameters(3);
    scale_1 = distributionParameters(4);
    mu2_1 = distributionParameters(5);
    B_1 = distributionParameters(6);
    beta_1 = distributionParameters(7);
    rho_1 = distributionParameters(8);
    amplitudeArrayExt = 0: amplitudeArray(2)-amplitudeArray(1) :amplitudeArray(end);
    riseTimeArrayExt = 0: riseTimeArray(2)-riseTimeArray(1) :riseTimeArray(end);
    nEvents_1 = logNormalBiv(mu1_1,A_1,alpha_1,mu2_1,B_1,beta_1,rho_1,scale_1,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Second Peak:
    mu1_2 = distributionParameters(9);
    A_2 = distributionParameters(10);
    alpha_2 = distributionParameters(11);
    scale_2 = distributionParameters(12);
    mu2_2 = distributionParameters(13);
    B_2 = distributionParameters(14);
    beta_2 = distributionParameters(15);
    rho_2 = distributionParameters(16);
    nEvents_2 = logNormalBiv(mu1_2,A_2,alpha_2,mu2_2,B_2,beta_2,rho_2,scale_2,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    nEvents = nEvents_1 + nEvents_2;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2);
elseif strcmpi(distributionType, 'Trimodal Log Normal')
    % First Peak:
    mu1_1 = distributionParameters(1);
    A_1 = distributionParameters(2);
    alpha_1 = distributionParameters(3);
    scale_1 = distributionParameters(4);
    mu2_1 = distributionParameters(5);
    B_1 = distributionParameters(6);
    beta_1 = distributionParameters(7);
    rho_1 = distributionParameters(8);
    amplitudeArrayExt = 0: amplitudeArray(2)-amplitudeArray(1) :amplitudeArray(end);
    riseTimeArrayExt = 0: riseTimeArray(2)-riseTimeArray(1) :riseTimeArray(end);
    nEvents_1 = logNormalBiv(mu1_1,A_1,alpha_1,mu2_1,B_1,beta_1,rho_1,scale_1,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Second Peak:
    mu1_2 = distributionParameters(9);
    A_2 = distributionParameters(10);
    alpha_2 = distributionParameters(11);
    scale_2 = distributionParameters(12);
    mu2_2 = distributionParameters(13);
    B_2 = distributionParameters(14);
    beta_2 = distributionParameters(15);
    rho_2 = distributionParameters(16);
    nEvents_2 = logNormalBiv(mu1_2,A_2,alpha_2,mu2_2,B_2,beta_2,rho_2,scale_2,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Third Peak:
    mu1_3 = distributionParameters(17);
    A_3 = distributionParameters(18);
    alpha_3 = distributionParameters(19);
    scale_3 = distributionParameters(20);
    mu2_3 = distributionParameters(21);
    B_3 = distributionParameters(22);
    beta_3 = distributionParameters(23);
    rho_3 = distributionParameters(24);
    nEvents_3 = logNormalBiv(mu1_3,A_3,alpha_3,mu2_3,B_3,beta_3,rho_3,scale_3,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    nEvents = nEvents_1 + nEvents_2 + nEvents_3;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2, nEvents_3);
elseif strcmpi(distributionType, 'Gaussian')
    mu1 = distributionParameters(1);
    sigma1 = distributionParameters(2);
    scale = distributionParameters(3);
    mu2 = distributionParameters(4);
    sigma2 = distributionParameters(5);
    theta = distributionParameters(6);
    nEvents = gaussianBiv(mu1,sigma1,scale,mu2,sigma2,theta,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    constMat = nEvents;
elseif strcmpi(distributionType, 'Bimodal Gaussian')
    % First peak:
    mu1_1 = distributionParameters(1);
    sigma1_1 = distributionParameters(2);
    scale_1 = distributionParameters(3);
    mu2_1 = distributionParameters(4);
    sigma2_1 = distributionParameters(5);
    theta_1 = distributionParameters(6);
    nEvents_1 = gaussianBiv(mu1_1,sigma1_1,scale_1,mu2_1,sigma2_1,theta_1,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Second peak:
    mu1_2 = distributionParameters(7);
    sigma1_2 = distributionParameters(8);
    scale_2 = distributionParameters(9);
    mu2_2 = distributionParameters(10);
    sigma2_2 = distributionParameters(11);
    theta_2 = distributionParameters(12);
    nEvents_2 = gaussianBiv(mu1_2,sigma1_2,scale_2,mu2_2,sigma2_2,theta_2,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    nEvents = nEvents_1 + nEvents_2;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2);
elseif strcmpi(distributionType, 'Trimodal Gaussian')
    % First peak:
    mu1_1 = distributionParameters(1);
    sigma1_1 = distributionParameters(2);
    scale_1 = distributionParameters(3);
    mu2_1 = distributionParameters(4);
    sigma2_1 = distributionParameters(5);
    theta_1 = distributionParameters(6);
    nEvents_1 = gaussianBiv(mu1_1,sigma1_1,scale_1,mu2_1,sigma2_1,theta_1,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Second peak:
    mu1_2 = distributionParameters(7);
    sigma1_2 = distributionParameters(8);
    scale_2 = distributionParameters(9);
    mu2_2 = distributionParameters(10);
    sigma2_2 = distributionParameters(11);
    theta_2 = distributionParameters(12);
    nEvents_2 = gaussianBiv(mu1_2,sigma1_2,scale_2,mu2_2,sigma2_2,theta_2,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Third peak:
    mu1_3 = distributionParameters(13);
    sigma1_3 = distributionParameters(14);
    scale_3 = distributionParameters(15);
    mu2_3 = distributionParameters(16);
    sigma2_3 = distributionParameters(17);
    theta_3 = distributionParameters(18);
    nEvents_3 = gaussianBiv(mu1_3,sigma1_3,scale_3,mu2_3,sigma2_3,theta_3,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    nEvents = nEvents_1 + nEvents_2 + nEvents_3;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2, nEvents_3);
elseif strcmpi(distributionType, 'Quadrimodal Gaussian')
    % First peak:
    mu1_1 = distributionParameters(1);
    sigma1_1 = distributionParameters(2);
    scale_1 = distributionParameters(3);
    mu2_1 = distributionParameters(4);
    sigma2_1 = distributionParameters(5);
    theta_1 = distributionParameters(6);
    nEvents_1 = gaussianBiv(mu1_1,sigma1_1,scale_1,mu2_1,sigma2_1,theta_1,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Second peak:
    mu1_2 = distributionParameters(7);
    sigma1_2 = distributionParameters(8);
    scale_2 = distributionParameters(9);
    mu2_2 = distributionParameters(10);
    sigma2_2 = distributionParameters(11);
    theta_2 = distributionParameters(12);
    nEvents_2 = gaussianBiv(mu1_2,sigma1_2,scale_2,mu2_2,sigma2_2,theta_2,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Third peak:
    mu1_3 = distributionParameters(13);
    sigma1_3 = distributionParameters(14);
    scale_3 = distributionParameters(15);
    mu2_3 = distributionParameters(16);
    sigma2_3 = distributionParameters(17);
    theta_3 = distributionParameters(18);
    nEvents_3 = gaussianBiv(mu1_3,sigma1_3,scale_3,mu2_3,sigma2_3,theta_3,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    % Fourth peak:
    mu1_4 = distributionParameters(19);
    sigma1_4 = distributionParameters(20);
    scale_4 = distributionParameters(21);
    mu2_4 = distributionParameters(22);
    sigma2_4 = distributionParameters(23);
    theta_4 = distributionParameters(24);
    nEvents_4 = gaussianBiv(mu1_4,sigma1_4,scale_4,mu2_4,sigma2_4,theta_4,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray);
    nEvents = nEvents_1 + nEvents_2 + nEvents_3 + nEvents_4;
    nEvents(nEvents < 0) = 0;
    constMat = cat(3, nEvents_1, nEvents_2, nEvents_3, nEvents_4);
elseif strcmpi(distributionType, 'Skew-normal')
    nEvents = skewNormalBiv(distributionParameters, amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt);
    nEvents(nEvents < 0) = 0;
    if sum(sum(isnan(nEvents)))
        nEvents(isnan(nEvents)) = 0;
    end
    if sum(sum(isinf(nEvents)))
        if sum(sum(isinf(nEvents))) == length(amplitudeArray)*length(riseTimeArray)
            nEvents(isinf(nEvents)) = 0;
        else
            nEvents(isinf(nEvents)) = max(nEvents(~isinf(nEvents)));
        end
    end
    constMat = nEvents;
elseif strcmpi(distributionType, 'Bimodal Skew-normal')
    nEvents1 = skewNormalBiv(distributionParameters([1:6 end-3:end-2]), amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt);
    nEvents2 = skewNormalBiv(distributionParameters([7:12 end-1:end]), amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt);
    nEvents = nEvents1 + nEvents2;
    nEvents(nEvents < 0) = 0;
    if sum(sum(isnan(nEvents)))
        nEvents(isnan(nEvents)) = 0;
    end
    if sum(sum(isinf(nEvents)))
        if sum(sum(isinf(nEvents))) == length(amplitudeArray)*length(riseTimeArray)
            nEvents(isinf(nEvents)) = 0;
        else
            nEvents(isinf(nEvents)) = max(nEvents(~isinf(nEvents)));
        end
    end
    constMat = cat(3, nEvents1, nEvents2);
elseif strcmpi(distributionType, 'Trimodal Skew-normal')
    nEvents1 = skewNormalBiv(distributionParameters([1:6 end-5:end-4]), amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt);
    nEvents2 = skewNormalBiv(distributionParameters([7:12 end-3:end-2]), amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt);
    nEvents3 = skewNormalBiv(distributionParameters([13:18 end-1:end]), amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt);
    nEvents = nEvents1 + nEvents2 + nEvents3;
    nEvents(nEvents < 0) = 0;
    if sum(sum(isnan(nEvents)))
        nEvents(isnan(nEvents)) = 0;
    end
    if sum(sum(isinf(nEvents)))
        if sum(sum(isinf(nEvents))) == length(amplitudeArray)*length(riseTimeArray)
            nEvents(isinf(nEvents)) = 0;
        else
            nEvents(isinf(nEvents)) = max(nEvents(~isinf(nEvents)));
        end
    end
    constMat = cat(3, nEvents1, nEvents2, nEvents3);
end

draw = 0;
if draw == 1
    if dimension == 1
        f2 = figure('position', [640 300 600 300]);
        figure(f2);
        plot(amplitudeArray,nEvents,'k');
        title('Distribution of simulated minis');
        xlabel('Amplitude(mV)');
        ylabel('Number of Events');
    end
end
end










function nEvents = normalBiv(mu1,sigma1,scale,mu2,sigma2,rho,amplitudeArray,riseTimeArray,ampLim2D,RTLim2D)

mu = [mu1 mu2];
covarMatrix = [sigma1^2 rho*sigma1*sigma2; rho*sigma1*sigma2 sigma2^2];
minorDet = det(covarMatrix(1,1));
majorDet = det(covarMatrix);
while minorDet <= 0 || majorDet <= 0
    %disp('error: covariance matrix is not positive definite. Scalling down the off-diagonal elements.');
    if covarMatrix(1,2) > 0
        covarMatrix = [covarMatrix(1,1) covarMatrix(1,2)-.01^2; covarMatrix(2,1)-.01^2 covarMatrix(2,2)];
    elseif covarMatrix(1,2) < 0
        covarMatrix = [covarMatrix(1,1) covarMatrix(1,2)+.01^2; covarMatrix(2,1)+.01^2 covarMatrix(2,2)];
    end
    minorDet = det(covarMatrix(1,1));
    majorDet = det(covarMatrix);
end
R = chol(covarMatrix);
events = repmat(mu,round(abs(scale)),1) + randn(round(abs(scale)),2)*R;
if ~isempty(events) && (min(events(:,1)) < amplitudeArray(1)-ampLim2D || max(events(:,1)) > amplitudeArray(end)+ampLim2D...
        || min(events(:,2)) < riseTimeArray(1)-RTLim2D || max(events(:,2)) > riseTimeArray(end)+RTLim2D)
    filtIndicesAmp = find(events(:,1) > amplitudeArray(1)-ampLim2D & events(:,1) < amplitudeArray(end)+ampLim2D);
    filtIndicesAmpLogical = zeros(length(events),1);
    filtIndicesAmpLogical(filtIndicesAmp) = 1; %#ok<FNDSB>
    filtIndicesRT = find(events(:,2) > riseTimeArray(1)-RTLim2D & events(:,2) < riseTimeArray(end)+RTLim2D);
    filtIndicesRTLogical = zeros(length(events),1);
    filtIndicesRTLogical(filtIndicesRT) = 1; %#ok<FNDSB>
    events = events(filtIndicesAmpLogical & filtIndicesRTLogical, :);
end
if isempty(events)
    nEvents = zeros(length(riseTimeArray),length(amplitudeArray));
else
    nEvents = hist2d(events(:,1),events(:,2),amplitudeArray,riseTimeArray);
end
end



function nEvents = logNormalBiv(mu1,A,alpha,mu2,B,beta,rho,scale,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray)

[X,Y] = meshgrid(amplitudeArrayExt,riseTimeArrayExt);
X = X - mu1;
X(X < 0) = 0;
Y = Y - mu2;
Y(Y < 0) = 0;
logNormal = ((2*pi*(X).*(Y)*alpha*beta*((1.000001-rho^2)^1/2)).^(-1)).*exp((-1/2*(1.000001-rho^2))...
    *(((log((X)/A)/alpha).^2) - 2*rho*(log((X)/A)/alpha).*(log((Y)/B)/beta) + (log((Y)/B)/beta).^2));
if sum(sum(isnan(logNormal)))
    logNormal(isnan(logNormal)) = 0;
end
if sum(sum(logNormal))
    logNormal = abs(scale*logNormal/sum(sum(logNormal)));
end
nEvents = logNormal(end-length(riseTimeArray)+1:end,end-length(amplitudeArray)+1:end);
end



function nEvents = gaussianBiv(mu1,sigma1,scale,mu2,sigma2,theta,amplitudeArrayExt,amplitudeArray,riseTimeArrayExt,riseTimeArray)

a = cos(theta)^2/2*sigma1^2 + sin(theta)^2/2*sigma2^2;
b = -sin(2*theta)/4*sigma1^2 + sin(2*theta)/4*sigma2^2;
c = sin(theta)^2/2*sigma1^2 + cos(theta)^2/2*sigma2^2;
[X,Y] = meshgrid(amplitudeArrayExt,fliplr(riseTimeArrayExt));
gaussian = exp(-(a*(X-mu1).^2 - 2*b*(X-mu1).*(Y-mu2) + c*(Y-mu2).^2));
gaussian = abs(scale*gaussian/sum(sum(gaussian)));
nEvents = gaussian(end-length(riseTimeArray)+1:end,end-length(amplitudeArray)+1:end);
end



function nEvents = skewNormalBiv(distributionParameters, amplitudeArray, riseTimeArray, amplitudeArrayExt, riseTimeArrayExt)

mu1 = distributionParameters(1);
sigma1 = distributionParameters(2);
scale = distributionParameters(3);
mu2 = distributionParameters(4);
sigma2 = distributionParameters(5);
rho = distributionParameters(6);
skew1 = distributionParameters(7);
skew2 = distributionParameters(8);
covarMatrix = [sigma1^2 rho*sigma1*sigma2; rho*sigma1*sigma2 sigma2^2];

minorDet = det(covarMatrix(1,1));
majorDet = det(covarMatrix);
while minorDet <= 0 || majorDet <= 0
    %disp('error: covariance matrix is not positive definite. Scalling down the off-diagonal elements.');
    if covarMatrix(1,2) > 0
        covarMatrix = [covarMatrix(1,1) covarMatrix(1,2)-.01^2; covarMatrix(2,1)-.01^2 covarMatrix(2,2)];
    elseif covarMatrix(1,2) < 0
        covarMatrix = [covarMatrix(1,1) covarMatrix(1,2)+.01^2; covarMatrix(2,1)+.01^2 covarMatrix(2,2)];
    end
    minorDet = det(covarMatrix(1,1));
    majorDet = det(covarMatrix);
end

[x,y] = meshgrid(amplitudeArrayExt,riseTimeArrayExt); x = x(:); y = y(:);

% Normal CDF with skew:
skew = [skew1 skew2]*(covarMatrix^-.5);
normCDF = .5*(1 + erf(((x-mu1)*skew(1) + (y-mu2)*skew(2))/sqrt(2)));
normCDF = reshape(normCDF, length(riseTimeArrayExt), length(amplitudeArrayExt));

% Normal PDF:
normPDF = reshape((1/(2*pi*sigma1*sigma2*((1-rho^2)^.5)))*exp((-1/(2*(1-rho^2)))*(((x-mu1)/sigma1).^2 - 2*rho*((x-mu1)/sigma1).*((y-mu2)/sigma2)...
    + ((y-mu2)/sigma2).^2)), length(riseTimeArrayExt), length(amplitudeArrayExt));

% Skew-normal:
skewNorm = 2*normPDF.*normCDF;
nEventsExt = scale*skewNorm/sum(sum(skewNorm));
nEvents = nEventsExt(end-length(riseTimeArray)+1:end,end-length(amplitudeArray)+1:end);
end