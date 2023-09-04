function cdat = gsmooth(dat,win,f,filtfunc)
% Over-rides in-built function smooth
% for gaussian, f=3,
%%[cut-off freq -3dB = .3748*sample_freq/(win)  where win is no. of data pts. in filter function
% so 30 pts = 250 Hz ~ 4ms
% 40 pts ~ sample_freq/100
% 20 pts ~ 0.019*sample_freq ~ sample_freq/50
% 10 pts ~ .038*sample_freq
% 4 pts ~ sample_freq/10  all error in Matlab manual! ]
%% win indicates +/- win # of points about a given point to weight over.
%% total smoothing window will be 2*win+1
%% supplied filtfunc should be normalized
%% the first and last win points will not be more lightly weighted than the rest.
% modified to include gaussian GM 2/24/02
% half width = 2* sqrt(2*ln(2)) sigma since exp(-x^2/2) = .5 when x= sqrt(2ln(2))
% =  2* sqrt(2*ln(2))  * win/(2*sqrt(2)) =  win*sqrt(ln(2)) =0.8326*win
% so new effective sample interval at half height = 0.8326*win*sample_interval
% so new effective frequency roll off (half height = -3dB) =
% 1/sqrt(log(2))*(sample_frequency/win) = 1.2*sample_frequency/win <-------------
if ~exist('filtfunc','var')
    if f == 1 % box window
        filtfunc = ones(1,2*win+1);
    end
    if f == 2 % triangular window
        filtfunc = [(0:win-1) win (win-1:-1:0)];
    end
    if f == 3 % gaussian with gsd=win/2  so real sd = win/(2*sqrt(2)) = win/2.828 = .3536*win
        gsd = win/2; % gsd is sqrt(2)* real sd
        for i = 1:win-1 % should be exp(-.5*(i/sd)^2);
            filtfunc(win-i) = exp(-(i/gsd)^2);  % 1.8% at edges should be e-(x^2/2sd^2)
            filtfunc(win+i) = filtfunc(win-i);
        end
        filtfunc(win) = 1;
    end
end
filtfunc = filtfunc/sum(filtfunc);
lrmv = length(filtfunc)-1;

cdat = conv(dat,filtfunc);
lc = length(cdat);
cdat = cdat(1+round(lrmv/2):lc-round(lrmv/2));
cdat(1:win) = dat(1:win);
cdat(end-win:end) = dat(end-win:end); % to avoid end effects due to asymmetrical windowing
if length(cdat) >length(dat)
    cdat = cdat(1:length(dat));
elseif ~isempty(cdat<length(dat))
    %dl = length(dat) - length(cdat);
    %cat= [cdat dat(end-dl+1:end)];
end
end