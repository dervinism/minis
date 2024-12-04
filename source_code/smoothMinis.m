function minis2D = smoothMinis(minis2D, smoothWindow, f)

for iRT = 1:size(minis2D,1)
    minis2D(iRT,:) = gsmooth(minis2D(iRT,:), smoothWindow, f);
end