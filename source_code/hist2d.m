function [counto,xco,yco] = hist2d(x,y,nx,ny)
% HIST2D calculates a 2-dimensional histogram
%    N = HIST2D(X,Y) bins the X and Y data into 10 equally spaced bins in
%    both dimension 
%
%    N = HIST2D(X,Y,M), where M is a scalar, bins the data into M equally
%    spaced bins in both dimensions
%
%    N = HIST2D(X,Y,B), where B is a vector, bins the data with centers
%    specified by B 
%
%    The number of bins or centers can be specified individually for either
%    dimension using N = HIST2D(X,Y,NX,NY) or N = HIST2D(X,Y,BX,BY)
%
%    [N,BX,BY] = HIST2D(...) also returns the positions of the bin centers
%    as two matrices in BX and BY
%
%    HIST2D(...) without output arguments produces a colormapped image plot
%    of the 2d histogram
%
% EXAMPLE
%   yin = randn(1,1000);
%   xin = randn(1,1000);
%   [n,x,y] = hist2d(xin,yin,11);
%   imagesc(x(1,:),y(:,1),n); hold on; plot(xin,yin,'y.'); colorbar

if ~exist('nx','var')
   nx = 10;
end

if ~exist('ny','var')
   ny = nx;
end
   
[~,xc] = hist(x,nx); %#ok<*HIST>
[~,yc] = hist(y,ny);

count = zeros(length(yc),length(xc));

for i = 1:length(yc)
   if i == 1
      lbound = -Inf;
   else
      lbound = (yc(i-1) + yc(i)) / 2;
   end
   if i == length(yc)
      ubound = inf;
   else
      ubound = (yc(i) + yc(i+1)) /2;
   end
   count(i,:) = hist(x((y >= lbound) & (y < ubound)),xc);
   %count(length(yc)-i+1,:) = hist(x((y >= lbound) & (y < ubound)),xc);
end

if nargout == 0
   [xc, yc] = meshgrid(xc, yc);
   imagesc(xc(1,:),yc(:,1),count);
else
   counto = count;
   xco = xc;
   yco = yc;
end