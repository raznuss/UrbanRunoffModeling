function [FSE,BIAS,CORR,RMSE,NSD,NASH]=statCalc(X,Y,varargin)
% 
%   F. Marra - 2013
% 
% calculates the comparison statistical parameters between two sets of
% values
% 
% Input:
% - X is the independent variable (gauge)
% - Y is the dependent variable (radar)
% X and Y must be of the same size
% 
% Output:
% - FSE: Fractional Standard Error
% - BIAS: Bias (in ratio form: Y/X)
% - CORR: Correlation index
% - RMSE: Root Mean Square Error
% - NSD: Normalized Standard Difference (equals the FSE of the unbiased values)
% - NASH: Nash-Sutcliff index

FSE=NaN; BIAS=NaN; CORR=NaN; RMSE=NaN; NSD=NaN; NASH=NaN;
if size(X)~=size(Y); warning('input error: size(X)~=size(Y);'); return; end

if nargin>2
    minRain=varargin{1};
else
    minRain=0.1;
end

vecX=reshape(X,numel(X),1);
vecY=reshape(Y,numel(Y),1);


% removes NaNs and uses only pairs of values where at least one value is
% higher than the minRain=0.1 threhsold
if minRain~=0
    condition=~isnan(vecX) & ~isnan(vecY) & (vecX>minRain | vecY>minRain);
elseif minRain==0
    condition=~isnan(vecX) & ~isnan(vecY);
end    

% % removes NaNs - old version
% condition=~isnan(vecX) & ~isnan(vecY);

x=vecX(condition);
y=vecY(condition);


if numel(x)<=1
    BIAS = sum(y)/sum(x);
    return
end

sx=sum(x);
nx=numel(x);
sqdx=sum((x-y).^2);

FSE=sqrt(sqdx*nx)/sx;
% FSE=sqrt(sum((x-y).^2)*numel(x))/sx;

BIAS=sum(y)/sx;

corr=corrcoef(x,y);
CORR=corr(1,2);

% 

RMSE = sqrt(sqdx/nx);
% RMSE = sqrt(sum((x-y).^2)/numel(x));

if nargout>4
    NSD = statCalc(X,Y./BIAS,minRain);
        
    if nargout>5
        NASH=1-sum((x-y).^2)/sum((x-nanmean(x)).^2);
    end    
end
    



end

% function [FSE,BIAS,CORR,RMSE]=statCalc(X,Y,varargin)
% % 
% %   F. Marra - 2013
% % 
% % calculates the comparison statistical parameters between two sets of
% % values
% % 
% % Input:
% % - X is the independent variable (gauge)
% % - Y is the dependent variable (radar)
% % X and Y must be of the same size
% % 
% % Output:
% % - FSE: Fractional Standard Error
% % - BIAS: Bias (in ratio form: Y/X)
% % - CORR: Correlation index
% % - RMSE: Root Mean Square Error
% % 
% 
% FSE=NaN; BIAS=NaN; CORR=NaN; RMSE=NaN;
% if size(X)~=size(Y); warning('input error: size(X)~=size(Y);'); return; end
% 
% if nargin>2
%     minRain=varargin{1};
% else
%     minRain=0.1;
% end
% 
% vecX=reshape(X,numel(X),1);
% vecY=reshape(Y,numel(Y),1);
% 
% 
% % removes NaNs and uses only pairs of values where at least one value is
% % higher than the minRain=0.1 threhsold
% if minRain~=0
%     condition=~isnan(vecX) & ~isnan(vecY) & (vecX>minRain | vecY>minRain);
% else
%     condition=~isnan(vecX) & ~isnan(vecY);
% end    
% 
% % % removes NaNs - old version
% % condition=~isnan(vecX) & ~isnan(vecY);
% 
% x=vecX(condition);
% y=vecY(condition);
% 
% 
% if numel(x)<=1
%     BIAS = sum(y)/sum(x);
%     return
% end
% 
% sx=sum(x);
% nx=numel(x);
% sqdx=sum((x-y).^2);
% 
% FSE=sqrt(sqdx*nx)/sx;
% % FSE=sqrt(sum((x-y).^2)*numel(x))/sx;
% 
% BIAS=sum(y)/sx;
% 
% corr=corrcoef(x,y);
% CORR=corr(1,2);
% 
% % NASH=1-sum((x-y).^2)/sum((x-nanmean(x)).^2);
% 
% RMSE = sqrt(sqdx/nx);
% % RMSE = sqrt(sum((x-y).^2)/numel(x));
% 
% 
% end