function y=kurt(xm,PM)
% KURT calculate kurtosis of distribution
%   y = ADDME(A) adds A to itself.
%   
%   xm - distribution cell centers 
%   PM - probability at cell center

y=(sum(PM.*(xm-sum(PM.*xm)).^4))/((sum(PM.*(xm-sum(PM.*xm)).^2)).^2);

end
