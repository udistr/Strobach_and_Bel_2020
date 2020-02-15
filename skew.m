function y=skew(xm,PM)
% SKEW calculate skewness of distribution
%   y = ADDME(A) adds A to itself.
%   
%   xm - distribution cell centers 
%   PM - probability at cell center

y=(sum(PM.*(xm-sum(PM.*xm)).^3))/((sum(PM.*(xm-sum(PM.*xm)).^2)).^(3/2));

end
