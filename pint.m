function x=pint(X,PDF,p)

prb=0;
i=1;
DX=diff(X);
while prb<=p
  prb=prb+PDF(i);
  if prb>p
    IDX=DX(i)/PDF(i)*(prb-p);
    x=X(i+1)-IDX;
  end
  i=i+1;
end


