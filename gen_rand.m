function y=gen_rand(n,XI,PI)

x=rand([1 n]);
c=[0 cumsum(PI)];
y=zeros(size(x));
for i=1:length(PI)
  ind=find(x>c(i)&x<=c(i+1));
  y(ind)=XI(i)+rand([1 length(ind)]).*(XI(i+1)-XI(i));
end
  
