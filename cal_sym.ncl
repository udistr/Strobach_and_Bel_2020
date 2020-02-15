function cal_sym(f,w,y,treshold)
local p,sigma,rmse,lim,n0,iu,id,pru,prd

begin

p=dim_sum_n(w*f,1)

sigma=sqrt(dim_sum_n(w*(f-conform(f,dim_sum_n(w*f,1),0))^2,1))
rmse=sqrt(avg((p-y)^2))

lim=max(abs(p-y)/sigma)
print(lim)
delta=10^floor(log10(lim))*10^-2
print(delta)
n0=floattointeger(ceil(lim/delta))
print(n0)
lrt=dimsizes(y)
g=conform_dims((/lrt,n0/),fspan(delta,n0*delta,n0),1);

;fraction of observations inside
pru=dim_avg_n(heavyside((conform(g,p,0)+g*conform(g,sigma,0))-conform(g,y,0)),0)
prd=dim_avg_n(heavyside(conform(g,y,0)-(conform(g,p,0)-g*conform(g,sigma,0))),0)

iu=min(ind(pru.ge.(1-treshold/2)))
id=min(ind(prd.ge.(1-treshold/2)))

cal=new(2,float)
;cal(0) is the upper limit gu
;cal(1) is the lower limit gd
cal=(/g(0,iu),g(0,id)/)
return(cal)

end


