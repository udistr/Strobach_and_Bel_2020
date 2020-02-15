;==================================================================================  
; Heavyside function
;==================================================================================  

function heavyside(x)
local X,s,y
begin
X=ndtooned(x)
s=dimsizes(X)
y=new(s,float)
do i=0,s-1
	if X(i).gt.0 then
		y(i)=1
	else if X(i).lt.0
		y(i)=0
	else
		y(i)=0.5	
	end if
	end if
end do

return(reshape(y,dimsizes(x)))

end




function signum(x)
local X,s,y
begin
X=ndtooned(x)
s=dimsizes(X)
y=new(s,float)
do i=0,s-1
	if X(i).gt.0 then
		y(i)=1
	else if X(i).lt.0
		y(i)=-1
	else
		y(i)=0.5	
	end if
	end if
end do

return(reshape(y,dimsizes(x)))

end




