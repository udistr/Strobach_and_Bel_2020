;==================================================================================  
; the gradient of the square loss function used by EGA
;==================================================================================
function loss_ega(f,y,p)
local f,y,p
begin

	l=2*(p-y)*f
	return(l)

end

;==================================================================================  
;The weight of the Exponential Gradient Average (EGA) Forecaster
;y(n,ltlon) is the array of observation for n time steps and ltln grid cells
;f(n,N,ltln) is the array of the forecasts by N models.
;e is a scalar of the learning rate
;The return value is the weight of each expert.
;==================================================================================  
function exp_grad_weight(f,y,e)
local dim,n,N,ltln,los,w,Y,m,eta,w_start,los,p
begin

	;************************
	; determining dimentions
	;************************
	dim=dimsizes(f)
	
	if dimsizes(dim).eq.3 then
	  lrt=dim(0)
	  N=dim(1)
	  ltln=dim(2)
	else if dimsizes(dim).eq.2
	  lrt=dim(0)
	  N=dim(1)
	  ltln=1
	else
	  print("dimention number not supported")
	  exit
	end if
	end if

	;************************
	; defining weights array
	;************************
	w=new((/lrt,N,ltln/),double)
	w!0="time"
	w!1="model"
	w!2="latlon"

	;************************
	;bound for the loss	
	;************************
	Y=conform(f,y,(/0,2/))
	M=conform(f(0,:,:),dim_max_n(dim_max_n(abs(loss_ega(f,Y,f)),0),0),(/1/))

	;************************************
	;eta and starting weight definitions	
	;************************************
	eta=new(dimsizes(e),double)
	eta=e;*sqrt(2*log(N)/n)
	w_start=new((/N,ltln/),double)
	w_start=1./N

	;************************************
	;first time step weight	
	;************************************
	p = new((/lrt+1,ltln/),double)
	p(0,:)=dim_sum_n(w_start*f(0,:,:),0)
	los=loss_ega(f(0,:,:),Y(0,:,:),conform(Y(0,:,:),p(0,:),(/1/)))/M
	w(0,:,:)=w_start*exp(-conform(los,eta,(/1/))*los)/conform(w_start,dim_sum_n(w_start*exp(-conform(los,eta,(/1/))*los),0),(/1/)) ; equation 3 in strobach and bel (2016)

	;************************************
	;Updating the weights	
	;************************************
	do t=1,lrt-1
		p(t,:)=dim_sum_n(w(t-1,:,:)*f(t,:,:),0)
		los=(loss_ega(f(t,:,:),Y(t,:,:),conform(Y(t,:,:),p(t,:),(/1/))))/M
		w(t,:,:)=w(t-1,:,:)*exp(-conform(los,eta,(/1/))*los)/conform(w(t-1,:,:),dim_sum_n(w(t-1,:,:)*exp(-conform(los,eta,(/1/))*los),0),(/1/)) ; equation 3 in strobach and bel (2016)
		; avoid 0 weight that may be caused due to computer accuracy
		w(t,:,:)=w(t,:,:)>10^-6 
		w(t,:,:)=w(t,:,:)/conform(w(t,:,:),dim_sum_n(w(t,:,:),0),(/1/))
	end do

	;************************************
	;Prepering the output	
	;************************************
	w1=new((/lrt+1,N,ltln/),double)
	w1(0,:,:)=w_start
	w1(1:,:,:)=w
	w1@M=M(0,:)
	return(dble2flt(w1))

end


;==================================================================================  
; This is our outer loop of the EGA forcaster. Here ww choose the learning
; rate the pruduce the minimum error during the larning period
;==================================================================================  
function EGA(f,y)
local start,mind,eta0,expf0,weight0,eta1,expf1,weight1,deta,metric,lp
begin

	;************************
	; determining dimentions
	;************************
	dim=dimsizes(f)
	
	if dimsizes(dim).eq.3 then
	  lrt=dim(0)
	  N=dim(1)
	  ltln=dim(2)
	  f1=f
	  y1=y
	else if dimsizes(dim).eq.2
	  lrt=dim(0)
	  N=dim(1)
	  ltln=1
	  f1=new((/lrt,N,ltln/),float)
	  f1(:,:,0)=f
	  y1=new((/lrt,ltln/),float)
	  y1(:,0)=y
	else
	  print("dimention number not supported")
	end if
	end if
;*************************************************************
; Defining variables
;*************************************************************

; vector that contains the best eta for each grid cell and updates after every iteration of the resolution
start=new(ltln,float) 
; index of minimum metric value
mind=new(ltln,integer)
; optimal learning parameter
eta1=new(ltln,float)
; predictions of the algorithm
expf1=new((/lrt,ltln/),float)
; weights of the models
weight1=new((/lrt+1,N,ltln/),float)
	
; first iteration
start=0. ; at the beginning eta=0 everywhere
n=69 ; we increase eta n times
deta=10. ; each time by deta

do res=0,3 ; 4 iterations
;("res=" + (/res/))
	
	; array that contains the metric for each eta
	metric=new((/ltln,n/),float)
	; n etas
	eta0=new((/ltln,n/),float)
	; n predictions for each eta
	expf0=new((/lrt,ltln,n/),float)
	; n weights for each eta and model
	weight0=new((/lrt+1,N,ltln,n/),float)

	do lp=0,n-1	; n predictions for each eta
		;print("lp=" + (/lp/))
		eta0(:,lp)=start+lp*deta ; adding lp*deta every loop

		;*************************************************************
		; Running the algorithm for a given eta
		;*************************************************************
		
	  temp_weight=exp_grad_weight(f1(0:(lrt-1),:,:),y1(0:(lrt-1),:),eta0(:,lp))

		; output weight for a given eta
		weight0(0:(lrt),:,:,lp)=(/temp_weight/)
		; output forecast for a given eta
		expf0(0:(lrt-1),:,lp)=dim_sum_n(weight0(0:(lrt-1),:,:,lp)*f1(0:(lrt-1),:,:),1) 

		rmse=dim_rmsd_n(expf0(0:(lrt-1),:,lp),y1(0:(lrt-1),:),0)
		if ((lp.eq.0).and.(res.eq.0)) then
			rmse0=dim_rmsd_n(dim_avg_n(f1(0:(lrt-1),:,:),1),y1(0:(lrt-1),:),0)
		end if
		rmse=rmse/rmse0 ; RMSE relative to the RMSE of simple average

		; calculation of the stability creteria for the eta optimization metric
		; see strobach and bel (2016) equation 6, the term inside the square brackets.
		stab=floor(dim_max_n(dim_max_n(abs(weight0(1:(lrt),:,:,lp)-weight0(0:(lrt-1),:,:,lp)),0),0)*N)

		;calculation of the matric
		metric(:,lp)=round(rmse*(1+stab)*10^5,0)/10^5

	end do
	
	do i=0,(ltln-1)
		mind(i)=minind(metric(i,:)) ; index of minimum metric value
		eta1(i)=eta0(i,mind(i)) ; choosing best eta
		expf1(:,i)=expf0(:,i,mind(i)) ; predictions of the best eta
		weight1(:,:,i)=weight0(:,:,i,mind(i)) ; weight of best eta
	end do

	deta=0.1*deta ; increasing eta resolution for the next iteration
	start=(eta1-9*deta)>deta ; defining new start points for the eta scanning
	n=19 ; running between eta-9*deta to eta+9*deta - 19 scanns

	; cleanning variables for the next iteration
	delete(expf0)
	delete(weight0)
	delete(metric)
	delete(eta0)

end do

if dimsizes(dim).eq.3 then
  return(weight1(lrt,:,:))
else if dimsizes(dim).eq.2
  return(weight1(lrt,:,0))
end if
end if

end

