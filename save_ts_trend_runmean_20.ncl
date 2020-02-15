load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "Heavyside.ncl"
load "cal_sym.ncl"
load "EGA.ncl"
begin

a=addfile("tas.nc","r")
;rcp=1
R="NCEP"
;SLA="EGA"

print(R)
print(rcp)

if rcp.eq.0
  RCP=26
else if rcp.eq.1
  RCP=45
else if rcp.eq.2
  RCP=60
else if rcp.eq.3
  RCP=85
end if
end if
end if
end if
  
nt=152  

tmp=a->tas(:,rcp,nt-1)
aind=ind(.not.ismissing(tmp))
tas=a->tas(aind,rcp,:)
delete(tmp)

lrt=69

b=addfile("../NCEP/air.mon.mean.nc","r")
nlat=73
nlon=144
pi=3.14159265359
L=new((nlat+1),float)
L(1:(nlat-1))=fspan((90-1.25),(-90+1.25),72)
L(0)=90
L(73)=-90
L1=abs(L(0:(nlat-1))-90)/180*pi
L2=abs(L(1:nlat)-90)/180*pi
W=(conform_dims((/nlat,nlon/),(cos(L1)-cos(L2))/2,0))
W=W/sum(W)

ncep=tofloat(dim_avg_n(reshape(dim_sum_n(dim_sum_n(conform(b->air(0:lrt*12-1,:,:),W,(/1,2/))*b->air(0:lrt*12-1,:,:),1),1),(/lrt,12/)),1))
y=(/ncep/)
  
f=transpose(tas)

nmean=20

tmp=new(dimsizes(f),float)
do i=nmean-1,nt-1
 tmp(i,:)=(/dim_avg_n(f(i-nmean+1:i,:),0)/)
end do
f=(/tmp/)
delete(tmp)

tmp=new(dimsizes(y),float)
do i=nmean-1,lrt-1
 tmp(i)=(/dim_avg_n(y(i-nmean+1:i),0)/)
end do
y=(/tmp/)
delete(tmp)

f(0:nmean-2,:)=f@_FillValue
y(0:nmean-2)=f@_FillValue

f=f-conform(f,dim_avg_n(f(0:lrt-1,:),0),1)+dim_avg_n(y,0)

f1=f(nmean-1:lrt-1,:)
y1=y(nmean-1:lrt-1)

w=new(dimsizes(f1),float)
if SLA.eq."EGA" then
  w=conform(f1,EGA(f1,y1),1)
  print(EGA(f1,y1))
else if SLA.eq."AVG"
  w=1./dimsizes(aind)
end if
end if

sig=fspan(0.,0.9,10)
cal=new((/dimsizes(sig),2,3/),float)

cal0=new((/2,3/),float)

do i=0,dimsizes(sig)-1
  cal(i,:,0)=cal_sym(f1,w,y1,sig(i))
  cal(i,:,1)=0
  cal(i,:,2)=sig(i)/2 
end do
print(cal)

weight=w(0,:)
nm=dimsizes(weight)

model=aind

time=ispan(0,nt-1,1)

;===================================================================  
; Assume variables T, PS and ORO exist and that they have 
; associated meta data: (a) coordinate variables time, lev, lat, lon       
; and (b) attributes
;===================================================================  

var="cal"

diro = "./CAL/"                     ; Output directory
filo = var + "_RCP" + RCP + "_" + R + "_" +  SLA + "_runmean_20.nc"             ; Output file
system("/bin/rm -f " + diro + filo)    ; remove if exists
fout  = addfile (diro + filo, "c")  ; open output file

;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
setfileoption(fout,"DefineMode",True)

;===================================================================
; create global attributes of the file
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "Calibration factors"  
fAtt@source_file   =  "original-file.nc"
fAtt@Conventions   = "None"   
fAtt@creation_date = systemfunc ("date")        
fileattdef( fout, fAtt )            ; copy file attributes    

;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================
dimNames = (/ "time", "sig", "side" , "val","model"/)  
dimSizes = (/ nt,  dimsizes(sig)   , 2  ,  3 , nm /) 
dimUnlim = (/ True, False, False, False, False /)   
filedimdef(fout,dimNames,dimSizes,dimUnlim)

;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
; Here we are using NCL functions to facilitate defining 
; each variable's dimension name(s) and type. 
; The following could be replaced with explicit, user defined dimension 
; names different from those associated with the variable in memory. 
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/)) 
;===================================================================
filevardef(fout, "sig"   ,"float"  ,(/"sig"/))    
filevardef(fout, "side" ,"integer",(/"side"/))                          
filevardef(fout, "val"    ,"integer"  ,(/"val"/))                          
filevardef(fout, "cal"   ,"float"  ,(/"sig","side","val"/))        
filevardef(fout, "model"    ,"float"  ,(/"model"/))  
filevardef(fout, "weight"    ,"float"  ,(/"model"/))
filevardef(fout, "time"    ,"float"  ,(/"time"/))
filevardef(fout, "f"    ,"float"  ,(/"time","model"/))
filevardef(fout, "y"    ,"float"  ,(/"time"/))

;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
setfileoption(fout,"DefineMode",False)

;===================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;====================================================================
fout->cal   = (/cal/)     
fout->sig    = (/sig/)
fout->side    = (/1,2/)
fout->model = (/model/)
fout->time = (/time/)
fout->val    = (/1,2,3/) 
fout->weight = (/weight/)
fout->f = (/f/)
fout->y(0:lrt-1) = (/y/)
end
