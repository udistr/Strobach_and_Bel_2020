load "Heavyside.ncl"
;load "cal_trend_lim.ncl"
load "EGA.ncl"
begin

pi=atan(1.)*4

R="NCEP"
SLA="AVG"

if R.eq."NCEP" then
  nt=152
else
  nt=142
end if
;-----------------------------------------------------------------------------
; get clibration factors
;----------------------------------------------------------------------------- 
diro="./CAL/"

c1=addfile(diro+"cal_RCP26_" + R + "_" + SLA + "_runmean_gaus_20.nc","r")
c2=addfile(diro+"cal_RCP45_" + R + "_" + SLA + "_runmean_gaus_20.nc","r")
c3=addfile(diro+"cal_RCP60_" + R + "_" + SLA + "_runmean_gaus_20.nc","r")
c4=addfile(diro+"cal_RCP85_" + R + "_" + SLA + "_runmean_gaus_20.nc","r")

print(c1->model)

f1=c1->f
f2=c2->f
f3=c3->f
f4=c4->f

w1=conform(f1,c1->weight,1)
w2=conform(f2,c2->weight,1)
w3=conform(f3,c3->weight,1)
w4=conform(f4,c4->weight,1)

cal1=c1->cal(1:,:,:)
cal2=c2->cal(1:,:,:)
cal3=c3->cal(1:,:,:)
cal4=c4->cal(1:,:,:)

y=c1->y

clm=dim_avg_n(y(57),0)

;-----------------------------------------------------------------------------
; calculate mean and STD
;----------------------------------------------------------------------------- 

p1=dim_sum_n(w1*f1,1)
p2=dim_sum_n(w2*f2,1)
p3=dim_sum_n(w3*f3,1)
p4=dim_sum_n(w4*f4,1)

p1r=dim_sum_n(c1->weight*f1(57,:),0)
p2r=dim_sum_n(c2->weight*f2(57,:),0)
p3r=dim_sum_n(c3->weight*f3(57,:),0)
p4r=dim_sum_n(c4->weight*f4(57,:),0)

sigma1=sqrt(dim_sum_n(w1*(f1-conform(f1,dim_sum_n(w1*f1,1),0))^2,1))
sigma2=sqrt(dim_sum_n(w2*(f2-conform(f2,dim_sum_n(w2*f2,1),0))^2,1))
sigma3=sqrt(dim_sum_n(w3*(f3-conform(f3,dim_sum_n(w3*f3,1),0))^2,1))
sigma4=sqrt(dim_sum_n(w4*(f4-conform(f4,dim_sum_n(w4*f4,1),0))^2,1))

;-----------------------------------------------------------------------------
; prepare for plot
;----------------------------------------------------------------------------- 

var1=new((/19,nt/),float)
var1(0:8,:)=conform_dims((/9,nt/),p1,1)- \
           (conform_dims((/9,nt/),cal1(:,1,0),0)+ \
           tan(conform_dims((/9,nt/),cal1(:,1,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
           conform_dims((/9,nt/),sigma1,1)
var1(9,:)=p1
var1(10:18,:)=conform_dims((/9,nt/),p1,1)+ \
             (conform_dims((/9,nt/),cal1(::-1,0,0),0)+ \
             tan(conform_dims((/9,nt/),cal1(::-1,0,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
             conform_dims((/9,nt/),sigma1,1)

var2=new((/19,nt/),float)
var2(0:8,:)=conform_dims((/9,nt/),p2,1)- \
           (conform_dims((/9,nt/),cal2(:,1,0),0)+ \
           tan(conform_dims((/9,nt/),cal2(:,1,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
           conform_dims((/9,nt/),sigma2,1)
var2(9,:)=p2
var2(10:18,:)=conform_dims((/9,nt/),p2,1)+ \
             (conform_dims((/9,nt/),cal2(::-1,0,0),0)+ \
             tan(conform_dims((/9,nt/),cal2(::-1,0,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
             conform_dims((/9,nt/),sigma2,1)

var3=new((/19,nt/),float)
var3(0:8,:)=conform_dims((/9,nt/),p3,1)- \
           (conform_dims((/9,nt/),cal3(:,1,0),0)+ \
           tan(conform_dims((/9,nt/),cal3(:,1,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
           conform_dims((/9,nt/),sigma3,1)
var3(9,:)=p3
var3(10:18,:)=conform_dims((/9,nt/),p3,1)+ \
             (conform_dims((/9,nt/),cal3(::-1,0,0),0)+ \
             tan(conform_dims((/9,nt/),cal3(::-1,0,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
             conform_dims((/9,nt/),sigma3,1)


var4=new((/19,nt/),float)
var4(0:8,:)=conform_dims((/9,nt/),p4,1)- \
           (conform_dims((/9,nt/),cal4(:,1,0),0)+ \
           tan(conform_dims((/9,nt/),cal4(:,1,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
           conform_dims((/9,nt/),sigma4,1)
var4(9,:)=p4
var4(10:18,:)=conform_dims((/9,nt/),p4,1)+ \
             (conform_dims((/9,nt/),cal4(::-1,0,0),0)+ \
             tan(conform_dims((/9,nt/),cal4(::-1,0,1),0)*pi/180)*conform_dims((/9,nt/),ispan(0,nt-1,1),1))* \
             conform_dims((/9,nt/),sigma4,1)


if R.eq."NCEP" then
  syear=1948
else
  syear=1958
end if

byear=2046;
eyear=2065;
print(avg((var1(0,(byear-syear):(eyear-syear))-clm)))
print(avg((var2(0,(byear-syear):(eyear-syear))-clm)))
print(avg((var3(0,(byear-syear):(eyear-syear))-clm)))
print(avg((var4(0,(byear-syear):(eyear-syear))-clm)))

print(avg((var1(9,(byear-syear):(eyear-syear))-clm)))
print(avg((var2(9,(byear-syear):(eyear-syear))-clm)))
print(avg((var3(9,(byear-syear):(eyear-syear))-clm)))
print(avg((var4(9,(byear-syear):(eyear-syear))-clm)))

print(avg((var1(18,(byear-syear):(eyear-syear))-clm)))
print(avg((var2(18,(byear-syear):(eyear-syear))-clm)))
print(avg((var3(18,(byear-syear):(eyear-syear))-clm)))
print(avg((var4(18,(byear-syear):(eyear-syear))-clm)))

byear=2080;
eyear=2099;
print(avg((var1(0,(byear-syear):(eyear-syear))-clm)))
print(avg((var2(0,(byear-syear):(eyear-syear))-clm)))
print(avg((var3(0,(byear-syear):(eyear-syear))-clm)))
print(avg((var4(0,(byear-syear):(eyear-syear))-clm)))

print(avg((var1(9,(byear-syear):(eyear-syear))-clm)))
print(avg((var2(9,(byear-syear):(eyear-syear))-clm)))
print(avg((var3(9,(byear-syear):(eyear-syear))-clm)))
print(avg((var4(9,(byear-syear):(eyear-syear))-clm)))

print(avg((var1(18,(byear-syear):(eyear-syear))-clm)))
print(avg((var2(18,(byear-syear):(eyear-syear))-clm)))
print(avg((var3(18,(byear-syear):(eyear-syear))-clm)))
print(avg((var4(18,(byear-syear):(eyear-syear))-clm)))

; Colors for area between curves
below_colors = (/"skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue1","slateblue2","slateblue3","slateblue4" \ 
                ,"slateblue4","slateblue3","slateblue2","slateblue1","skyblue4","skyblue3","skyblue2","skyblue1","skyblue"/)
; Colors for plot lines
line_colors  = (/"skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue1","slateblue2","slateblue3","slateblue4","black" \ 
                ,"slateblue4","slateblue3","slateblue2","slateblue1","skyblue4","skyblue3","skyblue2","skyblue1","skyblue"/)
line_colors  = (/"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","black" \ 
 ,"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent"/)
; Line thicknesses
line_thick = (/1,1,1,1,1,1,1,1,1,0.2,1,1,1,1,1,1,1,1,1/)

; ###########################################################
above_colors = new(dimsizes(below_colors),string)       
above_colors = "transparent"

; ###########################################################
; Plot data
; ###########################################################
; Open file
plot=new(4,graphic)
wks = gsn_open_wks("pdf","PLOTS/plot_ts_trend_lim_rcp_" + R + "_" + SLA + "_runmean_gaus_20")	; Open an pdf workstation

res                      = True
res@gsnDraw              = False                ; Don't draw plot 
res@gsnFrame             = False                ; Don't advance frame
res@xyMonoDashPattern    = True                 ; Solid lines for all curves
res@xyLineThicknesses    = line_thick           ; Plot line thickness
res@xyLineColors         = line_colors          ; Plot line colors
res@gsnXYFillColors = below_colors         ; Plot line color fill
res@gsnXYAboveFillColors = above_colors         ; A fix for V6.2.0. Set to transparent.
;res@tmXBMinorOn          = False                ; Turn off x-axis (bottom) minor ticks
res@tmXTMinorOn          = False                ; Turn off x-axis (bottom)minor ticks
res@tmXTOn               = False                ; Turn off x-axis (top) major ticks
;res@tmYRMinorOn          = False                ; Turn off y-axis (right) minor ticks
;res@tmYLMinorOn          = False                ; Turn off y-axis (left) minor ticks
;res@tmYROn               = False                ; Turn off y-axis (right) major ticks
res@tiYAxisString        = "Global surface warming [DegC]"        ; Y-axis title
res@tiYAxisFontHeightF   = 0.0175               ; Y-axis string size (X/0.025 for scaling)
res@tiMainFontHeightF    = 0.0175               ; Title string size (X/0.025 for scaling)
res@tmXBLabelFontHeightF = 0.015                ; Label font height (X/0.020 for scaling)
res@tmYLLabelFontHeightF = 0.015                ; Label font height (X/0.020 for scaling)

if R.eq."NCEP" then
  res@trXMinF = 1967
else
  res@trXMinF = 1957
end if
res@trXMaxF = 2100
res@trYMinF = -1.5
res@trYMaxF = 5.5

res@tmYLMode        = "Explicit"	
res@tmYLValues      = ispan(-2,6,1)
res@tmYLLabels      = "" + res@tmYLValues 
res@tmYLMinorValues = fspan(-2,6,17)
;res@gsnLeftString   = "tmXBMode = Explicit"

if R.eq."NCEP" then
  time=ispan(1948,2099,1)
else
  time=ispan(1958,2099,1)
end if
res@tiMainString = "RCP2.6"
plot(0) = gsn_csm_xy(wks,time,var1-clm,res)
res@tiMainString = "RCP4.5"
plot(1) = gsn_csm_xy(wks,time,var2-clm,res)
res@tiMainString = "RCP6.0"
plot(2) = gsn_csm_xy(wks,time,var3-clm,res)
res@tiMainString = "RCP8.5"
plot(3) = gsn_csm_xy(wks,time,var4-clm,res)

; Create legend
lgres                    = True
lgres@lgLineColors       = (/"skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue1","slateblue2","slateblue3","slateblue4","black"/)
lgres@lgLineThicknessF   = (/5,5,5,5,5,5,5,5,5,5/) 
lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
lgres@vpWidthF           = 0.15           ; width of legend (NDC)
lgres@vpHeightF          = 0.3            ; height of legend (NDC)
lgres@lgMonoDashIndex    = True   
lgres@lgPerimOn          = False          ; draw the box perimeter in orange
labels = " " + (/"0.9","0.8","0.7","0.6","0.5","0.4","0.3","0.2","0.1","Mean"/)
lbid = gsn_create_legend(wks,10,labels,lgres)

; Add legend to bottom right corner of last plot
amres                  = True
amres@amJust           = "TopLeft"  ; Corner of plot for positioning
amres@amOrthogonalPosF =  -0.48   ;  0.5 is the bottom edge of the plot.
amres@amParallelPosF   =  -0.47   ;  0.5 is the right edge of the plot.
annoid = gsn_add_annotation(plot(3), lbid, amres)


;----------------------------------------------------------------------------------------
; create panel
;----------------------------------------------------------------------------------------

pres=True
;pres@gsnFrame                     = False
pres@gsnPanelBottom               = 0.16
;pres@gsnPanelTop                  = 0.98
;pres@gsnPanelLabelBar = True
;pres@lbTitleString    = units
;pres@lbTitleFontHeightF= .015
;pres@pmLabelBarOrthogonalPosF = -0.03
;pres@gsnPanelRowSpec = True
gsn_panel(wks,plot,(/2,2/),pres)             ; now draw as one plot

















end
