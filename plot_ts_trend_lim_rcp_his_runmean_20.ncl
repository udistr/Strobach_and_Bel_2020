load "Heavyside.ncl"
load "cal_sym.ncl"
load "EGA.ncl"
begin

R="NCEP"
SLA="EGA"

lrt=69

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

n=15
pi=atan(1)/4

;-----------------------------------------------------------------------------
; get clibration factors
;----------------------------------------------------------------------------- 

c1=addfile("CAL/cal_RCP"+RCP+"_" + R + "_" + "AVG" + "_" + n + "_his_runmean_20.nc","r")
c2=addfile("CAL/cal_RCP"+RCP+"_" + R + "_" + SLA + "_" + n + "_his_runmean_20.nc","r")

cal1=c1->cal(1:,:,:)
cal2=c2->cal(1:,:,:)

f1=c1->f
f2=c2->f

w1=conform(f1,c1->weight,1)
w2=conform(f2,c2->weight,1)

ncep=c1->y
clm=ncep(57)

;-----------------------------------------------------------------------------
; calculate mean and STD
;----------------------------------------------------------------------------- 
p1=dim_sum_n(w1*f1,1)
p2=dim_sum_n(w2*f2,1)

sigma1=sqrt(dim_sum_n(w1*(f1-conform(f1,dim_sum_n(w1*f1,1),0))^2,1))
sigma2=sqrt(dim_sum_n(w2*(f2-conform(f2,dim_sum_n(w2*f2,1),0))^2,1))

;-----------------------------------------------------------------------------
; prepare for plot
;----------------------------------------------------------------------------- 
nt=lrt

var1=new((/19,nt/),float)
var1(0:8,:)=conform_dims((/9,nt/),p1,1)+ \
           (conform_dims((/9,nt/),cdfnor_x(fspan(0.05,0.45,9),conform_dims(9,0,0),conform_dims(9,1,0)),0))* \
           conform_dims((/9,nt/),sigma1,1)
var1(9,:)=p1
var1(10:18,:)=conform_dims((/9,nt/),p1,1)+ \
             (conform_dims((/9,nt/),cdfnor_x(fspan(0.55,0.95,9),conform_dims(9,0,0),conform_dims(9,1,0)),0))* \
             conform_dims((/9,nt/),sigma1,1)
;print(cdfnor_x(fspan(0.05,0.45,9),conform_dims(9,0,0),conform_dims(9,1,0)))
;print(cdfnor_x(fspan(0.55,0.95,9),conform_dims(9,0,0),conform_dims(9,1,0)))

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

;print(p2-clm)
; Colors for area between curves
below_colors = (/"skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue1","slateblue2","slateblue3","slateblue4" \ 
                ,"slateblue4","slateblue3","slateblue2","slateblue1","skyblue4","skyblue3","skyblue2","skyblue1","skyblue"/)
; Colors for plot lines
line_colors  = (/"skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue1","slateblue2","slateblue3","slateblue4","black" \ 
                ,"slateblue4","slateblue3","slateblue2","slateblue1","skyblue4","skyblue3","skyblue2","skyblue1","skyblue","red"/)
line_colors  = (/"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","black" \ 
 ,"transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","transparent","red"/)
; Line thicknesses
line_thick = (/1,1,1,1,1,1,1,1,1,0.2,1,1,1,1,1,1,1,1,1/)

; ###########################################################
above_colors = new(dimsizes(below_colors),string)       
above_colors = "transparent"

; ###########################################################
; Plot data
; ###########################################################
; Open file

wks = gsn_open_wks("pdf","PLOTS/plot_ts_trend_lim_rcp"+RCP+"_" + R + "_" + SLA + "_" + n + "_his_runmean_20"); Open an pdf workstation

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
res@trXMaxF = 2017
res@trYMinF = -0.5
res@trYMaxF = 0.5

res@tmYLMode        = "Explicit"	
res@tmYLValues      = fspan(-0.5,0.5,11)
res@tmYLLabels      = "" + res@tmYLValues 
res@tmYLMinorValues = fspan(-2,6,17)
;res@gsnLeftString   = "tmXBMode = Explicit"

if R.eq."NCEP" then
  time=ispan(1948,1948+lrt-1,1)
else
  time=ispan(1958,2099,1)
end if
res@tiMainString = "AVG+Gaussian"
plot1 = gsn_csm_xy(wks,time,var1-clm,res)
plres=True
plres@gsLineDashPattern  = 2
dum1=gsn_add_polyline(wks, plot1, (/1948+lrt-n,1948+lrt-n/), (/-2,6/), plres)
res@tiMainString = "EGA+AR"
plot2 = gsn_csm_xy(wks,time,var2-clm,res)
dum2=gsn_add_polyline(wks, plot2, (/1948+lrt-n,1948+lrt-n/), (/-2,6/), plres)

res1=True
res1@gsnDraw              = False                ; Don't draw plot 
res1@gsnFrame             = False                ; Don't advance frame
if R.eq."NCEP" then
  res1@trXMinF = 1967
else
  res1@trXMinF = 1957
end if
res1@trXMaxF = 2007
res1@trYMinF = -0.5
res1@trYMaxF = 0.5
res1@xyLineColor      = "red"
res1@xyLineThicknessF = 2.0
plot1b=gsn_csm_xy(wks,time,ncep-clm,res1)
overlay(plot1,plot1b)
plot2b=gsn_csm_xy(wks,time,ncep-clm,res1)
overlay(plot2,plot2b)

; Create legend
lgres                    = True
lgres@lgLineColors       = (/"skyblue","skyblue1","skyblue2","skyblue3","skyblue4","slateblue1","slateblue2","slateblue3","slateblue4","black","red"/)
lgres@lgLineThicknessF   = (/5,5,5,5,5,5,5,5,5,5,5/) 
lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
lgres@vpWidthF           = 0.15           ; width of legend (NDC)
lgres@vpHeightF          = 0.3            ; height of legend (NDC)
lgres@lgMonoDashIndex    = True   
lgres@lgPerimOn          = False          ; draw the box perimeter in orange
labels = " " + (/"0.9","0.8","0.7","0.6","0.5","0.4","0.3","0.2","0.1","Mean","NCEP"/)
lbid = gsn_create_legend(wks,10,labels,lgres)

; Add legend to bottom right corner of last plot
amres                  = True
amres@amJust           = "TopLeft"  ; Corner of plot for positioning
amres@amOrthogonalPosF =  -0.48   ;  0.5 is the bottom edge of the plot.
amres@amParallelPosF   =  -0.47   ;  0.5 is the right edge of the plot.
annoid = gsn_add_annotation(plot2, lbid, amres)


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
gsn_panel(wks,(/plot1,plot2/),(/1,2/),pres)             ; now draw as one plot

;print(sum(where(dim_max_n(var2(:,0:lrt-n-1),0)-ncep(0:lrt-n-1).le.0,1,0)))
;print(sum(where(-dim_min_n(var2(:,0:lrt-n-1),0)+ncep(0:lrt-n-1).le.0,1,0)))

print("rmse of the validation pariod (n years)")
print("AVG-"+sqrt(avg((p1(lrt-n:lrt-1)-ncep(lrt-n:lrt-1)))^2))
print("EGA-"+sqrt(avg((p2(lrt-n:lrt-1)-ncep(lrt-n:lrt-1)))^2))
print("EGA/AVG-"+sqrt(avg((p2(lrt-n:lrt-1)-ncep(lrt-n:lrt-1)))^2)/sqrt(avg((p1(lrt-n:lrt-1)-ncep(lrt-n:lrt-1)))^2))

;print(sqrt(avg((p1(lrt-n:lrt-1)-ncep(lrt-n:lrt-1)))^2))
;print(sqrt(avg((p2(lrt-n:lrt-1)-ncep(lrt-n:lrt-1)))^2))

; number of times that the predctions are higher than obs
print(sum(where(p1(0:lrt-1-n)-ncep(0:lrt-1-n).ge.0,1,0)))
print(sum(where(p2(0:lrt-1-n)-ncep(0:lrt-1-n).ge.0,1,0)))


print("90 percentile range of the validation pariod (n years)")
print("AVG")
;print(dim_sum_n(heavyside(var1(:,lrt-n:lrt-1)-conform(var1(:,lrt-n:lrt-1),ncep(lrt-n:lrt-1),1)),1))
print(dim_avg_n(var1(18,lrt-n:lrt-1),0)-dim_avg_n(var1(0,lrt-n:lrt-1),0))
print("EGA")
;print(dim_sum_n(heavyside(var2(:,lrt-n:lrt-1)-conform(var2(:,lrt-n:lrt-1),ncep(lrt-n:lrt-1),1)),1))
print(dim_avg_n(var2(18,lrt-n:lrt-1),0)-dim_avg_n(var2(0,lrt-n:lrt-1),0))
print("EGA/AVG")
print(1/((dim_avg_n(var1(18,lrt-n:lrt-1),0)-dim_avg_n(var1(0,lrt-n:lrt-1),0))/(dim_avg_n(var2(18,lrt-n:lrt-1),0)-dim_avg_n(var2(0,lrt-n:lrt-1),0))))



print("AVG")
r1=(dim_sum_n(heavyside(var1(18:10,lrt-n:lrt-1)-conform(var1(18:10,lrt-n:lrt-1),ncep(lrt-n:lrt-1),1))* \
                heavyside(conform(var1(ispan(0,8,1),lrt-n:lrt-1),ncep(lrt-n:lrt-1),1)-var1(ispan(0,8,1),lrt-n:lrt-1)),1)/(n))
print(r1)
print("EGA")
r2=(dim_sum_n(heavyside(var2(18:10,lrt-n:lrt-1)-conform(var2(0:8,lrt-n:lrt-1),ncep(lrt-n:lrt-1),1))* \
                heavyside(conform(var2(ispan(18,10,1),lrt-n:lrt-1),ncep(lrt-n:lrt-1),1)-var2(ispan(0,8,1),lrt-n:lrt-1)),1)/(n))

print(r2)

x=fspan(0.9,0.1,9)

data=new((/3,9/),float)
data(0,:)=r1
data(1,:)=r2
data(2,:)=x

 wks   = gsn_open_wks ("pdf","PLOTS/relibility_rcp"+RCP+"_" + n)    ; send graphics to PNG file

;---Set plotting parameters

 resr                   = True                      ; plot mods desired
; res@tiMainString      = "Two curve XY plot"       ; add title
;
; Similiar resources are xyLineThicknessF and xyLineColor,
; which will affect all lines in the array.
;

 resr@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
 resr@xyMarkers      = (/1,1,1/)               ; 3 different markers
 resr@xyMarkerThicknesses = (/  3.0,   3.0, 3.0/)          ; make second line thicker
 resr@xyMarkerColors      = (/"blue","red","black"/)          ; change line color
 resr@xyLineThicknesses = (/  2.0,   2.0, 2.0/)          ; make second line thicker
 resr@xyLineColors      = (/"blue","red","black"/)          ; change line color
 resr@xyDashPattern  = 0                         ; Make curves all solid

 resr@tiMainString  = "Fraction of observations inside range"          ; title
 resr@tiXAxisString = "Predicted"   ; axis string
 resr@tiYAxisString = "Observed"   ; axis string


  resr@pmLegendDisplayMode    = "Always"            ; turn on legend
  resr@pmLegendSide           = "Top"               ; Change location of 
  resr@pmLegendParallelPosF   = .2                  ; move units right
  resr@pmLegendOrthogonalPosF = -0.4                ; move units down
  resr@pmLegendWidthF         = 0.15                ; Change width and
  resr@pmLegendHeightF        = 0.15                ; height of legend.
  resr@lgPerimOn              = False               ; turn off box around
  resr@lgLabelFontHeightF     = .02                ; label font height
  resr@xyExplicitLegendLabels = (/"  AVG","  EGA","  reliable"/)

 plot  = gsn_csm_xy (wks,x,data,resr) ; create plot


end
