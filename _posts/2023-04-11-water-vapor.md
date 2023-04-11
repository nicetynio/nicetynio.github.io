---
layout: post
title: 水汽收支方程
categories: ncl
tags: 水汽
author: nicetynio
---

* content
{:toc}

##ncl

ncl官网提供了两种计算水汽收支的方法，现介绍第二种

```
;*************************************************
; mfc_div_2.ncl
;
; Similar to mfc_div_1.ncl except a different approach is used.
;
; Rather tha using DIV(U*Q) directly, this script expands this into two separate components.
;
; Concepts illustrated:
;   MFC = Moisture Flux Convergence
;
;   MFC_advect = -(u*(dq/dx)+v*(dq/dy) )    ; advection term
;   MFC_conv   = -q*((du/dx)+(dv/dy) )      ; con(div)-vergence
;
;   MFC = MFC_advect + MFC_convection         
;
;   - Plot a number of quantities
;*************************************************
;---Calculate the Horizontal Moisture Flux Convergence [MFC]
;*************************************************
;---High frequency source data: hourly/3hr/6hr/12hr/daily .... NOT monthly values
;---References:
;---http://www.cgd.ucar.edu/cas/catalog/newbudgets/
;---http://tornado.sfsu.edu/geosciences/classes/e260/AtmosphericRivers/Moisture%20Flux.pdf
;---https://www.spc.noaa.gov/publications/banacos/mfc-sls.pdf
;===================================================================
;   Data Source: ESRL Physical Sciences Division
;        https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;   NCEP Reanalysis data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, 
;   from their Web site at https://www.esrl.noaa.gov/psd/
;===================================================================
```

```
  ptop = 300              ; 'shum' upper level                  
  ptop@units = "hPa"
  g    = 9.80665          ; m/s2

  date  = 20080715        ; NH summer

;---ESRL: CDC data

  diri = "./"
  filq = "shum.2008.nc"   ; daily data for current year [366 days]
  filu = "uwnd.2008.nc"
  filv = "vwnd.2008.nc"
  filps= "pres.sfc.2008.nc"

  pthu = diri+filu      
  pthv = diri+filv
  pthq = diri+filq     
  pthps= diri+filps   

  fu   = addfile(pthu ,"r")
  fv   = addfile(pthv ,"r")
  fq   = addfile(pthq ,"r")
  fps  = addfile(pthps,"r")

;---Time

  ymd  = cd_calendar(fu->time, -2)    ; ymd[*]: human readable
  nt   = ind(ymd.eq.date)             ; date for plotting and testing
 
TEST = True        
if (.not.TEST) then                   ; all times                   
  u    = fu->uwnd(:,{1000:ptop},:,:)  ; m/s, (time,level,lat,lon)
  v    = fv->vwnd(:,{1000:ptop},:,:)
  q    = fq->shum                     ; [kg/kg], 1000-300 levels only
  ps   = fps->pres                    ; Pa=>[kg/(m-s2)], (time,lat,lon)

else                                  ; one time step; keep time dimension [ nt:nt: ]                     
  u    = fu->uwnd(nt:nt,{1000:ptop},:,:); m/s, (time,level,lat,lon)
  v    = fv->vwnd(nt:nt,{1000:ptop},:,:)
  q    = fq->shum(nt:nt,:,:,:)        ; [kg/kg], 1000-300 levels only
  ps   = fps->pres(nt:nt,:,:)         ; Pa=>[kg/(m-s2)], (time,lat,lon)
  nt   = 0                            ; only one time step
end if

;---Vertical levels

  ptop = ptop*100
  ptop@units = "Pa"

  plev = q&level                      ; hPa  
  plev = plev*100                     ; [100000,...,30000] Pa [kg/(m-s2)]
  plev@units = "Pa"

;---Change [kg/kg] to [g/kg]; not necessary: but common units for q

  q    = q*1000            
  q@units = "g/kg"

;---Divergence function [used later] requires S->N grid order

  u    = u(:,:,::-1,:)  
  v    = v(:,:,::-1,:)
  q    = q(:,:,::-1,:)     
  ps   =ps(:,  ::-1,:)       
```
经过处理后各变量单位

u  m/s

v  m/s

q  g/kg

ps pa => N/m2 => kg m/s2/m2 => kg/(m s2)

```
;---Layer thickness: ; Pa=>[kg/(m-s2)], (time,level,lat,lon) 
;---Mass weighting:  (dp/g) => [Pa/(m/s2)] => (Pa-s2)/m => [kg/(m-s2)][s2/m] =>  (kg/m2)
;---Reference: http://www.cgd.ucar.edu/cas/catalog/newbudgets/
```
dpres_plevel_Wrap:  Calculates the pressure layer thicknesses of a constant pressure level coordinate system.

计算气压层的厚度

```
  dp   = dpres_plevel_Wrap(plev, ps, ptop, 0) ; Pa; layar thickness 

  dpg  = dp/g    
  dpg@long_name = "Layer Mass Weighting"
  dpg@units     = "kg/m2"      ; dp/g, Pa/(m s-2), reduce to kg m-2
```

g    = 9.80665          ; m/s2

dpg的单位是   kg/(m s2) s2/m =>  kg/m2

```
;************************************************
; Calculate the MFC_advection term
;   MFC_advect = -(u*(dq/dx)+v*(dq/dy) ) 
; Internally, gradients are calculated via spherical harmonics
;*************************************************    


  long_name = "MFC_advection"
  units     = "g/(kg-s)"       ; (m/s)*(g/kg)*(1/m) => (m/s)*(g/kg-m) => g/(kg-s)

  gridType  = 1   ; global fixed grid ordered S->N
  opt_adv   = 0   ; return only the advected variable; no gradients
  mfc_adv   = advect_variable(u,v,q, gridType, long_name, units, opt_adv)
  mfc_adv   = -mfc_adv

  printVarSummary(mfc_adv)  
  printMinMax(mfc_adv, 0)  
  print("--------")
```
mfc_adv还没有进行整层积分

```
;************************************************
; Calculate the MFC_convergence term
;   MFC_conv   = -q*((du/dx)+(dv/dy) )      ; con(div)-vergence
;*************************************************    

  duv  = uv2dvF_Wrap(u, v)        ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)

  mfc_con   = -q*duv           
  mfc_con@long_name = "MFC_convergence"
  mfc_con@units     = "g/(kg-s)"  ; (g/kg)(1/s) => g/(kg-s)
  copy_VarCoords(duv,mfc_con)
  delete(duv)
```
mfc_con同样还没有进行整层积分

```
;************************************************
; Calculate the total MFC
;*************************************************    

  mfc = mfc_adv + mfc_con
  mfc@long_name = "Moisture Flux Convergence"
  mfc@units     = "g/(kg-s)"  ; (g/kg)(1/s) => g/(kg-s)
  

PRINT_RAW = True
if (PRINT_RAW) then
  printVarSummary(mfc_adv)                          ; (time,level,lat,lon)
  printMinMax(mfc_adv,0)
  printVarSummary(mfc_con)                          
  printMinMax(mfc_con,0)
  print("-----")
  printVarSummary(mfc)                          
  printMinMax(mfc,0)
  print("-----")

  printVarSummary(ps)                         ; (time,lat,lon); Pa => kg/(m-s2) 
  printMinMax(ps,0)
  print("-----")
  printVarSummary(dp)                         ; (time,level,lat,lon); Pa => kg/(m-s2)
  printMinMax(dp,0)
  print("-----")
                                ; examine layer thickness at selected locations
  print(dp(nt,:,{40},{180}))    ; mid-Pacific
  print(dp(nt,:,{40},{255}))    ; Boulder, CO 
  print("-----")
end if
```

下面的单位换算涉及积分

dim_sum_n: Computes the arithmetic sum of a variable's given dimension(s) at all other dimensions.

1代表 第二维  lev    也就是一个气压层高的加权和

```
;---Integrated mass weighted moisture flux components

  mfc_adv_dpg = mfc_adv*dpg                ; mass weighted 'uq'; [m/s][g/kg][kg/m2]=>[m/s][g/kg]
  imfc_adv    =  dim_sum_n(mfc_adv_dpg, 1)
  imfc_adv@long_name = "Integrated Mass Flux Advection" 
  imfc_adv@LONG_NAME = "Sum: Mass Weighted Integrated Mass Flux Advection: mfc_adv*dpg" 
  imfc_adv@units     = "[m/s][g/kg]"
  copy_VarCoords(u(:,0,:,:), imfc_adv); (time,lat,lon)
  delete(mfc_adv_dpg)

  mfc_con_dpg = mfc_con*dpg                ; mass weighted 'mfc_con'; [m/s][g/kg][kg/m2]=>[m/s][g/kg] 
  imfc_con    = dim_sum_n(mfc_con_dpg, 1)
  imfc_con@long_name = "Integrated  Mass Flux Convergence"
  imfc_con@LONG_NAME = "Sum: Mass Weighted Integrated Mass Flux Convergence [mfc_con*dpg]" 
  imfc_con@units     = "[m/s][g/kg]"
  copy_VarCoords(v(:,0,:,:), imfc_con); (time,lat,lon)
  delete(mfc_con_dpg)

  VIMFC =  imfc_adv +  imfc_con    
  VIMFC@long_name = "VIMFC"
  VIMFC@LONG_NAME = "VIMFC: [imfc_adv+imfc_con]"
  copy_VarCoords(q(:,0,:,:),VIMFC)            ; (time,lat,lon)

PRINT_RESULT = True
if (PRINT_RESULT) then
    printVarSummary(imfc_adv)                 ; (time,lat,lon)
    printMinMax(imfc_adv,0)
    print("-----")
    printVarSummary(imfc_con)                 ; (time,lat,lon)
    printMinMax(imfc_con,0)
    print("-----")
    printVarSummary(VIMFC)               ; (time,lat,lon)
    printMinMax(VIMFC,0)
    print("-----")
end if

;*************************************************
; plot results
;*************************************************    

  scl5  = 1e5                                  ; arbitrary: used for nicer plot values
  sclab5= "(10~S~-5~N~)"                       ; used later   
  SCLAB5= "(10~S~5~N~)"           

  scl6  = 1e6  
  sclab6= "(10~S~-6~N~)"         
  SCLAB6= "(10~S~6~N~)"         

  plot := new(2,graphic)

  wks   = gsn_open_wks("png","mfc_div_2")        ; send graphics to PNG file
        
;--- mfc_adv and mfc_con at a specified pressure level

  res                   = True             ; plot mods desired
  res@gsnDraw           = False            ; don't draw yet
  res@gsnFrame          = False            ; don't advance frame yet

  res@cnFillOn          = True             ; turn on color
  res@cnLinesOn         = False            ; turn off contour lines
  res@cnLineLabelsOn    = False            ; turn off contour lines
  res@cnFillPalette     = "ViBlGrWhYeOrRe" ; set White-in-Middle color map
  res@mpFillOn          = False            ; turn off map fill
  res@lbLabelBarOn      = False            ; turn off individual cb's
                                           ; Use a common scale
  res@cnLevelSelectionMode = "ManualLevels"; manual set levels so lb consistent
  res@cnMaxLevelValF       =   12.0        ; max level
  res@cnMinLevelValF       = -res@cnMaxLevelValF     ; min level
  res@cnLevelSpacingF      =    0.5        ; contour interval

  LEVP    = 700
  res@gsnCenterString      = LEVP+"hPa"

  MFC_ADV = mfc_adv(nt,{LEVP},:,:)         ; keep meta data
  MFC_ADV = MFC_ADV*scl5
  res@gsnRightString  = sclab5+" "+mfc_adv@units
  plot(0) = gsn_csm_contour_map(wks,MFC_ADV,res)

  MFC_CON = mfc_con(nt,{LEVP},:,:)
  MFC_CON = MFC_CON*scl5
  res@gsnRightString  = sclab5+" "+mfc_con@units
  plot(1) = gsn_csm_contour_map(wks,MFC_CON,res)

  resP                     = True                ; modify the panel plot
  resP@gsnPanelMainString  = date+": Unweighted MFC_ADV, MFC_CON"
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot

;--- Integrated Moisture Transport [iuq, ivq]

  delete([/res@gsnRightString, res@gsnCenterString/]) ; not used for this plot
    
  res@cnMaxLevelValF       =  0.50                   ; min level
  res@cnMinLevelValF       = -res@cnMaxLevelValF     ; min level
  res@cnLevelSpacingF      =  0.05                   ; contour interval

  IMFC_ADV = imfc_adv(nt,:,:)                    ; local array: keep meta data
  plot(0)  = gsn_csm_contour_map(wks,IMFC_ADV,res)

  IMFC_CON = imfc_con(nt,:,:)                    ; local array: keep meta data
  plot(1) = gsn_csm_contour_map(wks,IMFC_CON,res)

  resP@gsnPanelMainString  = date+": Integrated Moisture Flux: Advect, Convergence"
  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot

  delete( [/IMFC_ADV, IMFC_CON/] )                   ; no longer needed

  res@lbLabelBarOn      = True   
  res@gsnDraw = True
  res@gsnFrame= True

;---Integrated Divergence of Moisture Flux Convergence [no scaling]

 ;res@cnFillPalette        = "cmp_flux"
  res@cnLevelSelectionMode = "ManualLevels"; manual set levels so lb consistent
  res@cnMaxLevelValF       =  0.50                ; min level
  res@cnMinLevelValF       = -res@cnMaxLevelValF  ; min level
  res@cnLevelSpacingF      =  0.050               ; contour interval
  res@tiMainString         = date+": VIMFC: [IMFC_ADV+IMFC_CON]"

  plt = gsn_csm_contour_map(wks,VIMFC(nt,:,:) ,res)
```

shapefile，是美国环境系统研究所公司（ESRI）开发的空间数据开放格式，用于描述几何体对象：点、折线与多边形。可以存储井、河流、湖泊等空间对象的几何位置。除了几何位置，shp文件也可以存储这些空间对象的属性，例如河流的名字、城市的温度等等。

Shapefile文件指的是一种文件存储的方法，实际上该种文件格式是由多个文件组成的。其中三个文件必不可少（".shp", ".shx"与 ".dbf"）。表示同一数据的一组文件其文件名前缀应该相同。其中“真正”的Shapefile的后缀为shp，然而仅有这个文件数据是不完整的，必须要把其他两个附带上才能构成一组完整的地理数据。所有文件必须位于同一目录之中。

- .shp — **图形格式**，用于保存元素的几何实体。  
- .shx — **图形索引格式**。几何体位置索引，记录每一个几何体在shp文件之中的位置，能够加快向前或向后搜索一个几何体的效率。  
- .dbf — **属性数据格式**，以dBase III+ 的数据表格式存储每个几何形状的属性数据。  

## ncl
ncl可以用addfile读取shp文件，示例如下
```
	river                          = True
    river@gsLineThicknessF         = mp_thick+1.0
    river@gsLineColor              = "black"
	plotrv = gsn_add_shapefile_polylines(wks,plot,"/home/river.nc",river)
```

```
load "./shapefile_utils.ncl"
begin
  print_shapefile_info("gadm36_CHN_2.shp") #
  plot_shapefile("gadm36_CHN_2.shp") #

  f = addfile("gadm36_CHN_2.shp", "r")   ; Open shapefile
  segments = f->segments
  geometry = f->geometry
  segsDims = dimsizes(segments) ;2D
  geomDims = dimsizes(geometry) ;2D

  geom_segIndex = f@geom_segIndex ;=0
  geom_numSegs  = f@geom_numSegs  ;=1
  segs_xyzIndex = f@segs_xyzIndex ;=0
  segs_numPnts  = f@segs_numPnts  ;=1

  lines       = new(segsDims(0),graphic)   ; Array to hold polygons
  numFeatures = geomDims(0) ;=dimsizes(NAME_0)

  lon   = f->x
  lat   = f->y
  name  = f->NAME_2

  segNum = 0
  do i=0, numFeatures-1  
     startSegment = geometry(i, geom_segIndex)
     numSegments  = geometry(i, geom_numSegs)
     print(name(i)+" startSegment:"+startSegment+" numSegments:"+numSegments)
     do seg=startSegment, startSegment+numSegments-1
        startPT = segments(seg, segs_xyzIndex)
        endPT = startPT + segments(seg, segs_numPnts) - 1
        print("startPT:"+startPT+" endPT:"+endPT) 
        ;+" lon(startPT:endPT):"+lon(startPT:endPT)+"lat(startPT:endPT):"+lat(startPT:endPT))
        segNum = segNum + 1
     end do
  end do
end

