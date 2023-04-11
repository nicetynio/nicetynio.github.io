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

现在介绍另一种，是在气象家园找到的

```
;-------------------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
;-------------------------------------------------------------------------------------------------
begin
;-------------------------------------------------------------------------------------------------
g=9.80665  ;m/s2
;-------------------------------------------------------------------------------------------------
```

```
file_path="/mnt/d/JRA55_grib/processed/"
prefix = "anl_p125."
suffix = ".nc"
start_year=1958
end_year=2021
;-------------------------------------------------------------------------------------------------

do year=start_year,end_year

u_file_name = prefix +"033_ugrd."+ year+"01_" +year+"12"+ suffix
v_file_name = prefix +"034_vgrd."+ year+"01_" +year+"12"+ suffix
q_file_name = prefix +"051_spfh."+ year+"01_" +year+"12"+ suffix
ps_file_name ="anl_surf125.001_pres."+ year+"01_" +year+"12"+ suffix

f_u=addfile(file_path+"u/"+u_file_name,"r")
f_v=addfile(file_path+"v/"+v_file_name,"r")
f_q=addfile(file_path+"q/"+q_file_name,"r")
f_p=addfile(file_path+"ps/"+ps_file_name,"r")

ua=f_u->var33(:,::-1,1:143,:)
va=f_v->var34(:,::-1,1:143,:)
q=f_q->var51(:,::-1,1:143,:)
p=f_p->var1(:,1:143,:)

p@units = "Pa"

level=f_u->plev(::-1)
lat=f_u->lat(1:143)
lon=f_u->lon(:)
time=f_u->time(:)

level@positive="up"

q=q*1000.0;
q@units = "g/kg"
ndim=dimsizes(q)
```
ps pa => N/m2 => kg m/s2/m2 => kg/(m s2)
advect_variable_cfd:  adv_X =  U*(dX/dlon) + V*(dX/dlat)


m/s  g/kg  1/m  

```
long_name = "Horizontal_MFC_advection_term"
units     = "g/(kg-s)"
cyclic    = True
opt_adv   = 0
mfc_adv=advect_variable_cfd(ua,va,q,lat,lon, cyclic, long_name, units, opt_adv)


printVarSummary(level)
printVarSummary(mfc_adv)
printVarSummary(p)
```

下面是求整层积分

就是mfc_adv的量纲 乘以 dp/g 的量纲

ps pa => N/m2 => kg m/s2/m2 => kg/(m s2)
g=9.80665  ;m/s2


g/(kg-s)   kg/m2
adv  "g/(m2-s)"

要得到进一步的结果    adv要除以$/rho$



```
adv=vibeta(level,mfc_adv(time|:,lat|:,lon|:,plev|:),1,p,100000.0,10000.0)/g

adv@units     = "g/(m2-s)"
adv          := tofloat(24.*3600.*adv/1000.) ;double
adv@units     = "mm/day"
adv@long_name = "Horizontal_MFC_advection_term"
copy_VarCoords(ua(:,0,:,:), adv)

ncdf1 = addfile("/mnt/d/JRA55_grib/output/HMA/HMA_"+year+".nc","c")
setfileoption(ncdf1, "prefill", False)
setfileoption(ncdf1, "suppressclose", True)
setfileoption(ncdf1, "definemode", True)

dimNames = (/"time", "lat", "lon" /)
dimSizes = (/12,143,288/)         
dimUnlim = (/True,False,False/)      
filedimdef(ncdf1, dimNames,dimSizes,dimUnlim )
filevardef(ncdf1, "lat", typeof(lat), getvardims(lat))
filevarattdef(ncdf1, "lat", lat)
filevardef(ncdf1, "lon", typeof(lon), getvardims(lon) )
filevarattdef(ncdf1, "lon", lon)      
filevardef(ncdf1, "time", typeof(time), getvardims(time) )
filevarattdef(ncdf1, "time", time) 
filevardef(ncdf1,"adv",typeof(adv), getvardims(adv) ) 
filevarattdef(ncdf1,"adv", adv) 
ncdf1->lat    = (/ lat /)                                              
ncdf1->lon    = (/ lon /)
ncdf1->time   = (/ time /)                                              
ncdf1->adv   =  adv

;MFC_div   = q*((du/dx)+(dv/dy) )      ; con(div)-vergence
duv    = uv2dv_cfd(ua,va, lat, lon,3)  
mfc_div   = q*duv
mfc_div@long_name = "mfc_divergence"
mfc_div@units     = "g/(kg-s)"  
copy_VarCoords(ua,mfc_div)
delete(duv)

div             = vibeta(level,mfc_div(time|:,lat|:,lon|:,plev|:),1,p,100000.0,10000.0)/g
div@units  = "g/(m2-s)"
div             = 24.*3600.*div/1000.
div@units     = "mm/day"
div@long_name = "Horizontal_MFC_convergence_term"
copy_VarCoords(ua(:,0,:,:), div)

HMFC   = adv+div
HMFC@units = "mm/day"
HMFC@long_name = "Horizontal_MFC"
copy_VarCoords(adv,HMFC)


ncdf1 = addfile("/mnt/d/JRA55_grib/output/HMC/HMC_"+year+".nc","c")
setfileoption(ncdf1, "prefill", False)
setfileoption(ncdf1, "suppressclose", True)
setfileoption(ncdf1, "definemode", True)

dimNames = (/"time", "lat", "lon" /)
dimSizes = (/12,143,288/)
dimUnlim = (/True,False,False/)
filedimdef(ncdf1, dimNames,dimSizes,dimUnlim )
filevardef(ncdf1, "lat", typeof(lat), getvardims(lat))
filevarattdef(ncdf1, "lat", lat)
filevardef(ncdf1, "lon", typeof(lon), getvardims(lon) )
filevarattdef(ncdf1, "lon", lon)
filevardef(ncdf1, "time", typeof(time), getvardims(time) )
filevarattdef(ncdf1, "time", time)
filevardef(ncdf1,"div",typeof(div), getvardims(div))
filevarattdef(ncdf1,"div", div)
ncdf1->lat    = (/ lat /)                                              
ncdf1->lon    = (/ lon /)
ncdf1->time   = (/ time /)                                              
ncdf1->div   =  div



ncdf1 = addfile("/mnt/d/JRA55_grib/output/HMFC/HMFC_"+year+".nc","c")

setfileoption(ncdf1, "prefill", False)
setfileoption(ncdf1, "suppressclose", True)
setfileoption(ncdf1, "definemode", True)

dimNames = (/"time", "lat", "lon" /)
dimSizes = (/12,143,288/)       
dimUnlim = (/True,False,False/)
filedimdef(ncdf1, dimNames,dimSizes,dimUnlim )
filevardef(ncdf1, "lat", typeof(lat), getvardims(lat))
filevarattdef(ncdf1, "lat", lat)
filevardef(ncdf1, "lon", typeof(lon), getvardims(lon) )
filevarattdef(ncdf1, "lon", lon)      
filevardef(ncdf1, "time", typeof(time), getvardims(time) )
filevarattdef(ncdf1, "time", time) 
filevardef(ncdf1,"HMFC",typeof(HMFC), getvardims(HMFC) )
filevarattdef(ncdf1,"HMFC", HMFC)
ncdf1->lat    = (/ lat /) 
ncdf1->lon    = (/ lon /)
ncdf1->time   = (/ time /)
ncdf1->HMFC   =  HMFC


delete(u_file_name)
delete(v_file_name)
delete(q_file_name)
delete(ps_file_name)
delete(f_u)
delete(f_v)
delete(f_q)
delete(f_p)
delete(ua)
delete(va)
delete(q)
delete(p)
delete(level)
delete(lat)
delete(lon)
delete(time)
delete(ndim)
delete(mfc_adv)
delete(adv)
delete(mfc_div)
delete(div)
delete(HMFC)

print("+")
end do
print("Finally!")

end
```
