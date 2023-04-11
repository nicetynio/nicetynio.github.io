---
layout: post
title: 根据行政区划mask数据
categories: python ncl
tags: mask
author: renql
---

* content
{:toc}

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
```
## python
有一个叫regionmask的程序包，可以实现该功能。 其函数介绍页面 https://regionmask.readthedocs.io/en/stable/api.html
```py
import cmaps, regionmask
import geopandas as gpd

ds = xr.open_dataset('~/data/t2m.nc')
var = ds['t2m'] # 三维数组

shp_file = '~/data/gadm36_CHN_0.shp'
shpf = gpd.read_file(shp_file) 
mask = regionmask.mask_geopandas(shpf, ds)
var = var.where(mask==0.0).mean(('lat','lon')) # 虽然var和mask的维数不同，但不影响
# ds也可以换成var，只要是包含经纬度信息的xarray数据就可以
# 且默认ds中经纬度名字为lat，lon，若不是，需要指定，此时语法是
mask = regionmask.mask_geopandas(shpf, ds, 'longitude', 'latitude')
# 也可以在shpf后直接跟lon和lat数组,且该数组可以是一维或二维，可以是xarray或array
mask = regionmask.mask_geopandas(shpf, ds['longitude'], ds['latitude'])
# 无论哪种方法，都只返回一个二维数组xarray.DataArray mask，其中不位于shp_file定义的区域内的格点为nan
# 由于此处shp_file只定义了一个区域，所以在区域中的格点为0
# 若shp_file有多个区域，则mask中的值从0-len(区域)

# regionmask中已包含了许多常见的区域划分，可直接调用，不需要使用shapefile
# natural_earth_v5_0_0.countries_110中包含177个区域，其中China的编号是139
region = regionmask.defined_regions.natural_earth_v5_0_0.countries_110
mask = region.mask(var) #输入经纬度信息的方法和返回数组信息同regionmask.mask_geopandas
var = var.where(mask==139)
```

## 参考
- https://zh.m.wikipedia.org/zh-hans/Shapefile  
- https://www.ncl.ucar.edu/Applications/shapefiles.shtml
- https://regionmask.readthedocs.io/en/stable/notebooks/geopandas.html
