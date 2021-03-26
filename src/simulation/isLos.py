#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 13:05:04 2020

@author: terrylines
"""

from osgeo import ogr, gdal, osr

# make errors raise exceptions.
gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()


####Some test attributes
observer_time="na"
svid="tbd"
LLA=(-0.13453,51.52465,80)
def getSvidLocation(time,svid):
    #returns WGS84 cartesian co-ordinates of the svid at the given time (time format?)
    #currently return position of G15 on 11/02/2020 at 12.30
    return (22893939.149,593414.779,13458380.448)


def isLOS(self,observer_time, svid, *LLA):
        svid_position = getSvidLocationBNG(observer_time,svid)
        if svid_position[2] < 0: 
            return False  #below horizon
        observer_position = getObserverLocationBNG(LLA)
        line = convertToLine(observer_position, svid_position)
        return line.Intersects(self.map)
    
def getSvidLocationBNG(time,svid):
    #returns SVID position (projected to a bounding box to ensure validity of Helmert transform) in British National Grid co-ordinates.
    EPSG_WGS84_CART=4978
    EPSG_BNG=27700
    svid_position = getSvidLocation(time,svid) 
    bounded_position = clip(svid_position)
    transformed_position = reproject(bounded_position,EPSG_WGS84_CART,EPSG_BNG)
    return transformed_position

def getObserverLocationBNG(position):
    EPSG_WGS84=4979
    EPSG_BNG=27700  
    return reproject(position,EPSG_WGS84,EPSG_BNG)

def convertToLine(start_position,end_position):
    line = ogr.Geometry(ogr.wkbLineString)
    line.AddPoint(*start_position)
    line.AddPoint(*end_position)
    return line
   
def isIntersecting(line,polygon):
    return   



def setMap(self):
    #currently hardcoded to UCL quad. BNG crs.
    with open("/Users/terrylines/Documents/GNSS/map/UCL.csv") as f:
        wkt=f.read()
        
    poly=ogr.CreateGeometryFromWkt(wkt)
    self.Map = poly.GetBoundary()

def clip(position):
    #clips position to a 3D bounding box of 100km surrounding London 
    ORIGIN = (3980000,-10000,4970000) #London (WGS84 cartesian co-ordinates)
    BBOX_SIDE_LENGTH = 100000
    delta = [p - q for p,q in zip(position, ORIGIN)]
    scale = BBOX_SIDE_LENGTH / max(abs(i) for i in delta)
    return tuple(p + q * scale for p,q in zip(ORIGIN, delta))

def reproject(coordinates,source,target):
    #reproject point coords given as EPSG refs
    crs_source = osr.SpatialReference()
    crs_source.ImportFromEPSG(source)
    crs_target = osr.SpatialReference()
    crs_target.ImportFromEPSG(target)
    transformer = osr.CoordinateTransformation(crs_source,crs_target)
    
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(*coordinates)
    point.Transform(transformer)
    
    return point.GetPoint()






    
    def __init__(self):
        self.position.AddPoint(529513,182286,33) #assumes BNG position

    
    