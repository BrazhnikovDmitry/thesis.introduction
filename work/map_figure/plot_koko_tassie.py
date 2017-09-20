from PyNGL import Ngl
import scipy
import scipy.io
from numpy import *
import sys
import os

path = "/home/dmitry/Work/datasets/TSUNAMI/ETOPO1_Ice_g_gmt4.grd"
nfile = scipy.io.netcdf_file(path)

i1 = 18601
i2 = 21601
j1 = 1801
j2 = 9001

idx = range(i1, i2, 25)
jdx = range(j1, j2, 25)

lat = nfile.variables['y']
lat = lat.data[jdx]
lon = nfile.variables['x']
lon = lon.data[idx]
H = -nfile.variables['z'][j1:j2, i1:i2]
H = H[::25, ::25]
H[where (H < 0) ] = 0

rlist            = Ngl.Resources()
rlist.wkColorMap = "WhiteBlue"#"BlueWhiteOrangeRed"#'posneg_1'
wks_type = "pdf"
rlist.wkPaperSize = "A4"
rlist.wkOrientation = "portrait"
wks = Ngl.open_wks(wks_type,"aaa",rlist)

ws_id = Ngl.get_workspace_id()
rlist = Ngl.Resources()
rlist.wsMaximumSize = 53554432
Ngl.set_values(ws_id,rlist)

cnres  = Ngl.Resources()
cnres2  = Ngl.Resources()
mapres  = Ngl.Resources()

mapres.nglDraw  = False
mapres.nglFrame = False
cnres.nglDraw  = False
cnres.nglFrame = False
cnres2.nglDraw  = False
cnres2.nglFrame  = False

	### Contour related resources ###
cnres.sfXArray        = lon[:]		#Ngl.add_cyclic
cnres.sfYArray        = lat[:]

cnres.cnFillDrawOrder       = "Predraw"
cnres.cnLineLabelsOn        = False
cnres.cnLinesOn             = False
cnres.cnFillOn             = True
cnres.cnLevelSelectionMode = "ManualLevels"
cnres.cnLevelSpacingF      = 500#nice_spc
#cnres.cnMinLevelValF       = nice_min
#cnres.cnMaxLevelValF       = nice_max
## Legend
#cnres.lbTitleString  = "Ocean depth, m"
cnres.lbTitleFontHeightF        = 0.018
cnres.lbLabelFontHeightF        = 0.012
cnres.lbTitleOffsetF            = -0.27
cnres.lbBoxMinorExtentF         = 0.15
cnres.pmLabelBarOrthogonalPosF  = -0.01
cnres.lbOrientation             = "Vertical"
cnres.pmLabelBarSide		= "Left"
cnres.lbLabelBarOn		= False

###############################################################
#line_contour_plot = Ngl.contour(wks,(H)[:,:],cnres)		## Ngl.add_cyclic
###############################################

red   = 0#245./255.
green = 0#222./255.
blue  = 0#179./255.
#	tan = Ngl.new_color(wks,red,green,blue)   # Add tan.
black = Ngl.get_named_color_index(wks, "Black")
cnres.mpProjection          = "CylindricalEquidistant"
cnres.mpFillOn              = True
cnres.mpFillColors          = [0,-1,black,-1]     # -1 is transparent, tan
cnres.mpDataBaseVersion     = "MediumRes"
cnres.mpLimitMode           = "LatLon"
cnres.mpMinLonF             = min(lon)-360
cnres.mpMaxLonF             = max(lon)-360
cnres.mpMinLatF             = -60
cnres.mpMaxLatF             = 60
#cnres.mpCenterLonF           = (LON0+LON1)/2
cnres.mpGridAndLimbOn       = False

map_plot          = Ngl.contour_map(wks,(H)[:,:], cnres)

cnres2  = Ngl.Resources()
cnres2.sfXArray        = lon[:]		#Ngl.add_cyclic
cnres2.sfYArray        = lat[:]

#cnres2.cnFillDrawOrder       = "Predraw"
cnres2.cnInfoLabelOn = False
cnres2.cnLineLabelsOn        = False
cnres2.cnLinesOn             = True
cnres2.cnFillOn             = False
cnres2.cnLevelSelectionMode = "ManualLevels"
cnres2.cnLevelSpacingF      = 2500#nice_spc
#cnres.cnMinLevelValF       = nice_min
#cnres.cnMaxLevelValF       = nice_max

lines = Ngl.contour(wks,(H)[:,:],cnres2)		## Ngl.add_cyclic

### Boxes
gsres                   = Ngl.Resources()

# Polyline resources.
gsres.gsLineColor       = 255
gsres.gsLineThicknessF  = 5.0      # thrice thickness

xx = [147.0,152.5,152.5,147.0,147.0]
yy = [-45.0,-45.0,-40.5,-40.5,-45.0]
box1 = Ngl.add_polyline(wks,map_plot,xx,yy,gsres)

xx = [170.0,173.0,173.0,170.0,170.0]
yy = [33.75,33.75,36.0,36.0,33.75]
box2 = Ngl.add_polyline(wks,map_plot,xx,yy,gsres)


#Ngl.overlay(map_plot,line_contour_plot)
Ngl.overlay(map_plot,lines)

#srlist = Ngl.Resources()
#srlist.tiMainString = title
#Ngl.set_values(map_plot,srlist)

Ngl.maximize_plot(wks,map_plot)    # Maximize size of plot in frame.
Ngl.draw(map_plot)
Ngl.frame(wks)
Ngl.destroy(wks)
