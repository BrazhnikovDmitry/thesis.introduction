from PyNGL import Ngl
import scipy
import scipy.io
from numpy import *
import sys
import os

path = "/home/dmitry/Work/datasets/TSUNAMI/topo30.grd"
nfile = scipy.io.netcdf_file(path)

i1 = 20400 
i2 = 20760
j1 = 14810
j2 = 15550

idx = range(i1, i2, 1)
jdx = range(j1, j2, 1)

lat = nfile.variables['lat']
lat = lat.data[jdx]
lon = nfile.variables['lon']
lon = lon.data[idx]
H = nfile.variables['z'][j1:j2, i1:i2]
H = -H[::1, ::1]
#H[where (H < 0) ] = 0

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
cnres.pmLabelBarSide		= "Right"
#cnres.lbLabelBarOn		= True

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
cnres.mpMinLatF             = 33.75
cnres.mpMaxLatF             = 36.5
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
cnres2.cnLevelSpacingF      = 500#nice_spc
#cnres.cnMinLevelValF       = nice_min
#cnres.cnMaxLevelValF       = nice_max

lines = Ngl.contour(wks,(H)[:,:],cnres2)		## Ngl.add_cyclic
#Ngl.overlay(map_plot,line_contour_plot)
Ngl.overlay(map_plot,lines)

#srlist = Ngl.Resources()
#srlist.tiMainString = title
#Ngl.set_values(map_plot,srlist)

Ngl.maximize_plot(wks,map_plot)    # Maximize size of plot in frame.
Ngl.draw(map_plot)
Ngl.frame(wks)
Ngl.destroy(wks)
