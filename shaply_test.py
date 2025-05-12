# -*- coding: utf-8 -*-
"""
Created on Mon May 12 13:38:45 2025

@author: nesseler
"""

import matplotlib.pyplot as plt
import shapely
import math
import numpy as np

from shapely import LineString
from shapely.wkt import loads

line = LineString([(0, 0), (2, 2)])

shapely.intersection(line, LineString([(1, 1), (3, 3)]))
# <LINESTRING (1 1, 2 2)>

box1 = shapely.box(0, 0, 2, 2)

box2 = shapely.box(1, 1, 3, 3)

shapely.intersection(box1, box2).normalize()
# <POLYGON ((1 1, 1 2, 2 2, 2 1, 1 1))>

box1 = shapely.box(0.1, 0.2, 2.1, 2.1)

shapely.intersection(box1, box2, grid_size=1)
# <POLYGON ((2 2, 2 1, 1 1, 1 2, 2 2))>

r = 1
n = 64
t = np.linspace(0, 2 * np.pi, n + 1)

x_circle = r * np.cos(t) + 1
y_circle = r * np.sin(t) + 1

# create the linestring of circle's perimeter
wktcode1 = "LINESTRING ("
for i,(x,y) in enumerate(zip(x_circle, y_circle)):
    if i!=len(x_circle)-1:
        wktcode1 += str(x)+" "+str(y)+", "
    else:
        wktcode1 += str(x)+" "+str(y)+")"
    #print(wktcode1)

circle_perim = loads(wktcode1)  #a geometry object

# create another geometry, for the line
wktcode2 = "LINESTRING ("
xs = range(4)
ys = np.array(range(4))*0.42
for i,(x,y) in enumerate(zip(xs, ys)):
    if i!=len(range(4))-1:
        wktcode2 += str(x)+" "+str(y)+", "
    else:
        wktcode2 += str(x)+" "+str(y)+")"
    #print(wktcode2)
    pass

line_str = loads(wktcode2)    #a geometry object

# check if they are intersect
# ixx = circle_perim.intersection(box1)


# check if they are intersect
ixx = circle_perim.intersection(line_str)

intersec = line_str.intersection(circle_perim)
print(intersec.wkt)

# visualization of the intersection
# plot circle
plt.scatter(x_circle, y_circle, s = 1)
    
# x1,y1 = 
plt.plot(*box1.exterior.xy)
plt.plot(*line_str.xy)
# plt.scatter(*intersec.xy)

# # if ixx:
# #     # plot intersection points
for p in list(intersec.geoms):
    plt.scatter(*p.xy, s = 10, marker = 'x', lw = 1, c = 'r')


plt.gca().set_aspect(1)
plt.show()