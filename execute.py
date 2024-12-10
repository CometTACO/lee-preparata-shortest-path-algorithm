from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import geometry. __init__
import draw.__init__
import gsp.lee_preparata #shortest path algorithms


size = 10
digits = 3
font_size = 8
ticks = 10
timeout = .5
dpi = 300

points = [(-5.5, 0.96), (-4.5, 2.46), (-3.56, 1.44), (-2.5, 3.26), (-1.92, 1.86), (-0.48, 4.32), (0.28, 2.64), (1.38, 4.9), (2.64, 2.9), (3.8, 1.9), (4.86, 0.98), (5, 0), (5.16, -1.16), (4.32, -2.38), (3, -2), (2.76, -0.86), (1.98, -2.52), (1, -1), (0.52, -2.62), (-1.2, -2.72), (-1.1, -1.2), (-2.82, -1.56), (-2.84, -2.7), (-4.46, -2.84), (-4.38, -1.46), (-6.48, -1.62), (-7.68, -0.76), (-7.6, 0.18)]

start = geometry.Point(-4.38, 1.8)
end = geometry.Point(3.8, -1.58)

p_list = []
for point in points:
    p_list.append(geometry.Point(point[0],point[1]))
polygon = geometry.Polygon(p_list)
triangulated = geometry.TriangulatedPolygon(polygon)
shortest_path = gsp.lee_preparata.shortest_path(triangulated,start,end)

#-------Original-------#
fg, ax = plt.subplots()
draw.draw_polygon(ax,polygon)
plt.plot([start.x, end.x],[start.y, end.y],'ro')
plt.show()
#-------Triangulated-------#
fg, ax = plt.subplots()
draw.draw_polygon(ax,polygon)
draw.draw_polygon(ax,triangulated)
plt.plot([start.x, end.x],[start.y, end.y],'ro')
plt.show()
#-------Shortest Path-------#
fg, ax = plt.subplots()
draw.draw_polygon(ax,polygon)
x_coords = []
y_coords = []
for point in shortest_path:
    x_coords.append(point.x)
    y_coords.append(point.y)
plt.plot(x_coords,y_coords,'ro',linestyle = 'solid')
plt.show()
