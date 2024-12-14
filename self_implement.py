"""This file presents the main part of the program. For the full source, refer to our repository: https://github.com/CometTACO/lee-preparata-shortest-path-algorithm
We expect you to already have these modules: math, typing, triangle, collections, math, matplotlib.
"""

import math
from typing import Iterable, List, Optional, Any, Tuple, Union
from triangle import triangulate
from collections import deque
from math import ceil, floor
from matplotlib import pyplot, ticker
from matplotlib.axes import Axes


#################################################################################################
#################################################################################################
####                        ######  ######  ########  ##   ##  ######                        ####
####                        ##      ##         ##     ##   ##  ##   ##                       ####
####                        ######  ######     ##     ##   ##  ######                        ####
####                            ##  ##         ##     ##   ##  ##                            ####
####                        ######  ######     ##     #######  ##                            ####
#################################################################################################
################################################################################################# 


#Point
"""Defines a point class with various helper methods."""

class Point(object):
    """Defines a point in R^2."""

    CCW_TURN = 1
    NO_TURN = 0
    CW_TURN = -1

    def __init__(self, x: Union[int, float], y: Union[int, float]):
        """Initialize a new point with x and y coordinate."""
        assert isinstance(x, (int, float))
        assert isinstance(y, (int, float))

        self.x = x
        self.y = y

    @staticmethod
    def turn(p1: 'Point', p2: 'Point', p3: 'Point') -> int:
        """Return the direction of the turn three points form.

        The return value is one of the following:
         1 = left (counterclockwise) turn,
         0 = no turn at all (points are colinear),
        -1 = right (clockwise) turn.
        """
        assert isinstance(p1, Point)
        assert isinstance(p2, Point)
        assert isinstance(p3, Point)

        res = (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)

        if res > 0:
            return 1
        if res < 0:
            return -1
        return 0
    
    def tuple(self) -> Tuple[Union[int, float], Union[int, float]]:
        """Return the coordinates as a tuple."""
        return (self.x, self.y)
    
#PolygonPoint
class PolygonPoint(Point):
    """Extend the point to include an index.

    The index defines the position in the polygon's list of points.
    """

    def __init__(self, x: Union[int, float, Point], y: int = None, index: int = None):
        """Initialize a new polygon point with additional index.

        It can be given either as a point and an index or two coordinates and an
        index.
        """
        if isinstance(x, Point):
            assert isinstance(y, int) or y is None

            self.x = x.x
            self.y = x.y
            self.index = y if y is not None else index
        else:
            assert isinstance(index, int) or index is None

            super(PolygonPoint, self).__init__(x, y)
            self.index = index


#Line Segment
"""Some polygon helper classes needed by different modules."""


class LineSegment(object):
    """Defines a line segment in R^2."""

    def __init__(self, a: Point, b: Point):
        """Initialize a new line segment given by two points."""
        assert isinstance(a, Point)
        assert isinstance(b, Point)

        self.a = a
        self.b = b

#Edge

class Edge(LineSegment):
    """Just another name for a line segment."""

    def __eq__(self, other: Any) -> bool:
        """Two edges are equal if the defining points are equal not looking at the direction."""
        if isinstance(other, Edge):
            return (self.a == other.a and self.b == other.b or
                    self.a == other.b and self.b == other.a)

        raise NotImplementedError()
    

#Polygon
"""Defines a polygon class with various helper classes and methods."""


class Polygon(object):
    """Defines a polygon in R^2."""

    def __init__(self, points: Union[Iterable[Point], 'Polygon']):
        """O(n): Initialize a new polygon with a list of points.

        The list can be anything which ist convertible to a list when passed to
        the builtin list function. It is furthermore assumed that all objects in
        this list are Point instances.

        We also assume that the points are given in counterclockwise order.
        Otherwise the names of some functions do not make sense.
        """
        if isinstance(points, Polygon):
            points = points.points

        if not isinstance(points, list):
            points = list(points)

        for ix, point in enumerate(points):
            assert isinstance(point, Point)
            points[ix] = PolygonPoint(point, ix)

        self.points = tuple(points)
        self.len = len(self.points)

    def points_as_tuples(self) -> Iterable[Tuple[Union[int, float], Union[int, float]]]:
        """O(n): Return a map object of all polygon points as tuples.

        Returns:
            [(int|float, int|float)]: An iterable of tuples (x, y) each consisting of the x- and y-coordinate of each
                point.
        """
        return map(lambda p: p.tuple(), self.points)

"""Some polygon helper classes neede by different modules."""


#Triangle with ordered points set
"""Defines a triangle class with various helper methods."""

class DelaunayTriangle(object):
    """A triangle class."""

    def __init__(self, a: PolygonPoint, b: PolygonPoint, c: PolygonPoint):
        """Create a new triangle consisting of 3 points.

        The points NEED NOT be given in any special order but they MUST NOT be collinear.
        After construction the points of the object are in counter-clockwise order.
        """
        assert isinstance(a, PolygonPoint)
        assert isinstance(b, PolygonPoint)
        assert isinstance(c, PolygonPoint)


        if Point.turn(a, b, c) == Point.CCW_TURN:
            self.a = a
            self.b = b
            self.c = c
        else:
            self.a = a
            self.b = c
            self.c = b

        # Make sure the first point is the one with the lowest lexicographic (x, y) pair
        if self.b.tuple() < self.a.tuple() and self.b.tuple() < self.c.tuple():
            self.a, self.b, self.c = self.b, self.c, self.a
        elif self.c.tuple() < self.a.tuple() and self.c.tuple() < self.b.tuple():
            self.a, self.b, self.c = self.c, self.a, self.b
        assert self.a.tuple() < self.b.tuple() and self.a.tuple() < self.c.tuple()

        self.edges = (Edge(self.a, self.b), Edge(self.b, self.c), Edge(self.c, self.a))

    def contains(self, p: Point) -> bool:
        """O(1): Return whether the triangle contains the given point.

        Args:
            p (Point): The point to check.
        """
        return (Point.turn(self.a, self.b, p) == Point.CCW_TURN and
                Point.turn(self.b, self.c, p) == Point.CCW_TURN and
                Point.turn(self.c, self.a, p) == Point.CCW_TURN)
    
    def common_edge(self, other: 'DelaunayTriangle') -> Optional[Edge]:
        """O(1): Return the common edge of two triangles. Returns None if both are equal or no common edge exists."""
        if self == other:
            return None

        for edge in self.edges:
            if edge in other.edges:
                return edge

        return None


"""Defines a polygon class with various helper classes and methods."""


class TriangulatedPolygonTriangle(DelaunayTriangle):
    """A 2D triangle with neighbours."""

    def __init__(self, a: Point, b: Point, c: Point):
        """Create a new triangle."""
        super(TriangulatedPolygonTriangle, self).__init__(a, b, c)

        self.neighbour_indices = []

class TriangulatedPolygon(Polygon):
    """Defines a polygon in R^2 together with it's triangulation."""

    def __init__(self, points: Union[Iterable[Point], Polygon]):
        """O(n^3): Initialize a new polygon with a list of points.

        The list can be anything which ist convertible to a list when passed to
        the builtin list function. It is furthermore assumed that all objects in
        this list are Point instances.

        We also assume that the points are given in counterclockwise order.
        Otherwise the names of some functions do not make sense.
        """
        super(TriangulatedPolygon, self).__init__(points)

        self._triangulate()

    def _triangulate(self):
        """O(n^3): Triangulate the polygon."""
        data = dict(
            vertices=list(self.points_as_tuples()),
            segments=[(i, (i + 1) % self.len) for i in range(self.len)],
        )
        result = triangulate(data, 'p')

        mapping = dict()
        self.triangles = []

        for triple in result['triangles']:
            triangle = TriangulatedPolygonTriangle(*map(self.points.__getitem__, triple))
            self.triangles.append(triangle)
            for edge in triangle.edges:
                fst, snd = edge.a.index, edge.b.index
                if snd < fst:
                    fst, snd = snd, fst
                if (fst, snd) in mapping:
                    self.triangles[mapping[(fst, snd)]].neighbour_indices.append(len(self.triangles) - 1)
                    triangle.neighbour_indices.append(mapping[(fst, snd)])
                else:
                    mapping[(fst, snd)] = len(self.triangles) - 1

    def locate_point_in_triangle(self, p: Point) -> TriangulatedPolygonTriangle:
        """O(n): Find triangle in which p is located."""
        for triangle in self.triangles:
            if triangle.contains(p):
                return triangle

        raise ValueError()


#################################################################################################
#################################################################################################
####           ####   ##      ######  ######  ######   ##  ######  ##  ##  ##   ##           ####
####          ##  ##  ##      ##      ##  ##  ##   ##  ##    ##    ##  ##  #######           ####
####   	      ######  ##      ## ###  ##  ##  ######   ##    ##    ######  ## # ##           ####
####          ##  ##  ##      ##  ##  ##  ##  ##   ##  ##    ##    ##  ##  ##   ##           ####
####          ##  ##  ######  ######  ######  ##   ##  ##    ##    ##  ##  ##   ##           ####
#################################################################################################
################################################################################################# 


#LP_shortest_path_algorithm

def shortest_path_as_diagonals(polygon: TriangulatedPolygon, s_triangle: TriangulatedPolygonTriangle,
                               t_triangle: TriangulatedPolygonTriangle) -> List[Edge]:
    """O(n): Return all diagonals the shortest path crosses from s_triangle to t_triangle."""
    def recurse(triangle: TriangulatedPolygonTriangle, predecessor: TriangulatedPolygonTriangle = None
                ) -> List[Edge]:
        """
        Return a list of diagonals (last one first) from triangle to t_triangle ignoring predecessor.

        The diagonals are the dual edges of all the edges we need to visit in the tree towards t_triangle.
        """
        # Loop trough all the neighbours
        for neighbour_index in triangle.neighbour_indices:
            neighbour = polygon.triangles[neighbour_index]
            if neighbour == t_triangle:
                # We reached the goal t_triangle
                return [triangle.common_edge(neighbour)]
            if neighbour != predecessor:
                try:
                    # Try to look for t_triangle recursively. If t_triangle cannot be found an exception is thrown.
                    result = recurse(neighbour, triangle)
                except ValueError:
                    pass
                else:
                    # No exception thrown => t_triangle was found.
                    result.append(triangle.common_edge(neighbour))
                    return result

        raise ValueError()

    diagonals = recurse(s_triangle)
    diagonals.reverse()
    return diagonals

def shortest_path(polygon: TriangulatedPolygon, s: Point, t: Point) -> Iterable[Point]:
    """Find the shortest path from s to t inside polygon."""
    # ==========================================================================
    # Reset properties which can be accessed later on.
    # ==========================================================================
    shortest_path.properties = dict(iterations=0)

    # ==========================================================================
    # Type checking.
    # ==========================================================================
    assert isinstance(polygon, TriangulatedPolygon)
    assert isinstance(s, Point)
    assert isinstance(t, Point)

    # ==========================================================================
    # Trivial case: s == t.
    # ==========================================================================
    # In the very trivial case the start and end point are identical thus we can
    # just return without any calculation.
    if s == t:
        yield s
        return

    # ==========================================================================
    # Locate s and t. Trivial case: both in same triangle.
    # ==========================================================================
    # Locate start and end point inside our polygon.
    s_triangle = polygon.locate_point_in_triangle(s)
    t_triangle = polygon.locate_point_in_triangle(t)

    # If any point is not inside the polygon return
    if s_triangle is None or t_triangle is None:
        return

    # If both points are located inside the same triangle just return both in order
    if s_triangle == t_triangle:
        yield s
        yield t
        return

    # ==========================================================================
    # Find the shortest path in the dual tree.
    # ==========================================================================
    diagonals = shortest_path_as_diagonals(polygon, s_triangle, t_triangle)
    # Append one final edge from one end-point of the last diagonal to t. This is needed to make sure we ultimately
    # visit t.
    diagonals.append(Edge(diagonals[-1].a, t))

    # ==========================================================================
    # Preparation.
    # ==========================================================================
    apex = s
    funnel = deque([diagonals[0].a, apex, diagonals[0].b])

    if Point.turn(s, diagonals[0].a, diagonals[0].b) == Point.CCW_TURN:
        funnel.reverse()

    # ==========================================================================
    # Walking the triangles.
    # ==========================================================================
    for diagonal in diagonals[1:]:
        shortest_path.properties['iterations'] += 1
        # Save the new points as left and right
        left = diagonal.a
        right = diagonal.b

        # We know, that every new diagonal has exactly one end point common with the current funnel. We check whether
        # the common vertex is the "wrong" one and in this case swap the diagonal's "left" and "right".
        if funnel[0] == right or funnel[-1] == left:
            left, right = right, left

        if left == funnel[0]:
            # As long as the new point does not extend the existing funnel in a concave fashion we remove the last
            # funnel vertex. We also stop when reaching the apex since we need to think differently after reaching it.
            while funnel[-1] != apex and Point.turn(funnel[-2], funnel[-1], right) == Point.CCW_TURN:
                funnel.pop()
            if funnel[-1] == apex:
                # If we removed all funnel vertices of one side we might find the need to remove the apex and maybe some
                # more vertices from the other side. In the end the first not removed funnel vertex is the new apex.
                while len(funnel) > 1 and Point.turn(funnel[-1], funnel[-2], right) == Point.CCW_TURN:
                    yield funnel.pop()
                apex = funnel[-1]
            # Our new vertex definitely extends the (new) funnel.
            funnel.append(right)
        else:
            # This is exactly analogous to the left case.
            while funnel[0] != apex and Point.turn(funnel[1], funnel[0], left) == Point.CW_TURN:
                funnel.popleft()
            if funnel[0] == apex:
                while len(funnel) > 1 and Point.turn(funnel[0], funnel[1], left) == Point.CW_TURN:
                    yield funnel.popleft()
                apex = funnel[0]

            funnel.appendleft(left)

    # If either end of the funnel is our final point t we remove the other side until the apex and then yield all points
    # belonging to the apex because we need to visit all of them on our way to t.
    if funnel[0] == t:
        while funnel[-1] != apex:
            funnel.pop()
        while funnel:
            yield funnel.pop()
    elif funnel[-1] == t:
        while funnel[0] != apex:
            funnel.popleft()
        while funnel:
            yield funnel.popleft()
    else:
        yield apex
        yield t


##################################################################################################
##################################################################################################
####                ######  ##   ##  ######  ######  ##   ##  ########  ######                ####               
####                ##       ## ##   ##      ##      ##   ##     ##     ##                    ####
####   		        ######    ###    ######  ##      ##   ##     ##     ######                ####
####   		        ##       ## ##   ##      ##      ##   ##     ##     ##                    ####
####   		        ######  ##   ##  ######  ######  #######     ##     ######                ####
##################################################################################################
################################################################################################## 


#plotting
"""Collection of drawing functions."""

def draw_polygon(ax: Axes, polygon: Polygon, size: int, ticks: int, font_size: int) -> None:
    """Plot a polygon to the current matplotlib figure."""
    max_x, min_x = ceil(max(p.x for p in polygon.points)), floor(min(p.x for p in polygon.points))
    max_y, min_y = ceil(max(p.y for p in polygon.points)), floor(min(p.y for p in polygon.points))
    size_x = max_x - min_x
    size_y = max_y - min_y
    ax.set_aspect('equal')
    ax.set_xlim(min_x - size_x * 0.05, max_x + size_x * 0.05)
    ax.set_ylim(min_y - size_y * 0.05, max_y + size_y * 0.05)
    ticks -= 1
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=round(size_x / ticks)))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=round(size_y / ticks)))

    polygon_points = list(polygon.points_as_tuples())

    pyplot_polygon = pyplot.Polygon(polygon_points, fill=None, color='0.25', linewidth=.5)
    ax.add_patch(pyplot_polygon)
    pyplot.plot(list([x[0] for x in polygon_points]), list([x[1] for x in polygon_points]),
                color='black', linestyle='None', marker='.', markersize=3)

    if isinstance(polygon, TriangulatedPolygon):
        edges = set()
        for triangle in polygon.triangles:
            for edge in triangle.edges:
                if (edge.a.index, edge.b.index) in edges:
                    # The edge was already plotted
                    continue
                if edge.a.index in ((edge.b.index - 1) % polygon.len, (edge.b.index + 1) % polygon.len):
                    # The edge is a polygon edge, thus we ignore it
                    continue

                pyplot.plot([edge.a.x, edge.b.x], [edge.a.y, edge.b.y], color='0.5', linestyle='--', alpha=0.5)
                edges.add((edge.a.index, edge.b.index))
                edges.add((edge.b.index, edge.a.index))

    for i in range(0, polygon.len, 2):
        ax.annotate(str(i), xy=polygon_points[i],
                    xytext=(
                        polygon_points[i][0] + round(size_x / ticks) / 20,
                        polygon_points[i][1] + round(size_y / ticks) / 20),
                    fontsize=font_size, color='0.3', alpha=.8)


#execute

size = 10
digits = 3
font_size = 8
ticks = 10
timeout = .5
dpi = 300

points = [(-5.5, 0.96), (-4.5, 2.46), (-3.56, 1.44), (-2.5, 3.26), (-1.92, 1.86), (-0.48, 4.32), (0.28, 2.64), (1.38, 4.9), (2.64, 2.9), (3.8, 1.9), (4.86, 0.98), (5, 0), (5.16, -1.16), (4.32, -2.38), (3, -2), (2.76, -0.86), (1.98, -2.52), (1, -1), (0.52, -2.62), (-1.2, -2.72), (-1.1, -1.2), (-2.82, -1.56), (-2.84, -2.7), (-4.46, -2.84), (-4.38, -1.46), (-6.48, -1.62), (-7.68, -0.76), (-7.6, 0.18)]

start = Point(-4.38, 1.8)
end = Point(3.8, -1.58)

p_list = []
for point in points:
    p_list.append(Point(point[0],point[1]))
polygon = Polygon(p_list)
triangulated = TriangulatedPolygon(polygon)
s_path = shortest_path(triangulated,start,end)

#-------Original-------#
fg, ax = pyplot.subplots()
draw_polygon(ax,polygon,size,ticks,font_size)
pyplot.plot([start.x, end.x],[start.y, end.y],'ro')
pyplot.show()
#-------Triangulated-------#
fg, ax = pyplot.subplots()
draw_polygon(ax,polygon, size,ticks,font_size)
draw_polygon(ax,triangulated, size,ticks,font_size)
pyplot.plot([start.x, end.x],[start.y, end.y],'ro')
pyplot.show()
#-------Shortest Path-------#
fg, ax = pyplot.subplots()
draw_polygon(ax,polygon, size,ticks,font_size)
x_coords = []
y_coords = []
for point in s_path:
    x_coords.append(point.x)
    y_coords.append(point.y)
pyplot.plot(x_coords,y_coords,'ro',linestyle = 'solid')
pyplot.show()

