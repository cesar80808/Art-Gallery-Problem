import numpy as np
import matplotlib.pyplot as plt

EPS = 1e-8
# Color constants, used in coloring_BFS()
NO_CLR = -1
WHITE = 0
GRAY = 1
BLACK = 2

def cross2d(a, b):
    return (a[0] * b[1] - a[1] * b[0])


def CCW(A, B, C):
    """
    Determines if points A, B, C perform a counter clock-wise turn.
    
    Parameters
    ----------
    A, B, C : tuple[float, float]
        Coordinates of the query points.
    
    Returns
    ----------
    bool
        True if A, B, C perform a counter clock-wise turn, False otherwise. 
    """
    AB = (B[0] - A[0], B[1] - A[1])
    AC = (C[0] - A[0], C[1] - A[1])
    return (cross2d(AB, AC) > -EPS)


def is_point_in_triangle(P, T):
    """
    Check whether a point lies inside a triangle.

    Parameters
    ----------
    P : tuple[float, float]
        Coordinates of the query point.
    T: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
        Coordinates of triangle vertices.

    Returns
    ----------
    bool
        True if P lies on the inside of T, False otherwise.
    """
    ccw = [False, False, False]
    for i in range(3):
        ccw[i] = CCW(T[i], T[(i + 1) % 3], P)
    if ccw == [True, True, True] or ccw == [False, False, False]:
        return True
    return False
    

def shares_edge(T1, T2):
    """
    Determines if two triangles share exactly one edge.
    
    Parameters
    ----------
    T1, T2 : tuple[int, int, int]
        Triangles, represented by indices of polygon vertices.
        
    Returns
    ----------
    bool
        True if the two triangles share exactly one edge (i.e. have exactly
        two vertices in common), False otherwise.
    """
    matches = 0
    for i in T1:
        for j in T2:
            if i == j:
                matches += 1
    if matches == 2:
        return True
    return False


def ear_clipping_triangulation(polygon):
    """
    Triangulates a given simple polygon using the ear clipping algorithm. 

    Parameters
    ----------
    polygon : list[tuple[float, float]]
        Vertices of a simple polygon listed in counter clock-wise order.

    Returns
    ----------
    list[tuple[int, int, int]]
        The resulting triangulation, given as a list of triangles. Each
        triangle is given as indices of polygon vertices.
    """
    n = len(polygon)
    vertices = list(range(n))
    final_triangle_list = []

    i = 0
    while n > 2:
        pi = vertices[(i - 1 + n) % n]  
        ci = vertices[i]                
        ni = vertices[(i + 1) % n]    

        prev = polygon[pi]
        curr = polygon[ci]
        next = polygon[ni]

        # Check current vertex is not reflex
        if not CCW(prev, curr, next):
            i = (i + 1) % n
            continue

        # Check that triangle (prev, curr, next) doesn't contain any point
        # from 'vertices' in its interior
        T = (pi, ci, ni)
        point_in_triangle = False
        for j in range(n):
            idx = vertices[j]
            if idx in T:  # Exclude vertices of T itself
                continue
            if is_point_in_triangle(polygon[idx], [prev, curr, next]):
                point_in_triangle = True
                break
        if point_in_triangle:
            i = (i + 1) % n
            continue
        
        # At this point, current vertex is an ear
        final_triangle_list.append(T)
        del vertices[i]
        n -= 1
        i %= n
    
    return final_triangle_list


def dual_graph(triangles):
    """
    Returns the dual graph of a given triangulation.

    Each node of the graph corresponds to a triangle, and an undirected edge
    connects two nodes if the corresponding triangles share an edge.

    Parameters
    ----------
    triangles : list[tuple[int, int, int]]
        A triangulation of a polygon given as a list of triangles, where each
        triangle is given by indices of polygon vertices.

    Returns
    ----------
    list[list[int]]
        Adjacency list representation of the dual graph.
    """
    t = len(triangles)
    DG = [[] for _ in range(t)]
    for i in range(t):
        for j in range(i+1, t):
            if shares_edge(triangles[i], triangles[j]):
                DG[i].append(j)
                DG[j].append(i)
    return DG

def coloring_BFS(triangles, DG):
    """
    Compute a 3-coloring of polygon vertices induced by a triangulation.

    The algorithm assigns one of three colors to each polygon vertex such that
    every triangle in the triangulation has all three vertices with distinct
    colors. The coloring order is determined by a breadth-first traversal of
    the dual graph.

    Parameters
    ----------
    triangles : list[tuple[int, int, int]]
        Triangulation of a polygon, given as a list of triangles, where each
        triangles is represented by indices of polygon vertices.

    DG : list[list[int]]
        Dual graph of the triangulation, represented by adjacency lists.

    Returns
    ----------
    list[int]
        Indices of polygon vertices having least frequent color.
    """
    # List of vertices of each color.
    white = []
    gray = []
    black = []
    # color[i] is the color of vertex i.
    color = [NO_CLR for _ in range(len(triangles) + 2)]

    # We color white and gray vertices 0 and 1 of triangle[0].
    white.append(triangles[0][0])
    color[triangles[0][0]] = WHITE
    gray.append(triangles[0][1])
    color[triangles[0][1]] = GRAY

    # BFS on dual graph
    visited = [False for _ in range(len(triangles))]
    visited[0] = True
    queue = []
    queue.append(0)
    while len(queue) > 0:
        t_idx = queue.pop(0)
        T = triangles[t_idx]
        nc = 0      # Not colored index.
        sum = 0     # Sum of colored vertices.
        for i in range(3):
            if color[T[i]] == NO_CLR:
                nc = i
            else:
                sum += color[T[i]]
        # Due to the way we defined colors as values 0, 1, 2, every color is
        # equal to 3 minus the sum of the other two. 
        color[T[nc]] = 3 - sum
        # Update color list
        if color[T[nc]] == WHITE:
            white.append(T[nc])
        elif color[T[nc]] == GRAY:
            gray.append(T[nc])
        else:
            black.append(T[nc])
        # Expand BFS.
        for neighbour in DG[t_idx]:
            if not visited[neighbour]:
                visited[neighbour] = True
                queue.append(neighbour)
    
    bound = min(len(white), len(gray), len(black))

    if len(white) == bound:
        return white
    elif len(gray) == bound:
        return gray
    else:
        return black
    
def solve_art_gallery_problem(polygon):
    """
    Finds a set of solution points to the Art Gallery Problem of a simply 
    polygon.

    The returned points are vertices of the input polygon. For a polygon with
    n vertices, the algorithm returns at most floor(n/3) vertices, as 
    guaranteed by the Art Gallery Theorem.

    Parameters
    ----------
    polygon : list[tuple[float, float]]
        Vertices of the polygon, listed in counter clock-wise order.

    Returns
    ----------
    list[tuple[float, float]]
        Solution points to the Art Gallery Problem.
    """
    triangles = ear_clipping_triangulation(polygon)
    DG = dual_graph(triangles)
    indices = coloring_BFS(triangles, DG)
    solution_points = [polygon[i] for i in indices]
    return solution_points


def display_solution(polygon, solution_points):
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    # Plot polygon
    xs = [polygon[i][0] for i in range(len(polygon))]
    xs.append(polygon[0][0])
    ys = [polygon[i][1] for i in range(len(polygon))]
    ys.append(polygon[0][1])
    ax.plot(xs, ys)
    # Plot solution points
    for point in solution_points:
        ax.plot(point[0], point[1], "o", color = "red")
        plt.title("Solution to the Art Gallery Problem")
    plt.show()
    return