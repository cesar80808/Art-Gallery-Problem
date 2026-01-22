def read_polygon(filename):
    """Reads polygon from file, returning a list of tuples."""
    with open(filename, 'r') as file:
        points = []
        for line in file:
            x, y = map(float, line.split())
            points.append((x, y))
    return points

def output_solution(filename, solution_points):
    """Writes list of points to a file.
    
    args:
        - filename (str): name of the output file
        - soultion_points (list): a list of tuples representing the points to
            to be outputed.   
    """
    with open(filename, 'a') as file:
        for x, y in solution_points:
            file.write(f"{x} {y}\n")
    return
