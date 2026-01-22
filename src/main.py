import my_io 
from artgallprob import solve_art_gallery_problem, display_solution

input_filename = "input/points1.txt"
output_filename = "output/solution1.txt"

def main():
    polygon = my_io.read_polygon(input_filename)
    solution_points = solve_art_gallery_problem(polygon)
    my_io.output_solution(output_filename, solution_points)
    display_solution(polygon, solution_points)

if __name__ == "__main__":
    main()
