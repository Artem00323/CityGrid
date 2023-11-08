'''
Task 1: Grid Representation
Create a class CityGrid that can represent the city as an N x M grid. 
During the initialization of the class, obstructed blocks are 
randomly placed with coverage >30% (we can change this parameter).
'''
import numpy as np
import random
from collections import deque # for task4 (explained why in task)
import matplotlib.pyplot as plt # for task5
import matplotlib.colors as mcolors # for task5

class CityGrid:
    def __init__(self, n, m, coverage=30): # Parameters - N and M (rows x columns), coverage - obstructed percentage of city.
        self.n = n
        self.m = m
        # self.coverage = random.randint(coverage + 1, 100) # coverage > 30% (given number). coverage is random
        self.coverage = coverage
        self.grid = np.zeros((n, m), dtype=int)  # Create an N x M grid initialized to 0
        self.place_obstructions()

    def place_obstructions(self): # Randomly place obstructions in the grid based on the coverage percentage.
        total_blocks = self.n * self.m
        obstructed_blocks = total_blocks * self.coverage // 100

        while obstructed_blocks > 0:
            i = random.randint(0, self.n - 1)
            j = random.randint(0, self.m - 1)
            if self.grid[i][j] == 0:  # If the block is not already obstructed
                self.grid[i][j] = 1  # Mark the block as obstructed
                obstructed_blocks -= 1

    def place_tower(self, x, y, R): # Parameters - x and y (rows x columns), R - range of tower.
        self.R = R
        '''
        Task 2: Tower Coverage
        Each tower has a fixed range R (in blocks) within which it provides coverage. 
        This coverage is a square, with the tower in the center.
        Implement a method in the CityGrid class to place a tower and visualize its coverage.
        (Assuming coverage of the tower doesn't depend on obstructed blocks)
        '''
        # Bounds of the coverage area
        top = max(0, x - R)
        bottom = min(self.n - 1, x + R)
        left = max(0, y - R)
        right = min(self.m - 1, y + R)

        self.grid[x][y] = 8 # Initiate tower (with different symbol)

        # Coverage area on the grid
        for i in range(top, bottom + 1):
            for j in range(left, right + 1):
                if self.grid[i][j] not in [1, 8]:  # Only mark unobstructed blocks
                    self.grid[i][j] = 2  # Mark as covered

    def optimize_tower_placement(self, R): # Parameters - R - range of tower
        self.R = R
        '''
        Task 3: Optimization Problem
        Design an algorithm to place the minimum number of towers such 
        that all of non-obstructed blocks are within the coverage of 
        at least one tower. The algorithm cannot place towers on obstructed blocks.
        Implement a method in the CityGrid class to display the placement of towers.
        '''
        # Initialize towers list to keep track of tower positions
        self.towers = []

        # Create a copy of the grid to mark covered areas
        coverage_grid = np.copy(self.grid)

        # Function to mark coverage
        def mark_coverage(x, y, R):
            top = max(0, x - R)
            bottom = min(self.n - 1, x + R)
            left = max(0, y - R)
            right = min(self.m - 1, y + R)
            for i in range(top, bottom + 1):
                for j in range(left, right + 1):
                    if coverage_grid[i][j] not in [1, 8]:  # Only mark unobstructed blocks
                        coverage_grid[i][j] = 2  # Mark as tower coverege

        # Find the next best block for tower placement
        def find_next_block():
            candidates = np.argwhere(coverage_grid == 0) # maybe tower can be placed on other tower's coverage (but skiped)
            if len(candidates) > 0:
                # Heuristic: Choose the block closest to the center of the grid
                '''
                Can lead to lower number of towers to cover all unobstructed 
                blocks (as it starts from the closes to the center unobstructed block)
                => will probably cover more unobstructed blocks. "Just Assumption"
                '''
                center = np.array([self.n // 2, self.m // 2])
                closest_block = candidates[((candidates - center)**2).sum(axis=1).argmin()]
                return tuple(closest_block)
            return None

        # Place towers
        while True:
            block = find_next_block()
            if block is None:  # No more blocks to cover
                break
            self.place_tower(*block, R)
            mark_coverage(*block, R)
            self.towers.append(block)

    def find_most_reliable_path(self, start_tower, end_tower): # Parameters - start_tower and end_tower (row, column indices), R...
        '''
        Task 4: Path Reliability
        Imagine that data is transmitted between towers. For simplicity, 
        assume that each tower can directly communicate with any other tower within its range.
        Design an algorithm to find the most reliable path between two towers. 
        The reliability of a path decreases with the number of hops (tower-to-tower links). 
        So, a path with fewer hops is more reliable.
        (Used BFS algorithm with deque, which is O(1), whereas list in python is of higher complexity)
        '''
        graph = {}
        for tower in self.towers:
            x, y = tower
            for dx in range(-self.R, self.R + 1):
                for dy in range(-self.R, self.R + 1):
                    if 0 <= x + dx < self.n and 0 <= y + dy < self.m:
                        if self.grid[x + dx][y + dy] != 1:  # If the cell is not obstructed
                            graph.setdefault((x + dx, y + dy), set()).add(tower)
                            graph.setdefault(tower, set()).add((x + dx, y + dy))

        # BFS to find the shortest path
        queue = deque([[start_tower]])
        visited = set()
        while queue:
            path = queue.popleft()
            vertex = path[-1]
            if vertex == end_tower:
                self.most_reliable_path = path
                return path
            elif vertex not in visited:
                visited.add(vertex)
                for neighbor in graph[vertex]:
                    if neighbor not in visited:
                        new_path = list(path)
                        new_path.append(neighbor)
                        queue.append(new_path)
        return None


    def visualize_grid(self):
        '''
        Task 5: Visualization
        Implement functions to visualize the CityGrid, 
        including obstructed blocks, towers, coverage areas, and data paths.
        Use any Python plotting library of your choice, such as matplotlib or seaborn.
        '''
        # Colors for different elements in the grid
        cmap = mcolors.ListedColormap(['white', 'black', 'green', 'red', 'blue'])
        bounds = [0, 1, 2, 8, 9, 10]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        # Figure and axis for plottin the grid
        fig, ax = plt.subplots()
        ax.imshow(self.grid, cmap=cmap, norm=norm)

        # Towers on the grid
        for tower in self.towers:
            ax.text(tower[1], tower[0], 'T', va='center', ha='center', color='white', weight='bold')

        # The most reliable path if it exists
        if hasattr(self, 'most_reliable_path') and self.most_reliable_path:
            path_x, path_y = zip(*self.most_reliable_path)
            ax.plot(path_y, path_x, marker='o', markersize=5, linestyle='--', color='blue')

        # Grid labels and title
        ax.set_xticks(range(self.m))
        ax.set_yticks(range(self.n))
        ax.set_xticklabels(range(1, self.m + 1))
        ax.set_yticklabels(range(1, self.n + 1))
        ax.set_title('City Grid Visualization')

        plt.show()

    def display_towers(self):
        print(f"Number of towers = {len(self.towers)}")
        for x, y in self.towers:
            print(f"Tower placed at: ({x}, {y})")
        self.display_grid()

    def display_grid(self):
        for row in self.grid:
            print(' '.join(str(x) for x in row))



''' TESTING '''

city = CityGrid(10, 20, 30)
city.display_grid()

# city.place_tower(5, 5, 2)
# print("\nGrid after placing a tower at position (5,5) with range 2:")
# city.display_grid()

city.optimize_tower_placement(2)
print("\nOptimal tower placement:")
city.display_towers()
# city.display_grid()

start = city.towers[0]
end = city.towers[len(city.towers) - 1]
print(f'{start} -> {end}')
path = city.find_most_reliable_path(start, end)
if path:
    print("Most reliable path:", path)
else:
    print("No path found between the towers.")

city.visualize_grid()