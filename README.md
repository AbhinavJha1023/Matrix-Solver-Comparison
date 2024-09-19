# Matrix-Solver-Comparison


This repository compares various matrix-solving methods for linear systems of equations. It includes implementations of Gauss Elimination, Gauss-Jordan, Gauss-Seidel, Gauss-Jacobi, Crout Method, and Doolittle Method. Each method's performance is evaluated based on the time taken to solve randomly generated matrices of increasing size.

## Features

- **Matrix Methods Implemented:**
  - Gauss Elimination
  - Gauss-Jordan
  - Gauss-Seidel
  - Gauss-Jacobi
  - Crout Method
  - Doolittle Method

- **Performance Metrics:**
  - Average time for each method across multiple trials
  - Median time for each method
  - Fastest method for each matrix size
  - Bar graph representation of time vs matrix size

## Performance Evaluation

The script evaluates each matrix-solving method by solving random matrices of increasing size (from 2x2 to 20x20) and records the time taken by each method to find the solution.

- **Average Time**: The average time taken by each method across 50 trials for each matrix size.
- **Median Time**: The median time for each method to ensure robustness against outliers.
- **Fastest Method**: For each matrix size, the method that solves the system in the least amount of time is recorded.
- **Bar Graph**: A bar graph is generated to show the fastest times per matrix size.

## Methods Used

1. **Gauss Elimination**
2. **Gauss-Jordan**
3. **Gauss-Seidel**
4. **Gauss-Jacobi**
5. **Crout Method**
6. **Doolittle Method**

## How to Run

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/Matrix-Solver-Comparison.git
    ```

2. Run the MATLAB code to compare the performance:
    ```matlab
    % Execute the main MATLAB script that compares the methods
    ```

3. The script will display:
    - Mean, median, and fastest method for each matrix size
    - A bar graph showing the fastest method for each size

## Dependencies

- MATLAB (Version 2021 or later)

## Example Output

![Screenshot 2024-09-19 193331](https://github.com/user-attachments/assets/8e301dcb-6f92-464c-9108-fc966f6be17d)

![Screenshot 2024-09-19 193343](https://github.com/user-attachments/assets/db8fc7f0-b657-42ac-bc79-8e0b090c3e09)


