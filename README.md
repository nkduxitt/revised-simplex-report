# Revised Simplex Method for Solving Linear Programming Problems

This repository presents a scientific study and implementation of the **Revised Simplex Method**, an efficient variant of the classical Simplex Algorithm for solving **Linear Programming Problems (LPP)**. The project includes a complete LaTeX report, Python implementation, worked examples, and a comparison with SciPy’s built-in solver.

---

## Project Structure

```
project/
│
├─ main.tex                     # Main LaTeX file (entry point)
│
├─ 01_motivation.tex            # Motivation / Introduction
├─ 02_theory.tex                # Mathematical formulation
├─ 03_algorithm.tex             # Pseudocode and method description
├─ 04_implementation.tex        # Implementation notes and Python skeleton
├─ 05_examples.tex              # Solved examples and experiments
├─ 06_comparison.tex            # Comparison with SciPy solver
├─ 07_conclusion.tex            # Final discussion and summary
│
├─ code/
│   └─ revised_simplex.py       # Python implementation of the Revised Simplex
│
├─ results/
│   ├─ example_output.txt       # Example iteration logs
│   └─ comparison.csv           # Benchmark comparison with SciPy
│
├─ images/                      # Figures and diagrams
│   └─ flowchart.png
│
├─ refs.bib                     # BibTeX references
├─ requirements.txt             # Python dependencies
└─ README.md                    # Project documentation
```

---

## Report Overview

The LaTeX report is structured into seven chapters:

| Chapter        | Description                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Motivation     | Explains the need for the Revised Simplex Method and its advantages over the classical tableau approach.                               |
| Theory         | Formulates the linear programming problem, defines basis matrices, reduced costs, and optimality conditions.                           |
| Algorithm      | Presents detailed pseudocode for Phase I and Phase II of the Revised Simplex method, including pivot selection and feasibility checks. |
| Implementation | Discusses implementation choices, numerical safeguards, and data structures in Python.                                                 |
| Examples       | Provides worked examples for regular, two-phase, and special-case LP problems.                                                         |
| Comparison     | Compares the custom implementation with SciPy’s `linprog` in terms of runtime, accuracy, and feasibility.                              |
| Conclusion     | Summarizes the study and suggests extensions for large-scale or sparse LP problems.                                                    |

---

## Requirements

### LaTeX

To compile the report, ensure a LaTeX environment with the following packages:

* `amsmath`, `amssymb`, `graphicx`, `geometry`, `algorithm`, `algpseudocode`, `booktabs`, `hyperref`

### Python

The implementation requires:

* Python 3.9 or higher
* `numpy`
* `scipy`
* `pandas` (optional, for data comparison)

Install dependencies using:

```bash
pip install -r requirements.txt
```

---

## Compiling the Report

### Overleaf

1. Upload all `.tex` files, the `code/` folder, and `refs.bib`.
2. Open `main.tex`.
3. Click **Recompile** to produce `main.pdf`.

### Local Compilation

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

The output file will be `main.pdf`.

---

## Python Implementation 

The `revised_simplex.py` module contains a functional Python implementation of the Revised Simplex algorithm. Example usage:

```python
import numpy as np
from revised_simplex import RevisedSimplex

A = np.array([[1, 1, 1, 0],
              [1, 3, 0, 1]], dtype=float)
b = np.array([4, 6], dtype=float)
c = np.array([3, 2, 0, 0], dtype=float)
B_init = [2, 3]  # slack variables

solver = RevisedSimplex(A, b, c, B_init)

while True:
    status = solver.iterate()
    if status != 'continue':
        print("Status:", status)
        break
```

Expected output for the example:

```
optimal
```

---

## Comparison with SciPy

To verify results against SciPy:

```python
from scipy.optimize import linprog

res = linprog(-c, A_ub=A, b_ub=b, method='highs')
print(res.x, -res.fun)
```

| Problem   | Method          | Optimal Value | Feasibility | Notes                |
| --------- | --------------- | ------------- | ----------- | -------------------- |
| Example 1 | Revised Simplex | 10.0          | Feasible    | Matches SciPy result |
| Example 1 | SciPy (HiGHS)   | 10.0          | Feasible    | Faster runtime       |
| Example 2 | Revised Simplex | TBD           | Feasible    | Phase I test         |
| Example 2 | SciPy           | TBD           | Feasible    | Reference            |

---

## Extending the Project

Potential extensions include:

* Implementing full Phase I–Phase II procedure for infeasible LPs
* Integrating LU factorization updates for improved performance
* Testing large-scale or sparse LP problems
* Comparing with industrial solvers (Gurobi, CPLEX)

---

## Authors

* **Omar Hazem Ahmed** — Zewail City, Department of CSAI
* **Amr Yasser Anwar** — Zewail City, Department of CSAI

---

## References

1. Dantzig, G. B., *Linear Programming and Extensions*, Princeton University Press, 1963.
2. Chvátal, V., *Linear Programming*, W. H. Freeman, 1983.
3. Hall, J. A. J., “HiGHS: A high-performance linear optimization solver,” University of Edinburgh, 2020.

---

## License

This project is intended for **educational and academic purposes**. Proper citation is required for use in research or coursework.

