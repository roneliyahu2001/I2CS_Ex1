# Ex1 – Polynomials, Static Functions and JUnit

## Author
- Name: Ron Eliyahu
- ID: 323951038
- Course: Introduction to Computer Science 2026, Ariel University
- Exercise: Ex1 – arrays, static functions and JUnit

---

## Project Description

This exercise implements a set of static functions that work with polynomial functions,
represented as arrays of doubles.  
The project also includes a JUnit test class and a simple GUI for visualizing the polynomials.

A polynomial is represented as:  
`p[i]` = the coefficient of `x^i`.  
For example, the array `{2, 0, 3.1, -1.2}` represents the function `-1.2x^3 + 3.1x^2 + 2`.

---

## Implemented Functions in `Ex1.java`

- `f(double[] poly, double x)`  
  Computes the value of the polynomial at a given x.

- `poly(double[] poly)`  
  Returns a String representation of the polynomial (for example: `"-1.2x^3 +3.1x^2 +2.0"`).

- `getPolynomFromString(String p)`  
  Parses a String representing a polynomial and returns the matching coefficients array.

- `add(double[] p1, double[] p2)`  
  Computes the sum `p1 + p2` as a new polynomial.

- `mul(double[] p1, double[] p2)`  
  Computes the multiplication `p1 * p2` as a new polynomial.

- `derivative(double[] p)`  
  Computes the derivative of a polynomial.

- `equals(double[] p1, double[] p2)`  
  Checks if two polynomials represent the same function (up to EPS).

- `sameValue(double[] p1, double[] p2, double x1, double x2, double eps)`  
  Finds an x in `[x1, x2]` such that `|p1(x) - p2(x)| < eps`.

- `length(double[] p, double x1, double x2, int numberOfSegments)`  
  Approximates the length of the curve `y = p(x)` between `x1` and `x2`.

- `area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid)`  
  Approximates the area between two polynomials on `[x1, x2]` using trapezoids.

- `PolynomFromPoints(double[] xx, double[] yy)`  
  Computes a polynomial that passes through 1–3 given points.

---

## JUnit Tests – `Ex1Test.java`

The file `Ex1Test.java` includes:

- Tests for `f` (checking basic values of a simple polynomial).
- Tests for `add` (associativity with ZERO, commutativity, adding the negative polynomial).
- Tests for `mul` (multiplication with ZERO, commutativity, and checking `p1(x)*p2(x) == (p1*p2)(x)`).
- Tests for `derivative` (repeated derivatives until reaching ZERO).
- Tests for `equals` (arrays that should be considered equal up to EPS and arrays that should not).
- Tests for `poly` and `getPolynomFromString` (round–trip from array → String → array).
- Tests for `sameValue` (symmetry and a simple example of intersection).
- Tests for `area` (symmetry, simple case with `f(x)=0` and `f(x)=x`, and a more complex example).

Additional tests I added:

- `testLengthSimpleLine` – checks that the length of the line `y = x` between 0 and 1 is approximately √2.
- `testPolynomFromPointsLine` – checks that `PolynomFromPoints` returns the correct line through two points.
- `testSameValueSimple` – checks that `sameValue` finds the intersection between `y = x` and `y = 1` at `x = 1`.
- `testGetPolynomFromStringZero` – checks that parsing the String `"0"` returns the ZERO polynomial.
- `testDerivativeZeroAndConstant` – checks that the derivative of a constant polynomial is ZERO, and the derivative of ZERO is also ZERO.

---

## How to Run

### Run the GUI

1. Open `Ex1_GUI.java`.
2. Run the `main` function.
3. A window will open showing the graphs of the polynomials and allowing interaction.

### Run the JUnit tests

1. Open `Ex1Test.java`.
2. Click the green ▶ icon near the class name `Ex1Test`.
3. All tests in the class will run and show green (pass) or red (fail).

---

## Notes

- Polynomials are compared using an EPS value (`Ex1.EPS`) to allow small floating point differences.
- The implementation focuses on clear and simple code, written in the style of a first–semester student.# I2CS_Ex1
![תמונה 25 11 2025 ב-21 03](https://github.com/user-attachments/assets/f37d26d3-c691-4217-a206-227e3098f103)
