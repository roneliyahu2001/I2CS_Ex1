/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
    /**
     * Epsilon value for numerical computation, it serves as a "close enough" threshold.
     */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /**
     * The zero polynomial function is represented as an array with a single (0) entry.
     */
    public static final double[] ZERO = {0};

    /**
     * Computes the f(x) value of the polynomial function at x.
     *
     * @param poly - polynomial function
     * @param x
     * @return f(x) - the polynomial function value at x.
     */
    public static double f(double[] poly, double x) {
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * This function should be implemented recursively.
     *
     * @param p   - the polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    /**
     * This function computes a polynomial representation from a set of 2D points on the polynom.
     * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
     * Note: this function only works for a set of points containing up to 3 points, else returns null.
     *
     * @param xx
     * @param yy
     * @return an array of doubles representing the coefficients of the polynom.
     * <p>
     *
    check that both input arrays are valid
    check that the arrays have the same number of points
    check that the number of points is between 1 and 3
    handle the case of one point (return constant polynomial)
    handle the case of two points (compute line passing through the points)
    compute the slope using the two points
    compute the intercept using one of the points
    handle the case of three points (compute quadratic polynomial)
    compute differences between x values
    compute differences between y values
    compute sums needed for the quadratic formulas
    compute the coefficient a of the quadratic polynomial
    compute the coefficient b of the quadratic polynomial
    compute the constant term c of the quadratic polynomial
    return the polynomial coefficients in an array
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        if (xx == null || yy == null) return null;
        if (xx.length != yy.length) return null;
        if (xx.length < 1 || xx.length > 3) return null;

        if (xx.length == 1) {
            double c = yy[0];
            return new double[]{c};
        }

        if (xx.length == 2) {
            double x1 = xx[0], y1 = yy[0];
            double x2 = xx[1], y2 = yy[1];

            double a = (y2 - y1) / (x2 - x1);
            double b = y1 - a * x1;

            return new double[]{a, b};
        }

        if (xx.length == 3) {
            double x1 = xx[0], y1 = yy[0];
            double x2 = xx[1], y2 = yy[1];
            double x3 = xx[2], y3 = yy[2];

            double dx21 = x2 - x1;
            double dx31 = x3 - x1;

            double dy21 = y2 - y1;
            double dy31 = y3 - y1;

            double sum21 = x2 + x1;
            double sum31 = x3 + x1;

            double a = (dy21 * dx31 - dy31 * dx21) /
                    (dx21 * sum21 * dx31 - dx31 * sum31 * dx21);

            double b = (dy21 - a * dx21 * sum21) / dx21;

            double c = y1 - a * x1 * x1 - b * x1;

            return new double[]{a, b, c};
        }

        return null;
    }

    /**
     * Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
     *
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true iff p1 represents the same polynomial function as p2.
     *
     * check if both polynomials are null
     * check if only one of the polynomials is null
     * compute the degree of the first polynomial
     * compute the degree of the second polynomial
     * find the maximum degree between the two
     * loop over n+1 different x values
     * compute the value of the first polynomial at x
     * compute the value of the second polynomial at x
     * compare the two values using EPS
     * if the difference is too big, the polynomials are not equal
     * return the final result (true or false)
     */
    public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;

        if (p1 == null && p2 == null) return true;
        if (p1 == null || p2 == null) return false;

        int deg1 = p1.length - 1;
        int deg2 = p2.length - 1;
        int n = Math.max(deg1, deg2);

        for (int x = 0; x <= n && ans; x++) {
            double v1 = f(p1, x);
            double v2 = f(p2, x);

            if (Math.abs(v1 - v2) > EPS) {
                ans = false;
            }
        }

        return ans;
    }

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
     *
     * check if the polynomial is empty
     * start with an empty string for the result
     * this flag tells us if we are adding the first non-zero term
     * loop from the highest degree down to 0
     * get the coefficient of the current term
     * skip terms with coefficient zero
     * add a plus or space depending on the sign (not for the first term)
     * add the coefficient to the result
     * add x^i for powers greater than 1
     * add x for power 1
     * mark that we already added the first term
     * if no non-zero terms were added, the polynomial is zero
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            boolean first = true;

            for (int i = poly.length - 1; i >= 0; i--) {
                double c = poly[i];
                if (c == 0) continue;

                if (!first) {
                    if (c > 0) ans += " +";
                    else ans += " ";
                }

                ans += c;

                if (i > 1) ans += "x^" + i;
                else if (i == 1) ans += "x";

                first = false;
            }

            if (ans.equals("")) ans = "0";
		}
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     *
     * compute the difference at x1
     * compute the difference at x2
     * repeat while the interval is bigger than eps
     * take the midpoint between x1 and x2
     * compute the difference at the midpoint
     * check if the midpoint is close enough (difference < eps)
     * decide which side contains the zero crossing
     * update the interval to the correct half
     * save the current midpoint as the answer
	 */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double mid = (x1 + x2) / 2.0;
        double y1 = f(p1, x1) - f(p2, x1);
        double ymid = f(p1, mid) - f(p2, mid);

        while (Math.abs(ymid) > eps) {
            if (y1 * ymid <= 0) {
                x2 = mid;
            } else {
                x1 = mid;
                y1 = ymid;
            }
            mid = (x1 + x2) / 2.0;
            ymid = f(p1, mid) - f(p2, mid);
        }

        return mid;
    }
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
     *
     * compute the width of each segment (dx)
     * set the previous x value to the start of the interval
     * compute the previous y value using the polynomial
     * loop through all the segments
     * compute the current x position
     * compute the current y value
     * compute the horizontal difference between points
     * compute the vertical difference between points
     * compute the distance between the two points
     * add the distance to the total length
     * update the previous point to the current point
     * return the final length approximation
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        double dx = (x2 - x1) / numberOfSegments;
        double xPrev = x1;
        double yPrev = f(p, xPrev);

        for (int i = 1; i <= numberOfSegments; i++) {
            double xCurr = x1 + i * dx;
            double yCurr = f(p, xCurr);
            double dxSeg = xCurr - xPrev;
            double dySeg = yCurr - yPrev;
            ans += Math.sqrt(dxSeg * dxSeg + dySeg * dySeg);
            xPrev = xCurr;
            yPrev = yCurr;
        }

		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
     *
     * compute the width of each trapezoid (dx)
     * initialize a variable to store the total area
     * loop over all trapezoids between x1 and x2
     * compute the left x value of the current trapezoid
     * compute the right x value of the current trapezoid
     * compute the difference between the two polynomials at the left point
     * compute the difference between the two polynomials at the right point
     * take the absolute values to get positive heights
     * compute the trapezoid area using the average of the two heights times dx
     * add the trapezoid area to the total area
     * return the total approximated area

	 */
    public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0;
        double dx = (x2 - x1) / numberOfTrapezoid;

        for (int i = 0; i < numberOfTrapezoid; i++) {
            double xLeft = x1 + i * dx;
            double xRight = xLeft + dx;

            double y1 = f(p1, xLeft) - f(p2, xLeft);
            double y2 = f(p1, xRight) - f(p2, xRight);

            double trapArea = (Math.abs(y1) + Math.abs(y2)) * 0.5 * dx;
            ans += trapArea;
        }

        return ans;
    }
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
     *
     * check if the string is null or empty
     * remove extra spaces from the string
     * handle the case where the string represents zero
     * split the string into separate terms
     * find the highest power in the polynomial
     * create an array large enough to hold all coefficients
     * go over all terms again
     * extract the sign of the term
     * extract the coefficient value
     * determine the power of the term
     * place the coefficient in the correct index of the array
     * return the final polynomial coefficients array
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        if (p == null) {
            return ans;
        }

        p = p.trim();
        if (p.length() == 0 || p.equals("0")) {
            return ans;
        }

        String[] parts = p.split(" ");
        int maxPow = 0;

        for (int i = 0; i < parts.length; i++) {
            String term = parts[i];
            if (term.length() == 0) continue;

            int sign = 1;
            if (term.charAt(0) == '+') {
                term = term.substring(1);
            } else if (term.charAt(0) == '-') {
                sign = -1;
                term = term.substring(1);
            }

            int pow = 0;

            if (term.contains("x")) {
                int xIndex = term.indexOf('x');
                String coeffStr = term.substring(0, xIndex);
                if (coeffStr.length() == 0) {
                    coeffStr = "1.0";
                }
                if (term.contains("^")) {
                    int powIndex = term.indexOf('^');
                    String powStr = term.substring(powIndex + 1);
                    pow = Integer.parseInt(powStr);
                } else {
                    pow = 1;
                }
            } else {
                pow = 0;
            }

            if (pow > maxPow) {
                maxPow = pow;
            }
        }

        double[] coeffs = new double[maxPow + 1];

        for (int i = 0; i < parts.length; i++) {
            String term = parts[i];
            if (term.length() == 0) continue;

            int sign = 1;
            if (term.charAt(0) == '+') {
                term = term.substring(1);
            } else if (term.charAt(0) == '-') {
                sign = -1;
                term = term.substring(1);
            }

            double c;
            int pow;

            if (term.contains("x")) {
                int xIndex = term.indexOf('x');
                String coeffStr = term.substring(0, xIndex);
                if (coeffStr.length() == 0) {
                    coeffStr = "1.0";
                }
                c = Double.parseDouble(coeffStr) * sign;

                if (term.contains("^")) {
                    int powIndex = term.indexOf('^');
                    String powStr = term.substring(powIndex + 1);
                    pow = Integer.parseInt(powStr);
                } else {
                    pow = 1;
                }
            } else {
                c = Double.parseDouble(term) * sign;
                pow = 0;
            }

            coeffs[pow] += c;
        }

        ans = coeffs;

		return ans;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
     *
     * check if both polynomials are null
     * check if only one of the polynomials is null
     * find the length of the first polynomial
     * find the length of the second polynomial
     * compute the maximum length between the two
     * create a new array for the sum polynomial
     * loop over all positions up to the maximum length
     * take the coefficient from p1 if it exists
     * take the coefficient from p2 if it exists
     * add the two coefficients together
     * store the result in the sum array
     * return the final polynomial array
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        if (p1 == null && p2 == null) return ans;
        if (p1 == null) return p2;
        if (p2 == null) return p1;

        int n1 = p1.length;
        int n2 = p2.length;
        int n = Math.max(n1, n2);

        double[] sum = new double[n];

        for (int i = 0; i < n; i++) {
            double a = 0;
            double b = 0;

            if (i < n1) a = p1[i];
            if (i < n2) b = p2[i];

            sum[i] = a + b;
        }

        ans = sum;

		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
     *
     * check if one of the polynomials is null
     * get the length of the first polynomial
     * get the length of the second polynomial
     * create the result array with the correct size
     * loop over all coefficients of the first polynomial
     * loop over all coefficients of the second polynomial
     * multiply the two coefficients
     * add the product to the correct position in the result array
     * store the result as the answer
     * return the final polynomial
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//

        if (p1 == null || p2 == null) {
            return ans;
        }

        int n1 = p1.length;
        int n2 = p2.length;

        double[] res = new double[n1 + n2 - 1];

        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                res[i + j] += p1[i] * p2[j];
            }
        }

        ans = res;

		return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
     *
     * check if the polynomial is null or empty
     * check if the polynomial is constant (length 1)
     * create a new array for the derivative
     * loop from index 1 to the end
     * multiply each coefficient by its power
     * store the result in the new array
     * return the derivative polynomial
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//

        if (po == null || po.length == 0) {
            return ans;
        }

        if (po.length == 1) {
            return ZERO;
        }

        double[] d = new double[po.length - 1];

        for (int i = 1; i < po.length; i++) {
            d[i - 1] = po[i] * i;
        }

        ans = d;

		return ans;
	}
}
