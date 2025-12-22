# Table of Contents

### Table of Contents



- [Solution of Interpolation](#solution-of-interpolation)
  - [Newton's Forward Interpolation Method](#newtons-forward-interpolation-method)
    - [Theory](#newtons-forward-interpolation-theory)
      - [Introduction](#newtons-forward-interpolation-introduction)
      - [Formula](#newtons-forward-interpolation-formula)
      - [Algorithm Steps](#newtons-forward-interpolation-algorithm-steps)
      - [Application](#newtons-forward-interpolation-application)
    - [Code](#newtons-forward-interpolation-code)
    - [Input](#newtons-forward-interpolation-input)
    - [Output](#newtons-forward-interpolation-output)
  - [Newton's Backward Interpolation Method](#newtons-backward-interpolation-method)
    - [Theory](#newtons-backward-interpolation-theory)
      - [Introduction](#newtons-backward-interpolation-introduction)
      - [Formula](#newtons-backward-interpolation-formula)
      - [Algorithm Steps](#newtons-backward-interpolation-algorithm-steps)
      - [Application](#newtons-backward-interpolation-application)
    - [Code](#newtons-backward-interpolation-code)
    - [Input](#newtons-backward-interpolation-input)
    - [Output](#newtons-backward-interpolation-output)
  - [Divided Difference Method](#divided-difference-method)
    - [Theory](#divided-difference-theory)
      - [Introduction](#divided-difference-introduction)
      - [Formula](#divided-difference-formula)
      - [Algorithm Steps](#divided-difference-steps)
      - [Application](#divided-difference-application)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)

- [Solution of Curve Fitting Model](#solution-of-curve-fitting-model)
  - [Least Square Regression Method For Linear Equations](#least-square-regression-method-for-linear-equations)
    - [Theory](#least-square-regression-method-for-linear-equations-theory)
      - [Introduction](#least-square-regression-method-for-linear-equations-introduction)
      - [Formula](#least-square-regression-method-for-linear-equations-formula)
      - [Algorithm Steps](#least-square-regression-method-for-linear-equations-steps)
      - [Application](#least-square-regression-method-for-linear-equations-application)
    - [Code](#least-square-regression-method-for-linear-equations-code)
    - [Input](#least-square-regression-method-for-linear-equations-input)
    - [Output](#least-square-regression-method-for-linear-equations-output)
  - [Least Square Regression Method For Transcendental Equations](#least-square-regression-method-for-transcendental-equations)
    - [Theory](#least-square-regression-method-for-transcendental-equations-theory)
      - [Introduction](#least-square-regression-method-for-transcendental-equations-introduction)
      - [Formula](#least-square-regression-method-for-transcendental-equations-formula)
      - [Algorithm Steps](#least-square-regression-method-for-transcendental-equations-steps)
      - [Application](#least-square-regression-method-for-transcendental-equations-application)
    - [Code](#least-square-regression-method-for-transcendental-equations-code)
    - [Input](#least-square-regression-method-for-transcendental-equations-input)
    - [Output](#least-square-regression-method-for-transcendental-equations-output)
  - [Least Square Regression Method For Polynomial Equations](#least-square-regression-method-for-polynomial-equations)
    - [Theory](#least-square-regression-method-for-polynomial-equations-theory)
      - [Introduction](#least-square-regression-method-for-polynomial-equations-introduction)
      - [Formula](#least-square-regression-method-for-polynomial-equations-formula)
      - [Algorithm Steps](#least-square-regression-method-for-polynomial-equations-steps)
      - [Application](#least-square-regression-method-for-polynomial-equations-application)
    - [Code](#least-square-regression-method-for-polynomial-equations-code)
    - [Input](#least-square-regression-method-for-polynomial-equations-input)
    - [Output](#least-square-regression-method-for-polynomial-equations-output)


---




### Solution of Interpolation

### Newton's Forward Interpolation Method

#### Newton's Forward Interpolation Theory
##### Newton's Forward Interpolation Introduction

Newton’s Forward Interpolation method is a numerical technique used to estimate the value of a function at a point when the independent variable values are equally spaced. It is particularly effective when the required value lies near the beginning of the data table.
The method uses forward differences of the function values to construct the interpolation polynomial.


##### Newton's Forward Interpolation Formula
```
Interpolation formula:

y = y0 + p*Δy0 + (p*(p-1)/2!)*Δ²y0 + (p*(p-1)*(p-2)/3!)*Δ³y0 + ...

where:

p = (x - x0) / h

Explanation of terms:

- y0   : The value of the function at the first data point.
- Δy0, Δ²y0, Δ³y0, ... : Forward differences of the function values.
- h    : Spacing between consecutive values of x.
- x    : The point at which interpolation is required.
- p    : Ratio of distance from the first data point in units of h.

This formula estimates y at x by successively adding terms involving forward differences,
where each term accounts for higher-order variations of the data.
```


##### Newton's Forward Interpolation Algorithm Steps

1. Arrange the data:
   Ensure that the independent variable values x0, x1, ..., xn are equally spaced.

2. Construct a forward difference table:
   Compute the forward differences Δy0, Δ²y0, Δ³y0, ... from the given data.

3. Calculate p:
   p = (x - x0) / h
   where x is the value at which interpolation is required
   and h is the spacing between consecutive x values.

4. Apply the interpolation formula:
   y = y0 + p*Δy0 + (p*(p-1)/2!)*Δ²y0 + (p*(p-1)*(p-2)/3!)*Δ³y0 + ...
   Substitute y0, Δy0, Δ²y0, ... and p to compute the interpolated value of y.

5. Evaluate higher-order terms if necessary:
   Include as many terms as needed for the desired accuracy.



##### Newton's Forward Interpolation Application
```
[Add your output format here]
```

#### Newton's Forward Interpolation Code

```python
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    vector<vector<double>> diff(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) diff[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            diff[i][j] = diff[i+1][j-1] - diff[i][j-1];

    double value;
    fin >> value;
    double h = x[1] - x[0];
    double u = (value - x[0]) / h;
    double result = y[0];
    double term = 1.0;

    for (int i = 1; i < n; i++) {
        term *= (u - (i - 1));
        result += (term * diff[0][i]) / tgamma(i + 1);
    }

    fout << fixed << setprecision(6) << result << endl;
    return 0;
}
```
#### Newton's Forward Interpolation Input
```
[Add your output format here]
```
#### Newton's Forward Interpolation Output
```
[Add your output format here]
```




### Newton's Backward Interpolation Method

#### Newton's Backward Interpolation Theory
##### Newton's Backward Interpolation Introduction

Newton’s Backward Interpolation method is used when the data points are equally spaced
and the required value lies near the end of the data table.
This method uses backward differences to estimate the value of the function.



##### Newton's Backward Interpolation Formula
```
y = yn + p*∇y_n + (p*(p+1)/2!)*∇²y_n + (p*(p+1)*(p+2)/3!)*∇³y_n + ...

where:

p = (x - xn) / h

Explanation of terms:

- yn    : The value of the function at the last data point.
- ∇y_n, ∇²y_n, ∇³y_n, ... : Backward differences of the function values.
- h     : Spacing between consecutive values of x.
- x     : The point at which interpolation is required.
- p     : Ratio of distance from the last data point in units of h.

```

##### Newton's Backward Interpolation Algorithm Steps
```
1. Arrange the data:
   Ensure that the independent variable values x0, x1, ..., xn are equally spaced.

2. Construct a backward difference table:
   Compute the backward differences ∇y_n, ∇²y_n, ∇³y_n, ... from the given data.

3. Calculate p:
   p = (x - xn) / h
   where x is the value at which interpolation is required
   and h is the spacing between consecutive x values.

4. Apply the interpolation formula:
   y = yn + p*∇y_n + (p*(p+1)/2!)*∇²y_n + (p*(p+1)*(p+2)/3!)*∇³y_n + ...
   Substitute yn, ∇y_n, ∇²y_n, ... and p to compute the interpolated value of y.

5. Evaluate higher-order terms if necessary:
   Include as many terms as needed for the desired accuracy.

```
##### Newton's Backward Interpolation Application
```
[Add your output format here]
```

#### Newton's Backward Interpolation Code
```python
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    vector<vector<double>> diff(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) diff[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = j; i < n; i++)
            diff[i][j] = diff[i][j-1] - diff[i-1][j-1];

    double value;
    fin >> value;
    double h = x[1] - x[0];
    double u = (value - x[n-1]) / h;
    double result = y[n-1];
    double term = 1.0;

    for (int i = 1; i < n; i++) {
        term *= (u + (i - 1));
        result += (term * diff[n-1][i]) / tgamma(i + 1);
    }

    fout << fixed << setprecision(6) << result << endl;
    return 0;
}
```

#### Newton's Backward Interpolation Input
```
[Add your output format here]
```
#### Newton's Backward Interpolation Output
```
[Add your output format here]
```




### divided difference  Method

#### divided difference Theory
##### divided difference Introduction
Newton’s Divided Difference Interpolation is used when the data points are unequally spaced.
It constructs an interpolation polynomial that passes through all the given data points.
This method uses divided differences to compute the coefficients of the polynomial.


##### divided difference Formula
```
P(x) = y0 
       + (x - x0) * f[x0, x1] 
       + (x - x0)*(x - x1) * f[x0, x1, x2] 
       + (x - x0)*(x - x1)*(x - x2) * f[x0, x1, x2, x3] 
       + ...

where:

- y0                  : The value of the function at the first data point.
- f[x0, x1], f[x0,x1,x2], ... : First, second, and higher-order divided differences.
- x                    : The point at which interpolation is required.

```

##### divided difference Steps
```
1. Arrange the data:
   List the data points (x0, y0), (x1, y1), ..., (xn, yn).

2. Construct a divided difference table:
   - Compute first-order divided differences: f[xi, xi+1] = (yi+1 - yi) / (xi+1 - xi)
   - Compute second-order divided differences: f[xi, xi+1, xi+2] = (f[xi+1, xi+2] - f[xi, xi+1]) / (xi+2 - xi)
   - Continue calculating higher-order divided differences as needed.

3. Form the interpolation polynomial:
   P(x) = y0 + (x - x0)*f[x0,x1] + (x - x0)*(x - x1)*f[x0,x1,x2] + ...

4. Evaluate P(x) at the required value of x:
   Substitute the computed divided differences and the value of x to get the interpolated value of y.

```

##### divided difference Application
```
[Add your output format here]
```

#### divided difference Code
```python
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    vector<vector<double>> d(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) d[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            d[i][j] = (d[i+1][j-1] - d[i][j-1]) / (x[i+j] - x[i]);

    double value;
    fin >> value;
    double result = d[0][0];
    double term = 1.0;

    for (int i = 1; i < n; i++) {
        term *= (value - x[i-1]);
        result += d[0][i] * term;
    }

    fout << fixed << setprecision(6) << result << endl;
    return 0;
}
```

#### divided difference Input
```
[Add your output format here]
```
#### divided difference Output
```
[Add your output format here]
```



### solution-of-curve-fitting-model

### least square regression method for linear equations
#### least square regression method for linear equations theory

##### least square regression method for linear equations Introduction
Linear regression is used to model the relationship between a dependent variable \(y\) and an independent variable \(x\) by fitting a straight line to the observed data.  
```
The line is represented as:

\[
y = a + bx
\]

Where:  
- \(a\) = y-intercept  
- \(b\) = slope of the line  
```
The objective is to find the best-fit line that minimizes the sum of squared errors between observed and predicted values.

---

#####  least square regression method for linear equations Formula 
```
The coefficients \(a\) and \(b\) are obtained from the **normal equations**:

\[
\sum y = n a + b \sum x
\]

\[
\sum xy = a \sum x + b \sum x^2
\]

Where:  
- \(n\) = number of data points  
- \(\sum x\) = sum of x-values  
- \(\sum y\) = sum of y-values  
- \(\sum xy\) = sum of product of x and y  
- \(\sum x^2\) = sum of squares of x-values  

Solving these equations gives the values of \(a\) and \(b\).
```
---


#####  least square regression method for linear equations steps
1. Input data points \((x_i, y_i)\) for \(i = 1, 2, ..., n\).  
2. Compute the sums: \(\sum x\), \(\sum y\), \(\sum xy\), \(\sum x^2\).  
3. Form the normal equations:  
   - \(\sum y = n a + b \sum x\)  
   - \(\sum xy = a \sum x + b \sum x^2\)  
4. Solve the equations simultaneously to find \(a\) and \(b\).  
5. The best-fit line is \(y = a + bx\).  
6. Optionally, calculate predicted \(y_i\) for each \(x_i\) to compare with observed values.

---
#####  least square regression method for linear equations application
```
Used For

Fitting a straight line

𝑦=𝑎+𝑏𝑥
y=a+bx

to experimental or observed data
```
####  least square regression method for linear equations Code
```python
#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    double sumx = 0, sumy = 0, sumxy = 0, sumx2 = 0;

    for (int i = 0; i < n; i++) {
        sumx += x[i];
        sumy += y[i];
        sumxy += x[i] * y[i];
        sumx2 += x[i] * x[i];
    }

    double a = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
    double b = (sumy - a * sumx) / n;

    fout << fixed << setprecision(2);
    fout << "Linear Regression Line: y = " << a <<  " x + " << b << "\n";

    fin.close();
    fout.close();

    return 0;
}
```

####  least square regression method for linear equations Input
```
[Add your output format here]
```
####  least square regression method for linear equations Output
```
[Add your output format here]
```







### least square regression method for transcendental equations
#### least square regression method for transcendental equations theory

##### least square regression method for transcendental equations Introduction
Transcendental equations are equations involving transcendental functions such as exponential, logarithmic, trigonometric, or combinations of these, e.g.,  
```
\[
f(x) = e^x - 3x = 0
\]
```
These equations cannot be solved analytically in most cases. Numerical methods are used to approximate the roots of the equation. Common methods include:  
- Bisection Method  
- False Position (Regula Falsi) Method  
- Newton-Raphson Method  
- Secant Method  

---

#####  least square regression method for transcendental equations Formula 
```
y=ae^bx
y=ax^b
equations are similar to linear equaton
```
#####  least square regression method for transcendental equations steps
1. Select an initial guess or interval depending on the method.  
2. Evaluate the function \(f(x)\) at required points.  
3. Apply the iterative formula specific to the chosen method.  
4. Check for convergence:  
   \(|x_{n+1} - x_n| < \text{tolerance}\) or \(|f(x_{n+1})| < \text{tolerance}\).  
5. Repeat steps 2–4 until convergence.  
6. Report \(x_{root}\) as the approximate solution.

---

#####  least square regression method for transcendental equations application

####  least square regression method for transcendental equations Code
```python
#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n;
    fin >> n;

    vector<double> t(n), T(n);

    for (int i = 0; i < n; i++) fin >> t[i];
    for (int i = 0; i < n; i++) fin >> T[i];

    // Compute f(t_i) = exp(t_i / 4)
    vector<double> f(n);

    double sumf = 0, sumf2 = 0, sumT = 0, sumfT = 0;

    for (int i = 0; i < n; i++) {
        f[i] = exp(t[i] / 4.0);
        sumf += f[i];
        sumf2 += f[i] * f[i];
        sumT += T[i];
        sumfT += f[i] * T[i];
    }

    double denom = n * sumf2 - sumf * sumf;

    if (fabs(denom) < 1e-15) {
        fout << "Error: Singular matrix. Cannot compute regression.\n";
        return 0;
    }

    // Least squares formulas
    double b = (n * sumfT - sumf * sumT) / denom;
    double a = (sumT - b * sumf) / n;

    fout << fixed << setprecision(10);
    fout << "Computed parameters:\n";
    fout << "a = " << a << "\n";
    fout << "b = " << b << "\n";

    // Prediction
    double t_predict;
    fin >> t_predict;

    double f_predict = exp(t_predict / 4.0);
    double T_predict = a + b * f_predict;

    fout << "\nEstimated T(" << t_predict << ") = " << T_predict << "\n";

    fin.close();
    fout.close();

    return 0;
}
```

####  least square regression method for transcendental equations Input
```
[Add your output format here]
```
####  least square regression method for transcendental equations Output
```
[Add your output format here]
```
















### least square regression method for Polynomial equations
#### least square regression method for Polynomial equations theory

##### least square regression method for Polynomial equations – Introduction
Polynomial regression is a statistical method used to fit a higher-degree polynomial to a set of data points. Unlike linear regression, which fits a straight line, polynomial regression can capture curvature in the data.  

The assumed form of the polynomial is:

y = a0 + a1*x + a2*x^2 + a3*x^3 + ... + an*x^n

Where:  
- a0, a1, a2, ..., an are the coefficients of the polynomial  
- n is the degree of the polynomial  

The coefficients are determined such that the sum of the squares of the differences between the observed values and the predicted values is minimized.

---

##### least square regression method for Polynomial equations – Formula
```
For a polynomial of degree n, the **normal equations** are:

Σy = n*a0 + a1*Σx + a2*Σx^2 + ... + an*Σx^n

Σ(x*y) = a0*Σx + a1*Σx^2 + a2*Σx^3 + ... + an*Σx^(n+1)

Σ(x^2*y) = a0*Σx^2 + a1*Σx^3 + a2*Σx^4 + ... + an*Σx^(n+2)

...

Σ(x^n*y) = a0*Σx^n + a1*Σx^(n+1) + a2*Σx^(n+2) + ... + an*Σx^(2n)

Solve these equations simultaneously to determine the coefficients a0, a1, ..., an.
```
---

##### least square regression method for Polynomial equations – steps
1. Collect the data points (x, y).  
2. Decide the degree n of the polynomial.  
3. Compute the required sums: Σx, Σx^2, ..., Σx^(2n), Σy, Σ(x*y), Σ(x^2*y), ..., Σ(x^n*y).  
4. Form the normal equations:

Σy = n*a0 + a1*Σx + a2*Σx^2 + ... + an*Σx^n  
Σ(x*y) = a0*Σx + a1*Σx^2 + a2*Σx^3 + ... + an*Σx^(n+1)  
Σ(x^2*y) = a0*Σx^2 + a1*Σx^3 + a2*Σx^4 + ... + an*Σx^(n+2)  
...  
Σ(x^n*y) = a0*Σx^n + a1*Σx^(n+1) + a2*Σx^(n+2) + ... + an*Σx^(2n)  

5. Solve the normal equations simultaneously to find the coefficients a0, a1, ..., an.  
6. Construct the regression polynomial:

y = a0 + a1*x + a2*x^2 + ... + an*x^n  

7. Use the polynomial to predict y for any given x.
##### least-square-regression-method-for-polynomial-equations-application

 #### least square regression method for Polynomial equations Code
```python
#include <bits/stdc++.h>
using namespace std;

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++)
            if (abs(A[j][i]) > abs(A[pivot][i]))
                pivot = j;

        swap(A[i], A[pivot]);
        swap(b[i], b[pivot]);

        double factor = A[i][i];
        for (int j = i; j < n; j++) A[i][j] /= factor;
        b[i] /= factor;

        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double f = A[j][i];
            for (int k = i; k < n; k++)
                A[j][k] -= f * A[i][k];
            b[j] -= f * b[i];
        }
    }
    return b;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int n, degree;
    fin >> n >> degree;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
        fin >> x[i];
    for (int i = 0; i < n; i++)
         fin >> y[i];

    vector<vector<double>> A(degree + 1, vector<double>(degree + 1, 0));
    vector<double> B(degree + 1, 0);

    for (int i = 0; i <= degree; i++) {
        for (int j = 0; j <= degree; j++)
            for (int k = 0; k < n; k++)
                A[i][j] += pow(x[k], i + j);

        for (int k = 0; k < n; k++)
            B[i] += y[k] * pow(x[k], i);
    }

    vector<double> coeff = gaussianElimination(A, B);

    fout << fixed << setprecision(6);
    for (int i = 0; i <= degree; i++)
        fout << "a" << i << " = " << coeff[i] << "\n";

    return 0;
}
```

#### least square regression method for Polynomial equations Input
```
[Add your output format here]
```
#### least square regression method for Polynomial equations Output
```
[Add your output format here]
```
