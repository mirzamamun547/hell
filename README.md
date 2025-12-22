# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)

---

### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory
this is gauss eli

#### Gauss Elimination Code
-[code](Solution Of Linear equations\Gauss Elimination\)

#### Gauss Elimination Input
```
[Add your input format here]
```

#### Gauss Elimination Output
```
[Add your output format here]
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
# Polynomial Regression

## Polynomial Regression – Introduction
Polynomial regression is a statistical method used to fit a higher-degree polynomial to a set of data points. Unlike linear regression, which fits a straight line, polynomial regression can capture curvature in the data.  

The assumed form of the polynomial is:

y = a0 + a1*x + a2*x^2 + a3*x^3 + ... + an*x^n

Where:  
- a0, a1, a2, ..., an are the coefficients of the polynomial  
- n is the degree of the polynomial  

The coefficients are determined such that the sum of the squares of the differences between the observed values and the predicted values is minimized.

---

## Polynomial Regression – Formula / Concept
For a polynomial of degree n, the **normal equations** are:

Σy = n*a0 + a1*Σx + a2*Σx^2 + ... + an*Σx^n

Σ(x*y) = a0*Σx + a1*Σx^2 + a2*Σx^3 + ... + an*Σx^(n+1)

Σ(x^2*y) = a0*Σx^2 + a1*Σx^3 + a2*Σx^4 + ... + an*Σx^(n+2)

...

Σ(x^n*y) = a0*Σx^n + a1*Σx^(n+1) + a2*Σx^(n+2) + ... + an*Σx^(2n)

Solve these equations simultaneously to determine the coefficients a0, a1, ..., an.

---

## Polynomial Regression – Procedure
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


#### Gauss Jordan Code
```python
# Add your code here
```

#### Gauss Jordan Input
```
[Add your input format here]
```

#### Gauss Jordan Output
```
[Add your output format here]
```

---

### LU Decomposition Method

#### LU Decomposition Theory
[Add your theory content here]

#### LU Decomposition Code
```python
# Add your code here
```

#### LU Decomposition Input
```
[Add your input format here]
```

#### LU Decomposition Output
```
[Add your output format here]
```

---

### Matrix Inversion

#### Matrix Inversion Theory
[Add your theory content here]

#### Matrix Inversion Code
```python
# Add your code here
```

#### Matrix Inversion Input
```
[Add your input format here]
```

#### Matrix Inversion Output
```
[Add your output format here]
```

---

### Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory
[Add your theory content here]

#### Bisection Code
```python
# Add your code here
```

#### Bisection Input
```
[Add your input format here]
```

#### Bisection Output
```
[Add your output format here]
```

---

### False Position Method

#### False Position Theory
[Add your theory content here]

#### False Position Code
```python
# Add your code here
```

#### False Position Input
```
[Add your input format here]
```

#### False Position Output
```
[Add your output format here]
```

---

