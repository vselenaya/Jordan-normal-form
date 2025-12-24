# Jordan Normal Form (JNF) Computation

> [!TIP]
> 🇷🇺 **Русская версия:** [Читать на русском](README.ru.md)

This repository contains an implementation of the algorithm for constructing the **Jordan Normal Form** (JNF) of a matrix and the corresponding **Jordan basis** (transition matrix).

The project offers two implementation variants: for real numbers (with detailed comments) and for finite fields (exact arithmetic).

## 📂 Repository Structure

There are two main files in the project:

### 1. `with-comments.cpp` (Base Version)

* **Data Type:** `long double` (real numbers).
* **Features:**
    * Contains detailed theoretical notes and comments for every step of the algorithm.
    * Ideal for educational purposes and understanding the mechanics of root subspaces and Jordan chains.
    * Uses `EPSILON` to compensate for floating-point errors.

### 2. `with-check-and-module.cpp` (Extended Version)

* **Data Type:** Supports operations over residue fields (modulo a prime number, configurable at the beginning of the file).
* **Features:**
    * **Exact Arithmetic:** Completely eliminates rounding errors inherent to `float/double`. Operations (division, multiplication) are performed using modular inverse.
    * **Automatic Verification:** Includes a result verification block. The program multiplies the resulting matrices to verify the equality $A \cdot S = S \cdot J$.

## 🚀 Usage

The algorithm requires two types of input:

1.  **Square Matrix** of size $n \times n$.
2.  **Matrix Spectrum** (list of eigenvalues), including their algebraic multiplicities.

> **Note:** The algorithm focuses on linear algebra (constructing chains), so finding the roots of the characteristic polynomial (eigenvalues) must be done beforehand.

### Input Example

Consider the matrix $A$:

$$
A = \begin{pmatrix}
3 & 1 & -4 & -7 \\
-1 & 1 & 5 & 9 \\
0 & 0 & 4 & 4 \\
0 & 0 & -1 & 0
\end{pmatrix}
$$

The characteristic polynomial is $(\lambda - 2)^4$.
Therefore, the eigenvalue $\lambda = 2$ has an algebraic multiplicity of 4.

**Program Input:**

```text
4
3 1 -4 -7
-1 1 5 9
0 0 4 4
0 0 -1 0
2 2 2 2
```

## Build and Run
Use g++ for compilation:

```Bash
# For the commented version (float)
g++ with-comments.cpp -o jnf_float
./jnf_float

# For the modular arithmetic version
g++ with-check-and-module.cpp -o jnf_mod
./jnf_mod
```
