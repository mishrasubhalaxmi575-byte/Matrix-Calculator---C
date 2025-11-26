# 🧮 Advanced Matrix Calculator in C

A *unique, feature-rich, and educational matrix calculator* built entirely in C. This project is far more advanced than typical matrix calculators found online. It includes:

* Color-coded heatmap matrix printing
* Detailed step-by-step determinant solving (with minors)
* Gauss-Jordan inverse with full step-by-step augmented matrix
* Sparse matrix detection
* File input/output
* Menu-driven UI with stored matrices

---

## ⭐ Features

### 🔵 1. Color-Coded Heatmap Output

Every value in the matrix is printed with ANSI colors:

* 🔴 *Red* – very large values
* 🟡 *Yellow* – moderate values
* 🟢 *Green* – mid-range
* 🔵 *Cyan* – small values

This makes matrix patterns easy to understand visually.

---

### 🔷 2. Step-by-Step Determinant Calculation

This calculator prints:

* Minor matrices
* Cofactor signs
* Cofactor values
* Recursive expansion along the first row

Perfect for learning mathematics.

---

### 🔶 3. Step-by-Step Matrix Inverse (Gauss-Jordan)

The inverse is computed using *Gauss–Jordan elimination*, and after each step you will see:

* Row swaps
* Normalized pivot rows
* Row elimination steps
* Full augmented matrix after every operation

This is extremely helpful for students.

---

### 🔶 4. Sparse Matrix Detection

If a matrix has *>60% zeros*, the program warns:


[Info] A seems sparse (zeros > 60%).


---

### 📂 5. File Loading and Saving

* Load matrix from matrix_input.txt
* Save results to result.txt
* Useful for large matrices

---

### 🟣 6. Menu-Driven Design

* Store matrices A and B
* Perform operations anytime
* Keep last result in memory

---

## 📌 Available Operations

* Addition
* Subtraction
* Multiplication
* Transpose
* Determinant
* Inverse (Gauss-Jordan)
* File load & save

---

## 📁 Project Structure


📦 Matrix-Calculator

 ┣ 📜 matrix_calculator.c
 
 ┣ 📜 README.md
 
 ┣ 📜 output.png  
 
 ┗ 📜 matrix_input.txt 


---

## ▶ How to Compile and Run

### *Linux / MacOS*

bash
gcc matrix_calculator.c -o matrix_calc
./matrix_calc


### *Windows (MinGW)*

bash
gcc matrix_calculator.c -o matrix_calc.exe
matrix_calc.exe


---

## 💾 Example: matrix_input.txt


3 3 

1 2 3

4 5 6

7 8 9


---

## 📸 Adding Screenshot to GitHub

1. Take screenshot of program output
2. Save as output.png
3. Upload to GitHub
4. Add this line to README:


![Program Output](output.png)


---

## 🧠 Future Enhancements (Ideas)

* LU Decomposition
* Rank of matrix
* Eigenvalues & eigenvectors
* Save session history
* Export result to CSV

---

## 🤝 Contributions

Pull requests are welcome! Improve features, add operations, or optimize algorithms.

---

## 📜 License

This project is open-source. You may modify and use it freely.

---

### ⭐ Special Note

This project is designed to be *unique, educational, and visually informative*. Perfect for GitHub portfolios and academic submissions.

---

