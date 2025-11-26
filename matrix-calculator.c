#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double dtype;

/* ANSI color helpers for terminal "heatmap" */
#define COLOR_RESET  "\033[0m"
#define COLOR_RED    "\033[31m"
#define COLOR_YELLOW "\033[33m"
#define COLOR_GREEN  "\033[32m"
#define COLOR_CYAN   "\033[36m"

/* Dynamic matrix allocation */
dtype** alloc_matrix(int r, int c) {
    dtype *m = malloc(r * sizeof(dtype));
    if (!m) { perror("malloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < r; ++i) {
        m[i] = malloc(c * sizeof(dtype));
        if (!m[i]) { perror("malloc"); exit(EXIT_FAILURE); }
    }
    return m;
}

void free_matrix(dtype **m, int r) {
    if (!m) return;
    for (int i = 0; i < r; ++i) free(m[i]);
    free(m);
}

void input_matrix(dtype **m, int r, int c) {
    printf("Enter %d x %d elements row-wise:\n", r, c);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            scanf("%lf", &m[i][j]);
}

/* Print single value with color based on magnitude */
void print_val_colored(dtype v) {
    dtype av = fabs(v);
    if (av > 100.0) printf(COLOR_RED);
    else if (av > 10.0) printf(COLOR_YELLOW);
    else if (av > 1.0) printf(COLOR_GREEN);
    else printf(COLOR_CYAN);
    printf("%10.4g", v);
    printf(COLOR_RESET " ");
}

void print_matrix(dtype **m, int r, int c) {
    if (!m) { printf("NULL matrix\n"); return; }
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j)
            print_val_colored(m[i][j]);
        putchar('\n');
    }
}

/* Print matrix without colors (used for augmented prints optionally) */
void print_matrix_plain(dtype **m, int r, int c) {
    if (!m) { printf("NULL\n"); return; }
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j)
            printf("%10.4g ", m[i][j]);
        putchar('\n');
    }
}

/* Copy matrix */
dtype** copy_matrix(dtype **src, int r, int c) {
    dtype **dst = alloc_matrix(r, c);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            dst[i][j] = src[i][j];
    return dst;
}

/* Addition & subtraction (same dimensions) */
dtype** add_matrix(dtype **a, dtype **b, int r, int c) {
    dtype **res = alloc_matrix(r,c);
    for (int i=0;i<r;++i) for (int j=0;j<c;++j) res[i][j]=a[i][j]+b[i][j];
    return res;
}

dtype** sub_matrix(dtype **a, dtype **b, int r, int c) {
    dtype **res = alloc_matrix(r,c);
    for (int i=0;i<r;++i) for (int j=0;j<c;++j) res[i][j]=a[i][j]-b[i][j];
    return res;
}

/* Multiplication */
dtype** mul_matrix(dtype **a, int ar, int ac, dtype **b, int br, int bc) {
    if (ac != br) return NULL;
    dtype **res = alloc_matrix(ar, bc);
    for (int i=0;i<ar;++i) for (int j=0;j<bc;++j) {
        res[i][j] = 0.0;
        for (int k=0;k<ac;++k) res[i][j] += a[i][k]*b[k][j];
    }
    return res;
}

/* Scalar multiplication */
dtype** scalar_mul(dtype **a, int r, int c, dtype k) {
    dtype **res = alloc_matrix(r,c);
    for (int i=0;i<r;++i) for (int j=0;j<c;++j) res[i][j] = a[i][j]*k;
    return res;
}

/* Transpose */
dtype** transpose(dtype **a, int r, int c) {
    dtype **t = alloc_matrix(c, r);
    for (int i=0;i<r;++i) for (int j=0;j<c;++j) t[j][i] = a[i][j];
    return t;
}

/* Minor matrix (for determinant) */
dtype** minor_matrix(dtype **a, int n, int row, int col) {
    dtype **m = alloc_matrix(n-1, n-1);
    int mi = 0, mj;
    for (int i=0;i<n;++i) {
        if (i==row) continue;
        mj = 0;
        for (int j=0;j<n;++j) {
            if (j==col) continue;
            m[mi][mj] = a[i][j];
            mj++;
        }
        mi++;
    }
    return m;
}

/* --- Step-by-step determinant with Laplace expansion (prints minors for first-level expansion) --- */
dtype determinant(dtype **a, int n) {
    if (n==1) return a[0][0];
    if (n==2) return a[0][0]*a[1][1] - a[0][1]*a[1][0];

    dtype det = 0.0;
    for (int c = 0; c < n; ++c) {
        dtype **m = minor_matrix(a, n, 0, c);
        dtype sub = determinant(m, n-1);
        dtype cofactor_sign = ((c%2==0)?1:-1);
        dtype cofactor = cofactor_sign * a[0][c] * sub;

        /* Print step for this term (only top-level) */
        printf("\nMinor removing row 0 col %d:\n", c);
        print_matrix(m, n-1, n-1);
        printf("Cofactor term: (%g) * %g = %g  (sign=%g, minor_det=%g)\n",
               a[0][c], cofactor_sign, cofactor, cofactor_sign, sub);

        det += cofactor;
        free_matrix(m, n-1);
    }
    return det;
}

/* Utility to print augmented matrix during Gauss-Jordan */
void print_augmented(dtype **aug, int n, int cols) {
    for (int i=0;i<n;++i) {
        for (int j=0;j<cols;++j) {
            printf("%10.4g ", aug[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

/* Inverse by Gauss-Jordan elimination with step prints */
/* Returns NULL if singular or not square */
dtype** inverse(dtype **a, int n) {
    dtype **aug = alloc_matrix(n, 2*n);
    for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) aug[i][j] = a[i][j];
        for (int j=n;j<2*n;++j) aug[i][j] = (j-n==i)?1.0:0.0;
    }

    printf("\nInitial augmented matrix [A|I]:\n");
    print_augmented(aug, n, 2*n);

    for (int col = 0; col < n; ++col) {
        /* Find pivot: row with max absolute value in col at or below 'col' */
        int pivot = col;
        dtype maxv = fabs(aug[col][col]);
        for (int r = col+1; r < n; ++r) {
            dtype val = fabs(aug[r][col]);
            if (val > maxv) { maxv = val; pivot = r; }
        }
        if (fabs(aug[pivot][col]) < 1e-12) { free_matrix(aug, n); return NULL; } /* singular */

        /* swap rows if needed */
        if (pivot != col) {
            dtype *tmp = aug[pivot];
            aug[pivot] = aug[col];
            aug[col] = tmp;
            printf("Swapped row %d with row %d:\n", pivot, col);
            print_augmented(aug, n, 2*n);
        }

        /* Normalize pivot row */
        dtype piv = aug[col][col];
        for (int j = 0; j < 2*n; ++j) aug[col][j] /= piv;
        printf("Normalized row %d (pivot -> 1):\n", col);
        print_augmented(aug, n, 2*n);

        /* Eliminate other rows */
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            dtype factor = aug[r][col];
            if (fabs(factor) > 0) {
                for (int j = 0; j < 2*n; ++j)
                    aug[r][j] -= factor * aug[col][j];
            }
        }
        printf("After eliminating column %d:\n", col);
        print_augmented(aug, n, 2*n);
    }

    /* Extract inverse from augmented */
    dtype **inv = alloc_matrix(n, n);
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) inv[i][j] = aug[i][j+n];

    free_matrix(aug, n);
    return inv;
}

/* Utility: ask for matrix dimensions and elements, returns pointer and sets r,c */
dtype** create_and_input(int *r, int *c) {
    printf("Rows: "); if (scanf("%d", r) != 1) return NULL;
    printf("Columns: "); if (scanf("%d", c) != 1) return NULL;
    if (*r <= 0 || *c <= 0) { printf("Invalid dimensions.\n"); return NULL; }
    dtype **m = alloc_matrix(*r, *c);
    input_matrix(m, *r, *c);
    return m;
}

/* Load matrix from file: format first line "r c" then r*c numbers row-wise */
dtype** load_matrix_from_file(const char *filename, int *r, int *c) {
    FILE *f = fopen(filename, "r");
    if (!f) { perror("fopen"); return NULL; }
    if (fscanf(f, "%d %d", r, c) != 2) { fclose(f); return NULL; }
    if (*r <=0 || *c <=0) { fclose(f); return NULL; }
    dtype **m = alloc_matrix(*r, *c);
    for (int i=0;i<*r;++i) for (int j=0;j<*c;++j) {
        if (fscanf(f, "%lf", &m[i][j]) != 1) { free_matrix(m, *r); fclose(f); return NULL; }
    }
    fclose(f);
    return m;
}

/* Save matrix to file */
int save_matrix_to_file(const char *filename, dtype **m, int r, int c) {
    FILE *f = fopen(filename, "w");
    if (!f) { perror("fopen"); return -1; }
    fprintf(f, "%d %d\n", r, c);
    for (int i=0;i<r;++i) {
        for (int j=0;j<c;++j) fprintf(f, "%0.10g ", m[i][j]);
        fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}

/* Count zeros and detect sparse */
int is_sparse(dtype **m, int r, int c) {
    int zeros = 0;
    for (int i=0;i<r;++i) for (int j=0;j<c;++j) if (fabs(m[i][j]) < 1e-12) zeros++;
    int total = r*c;
    double frac = (double)zeros / (double)total;
    if (frac > 0.6) return 1;
    return 0;
}

void print_menu() {
    puts("\n========== Matrix Calculator (Unique Version) ==========");
    puts("1. Add (A + B)");
    puts("2. Subtract (A - B)");
    puts("3. Multiply (A * B)");
    puts("4. Scalar multiply (k * A)");
    puts("5. Transpose (A^T)");
    puts("6. Determinant (det(A))  -- prints steps");
    puts("7. Inverse (A^-1)       -- prints gauss-jordan steps");
    puts("8. Create / Input a matrix (store as A or B)");
    puts("9. Print stored matrices (with color heatmap)");
    puts("10. Load matrix from file (matrix_input.txt)");
    puts("11. Save last result to file (result.txt)");
    puts("0. Exit");
    puts("Note: Stored matrices: A and B (you can override them).");
    puts("========================================================");
}

int main() {
    dtype **A = NULL, **B = NULL, **last_result = NULL;
    int Ar=0, Ac=0, Br=0, Bc=0, Lr=0, Lc=0;
    int choice;
    while (1) {
        print_menu();
        printf("Choice: ");
        if (scanf("%d", &choice) != 1) { printf("Bad input. Exiting.\n"); break; }
        if (choice == 0) break;

        if (choice == 8) {
            printf("Which matrix to input? (1 for A, 2 for B): ");
            int which; scanf("%d", &which);
            if (which == 1) {
                if (A) { free_matrix(A, Ar); A = NULL; }
                A = create_and_input(&Ar, &Ac);
                if (A && is_sparse(A, Ar, Ac)) printf("[Info] A seems sparse (zeros > 60%%).\n");
            } else if (which == 2) {
                if (B) { free_matrix(B, Br); B = NULL; }
                B = create_and_input(&Br, &Bc);
                if (B && is_sparse(B, Br, Bc)) printf("[Info] B seems sparse (zeros > 60%%).\n");
            } else printf("Invalid choice.\n");
            continue;
        }

        if (choice == 10) {
            printf("Load into which matrix? (1 for A, 2 for B): ");
            int which; scanf("%d", &which);
            const char *fname = "matrix_input.txt";
            if (which == 1) {
                if (A) { free_matrix(A, Ar); A = NULL; }
                A = load_matrix_from_file(fname, &Ar, &Ac);
                if (!A) printf("Failed to load A from %s\n", fname);
                else {
                    printf("Loaded A from %s (%d x %d)\n", fname, Ar, Ac);
                    print_matrix(A, Ar, Ac);
                }
            } else if (which == 2) {
                if (B) { free_matrix(B, Br); B = NULL; }
                B = load_matrix_from_file(fname, &Br, &Bc);
                if (!B) printf("Failed to load B from %s\n", fname);
                else {
                    printf("Loaded B from %s (%d x %d)\n", fname, Br, Bc);
                    print_matrix(B, Br, Bc);
                }
            } else printf("Invalid choice.\n");
            continue;
        }

        if (choice == 11) {
            if (!last_result) { printf("No last result to save.\n"); continue; }
            if (save_matrix_to_file("result.txt", last_result, Lr, Lc) == 0)
                printf("Saved last result to result.txt\n");
            continue;
        }

        if (choice == 9) {
            printf("\nMatrix A (%d x %d):\n", Ar, Ac);
            if (A) { print_matrix(A, Ar, Ac); if (is_sparse(A,Ar,Ac)) printf("[Note] A is sparse.\n"); }
            else printf("A is empty.\n");
            printf("\nMatrix B (%d x %d):\n", Br, Bc);
            if (B) { print_matrix(B, Br, Bc); if (is_sparse(B,Br,Bc)) printf("[Note] B is sparse.\n"); }
            else printf("B is empty.\n");
            continue;
        }

        if (choice == 1) {
            if (!A || !B) { printf("Both A and B must be defined.\n"); continue; }
            if (Ar!=Br || Ac!=Bc) { printf("Matrices must have same dimensions.\n"); continue; }
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = add_matrix(A,B,Ar,Ac); Lr = Ar; Lc = Ac;
            printf("Result (A + B):\n"); print_matrix(last_result, Lr, Lc);
        } else if (choice == 2) {
            if (!A || !B) { printf("Both A and B must be defined.\n"); continue; }
            if (Ar!=Br || Ac!=Bc) { printf("Matrices must have same dimensions.\n"); continue; }
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = sub_matrix(A,B,Ar,Ac); Lr = Ar; Lc = Ac;
            printf("Result (A - B):\n"); print_matrix(last_result, Lr, Lc);
        } else if (choice == 3) {
            if (!A || !B) { printf("Both A and B must be defined.\n"); continue; }
            if (Ac != Br) { printf("A's columns must equal B's rows for multiplication.\n"); continue; }
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = mul_matrix(A, Ar, Ac, B, Br, Bc); Lr = Ar; Lc = Bc;
            printf("Result (A * B):\n"); print_matrix(last_result, Lr, Lc);
        } else if (choice == 4) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            printf("Enter scalar k: ");
            dtype k; scanf("%lf", &k);
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = scalar_mul(A, Ar, Ac, k); Lr = Ar; Lc = Ac;
            printf("Result (k * A):\n"); print_matrix(last_result, Lr, Lc);
        } else if (choice == 5) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = transpose(A, Ar, Ac); Lr = Ac; Lc = Ar;
            printf("Transpose of A:\n"); print_matrix(last_result, Lr, Lc);
        } else if (choice == 6) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            if (Ar != Ac) { printf("Determinant requires a square matrix.\n"); continue; }
            printf("Computing determinant with step-by-step expansion (top-level minors shown)...\n");
            dtype det = determinant(A, Ar);
            printf("\ndet(A) = %.10g\n", det);
            /* Save scalar result into last_result as 1x1 matrix for easy saving */
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = alloc_matrix(1,1); last_result[0][0] = det; Lr = 1; Lc = 1;
        } else if (choice == 7) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            if (Ar != Ac) { printf("Inverse requires a square matrix.\n"); continue; }
            dtype det = determinant(A, Ar);
            if (fabs(det) < 1e-12) { printf("Matrix is singular (det=0). No inverse.\n"); continue; }
            printf("Computing inverse using Gauss-Jordan with step outputs...\n");
            if (last_result) { free_matrix(last_result, Lr); last_result = NULL; }
            last_result = inverse(A, Ar);
            if (!last_result) { printf("Inverse computation failed (matrix may be singular).\n"); continue; }
            Lr = Ar; Lc = Ac;
            printf("Inverse of A:\n"); print_matrix(last_result, Lr, Lc);
        } else {
            printf("Invalid choice.\n");
        }
    }

    if (A) free_matrix(A, Ar);
    if (B) free_matrix(B, Br);
    if (last_result) free_matrix(last_result, Lr);
    printf("Goodbye!\n");
    return 0;
}
