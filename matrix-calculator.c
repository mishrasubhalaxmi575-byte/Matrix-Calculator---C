

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double dtype;

/* Dynamic matrix allocation */
dtype** alloc_matrix(int r, int c) {
    dtype **m = malloc(r * sizeof(dtype*));
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

void print_matrix(dtype **m, int r, int c) {
    if (!m) { printf("NULL matrix\n"); return; }
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

/* Determinant (recursive Laplace) */
dtype determinant(dtype **a, int n) {
    if (n==1) return a[0][0];
    if (n==2) return a[0][0]*a[1][1] - a[0][1]*a[1][0];
    dtype det = 0.0;
    for (int c = 0; c < n; ++c) {
        dtype **m = minor_matrix(a, n, 0, c);
        dtype sub = determinant(m, n-1);
        dtype cofactor = ((c%2==0)?1:-1) * a[0][c] * sub;
        det += cofactor;
        free_matrix(m, n-1);
    }
    return det;
}

/* Inverse by Gauss-Jordan elimination */
/* Returns NULL if singular or not square */
dtype** inverse(dtype **a, int n) {
    /* create augmented matrix [A | I] of size n x 2n */
    dtype **aug = alloc_matrix(n, 2*n);
    for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) aug[i][j] = a[i][j];
        for (int j=n;j<2*n;++j) aug[i][j] = (j-n==i)?1.0:0.0;
    }

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
        }

        /* Normalize pivot row */
        dtype piv = aug[col][col];
        for (int j = 0; j < 2*n; ++j) aug[col][j] /= piv;

        /* Eliminate other rows */
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            dtype factor = aug[r][col];
            if (fabs(factor) > 0) {
                for (int j = 0; j < 2*n; ++j)
                    aug[r][j] -= factor * aug[col][j];
            }
        }
    }

    /* Extract inverse from augmented */
    dtype **inv = alloc_matrix(n, n);
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) inv[i][j] = aug[i][j+n];

    free_matrix(aug, n);
    return inv;
}

/* Utility: ask for matrix dimensions and elements, returns pointer and sets r,c */
dtype** create_and_input(int *r, int *c) {
    printf("Rows: "); scanf("%d", r);
    printf("Columns: "); scanf("%d", c);
    if (*r <= 0 || *c <= 0) { printf("Invalid dimensions.\n"); return NULL; }
    dtype **m = alloc_matrix(*r, *c);
    input_matrix(m, *r, *c);
    return m;
}

void pause_and_clear_input() {
    printf("Press ENTER to continue...");
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF) {}
    getchar();
}

/* Menu and main */
void print_menu() {
    puts("\n========== Matrix Calculator ==========");
    puts("1. Add (A + B)");
    puts("2. Subtract (A - B)");
    puts("3. Multiply (A * B)");
    puts("4. Scalar multiply (k * A)");
    puts("5. Transpose (A^T)");
    puts("6. Determinant (det(A))");
    puts("7. Inverse (A^-1)");
    puts("8. Create / Input a matrix (store as A or B)");
    puts("9. Print stored matrices");
    puts("0. Exit");
    puts("Note: Stored matrices: A and B (you can override them).");
    puts("========================================");
}

int main() {
    dtype **A = NULL, **B = NULL;
    int Ar=0, Ac=0, Br=0, Bc=0;
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
            } else if (which == 2) {
                if (B) { free_matrix(B, Br); B = NULL; }
                B = create_and_input(&Br, &Bc);
            } else printf("Invalid choice.\n");
            continue;
        }

        if (choice == 9) {
            printf("\nMatrix A (%d x %d):\n", Ar, Ac);
            if (A) print_matrix(A, Ar, Ac); else printf("A is empty.\n");
            printf("\nMatrix B (%d x %d):\n", Br, Bc);
            if (B) print_matrix(B, Br, Bc); else printf("B is empty.\n");
            continue;
        }

        if (choice == 1) {
            if (!A || !B) { printf("Both A and B must be defined.\n"); continue; }
            if (Ar!=Br || Ac!=Bc) { printf("Matrices must have same dimensions.\n"); continue; }
            dtype **r = add_matrix(A,B,Ar,Ac);
            printf("Result (A + B):\n"); print_matrix(r, Ar, Ac);
            free_matrix(r, Ar);
        } else if (choice == 2) {
            if (!A || !B) { printf("Both A and B must be defined.\n"); continue; }
            if (Ar!=Br || Ac!=Bc) { printf("Matrices must have same dimensions.\n"); continue; }
            dtype **r = sub_matrix(A,B,Ar,Ac);
            printf("Result (A - B):\n"); print_matrix(r, Ar, Ac);
            free_matrix(r, Ar);
        } else if (choice == 3) {
            if (!A || !B) { printf("Both A and B must be defined.\n"); continue; }
            if (Ac != Br) { printf("A's columns must equal B's rows for multiplication.\n"); continue; }
            dtype **r = mul_matrix(A, Ar, Ac, B, Br, Bc);
            printf("Result (A * B):\n"); print_matrix(r, Ar, Bc);
            free_matrix(r, Ar);
        } else if (choice == 4) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            printf("Enter scalar k: ");
            dtype k; scanf("%lf", &k);
            dtype **r = scalar_mul(A, Ar, Ac, k);
            printf("Result (k * A):\n"); print_matrix(r, Ar, Ac);
            free_matrix(r, Ar);
        } else if (choice == 5) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            dtype **t = transpose(A, Ar, Ac);
            printf("Transpose of A:\n"); print_matrix(t, Ac, Ar);
            free_matrix(t, Ac);
        } else if (choice == 6) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            if (Ar != Ac) { printf("Determinant requires a square matrix.\n"); continue; }
            dtype det = determinant(A, Ar);
            printf("det(A) = %.10g\n", det);
        } else if (choice == 7) {
            if (!A) { printf("Matrix A must be defined.\n"); continue; }
            if (Ar != Ac) { printf("Inverse requires a square matrix.\n"); continue; }
            dtype det = determinant(A, Ar);
            if (fabs(det) < 1e-12) { printf("Matrix is singular (det=0). No inverse.\n"); continue; }
            dtype **inv = inverse(A, Ar);
            if (!inv) { printf("Inverse computation failed (matrix may be singular).\n"); continue; }
            printf("Inverse of A:\n"); print_matrix(inv, Ar, Ac);
            free_matrix(inv, Ar);
        } else {
            printf("Invalid choice.\n");
        }
    }

    if (A) free_matrix(A, Ar);
    if (B) free_matrix(B, Br);
    printf("Goodbye!\n");
    return 0;
}
