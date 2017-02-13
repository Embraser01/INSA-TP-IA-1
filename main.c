/* PageRank */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* allocate one object of given type */
#define    NEW(type)    ((type*)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define    NEW_A(num, type)    ((type*)calloc((size_t)(num),(size_t)sizeof(type)))

typedef unsigned int u_int;
typedef double Real;

/* vector definition */
typedef struct
{
    u_int dim;
    Real *ve;
} VEC;

/* matrix definition */
typedef struct
{
    u_int m, n;
    Real **me;
} MAT;

/* sparse matrix definition */
typedef struct graphnode
{
    u_int col;
    Real val;
    struct graphnode *next;
} NODE;

typedef struct
{
    u_int m, n;
    NODE *rows;
} SMAT;

/* v_get -- gets a VEC of dimension 'dim'
   Precondition: size >= 0
   Postcondition: initialized to zero */
VEC *v_get(u_int size)
{
    VEC *v;

    if ((v = NEW(VEC)) == (VEC *) NULL)
    {
        fprintf(stderr, "v_get memory error");
        exit(-1);
    }

    v->dim = size;
    if ((v->ve = NEW_A(size, Real)) == (Real *) NULL)
    {
        free(v);
        fprintf(stderr, "v_get memory error");
        exit(-1);
    }

    return (v);
}

/* v_free -- returns VEC & associated memory back to memory heap */
int v_free(VEC *vec)
{
    if (vec == (VEC *) NULL)
        return (-1);

    if (vec->ve == (Real *) NULL)
    {
        free(vec);
    } else
    {
        free(vec->ve);
        free(vec);
    }

    return (0);
}

SMAT *sm_get(u_int m, u_int n)
{
    SMAT *G;

    if ((G = NEW(SMAT)) == (SMAT *) NULL)
    {
        fprintf(stderr, "sm_get memory error");
        exit(-1);
    }

    G->m = m;
    G->n = n;

    if ((G->rows = NEW_A(m, NODE)) == (NODE *) NULL)
    {
        free(G);
        fprintf(stderr, "sm_get memory error");
        exit(-1);
    }

    for (u_int i = 0; i < G->m; i++)
        (G->rows[i]).val = -1;

    return (G);
}

int sm_free(SMAT *G)
{
    if (G == (SMAT *) NULL)
        return (-1);

    if (G->rows == (NODE *) NULL)
    {
        free(G);
    } else
    {
        NODE *n0;
        NODE *n1;
        for (u_int i = 0; i < G->m; i++)
        {
            n0 = &(G->rows[i]);
            if (n0->val < 0.0) break; /* empty line */
            n0 = n0->next;
            while (n0->val >= 0.0)
            {
                n1 = n0->next;
                free(n0);
                n0 = n1;
            }
            free(n0);
        }
        free(G->rows);
        free(G);
    }

    return (0);
}

NODE *sm_add(NODE *n0, u_int c, Real v)
{
    NODE *n1;
    n0->col = c;
    n0->val = v;
    if ((n1 = NEW(NODE)) == (NODE *) NULL)
    {
        fprintf(stderr, "sm_add memory error");
        exit(-1);
    }
    n1->val = -1;
    n0->next = n1;
    return (n1);
}

/* m_get -- gets an mxn matrix by dynamic memory allocation
   Precondition: m>=0 && n>=0
   Postcondition: initialized to zero */
MAT *m_get(u_int m, u_int n)
{
    MAT *g;

    if ((g = NEW(MAT)) == (MAT *) NULL)
    {
        fprintf(stderr, "m_get memory error");
        exit(-1);
    }

    g->m = m;
    g->n = n;

    if ((g->me = NEW_A(m, Real*)) == (Real **) NULL)
    {
        free(g);
        fprintf(stderr, "m_get memory error");
        exit(-1);
    }

    for (int i = 0; i < m; i++)
        if ((g->me[i] = NEW_A(n, Real)) == (Real *) NULL)
        {
            fprintf(stderr, "m_get memory error");
            exit(-1);
        }

    return (g);
}

/* m_free -- returns MAT & associated memory back to memory heap */
int m_free(MAT *mat)
{
    if (mat == (MAT *) NULL)
        return (-1);

    for (int i = 0; i < mat->m; i++)
        if (mat->me[i] != (Real *) NULL) free(mat->me[i]);

    if (mat->me != (Real **) NULL) free(mat->me);

    free(mat);

    return (0);
}

/* m_input -- file input of matrix */
MAT *m_input(FILE *fp)
{
    MAT *g;
    u_int m, n, val;

    /* get dimension */
    if (fscanf(fp, " Matrix: %u by %u", &m, &n) < 2)
    {
        fprintf(stderr, "m_input error reading dimensions");
        exit(-1);
    }

    /* allocate memory if necessary */
    g = m_get(m, n);

    /* get entries */
    for (u_int i = 0; i < m; i++)
    {
        if (fscanf(fp, " row %u:", &val) < 1)
        {
            fprintf(stderr, "m_input error reading line %u", i);
            exit(-1);
        }
        for (u_int j = 0; j < n; j++)
            if (fscanf(fp, "%lf", &g->me[i][j]) < 1)
            {
                fprintf(stderr, "m_input error reading line %u col %u", i, j);
                exit(-1);
            }
    }

    return (g);
}

/* sm_input -- file input of sparse matrix */
SMAT *sm_input(FILE *fp)
{
    SMAT *g;
    u_int m, n, row;
    Real col;
    NODE *n0;

    /* get dimension */
    if (fscanf(fp, " SparseMatrix: %u by %u", &m, &n) < 2)
    {
        fprintf(stderr, "sm_input error reading dimensions");
        exit(-1);
    }

    g = sm_get(m, n);

    /* get entries */
    for (u_int i = 0; i < m; i++)
    {
        if (fscanf(fp, " row %u:", &row) < 1)
        {
            fprintf(stderr, "sm_input error reading line %u", i);
            exit(-1);
        }
        n0 = &(g->rows[i]);
        for (;;)
        {
            if (fscanf(fp, "%lf", &col) < 1)
            {
                fprintf(stderr, "sm_input error reading line %u col x", i);
                exit(-1);
            }
            if (col < 0.0) break;
            n0 = sm_add(n0, (u_int) col, 1.0);
        }
    }

    return (g);
}

static char *format = "%1.5g ";

void sm_output(FILE *fp, SMAT *G)
{
    NODE *n0;

    fprintf(fp, "SparseMatrix: %d by %d\n", G->m, G->n);
    for (u_int i = 0; i < G->m; i++)
    {
        fprintf(fp, "row %u: ", i);
        n0 = &(G->rows[i]);
        while (n0->val >= 0.0)
        {
            fprintf(fp, format, (Real) n0->col);
            n0 = n0->next;
        }
        fprintf(fp, "-1\n");
    }
}

/* m_output -- file output of matrix
   Precondition: Memory already allocated for the matrix */
void m_output(FILE *fp, MAT *g)
{
    u_int tmp;

    fprintf(fp, "Matrix: %d by %d\n", g->m, g->n);
    for (u_int i = 0; i < g->m; i++)
    {
        fprintf(fp, "row %u: ", i);

        tmp = 2;
        for (u_int j = 0; j < g->n; j++, tmp++)
        {
            fprintf(fp, format, g->me[i][j]);
            if (!(tmp % 9)) putc('\n', fp);
        }
        if (tmp % 9 != 1) putc('\n', fp);
    }
}

/* v_output -- file output of vector */
void v_output(FILE *fp, VEC *v)
{
    fprintf(fp, "Vector: %d\n", v->dim);
    for (u_int i = 0; i < v->dim; i++) fprintf(fp, format, v->ve[i]);
    putc('\n', fp);
}

/* m_cp -- copy matrix M in OUT
   Precondition: memory is already allocated for M and OUT
   Precondition: sizes of M and OUT must match*/
MAT *m_cp(MAT *M, MAT *OUT)
{
    for (u_int i = 0; i < M->m; i++)
        memmove(&(OUT->me[i][0]), &(M->me[i][0]), (M->n) * sizeof(Real));

    return (OUT);
}

/* v_cp -- copy vector v in out
   Precondition: memory is already allocated for v and out*/
VEC *v_cp(VEC *v, VEC *out)
{
    memmove(&(out->ve[0]), &(v->ve[0]), (v->dim) * sizeof(Real));
    return (out);
}

VEC *vm_mult(VEC *vec, MAT *mat)
{
    if (vec->dim != mat->m)
    {
        fprintf(stderr, "vm_mult vector and mat height must be the same size");
        exit(-1);
    }

    VEC *res = v_get(mat->n);

    Real acc;
    for (u_int i = 0; i < res->dim; ++i)
    {
        acc = 0;
        for (u_int j = 0; j < vec->dim; ++j)
        {
            acc += vec->ve[j] * mat->me[j][i];
        }

        res->ve[i] = acc;
    }

    return res;
}

MAT *m_mult(MAT *mat1, MAT *mat2)
{
    if (mat1->n != mat2->m)
    {
        fprintf(stderr, "m_mult two matrix need to be compatible");
        exit(-1);
    }

    MAT *res = m_get(mat1->m, mat2->n);

    Real acc;

    for (u_int i = 0; i < res->m; ++i)
    {
        for (u_int j = 0; j < res->n; ++j)
        {
            acc = 0;

            for (u_int k = 0; k < mat1->n; ++k)
            {
                acc += mat1->me[i][k] * mat2->me[k][j];
            }

            res->me[i][j] = acc;
        }
    }

    return res;
}


MAT *m_pow(MAT *mat, u_int n)
{
    if (mat->n != mat->m)
    {
        fprintf(stderr, "m_pow matrix need to be square");
        exit(-1);
    }

    MAT *res = m_get(mat->m, mat->n);

    if (n == 0)
    {
        for (u_int j = 0; j < mat->m; ++j)
        {
            res->me[j][j] = 1;
        }

        return res;
    }

    m_cp(mat, res);

    for (u_int i = 1; i < n; ++i)
    {
        res = m_mult(res, mat);
    }

    return res;
}

void main()
{
    u_int i, j;
    u_int n = 200;
    Real alpha = 0.9;

    MAT *G;

    FILE *fp;
    fp = fopen("data/g.dat", "r");
    G = m_input(fp);
    fclose(fp);

    VEC *succ;
    succ = v_get(G->m);

    for (i = 0; i < succ->dim; ++i)
    {
        Real acc = 0;

        for (j = 0; j < G->n; ++j)
        {
            acc += G->me[i][j];
        }

        succ->ve[i] = acc;
    }

    MAT *H = m_get(G->m, G->n);
    m_cp(G, H);

    Real newVal;
    for (i = 0; i < H->m; ++i)
    {
        newVal = 1.0 / succ->ve[i];

        for (j = 0; j < H->n; ++j)
        {
            if (H->me[i][j])
            {
                H->me[i][j] = newVal;
            } else if (succ->ve[i] == 0)
            {
                H->me[i][j] = 1.0 / H->n;
            }
        }
    }


    MAT *E = m_get(H->m, H->n);
    for (i = 0; i < E->m; ++i)
    {
        for (j = 0; j < E->n; ++j)
        {
            E->me[i][j] = alpha * H->me[i][j] + (1 - alpha) * 1 / E->m;
        }
    }


    VEC *U0 = v_get(E->m);
    for (i = 0; i < U0->dim; ++i)
    {
        U0->ve[i] = 1.0 / E->m;
    }


    VEC *Un = vm_mult(U0, m_pow(E, n));

    m_output(stdout, E);
    v_output(stdout, Un);


    v_free(U0);
    v_free(Un);
    v_free(succ);
    m_free(H);
    m_free(E);
    m_free(G);

//    SMAT *SG;
//    fp = fopen("dataset/genetic.dat", "r");
//    SG = sm_input(fp);
//    fclose(fp);





//    fp = fopen("dataset/test.dat", "w");
//    sm_output(fp, SG);
//    sm_free(SG);

    exit(0);
}