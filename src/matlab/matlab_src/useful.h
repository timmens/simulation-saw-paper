void gaussj(double *a, int n, double *b, int m);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, const double *R, const int p, const double lam, const double *g,
	double (*func)(double,const double*, const int, const double, const double*));

int brent(double ax, double bx, double cx, const double *R, const int p1, const double lam, const double *g, double (*f)(double,const double*, const int, const double, const double*), double tol,
	double *xmin);

void vcopy(const double *y,int n, double *x);

void sumr(const double *x, int n, int p, int istart, int length, double *s);
    
void sumc(const double *x, int n, int p, int jstart, int length, double *s);

void summ(double *x, int n, int p, int dim, double *s);

void cumsum(double *x, int n, int p, double *ans);

void invshift(double *x, int n, int p, double *y);

void fillzero(double *x,int n);

void fillv(double *x,int n,double v);

void dotmult(double *A, int n, int p, double *B);

double dotprod(double *x, double n, double *y);

void getrow(double *R,int n,int p,int q,int row,double *A);

void insert_row(double *R,int n,int p,int q,int k,const double *A);

void getrow2(double *Rn,int n,int p,int row,double *A);

double norm2(double *x, int n);

void mtimesv(const double *a, int n, int p, const double *b, double *c);

void mtimes(const double *a, int n, int p, const double *b, int q, double *c);

void mtrans(const double *a, int n, int p, double *b);

void madd(double *a, int n, int p, const double *b);

void vadd(double *a, int n, const double *b);

void vsub(double *a, int n, const double *b);

void sandwitch(double *a, int n, int p, double *b, double *c, int q, double *d);

void getblock(double *R,int M,int N,int m,int n, int nrow, int ncol,double *A);

void putblock(double *R,int M,int N,int m,int n, int nrow, int ncol,double *A);


