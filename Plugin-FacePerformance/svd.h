#include <mkl.h>
template<class TMATA, class TMATU, class TMATV, class R>
int mklSvdcmp(const TMATA& a, int m, int n, TMATU& u, R *w, TMATV& vt)
{
    using namespace std;

    int lda = m;
    int ldu = m;
    int ldvt = n;

    vector<double> at(m*n);
    vector<double> ut(ldu*ldu);
    vector<double> vtt(ldvt*ldvt);
    int lwork = max(3*min(m, n) + max(m, n), 5*min(m, n));
    vector<double> work(lwork);
    vector<double> ww(m);

    int i;
    // Translate the row dominant matrix to column dominant
    for (i=0; i<m; i++)
        for (int j=0; j<n; j++)
            at[j*m+i] = a[i][j];
       
    char u_c = 'A';
    char v_c = 'A';
    int info;
   
    dgesvd(&u_c, &v_c, &m, &n, &at[0], &lda, &ww[0], &ut[0], &ldu, &vtt[0], &ldvt, &work[0], &lwork, &info);

    // Translate column dominant to row dominant matrix
    for (i=0; i<ldu; i++)
        for (int j=0; j<ldu; j++)
            u[j][i] = ut[i*ldu+j];

    for (i=0; i<ldvt; i++)
        for (int j=0; j<ldvt; j++)
            vt[i][j] = vtt[i*ldvt+j];

    for(i=0; i<m; i++) w[i] = ww[i];
       
    if(info<0)    fprintf(stderr, "\nError : the %dth parameter has an illegal value!");
    if(info>0)    fprintf(stderr, "\n%d superdiagonals of intermediate matrix do not converge to 0!");
    return info;
}
