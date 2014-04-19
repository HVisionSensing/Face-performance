//multiplication : c = a*b
//a	:	matrix with m*k
//b	:	matrix with k*n
//c	:	results

#include <mkl.h>
#include <vector>
using namespace std;
template<class TMATA, class TMATB, class TMATC>
int mklMulti(const TMATA& a, int m, int k, CBLAS_TRANSPOSE aTran, const TMATB& b, int n, CBLAS_TRANSPOSE bTran, TMATC& c, CBLAS_ORDER myMajor)
{

	vector<double>	at(m*k);
	vector<double>	bt(k*n);
	vector<double>	ct(m*n);

	for (int i=0; i<m; i++)
	{
		for (int j=0; j<k; j++)
		{
			at[i*k+j]	=	a[i][j]; 
		}
	}

	for (int i=0; i<k; i++)
	{
		for (int j=0; j<n; j++)
		{
			bt[i*n+j]	=	b[i][j]; 
		}
	}

	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			ct[i*n+j]	=	0.0; 
		}
	}

	cblas_dgemm(myMajor, aTran, bTran, m, n, k, 1, &at[0], k, &bt[0], n, 0, &ct[0], n);

	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			c[i][j]	=	ct[i*n+j]; 
		}
	}

	return 0;
}