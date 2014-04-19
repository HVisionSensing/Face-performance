#pragma   warning(disable:   4786)

#ifndef NUMC_CSRMATRIX_H
#define NUMC_CSRMATRIX_H

#include <cmath>
#include <map>
#include <vector>
#include <limits>
#include <algorithm>

#ifndef _NUMC_BEGIN
#define _NUMC_BEGIN	namespace numc {
#define _NUMC_END	}
#endif

_NUMC_BEGIN
using std::vector;
using std::map;
using std::swap;
using std::lower_bound;

template<typename R> class CSRMatrix;
template<typename R> struct RowMat;
template<typename R> void CreateCSRMatrixFromRowMap(CSRMatrix<R>&, const RowMat<R>&);
template<typename R> void CSRMatrixTranspose (const CSRMatrix<R>&, CSRMatrix<R>&);
template<typename R> bool Mul2MatricesSymmResult (const CSRMatrix<R>&, const CSRMatrix<R>&, RowMat<R>&);
template<typename R> bool Mul2MatricesSymmResult (const CSRMatrix<R>&, const CSRMatrix<R>&, CSRMatrix<R>&);
template<typename R> bool Mul2Matrices (const CSRMatrix<R>&, const CSRMatrix<R>&, RowMat<R>&);
template<typename R> bool Mul2Matrices (const CSRMatrix<R>&, const CSRMatrix<R>&, CSRMatrix<R>&);

template<typename R>
class CSRMatrix
{
public:
	enum MatType{ RealStrucSymm=1, RealSymmPosDef=2, RealSymmIndef=-2, 
		CompStrucSymm=3, CompHermPosDef=4, CompHermIndef=-4, 
		CompSymm=6, RealUnSymm=11, CompUnSymm=13};

	std::vector<R> mAv;
	std::vector<int> mAi, mAj;

	CSRMatrix():mNRow(0),mNCol(0),mMtype(RealUnSymm) {};
	CSRMatrix(const RowMat<R> &rm):mMtype(RealUnSymm)	{ CreateCSRMatrixFromRowMap(*this, rm); }
	MatType mMtype;

#if ( defined(_MSC_VER)&&(_MSC_VER>1300) )
	template<typename T> class CSRMatrix;
#define FFM <>	// FFM FRIEND_FUNCTION_MAGIC
#else
#define FFM
#endif

	friend void CreateCSRMatrixFromRowMap FFM(CSRMatrix<R>&, const RowMat<R>&);
	friend void CSRMatrixTranspose FFM(const CSRMatrix<R>&, CSRMatrix<R>&);
	friend bool Mul2MatricesSymmResult FFM(const CSRMatrix<R>&, const CSRMatrix<R>&, RowMat<R>&);
	friend bool Mul2MatricesSymmResult FFM(const CSRMatrix<R>&, const CSRMatrix<R>&, CSRMatrix<R>&);
	friend bool Mul2Matrices FFM(const CSRMatrix<R>&, const CSRMatrix<R>&, RowMat<R>&);
	friend bool Mul2Matrices FFM(const CSRMatrix<R>&, const CSRMatrix<R>&, CSRMatrix<R>&);

#undef FFM

private:
	int mNRow, mNCol;

public:
	inline int nRow() const { return mNRow; }
	inline int nCol() const { return mNCol; }
	inline int empty() const { return (0==nRow() || 0==nCol()) ;}
	inline bool onebase() const { return (!mAi.empty()) && (mAi[0]==1); }		//default: zerobase, including empty matrix!
	inline MatType mtype() const { return mMtype; }

	inline bool issymm() const
	{ return (RealSymmIndef==mMtype || RealSymmPosDef==mMtype || CompSymm==mMtype); }

	inline void clear() 
	{ mNRow = 0; mNCol = 0; mMtype = RealUnSymm; mAv.clear(); mAi.clear(); mAj.clear(); }

	R getElementV(int x, int y) const{
		double *p = const_cast<CSRMatrix*>(this)->getElementP(x, y);
		return (NULL==p)?0:*p;	
	}

	R* getElementP(int x, int y) {
		if (x<0||y<0)	return NULL;
		if ( issymm() && x>y) swap(x, y);
		const int off = onebase();
		std::vector<int>::const_iterator it = lower_bound(mAj.begin()+mAi[x]-off, mAj.begin()+(mAi[x+1]-off), y+off);
		if (it==mAj.begin()+(mAi[x+1]-off) || *it!=(y+off) ) return NULL;

		return &mAv.front() + ( it-mAj.begin() );
	}

	void MultiVect(R* in, R *out) const;

	bool ChangeBase(bool oneBase) {
		if( onebase() == oneBase )	return true;
		int ii;
		if (oneBase) {
			if (mAi[0] != 0){
				fprintf(stderr, "error matrix!");
				return false;
			}

			for (ii=0; ii<mAi.back(); ii++) mAj[ii]++;

			for (ii=0; ii<mNRow+1; ii++) mAi[ii]++;
		}
		else {
			if (mAi[0] != 1){
				fprintf(stderr, "error matrix!");
				return false;
			}

			for (ii=0; ii<mAi.back()-1; ii++) mAj[ii]--;

			for (ii=0; ii<mNRow+1; ii++) mAi[ii]--;
		}
		return true;
	}

	void GetSubMat(const std::vector<int> &xi, const std::vector<int> &yi, CSRMatrix& mat){
		mat.resize( yi.size() );
		for(int j=0; j<yi.size(); j++){
			for(int i=0; i<xi.size(); i++){

			}
		}

	}
};

template<typename R>
struct RowMat:public vector<map<int, R> >
{
	typedef map<int, R> RowType;
	typedef vector<RowType> BaseType;
	int mNCol;

	RowMat(int nrow=0, int ncol=0):BaseType(nrow),mNCol(ncol)	{ }
	~RowMat()	{ }
	RowMat(const CSRMatrix<R> &mat)
	{
		resize(mat.nCol(), mat.nRow());

		int off = mat.onebase()?1:0;
		for(unsigned int i=0; i+1<mat.mAi.size(); i++){
			for(int j=mat.mAi[i]; j<mat.mAi[i+1]; j++){
				(*this)(i, mat.mAj[j-off]-off) = mat.mAv[ j-off ];
			}
		}

	}

	virtual inline R operator()(int i, int j) const{
		const RowType &crow = (*this)[i];
		RowType::const_iterator it = crow.find(j);
		return ( it==crow.end() )?0:it->second;
	}

	virtual inline R& operator()(int i, int j){
		RowType &crow = (*this)[i];
		RowType::iterator it = crow.find(j);

		return ( it==crow.end() )?crow.insert( std::pair<int,R>(j, 0) ).first->second:it->second;
	}

	int clearZero() {
		int nZero;
		for (size_t i=0; i<size(); i++) {
			RowType &r = (*this)[i];
			for (RowType::const_iterator it=r.begin(); it!=r.end(); ) {
				RowType::const_iterator cit = it++;
				if(fabs(cit->second) < std::numeric_limits<R>::epsilon()){
					r.erase(cit);
					++nZero;
				}
			}
		}
		return nZero;
	}

	inline void resize(int nrow, int ncol=0) { BaseType::resize(nrow); mNCol = ncol>0?ncol:mNCol; }
	inline int nRow()	const { return BaseType::size(); }
	inline int nCol()	const { return mNCol; }

	int operator+=(const RowMat<R>& r) { return add(r, 1); }

	void operator *= (const double v) {
		for(int i=0; i<nRow(); i++){
			map<int, R>& crow = (*this)[i];
			for(map<int,R>::iterator it=crow.begin(); it!=crow.end(); ++it){
				it->second *= v;
			}
		}
	}

	int add(const RowMat<R>& r, R k){
		const int h = nRow();
		if( r.nRow() != h )	return -1;

		for(int i=0; i<h; i++){
			const map<int, R>& rrow = r[i];
			map<int, R>& lrow = (*this)[i];
			for(RowType::const_iterator it=rrow.begin(); it!=rrow.end(); ++it){
				int col = it->first;
				R val = it->second * k;

				map<int,R>::iterator f = lrow.find(col);
				if( f != lrow.end() )	f->second += val;
				lrow.insert( make_pair<int, R>(col, val) );
			}
		}

		return 0;
	}

};

template<typename R>
struct RowMatSym : public RowMat<R>
{
	typedef RowMat<R> Base;

	RowMatSym(int nrow=0, int ncol=0):Base(nrow, ncol)	{}
	~RowMatSym()	{ }

	RowMatSym(const RowMat &rm){
		resize(rm.nCol(), rm.nRow());

		for(unsigned int i=0; i<rm.nRow(); i++){
			const RowType &crow = rm[i];
			for(RowType::const_iterator it=crow.begin(); it!=crow.end(); ++it){
				(*this)(i, it->first) = it->second;
			}
		}
	}

	RowMatSym(const CSRMatrix<R> &mat)
	{
		resize(mat.nCol(), mat.nRow());

		int off = mat.onebase()?1:0;
		for(unsigned int i=0; i+1<mat.mAi.size(); i++){
			for(int j=mat.mAi[i]; j<mat.mAi[i+1]; j++){
				const int k = mat.mAj[j-off]-off;
				(*this)(i, k) = mat.mAv[ j-off ];
			}
		}
	}

	inline R operator()(int i, int j) const{
		if(i>j) std::swap(i, j);
		const RowType &crow = (*this)[i];
		RowType::const_iterator it = crow.find(j);
		return ( it==crow.end() )?0:it->second;
	}

	inline R& operator()(int i, int j){
		if(i>j) std::swap(i, j);
		RowType &crow = (*this)[i];
		RowType::iterator it = crow.find(j);

		return ( it==crow.end() )?crow.insert( std::pair<int,R>(j, 0) ).first->second
			:it->second;
	}

};


//////////////////////////////////////////////////////////////////////////
/// Get SubMatrix from rowmap
//////////////////////////////////////////////////////////////////////////
template <typename R>
void GetSubRowMap(const RowMat<R>& src,
				  const std::vector<int> &xi, const std::vector<int> &yi,
				  RowMat<R>& dst)
{
	assert( xi.size()>0 && yi.size()>0 );
	assert( yi.back() < (int)src.size() );
	if( xi.empty() || yi.empty() || yi.back() >= (int)src.size() ){
		fprintf(stderr, "\nerror matrix!");
		return;
	}


	typedef std::map<int, R> row;
	dst.resize( yi.size() );

	for(unsigned int j=0; j<yi.size(); j++){
		const row &srow = src[ yi[j] ];
		row &drow = dst[j];

		const int* pxi = &xi.front();

		for(row::const_iterator it=srow.begin(); it!=srow.end(); ++it ){
			pxi = std::lower_bound(pxi, &xi.back(), it->first);

			if(*pxi != it->first){			// done
				continue;
			}

			int col = pxi - &xi.front();
			drow[col] = it->second;
		}
	}
}

//////////////////////////////////////////////////////////////////////////
/// Create CSR matrix from vector of rowmap
//////////////////////////////////////////////////////////////////////////
template <typename R>
void CreateCSRMatrixFromRowMap(CSRMatrix<R> &matC, const RowMat<R>& rowC)
{
	if(rowC.mNCol <= 0){
		fprintf(stderr, "\nnumber of column not specified in the RowMat, exit!");
		return;
	}
	matC.clear();
	matC.mNRow = rowC.nRow();
	matC.mNCol = rowC.nCol();

	matC.mAi.reserve(matC.nRow()+1);
	int i,nnz=0;
	for (i=0; i<matC.nRow(); i++) {
		nnz += rowC[i].size();
	}
	matC.mAj.reserve(nnz);
	matC.mAv.reserve(nnz);

	// copy rows into matC
	matC.mAi.push_back(0);
	for (i=0; i<matC.nRow(); i++) {
		matC.mAi.push_back(matC.mAi.back());
		for (std::map<int, R>::const_iterator it=rowC[i].begin(); it!=rowC[i].end(); it++) {
			matC.mAi.back()++;
			matC.mAj.push_back(it->first);
			matC.mAv.push_back(it->second);
		}
	}
}


//////////////////////////////////////////////////////////////////////////
/// Computes the transpose of a matrix.
template <typename R>
void CSRMatrixTranspose(const CSRMatrix<R> &matA, CSRMatrix<R> &matAT)
{
	if (matA.issymm()) {		// symmetric - just copy the matrix
		matAT = matA;
		return;
	}

	matAT.mNRow = matA.nCol();
	matAT.mNCol = matA.nRow();
	matAT.mMtype = matA.mMtype;

	// non-symmetric matrix -> need to build data structure.
	// we'll go over the columns and build the rows
	int off = matA.onebase()?1:0;
	RowMat<R> rowC(matA.nCol(), matA.nRow());
	for (int i=0; i<matA.nRow(); i++) {
		for (int j=matA.mAi[i]; j<matA.mAi[i+1]; j++) {
			rowC[matA.mAj[j-off]-off][i] = matA.mAv[j-off];
		}
	}

	CreateCSRMatrixFromRowMap(matAT, rowC);
}


//////////////////////////////////////////////////////////////////////////
// multiplication of sparse matrix
// Assuming nothing about the result (the result is not stored symmetric)
//////////////////////////////////////////////////////////////////////////
template <typename R>
bool Mul2Matrices(const CSRMatrix<R> &matA, const CSRMatrix<R> &matB,
				  RowMat<R> &rowsC)
{
	if(matA.onebase() || matB.onebase() ){
		fprintf(stderr, "\nmatrix saved in 1-based format, pleased converted it to 0-based and try again!");
		return false;
	}
	// Compatibility of dimensions
	if (matA.nCol() != matB.nRow())
		return false;

	// (m x n)*(n x k) = (m x k)
	const int m=matA.nRow();
	const int n=matA.nCol();
	const int k=matB.nCol();

	rowsC = RowMat<R>(m, k);	// clean all coefficients

	R aiv, valB;
	int colInd, colB;

	for (int i=0; i<m; ++i) {					// creating row i of C
		std::map<int, R> &mapRow2Val = rowsC[i];
		for (int iAi = matA.mAi[i]; iAi < matA.mAi[i+1]; ++iAi) {			// travel on ai
			colInd = matA.mAj[iAi];
			aiv = matA.mAv[iAi];
			// make aiv*b_{rowInd} and insert into mapRow2Val
			for (int iB=matB.mAi[colInd]; iB<matB.mAi[colInd+1]; ++iB) {
				colB=matB.mAj[iB];
				valB=matB.mAv[iB];
				// insert valA*aiv into map
				std::map<int, R>::iterator it = mapRow2Val.find(colB);
				if (it == mapRow2Val.end()) {		// first time
					mapRow2Val[colB] = valB*aiv;
				}
				else {
					it->second = it->second + valB*aiv;
				}
			}
		}
	}// now column i is created

	return true;
}

template <typename R>
bool Mul2Matrices(const CSRMatrix<R> &matA, const CSRMatrix<R> &matB, CSRMatrix<R> &matC)
{
	RowMat<R> rowsC;
	if( !Mul2Matrices(matA, matB, rowsC) ) return false;

	const int k=matB.mNCol;
	rowsC.mNCol = k;
	CreateCSRMatrixFromRowMap(matC, rowsC);						// modified by jianwei hu @ 16/09/07
	return true;
}

//////////////////////////////////////////////////////////////////////////
// multiplication of sparse matrix
// The result is symmetric
//////////////////////////////////////////////////////////////////////////
template <typename R>
bool Mul2MatricesSymmResult(const CSRMatrix<R> &matA, const CSRMatrix<R> &matB,
							RowMat<R> &rowsC)
{
	if(matA.onebase() || matB.onebase()){
		fprintf(stderr, "\nmatrix saved in 1-based format, pleased converted it to 0-based and try again!");
		return false;
	}
	// Compatibility of dimensions
	if(matA.nCol() != matB.nRow() || matA.nRow() != matB.nCol())	return false;

	// (m x n)*(n x m) = (m x m)
	const int m=matA.nRow();
	const int n=matA.nCol();

	rowsC = RowMat<R>(m, m);	// clean all coefficients

	R aiv, valB;
	int colInd, colB;

	for(int i=0; i<m; ++i) {					// creating row i of C
		std::map<int, R> &mapRow2Val = rowsC[i];
		for (int iAi = matA.mAi[i]; iAi < matA.mAi[i+1]; ++iAi) {			// travel on ai
			colInd = matA.mAj[iAi];
			aiv = matA.mAv[iAi];
			// make aiv*b_{colInd} and insert into mapRow2Val
			for (int iB=matB.mAi[colInd]; iB<matB.mAi[colInd+1]; ++iB) {
				colB=matB.mAj[iB];
				if (colB >= i) {
					valB=matB.mAv[iB];
					// insert valA*aiv into map
					std::map<int, R>::iterator it = mapRow2Val.find(colB);
					if (it == mapRow2Val.end()) {		// first time
						mapRow2Val[colB] = valB*aiv;
					}
					else {
						it->second = it->second + valB*aiv;
					}
				}
			}
		}
	}// now column i is created

	return true;
}

template <typename R>
bool Mul2MatricesSymmResult(const CSRMatrix<R> &matA, const CSRMatrix<R> &matB, CSRMatrix<R> &matC)
{
	RowMat<R> rowsC;
	if( !Mul2MatricesSymmResult(matA, matB, rowsC) )	return false;

	rowsC.mNCol = matA.nCol();
	CreateCSRMatrixFromRowMap(matC, rowsC);
	matC.mMtype = CSRMatrix<R>::RealSymmIndef;
	return true;
}


template<typename R>
void DebugShowRowMatrix(const std::vector<std::map<int,R> > &rowA)
{
	printf("\n");
	for(unsigned int i=0; i<rowA.size(); i++){
		printf("\n%3d#", i);
		std::map<int,R>::const_iterator it=rowA[i].begin();
		int j=0;
		for(; it!=rowA[i].end(); ++it){
			for(int k=j; k<it->first; k++){
				printf("\t%3.2f", 0);
			}

			j=it->first+1;
			printf("\t%3.2f", it->second);
		}
	}
}

typedef CSRMatrix<float> CSRMatrixf;
typedef CSRMatrix<double> CSRMatrixd;


_NUMC_END

#endif