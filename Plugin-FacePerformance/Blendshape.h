#ifndef BLENDSHAPE_HH_INCLUDED
#define BLENDSHAPE_HH_INCLUDED
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Neutral_PCA.h"
#include "GraphLaplacian_PCA.h"
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
class Blendshape
{
public:
	Blendshape();
	~Blendshape();
	void refineOtherBlendshape(const Blendshape& neutralBlendshape, Eigen::SimplicialLDLT<SparseMatrixXd>& LDLT_solver, const SparseMatrixXd& GTHGF, const GraphLaplacian_PCA& gl_pca);
	void refineNeutralBlendshape(const Neutral_PCA& neutral_pca, const GraphLaplacian_PCA& gl_pca);

	bool isNeutral;
	Eigen::VectorXd y;
	Eigen::VectorXd z;
	
	Eigen::VectorXd blendshape;
//private:
	
};
#endif