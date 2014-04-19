#include "Blendshape.h"
void Blendshape::refineOtherBlendshape(const Blendshape& neutralBlendshape, Eigen::SimplicialLDLT<SparseMatrixXd>& LDLT_solver, const SparseMatrixXd& GTHGF, const GraphLaplacian_PCA& gl_pca)
{
	Eigen::VectorXd vec = GTHGF * neutralBlendshape.blendshape;
	blendshape = LDLT_solver.solve(vec) + gl_pca.E_eigenvectors*z;
}
void Blendshape::refineNeutralBlendshape(const Neutral_PCA& neutral_pca, const GraphLaplacian_PCA& gl_pca)
{
	blendshape = neutral_pca.mean_neutral + neutral_pca.P_eigenvectors * y + gl_pca.E_eigenvectors*z;
}