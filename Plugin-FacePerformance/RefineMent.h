#ifndef REFINEMENT_H_INCLUDED
#define REFINEMENT_H_INCLUDED
#include "stdafx.h"
#include <vector>
#include <string>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <TriMesh.h>
#include "TriMesh_algo.h"
#include <KDtree.h>
#include <XForm.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <algorithm>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
typedef Eigen::SparseMatrix<double> SparseMatrixXd;

class RefineMentParameters
{
public:
	//输入
	Eigen::VectorXd x;  //表情系数，不包括中性脸，46维
	std::vector<SparseMatrixXd> GTHGF;
	Eigen::SimplicialLDLT<SparseMatrixXd> LDLT_solver; //LDLT_solver.solve(GTHGF[i] * neutral_mesh)
	Eigen::MatrixXd P_eigenvectors; //一列一个特征向量，每列有3m维，m为网格的顶点数
	Eigen::VectorXd P_eigenvalues;  //特征值
	Eigen::MatrixXd E_eigenvectors;
	Eigen::VectorXd E_eigenvalues;

	//输出
	Eigen::MatrixXd delta_B; //3m*n ，m为网格顶点数(m=11510); n为表情数(n=46,不包括中性脸)
	Eigen::VectorXd b0;     //中性脸 3m*1
};

class RefineMent
{
	RefineMent(RefineMentParameters* paras)
	{
		this->paras = paras;
	}
	~RefineMent();
private:
	RefimentParameters * paras;
};
#endif