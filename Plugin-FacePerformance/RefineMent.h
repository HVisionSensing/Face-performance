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
	//����
	Eigen::VectorXd x;  //����ϵ������������������46ά
	std::vector<SparseMatrixXd> GTHGF;
	Eigen::SimplicialLDLT<SparseMatrixXd> LDLT_solver; //LDLT_solver.solve(GTHGF[i] * neutral_mesh)
	Eigen::MatrixXd P_eigenvectors; //һ��һ������������ÿ����3mά��mΪ����Ķ�����
	Eigen::VectorXd P_eigenvalues;  //����ֵ
	Eigen::MatrixXd E_eigenvectors;
	Eigen::VectorXd E_eigenvalues;

	//���
	Eigen::MatrixXd delta_B; //3m*n ��mΪ���񶥵���(m=11510); nΪ������(n=46,������������)
	Eigen::VectorXd b0;     //������ 3m*1
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