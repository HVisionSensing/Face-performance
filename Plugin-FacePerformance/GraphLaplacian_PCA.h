#ifndef GRAPHLAPLACIAN_PCA_H_INCLUDED
#define GRAPHLAPLACIAN_PCA_H_INCLUDED
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <QString>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <string>
#include <Eigen/Dense>
//#include<fstream>
class GraphLaplacian_PCA
{
public:
	GraphLaplacian_PCA();
	~GraphLaplacian_PCA();
	void Mesh2GLMat(cv::Mat& mat, TriMesh& mesh);
public:
	cv::PCA gl_pca;
	Eigen::MatrixXd E_eigenvectors; //ÐÐ
	Eigen::VectorXd E_eigenvalues;
	
};
#endif