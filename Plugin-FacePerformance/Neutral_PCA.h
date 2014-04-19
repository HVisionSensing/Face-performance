#ifndef NUETRAL_PCA_H_INCLUDED
#define NUETRAL_PCA_H_INCLUDED
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <QString>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <string>
#include <Eigen/Dense>
class Neutral_PCA
{
public:
	Neutral_PCA();
	~Neutral_PCA();
	static void Mesh2Mat(cv::Mat& mat, TriMesh* mesh);
	
	cv::PCA neutral_pca;
	int numNeutralMeshes;
	int numFeature_pca;

	Eigen::VectorXd mean_neutral;
	Eigen::VectorXd P_eigenvalues;
	Eigen::MatrixXd P_eigenvectors;//ÐÐ
};
#endif