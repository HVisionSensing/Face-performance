#ifndef TRACKINGMESH_HH_INCLUDED
#define TRACKINGMESH_HH_INCLUDED
/*
#include "stdafx.h"
#include <vector>
#include <string>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
//#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <TriMesh.h>
#include "TriMesh_algo.h"
#include <KDtree.h>
#include <XForm.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <algorithm>
#include "mkl_addon.h"
#include <math.h>
#include "lbfgs/HLBFGS.h"
#include "lbfgs/Lite_Sparse_Matrix.h"
#include "ICP.h"
#include "svd.h"
#include "multiplication.h"
//#include "Neutral_PCA.h"
//#include "GraphLaplacian_PCA.h"
#include <Eigen/Dense>
#include <fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"
typedef trimesh::XForm<double> xform;
using numc::RowMat;
using numc::CSRMatrix;
using numc::LeastSquareSparseSolver;

class TrackingMesh//: public SizedCostFunction<1 ,1 >
{
public:

	TrackingMesh(trimesh::TriMesh* pointCloudMesh);
	~TrackingMesh();
	//Neutral_PCA neutral_pca;
	//GraphLaplacian_PCA gl_pca;
	trimesh::TriMesh * pointCloudMesh;
	//trimesh::TriMesh * templateMesh;
	//trimesh::KDtree * treeTemplate;
	trimesh::KDtree * treeCloud;

	pcl::PointCloud<pcl::PointXYZ>::Ptr pclCloud;
	pcl::KdTreeFLANN<pcl::PointXYZ> pclCloudKdtree;

	std::vector<trimesh::TriMesh*> meshes;
	int numX;
	Eigen::VectorXd mean_neutral;
	Eigen::VectorXd P_eigenvalues;
	Eigen::MatrixXd P_eigenvectors;
	
	Eigen::VectorXd E_eigenvalues;
	Eigen::MatrixXd E_eigenvectors;

	std::vector<int> registerIndex;
	std::vector<int> nearestIndex;
	double * x;
	int startR;
	int startt;
	int starty;
	int startz;
	int numy;
	int numz;
	int numR;
	int numt;
	double beta1;
	double beta2;
	double beta3;
	Eigen::VectorXd neutral_blendshape;
	std::ofstream of;////////////Êä³ö

	double alpha;
	double mu;
	std::vector<double> x;
	Eigen::MatrixXd delta_B;
public:
	void readRegisterIndexFromFile(std::vector<int>& index, std::string filename);
	void function_f(double * x, double * f_sum);
	void function_df(double *x, double * g);
	void initRtyz();
	void update_blendshape();
	void update_all();
	void find_nearest();
	void trimesh2pcl(const trimesh::TriMesh *mesh, pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud);
	void blendshape2trimesh(trimesh::TriMesh * tmesh);
	void loadEigenMatrix(Eigen::MatrixXd& mat, string filename);
	void loadEigenVector(Eigen::VectorXd& vec, string filename);
	void loadEigenMatrix_gl(Eigen::MatrixXd& mat, string filename);
	void loadEigenVector_gl(Eigen::VectorXd& vec, string filename);

	update_x();
	update_d();
	update_w();
	
};*/
#endif