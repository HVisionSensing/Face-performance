#ifndef INITREGISTERMESH_HH_INCLUDED
#define INITREGISTERMESH_HH_INCLUDED
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
/*
#include "mkl_addon.h"

#include "lbfgs/HLBFGS.h"
#include "lbfgs/Lite_Sparse_Matrix.h"
#include "ICP.h"
#include "svd.h"
#include "multiplication.h"*/
#include <Eigen/Dense>
#include <fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "InitLM.h"
typedef trimesh::XForm<double> xform;
//using numc::RowMat;
//using numc::CSRMatrix;
//using numc::LeastSquareSparseSolver;


class InitRegisterMesh
{
public:
	InitRegisterMesh(InitLMParameters* paras);//trimesh::TriMesh* pointCloudMesh, std::string& meshIndexFilename, std::string& depthIndexFilename, cv::Mat& texture, int *colorCoordinate, int * marker_Coordinates);
	~InitRegisterMesh();
	trimesh::KDtree * treeCloud;

	pcl::PointCloud<pcl::PointXYZ>::Ptr pclCloud;
	pcl::KdTreeFLANN<pcl::PointXYZ> pclCloudKdtree;

	InitLMParameters* paras;

	std::vector<trimesh::TriMesh*> meshes;
	double * x;
	
	std::ofstream of;////////////Êä³ö
public:
	void readRegisterIndexFromFile(std::vector<int>& index, std::string filename);
	//void function_f(double * x, double * f_sum);
	//void function_df(double *x, double * g);
	void init_srtyz(std::string& meshIndexFilename, std::string& depthIndexFilename);
	void update_blendshape();
	void update_all();
	void find_nearest();
	void fit_mesh();
	void trimesh2pcl(const trimesh::TriMesh *mesh, pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud);
	void blendshape2trimesh(trimesh::TriMesh * tmesh);
	void getinitParas(Eigen::VectorXd& vec, double& s, Eigen::VectorXd& r, Eigen::VectorXd& t, Eigen::MatrixXd& _3d_to_2d_matrix);
	void loadEigenMatrix(Eigen::MatrixXd& mat, string filename);
	void loadEigenVector(Eigen::VectorXd& vec, string filename);
	void loadEigenMatrix_gl(Eigen::MatrixXd& mat, string filename);
	void loadEigenVector_gl(Eigen::VectorXd& vec, string filename);
	void get_3D_to_2D_Matrix(Eigen::MatrixXd & RT);
};
#endif