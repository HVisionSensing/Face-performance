#ifndef NONRIGIDMESH_H_INCLUDED
#define NONRIGIDMESH_H_INCLUDED
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
#include <fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "NonRigidLM.h"
class NonRigidMesh
{
public:
	NonRigidMesh(NonRigidLMParameters* nonRigidparas);
	~NonRigidMesh();
	trimesh::KDtree * treeCloud;

	pcl::PointCloud<pcl::PointXYZ>::Ptr pclCloud;
	pcl::KdTreeFLANN<pcl::PointXYZ> pclCloudKdtree;

	NonRigidLMParameters* paras;

	std::vector<trimesh::TriMesh*> meshes;
	Eigen::VectorXd x;
	Eigen::VectorXd pre_x;
	
	
	std::ofstream of;////////////Êä³ö
public:
	void readRegisterIndexFromFile(std::vector<int>& index, std::string filename);
	void read_configure();
	void init();
	void update_expression();
	void update_all();
	void find_nearest();
	void update_x();
	void update_w();
	void update_v();
	void update_dual_variable();
	bool check_convergence_and_update_penalty();
	void fit_mesh();
	void getExpression(Eigen::VectorXd& expression);
	void trimesh2pcl(const trimesh::TriMesh *mesh, pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud);
	void blendshape2trimesh(trimesh::TriMesh * tmesh);
	void get_3D_to_2D_Matrix(Eigen::MatrixXd & RT);
	double f_val();
};
#endif