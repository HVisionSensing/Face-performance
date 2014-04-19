#ifndef RIGIDMESH_H_INCLUDED
#define RIGIDMESH_H_INCLUDED
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
//#include "lbfgs/HLBFGS.h"
//#include "lbfgs/Lite_Sparse_Matrix.h"
//#include "multiplication.h"
#include <Eigen/Dense>
#include <fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "RigidLM.h"
class RigidMesh
{
public:
	RigidMesh(RigidLMParameters* rigidParas);
	~RigidMesh();
	trimesh::KDtree * treeCloud;

	pcl::PointCloud<pcl::PointXYZ>::Ptr pclCloud;
	pcl::KdTreeFLANN<pcl::PointXYZ> pclCloudKdtree;

	RigidLMParameters* paras;

	std::vector<trimesh::TriMesh*> meshes;
	double * x;
	
	std::ofstream of;////////////Êä³ö
public:
	void readRegisterIndexFromFile(std::vector<int>& index, std::string filename);
	void read_configure();
	void update_all();
	void find_nearest();
	void fit_mesh();
	void trimesh2pcl(const trimesh::TriMesh *mesh, pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud);
	void blendshape2trimesh(trimesh::TriMesh * tmesh);
};
#endif