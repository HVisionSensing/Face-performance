#ifndef REGISTERMESH_HH_INCLUDED
#define REGISTERMESH_HH_INCLUDED
/*
#include <vector>
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
//#include "mkl_addon.h"
#include <math.h>
//#include "lbfgs/HLBFGS.h"
//#include "lbfgs/Lite_Sparse_Matrix.h"
//#include "ICP.h"
//#include "svd.h"
//#include "multiplication.h"

typedef trimesh::XForm<double> xform;
//using numc::RowMat;
//using numc::CSRMatrix;
//using numc::LeastSquareSparseSolver;
class RegisterMesh
{
public:
	RegisterMesh(trimesh::TriMesh* pointCloudMesh);
	~RegisterMesh();
public:
	void DrawDeformCloud();
	void Transformation();
	void MySVD(std::vector<trimesh::point>& KPts, std::vector<trimesh::point>& TPts, float R[3][3], float T[3]);
	void Generate();
	void Non_Rigid_Reg();
	inline void WrapMesh(trimesh::TriMesh * , const xform& );
	void ICP(const xform& );

	void InitProject();
	void CalProjectMat();
	void CalColorIndex();

	std::vector<trimesh::point>		m_deformCloud;
	std::vector<int>         m_kinectcloudindex;
	std::vector<int>         m_templatecloudindex;

	std::vector<trimesh::TriMesh*>	meshes;
	trimesh::KDtree * m_treeTemplate;
	trimesh::KDtree * m_treeCloud;
	float y_min;
	std::vector<trimesh::point> m_pointArray;
	std::vector<trimesh::ivec2> m_imageArray;
	std::vector<trimesh::ivec2> m_imageArrayy;
	std::vector<cv::Point2d> pictureextrafeature;
	float m_fProjectMat[3][4];

	trimesh::TriMesh * pointCloudMesh;
	trimesh::TriMesh * templateMesh;
	trimesh::KDtree * treeTemplate;
	trimesh::KDtree * treeCloud;

	std::vector<int> m_templateDeformFeature;
	std::vector<int> m_cloudDeformFeature;
	std::vector<int> m_templateAdditionFeature;
};*/
#endif