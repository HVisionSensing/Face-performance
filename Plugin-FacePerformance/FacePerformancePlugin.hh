#ifndef FACEPERFORMANCEPLUGIN_HH_INCLUDED
#define FACEPERFORMANCEPLUGIN_HH_INCLUDED
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/common/Types.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <opencv/cv.h>
#include <opencv/highgui.h>
//#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <vector>
#include <string>
#include <NuiApi.h>
#include "FeatureExtraction.h"
#include <TriMesh.h>
#include "InitRegisterMesh.h"
#include "RigidMesh.h"
#include "NonRigidMesh.h"
#include "NonRigidLM.h"
#include "RigidLM.h"
//#include "InitRegister.h"
#include "Neutral_PCA.h"
#include "GraphLaplacian_PCA.h"
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
typedef Eigen::Triplet<double> T;
class FacePerformancePlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, LoadSaveInterface
{
	Q_OBJECT
	Q_INTERFACES(BaseInterface)
	Q_INTERFACES(LoggingInterface)
	Q_INTERFACES(ToolboxInterface)
	Q_INTERFACES(LoadSaveInterface)
signals:
	//BaseInterface
	void updateView();
	void updateObject(int _identifier, const UpdateType& _type);
	//LoggingInterface
	void log(Logtype _type, QString _message);
	void log(QString _message);
	//ToolboxInterface
	void addToolbox(QString _name, QWidget* _widget);
	//LoadSaveInterface
	void load(QString _filename, DataType _type, int& _id);
	void save(int _id, QString _filename);
public :
	FacePerformancePlugin();
	~FacePerformancePlugin(){};
	QString name() { return QString("FacePerformancePlugin"); };
	QString description() { return QString("FacePerformance work!"); };
	void compute_T_Star();
	void load_T_Star();
	void compute_S_face(TriMesh * trimesh, OpenMesh::FPropHandleT<Eigen::MatrixXd>& S, bool isInverse = false);
	void compute_H_Star(TriMesh * neutral_mesh, OpenMesh::FPropHandleT<Eigen::MatrixXd>& S_neutral, TriMesh * blendshape_mesh, OpenMesh::FPropHandleT<Eigen::MatrixXd>& S_blendshape, SparseMatrixXd & H_Star);
	void computeBlendShapeByNeutral(TriMesh* neutral_mesh, int id_mesh, Eigen::SimplicialLDLT<SparseMatrixXd>& LDLT_solver, const SparseMatrixXd& GTHGF);//SparseMatrixXd& T_Star);
	void computeDeltaBByNeutral(Eigen::VectorXd& neutral_mesh, Eigen::MatrixXd& blendshape);
	//kinect
	void getColorImage(HANDLE &colorEvent, HANDLE &colorStreamHandle, cv::Mat &colorImage); 
	void getDepthImage(HANDLE &depthEvent, HANDLE &depthStreamHandle, cv::Mat &depthImage); 
	
	
private:
	QSpinBox* iterationsSpinbox_;
	QPushButton* loadButton;
	QPushButton* openKinectButton;
	QPushButton* initRegisterButton;
	QPushButton* readFrameButton;
	QPushButton* trackingButton;

	QWidget* toolBox;
	QGridLayout* layout;
	const static int num_blendshape = 47;//47;
	int id_blendshape[num_blendshape];

	SparseMatrixXd GTGF;         
	vector<SparseMatrixXd> GTHGF;
	Eigen::SimplicialLDLT<SparseMatrixXd> LDLT_solver;

	//Kinect
	cv::Mat colorImage;
	cv::Mat depthImage;
	FeatureExtraction* featureExtraction;

	INuiSensor *pNuiSensor;
	HANDLE h1, h2, h3,h4;
	bool isDepthCollection;
	LONG *m_colorCoordinates;
	BYTE *m_colorRGBX;
	INuiSensor *m_pNuiSensor;
	float fx, fy;
	std::vector<int> m_DepthIndex;
	QTimer *timer;
	
	

	//trimesh
	trimesh::TriMesh * depthMesh;
	trimesh::TriMesh * templateMesh;
	std::ofstream of;

	bool isStopReadFrame;
	
private slots:
	//BaseInterface
	void initializePlugin();
	
public slots:
	void loadBlendshapes();
	void openKinect();
	void ReadFrame();
	void initRegisterMesh();
	void stopReadFrame();
	void trackingMesh();
	QString version(){return QString("1.0");}
};
#endif