#ifndef INITLM_H_INCLUDED
#define INITLM_H_INCLUDED
#include <vector>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include <TriMesh.h>
#include <fstream>
#include <opencv/cv.h>
#include <opencv/highgui.h>
using ceres::CostFunction;
using ceres::SizedCostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class InitLMParameters
{
public:
	double beta1;
	double beta2;
	double beta3;
	int outIters;
	int maxInIters;
	double markerWeight;
	int num_s;
	int num_r;
	int num_t;
	int num_y;
	int num_z;
	int num_all;
	int start_s;
	int start_r;
	int start_t;
	int start_y;
	int start_z;
	Eigen::VectorXd neutral_blendshape;
	std::vector<int> init_registerIndex;
	std::vector<int> init_nearestIndex;
	std::vector<double> init_nearestWeight;
	double registerRate;
	cv::Mat texture;//纹理
	double * colorCoordinate;
	double * marker_Coordinates;
	Eigen::MatrixXd _3d_to_2d_matrix;
	std::string meshIndexFilename;
	std::string depthIndexFilename;
	double maxDepth;

	std::vector<int> registerIndex;
	trimesh::TriMesh* pointCloudMesh;
	std::vector<int> nearestIndex;
	std::vector<double> nearestWeight;
	std::vector<int> marker_registerIndex;
	std::vector<int> marker_nearestIndex;
	std::vector<double> marker_nearestWeight;
	Eigen::VectorXd mean_neutral;
	Eigen::MatrixXd P_eigenvectors;
	Eigen::MatrixXd E_eigenvectors;
	Eigen::VectorXd P_eigenvalues;
	Eigen::VectorXd E_eigenvalues;
};
class InitLMCostFunction: public CostFunction
{
 public:
	 std::ofstream of;
	 InitLMParameters* paras;
	 double *s;
	 double *r;
	 double *t;
	 double *y;
	 double *z;
	 static Eigen::Vector2d f_vt(double x, double y, double z, cv::Mat& colorImage, double * colorCoordinate, Eigen::MatrixXd & _3d_to_2d_matrix)
	 {
		 Eigen::Vector2d vt;
		 double texture_x = x * _3d_to_2d_matrix(0,0) + y * _3d_to_2d_matrix(1,0) + z * _3d_to_2d_matrix(2,0) + _3d_to_2d_matrix(3,0);
		 double texture_y = x * _3d_to_2d_matrix(0,1) + y * _3d_to_2d_matrix(1,1) + z * _3d_to_2d_matrix(2,1) + _3d_to_2d_matrix(3,1);
		 texture_x *= 90.0;
		 texture_y *= 90.0;
		 if(texture_x<0)texture_x = 0;
		 else if(texture_x>=colorImage.cols)texture_x = colorImage.cols-1;
		 if(texture_y<0)texture_y = 0;
		 else if(texture_y>=colorImage.rows)texture_y = colorImage.rows-1;
		 vt(0) = texture_x * 1.0 / colorImage.cols;
		 vt(1) = 1.0-texture_y * 1.0 / colorImage.rows;
		 return vt;
		 /*
		 Eigen::Vector2d vt;
		 double fx = 1.2017214;
		 double fy = 0.9030345;
		 double maxDepth = 900;
		 Eigen::VectorXd color(3);//三通道
		 int m = (0.5-y/fy/(maxDepth-z))*colorImage.rows+0.5;
		 int n = (x/fx/(maxDepth-z)+0.5)*colorImage.cols+0.5;
		 if(m<0)m=0;
		 else if(m>=colorImage.rows)m=colorImage.rows-1;
		 if(n<0)n=0;
		 else if(n>colorImage.cols)n=colorImage.cols-1;

		 int k=m*colorImage.cols+n;
		 //vt(0) = n*1.0/colorImage.cols;
		 //vt(1) = m*1.0/colorImage.rows;
		 vt(0) = colorCoordinate[k*2]*1.0/colorImage.cols;
		 vt(1) = 1.0-colorCoordinate[k*2+1]*1.0/colorImage.rows;
		 return vt;
		 /*
		 for(int i=0; i<3; i++)
		 {
			// output[ll_l].at<Vec<unsigned char,3> >(ii,jj)[kk];
			 color(i) = colorImage.at<cv::Vec<unsigned char,3> >(m,n)[i];
		 }
		 return color;*/
	 }
	 static double getr1(double * r, Eigen::VectorXd& v, Eigen::VectorXd& n)
	 {
		 double alpha = r[0];
		 double beta = r[1];
		 double gamma = r[2];
		 double val = ( cos(alpha)*sin(beta)*cos(gamma)+sin(gamma)*sin(alpha) )*v(0)*n(1)
					+ (	sin(beta)*sin(gamma)*cos(alpha)-sin(alpha)*cos(gamma) )*v(1)*n(1)
					+  cos(alpha)*cos(beta) * v(2)*n(1)
					+ ( -sin(alpha)*sin(beta)*cos(gamma)+cos(alpha)*sin(gamma))*v(0)*n(2)
					+ ( -sin(alpha)*sin(beta)*sin(gamma)-cos(alpha)*cos(gamma))*v(1)*n(2)
					+ ( -sin(alpha)*cos(beta)) *v(2)*n(2);
		 return val;
	 }
	 static double getr2(double * r, Eigen::VectorXd& v, Eigen::VectorXd& n)
	 {
		 double alpha = r[0];
		 double beta = r[1];
		 double gamma = r[2];
		 double val = -sin(beta)*cos(gamma)*v(0)*n(0)
					  -sin(beta)*sin(gamma)*v(1)*n(0)
					  -cos(beta)*v(2)*n(1)
					  +sin(alpha)*cos(beta)*cos(gamma)*v(0)*n(1)
					  +sin(alpha)*cos(beta)*sin(gamma)*v(1)*n(1)
					  -sin(alpha)*sin(beta)*v(2)*n(1)
					  +cos(beta)*cos(alpha)*cos(gamma)*v(0)*n(2)
					  +cos(alpha)*cos(beta)*sin(gamma)*v(1)*n(2)
					  -cos(alpha)*sin(beta)*v(2)*n(2);
		 return val;
	 }
	 static double getr3(double * r, Eigen::VectorXd& v, Eigen::VectorXd& n)
	 {
		 double alpha = r[0];
		 double beta = r[1];
		 double gamma = r[2];
		 double val = -cos(beta)*sin(gamma)*v(0)*n(0)
					  +cos(beta)*cos(gamma)*v(1)*n(0)
					  +(-sin(alpha)*sin(beta)*sin(gamma)-cos(gamma)*cos(alpha))*v(0)*n(1)
					  +(sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma))*v(1)*n(1)
					  +(-cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma))*v(0)*n(2)
					  +(cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma))*v(1)*n(2);
		 return val;
	 }
	 static Eigen::MatrixXd getR(double* r)
	{
		double alpha = r[0];
		double beta = r[1];
		double gamma = r[2];
		/*Eigen::MatrixXd R(3,3);
		R(0,0) = cos(alpha);		R(0,1) = -sin(gamma);		R(0,2) = sin(beta);
		R(1,0) = sin(gamma);		R(1,1) = cos(beta);			R(1,2) = -sin(alpha);
		R(2,0) = -sin(beta);		R(2,1) = sin(alpha);		R(2,2) = cos(gamma);*/
		Eigen::MatrixXd R1(3,3);
		R1(0,0) = 1;		R1(0,1) = 0;				R1(0,2) = 0;
		R1(1,0) = 0;		R1(1,1) = cos(alpha);		R1(1,2) = sin(alpha);
		R1(2,0) = 0;		R1(2,1) = -sin(alpha);		R1(2,2) = cos(alpha);
		Eigen::MatrixXd R2(3,3);
		R2(0,0) = cos(beta);		R2(0,1) = 0;		R2(0,2) = -sin(beta);
		R2(1,0) = 0;				R2(1,1) = 1;		R2(1,2) = 0;
		R2(2,0) = sin(beta);		R2(2,1) = 0;		R2(2,2) = cos(beta);
		Eigen::MatrixXd R3(3,3);
		R3(0,0) = cos(gamma);		R3(0,1) = sin(gamma);		R3(0,2) = 0;
		R3(1,0) = -sin(gamma);		R3(1,1) = cos(gamma);		R3(1,2) = 0;
		R3(2,0) = 0;				R3(2,1) = 0;				R3(2,2) = 1;
		return R1*R2*R3;
	}
	 InitLMCostFunction(InitLMParameters* paras)//Eigen::VectorXd& neutral_blendshape, const std::vector<int>& registerIndex, trimesh::TriMesh* pointCloudMesh, const std::vector<int>& nearestIndex, int num_s, int num_r, int num_t, int num_y, int num_z)//s, r, t, y, z
	 {
		 of.open("error_InitM.txt");
		 of<<"init "<<' '<<paras->registerIndex.size() + paras->num_y + paras->num_z*2<<' '<<paras->num_all<<std::endl;
		this->paras = paras;
		set_num_residuals(paras->registerIndex.size() + paras->num_y + paras->num_z*2);
		mutable_parameter_block_sizes()->push_back(paras->num_all);//paras->num_s + paras->num_r + paras->num_t + paras->num_y + paras->num_z);
		
		s = new double[paras->num_s];
		r = new double[paras->num_r];
		t = new double[paras->num_t];
		y = new double[paras->num_y];
		z = new double[paras->num_z];
		
	 }
	virtual ~InitLMCostFunction() 
	{
		if(s!=NULL)delete[]s;
		if(r!=NULL)delete[]r;
		if(t!=NULL)delete[]t;
		if(y!=NULL)delete[]y;
		if(z!=NULL)delete[]z;
	}
	
	/*virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
		int ii=0;
		std::ofstream oof("error_evalue.txt",std::ios::app);
		
		oof<<"Evaluate"<<std::endl;
		for(int i=0; i<paras->num_s; i++,ii++)s[i] = parameters[0][ii];
		for(int i=0; i<paras->num_r; i++,ii++)r[i] = parameters[0][ii];
		for(int i=0; i<paras->num_t; i++,ii++)t[i] = parameters[0][ii];
		for(int i=0; i<paras->num_y; i++,ii++)y[i] = parameters[0][ii];
		for(int i=0; i<paras->num_z; i++,ii++)z[i] = parameters[0][ii];
		oof<<"copy finish"<<std::endl;
		Eigen::VectorXd T(3);
		Eigen::MatrixXd R = getR(r);
		oof<<R<<std::endl;
		oof<<T<<std::endl;
		int idx=0;
		// ||A0(Rb0+t)-c0||
		//oof<<"jac: "<<jacobians[2][0]<<std::endl;
		//int cnt=0;
		
		//oof<<"e: " <<jacobians[0][208]<<std::endl;
		for(int i=0; i<paras->num_all; i++)
		{
			jacobians[0][i] = 0.0;
		}
		for(int i=0; i<paras->registerIndex.size(); i++,idx++)
		{
			oof<<"i1: " << i<<std::endl;
			int register_idx = paras->registerIndex[i];oof<<"i2: " << i<<std::endl;
			int nearest_idx = paras->nearestIndex[i];oof<<"i3: " << i<<std::endl;
			Eigen::VectorXd v(3);oof<<"i4: " << i<<std::endl;
			Eigen::VectorXd c(3);oof<<"i5: " << i<<std::endl;
			Eigen::VectorXd n(3);oof<<"i6: " << i<<std::endl;
			for(int j=0; j<3; j++)
			{
				oof<<"i7: " << i<<std::endl;
				v(j) = paras->neutral_blendshape(register_idx*3+j);
				c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
				n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			}
			Eigen::VectorXd new_v = (R*v)*s[0]+T;
			oof<<"i8: " << i<<std::endl;
			residuals[idx] = new_v.dot(n);
			
			oof<<"i9: " << i<<' '<<idx<<std::endl;
			int ii=0;
			oof<<R<<std::endl;
			oof<<v<<std::endl;
			oof<<n<<std::endl;
			oof<<residuals[idx]<<std::endl;
			oof<<(R*v).dot(n)<<std::endl<<std::endl;
			oof<<idx<<' '<<ii<<' '<<std::endl;
			//oof<<jacobians[2][0]<<' '<<jacobians[2][1]<<std::endl;
			//oof<<jacobians[0][idx*paras->num_all+ii]<<std::endl;
			oof<<(R*v).dot(n) * residuals[idx]<<std::endl;
			jacobians[0][ii] += (R*v).dot(n) * residuals[idx];//s
			oof<<"i91: " << i<<' '<<idx<<std::endl;
			ii++;
			jacobians[0][ii] += (-v(2)*c(1)+v(1)*c(2))*s[0] * residuals[idx];  //r alpha
			oof<<"i92: " << i<<' '<<idx<<std::endl;
			oof<<"all: " << paras->num_all<<std::endl;
			ii++;
			jacobians[0][ii] += (-v(0)*c(2)+v(2)*c(0))*s[0] * residuals[idx];   //r beta
			oof<<"i93: " << i<<' '<<idx<<std::endl;
			ii++;
			jacobians[0][ii] += (-v(1)*c(0)+v(0)*c(1))*s[0] * residuals[idx];  //r gamma
			ii++;
			
			for(int j=0; j<3; j++, ii++) //t
			{oof<<"i10: " << i <<' '<<j<<' '<<ii<<std::endl;
				jacobians[0][ii] += n(j) * residuals[idx];
			}
			
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{oof<<"i11: " << i<<' '<<j<<' '<<ii<<std::endl;
				Eigen::VectorXd Pi(3);
				for(int k=0; k<3; k++)
				{
					Pi(k) = paras->P_eigenvectors.coeff(j, i*3+k);
				}
				jacobians[0][ii] += n.dot(R*Pi*s[0]) * residuals[idx];
			}
			
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{oof<<"i12: " << i<<' '<<j<<' '<<ii<<std::endl;
				Eigen::VectorXd Ei(3);
				for(int k=0; k<3; k++)
				{
					Ei(k) = paras->E_eigenvectors.coeff(j, i*3+k);
				}
				jacobians[0][ii] += n.dot(R*Ei*s[0]) * residuals[idx];
			}
		}

		//||Dpy||
		for(int i=0; i<paras->num_y; i++,idx++)
		{
			residuals[idx] = paras->P_eigenvalues[i]*y[i];
			
			int ii=0;
			jacobians[0][ii] += 0;//s
			ii++;
			jacobians[0][ii] += 0;  //r alpha
			ii++;
			jacobians[0][ii] += 0;   //r beta
			ii++;
			jacobians[0][ii] += 0;  //r gamma
			ii++;
			for(int j=0; j<3; j++, ii++) //t
			{
				jacobians[0][ii] += 0;
			}
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{
				if(j==i)jacobians[0][ii] += paras->P_eigenvalues[j] * residuals[idx] * paras->beta1;
				else jacobians[0][ii] += 0;
			}
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{
				jacobians[0][ii] += 0;
			}
		}
		// ||Ez||
		for(int i=0; i<paras->num_z; i++,idx++)
		{
			residuals[idx] = paras->E_eigenvalues[i]*z[i];
					
			int ii=0;
			jacobians[0][ii] += 0;//s
			ii++;
			jacobians[0][ii] += 0;  //r alpha
			ii++;
			jacobians[0][ii] += 0;   //r beta
			ii++;
			jacobians[0][ii] += 0;  //r gamma
			ii++;
			for(int j=0; j<3; j++, ii++) //t
			{
				jacobians[0][ii] += 0;
			}
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{
				jacobians[0][ii] += 0;
			}
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{
				if(j==i)jacobians[0][ii] += paras->E_eigenvalues[j] * residuals[idx] * paras->beta2;
				else jacobians[0][ii] += 0;
			}
		}
		//||z||
		for(int i=0; i<paras->num_z; i++,idx++)
		{
			residuals[idx] = z[i];
			
			int ii=0;
			jacobians[0][ii] += 0;//s
			ii++;
			jacobians[0][ii] += 0;  //r alpha
			ii++;
			jacobians[0][ii] += 0;   //r beta
			ii++;
			jacobians[0][ii] += 0;  //r gamma
			ii++;
			for(int j=0; j<3; j++, ii++) //t
			{
				jacobians[0][ii] += 0;
			}
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{
				jacobians[0][ii] += 0;
			}
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{
				if(j==i)jacobians[0][ii] += 1.0 * residuals[idx] * paras->beta3;
				else jacobians[0][ii] += 0;
			}
		}
		oof<<"iter finish"<<std::endl;
		oof.close();
    return true;
  }*/
	virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
		int ii=0;
		std::ofstream oof("error_evalue.txt",std::ios::app);
		//oof<<paras->num_all<<' '<<paras->
		//oof<<"Evaluate"<<std::endl;
		for(int i=0; i<paras->num_s; i++,ii++)s[i] = parameters[0][ii];
		for(int i=0; i<paras->num_r; i++,ii++)r[i] = parameters[0][ii];
		for(int i=0; i<paras->num_t; i++,ii++)t[i] = parameters[0][ii];
		for(int i=0; i<paras->num_y; i++,ii++)y[i] = parameters[0][ii];
		for(int i=0; i<paras->num_z; i++,ii++)z[i] = parameters[0][ii];
		//oof<<"copy finish"<<std::endl;
		Eigen::VectorXd T(3);
		Eigen::MatrixXd R = getR(r);
		//oof<<R<<std::endl;
		//oof<<T<<std::endl;
		int idx=0;
		for(int i=0; i<3; i++)T(i) = t[i];
		// ||A0(Rb0+t)-c0||
		////oof<<"jac: "<<jacobians[2][0]<<std::endl;
		//int cnt=0;

		////oof<<"e: " <<jacobians[0][208]<<std::endl;
		for(int i=0; i<paras->registerIndex.size(); i++,idx++)
		{
			//oof<<"i1: " << i<<std::endl;
			int register_idx = paras->registerIndex[i];//oof<<"i2: " << i<<std::endl;
			int nearest_idx = paras->nearestIndex[i];//oof<<"i3: " << i<<std::endl;
			Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
			Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
			Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
			for(int j=0; j<3; j++)
			{
				//oof<<"i7: " << i<<std::endl;
				v(j) = paras->mean_neutral(register_idx*3+j);
				for(int k=0; k<paras->num_y; k++)
				{
					v(j) = v(j) + paras->P_eigenvectors(k, register_idx*3+j)*y[k];
				}
				for(int k=0; k<paras->num_z; k++)
				{
					v(j) = v(j) + paras->E_eigenvectors(k, register_idx*3+j)*z[k];
				}
				//v(j) = paras->neutral_blendshape(register_idx*3+j);
				c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
				n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			}
			Eigen::VectorXd new_v = (R*v)*s[0]+T-c;
			//oof<<"i8: " << i<<std::endl;
			residuals[idx] = new_v.dot(n);// * paras->nearestWeight[i];
			
			//oof<<"i9: " << i<<' '<<idx<<std::endl;
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				int ii=0;
				//oof<<R<<std::endl;
				//oof<<v<<std::endl;
				//oof<<n<<std::endl;
				//oof<<residuals[idx]<<std::endl;
				//oof<<(R*v).dot(n)<<std::endl<<std::endl;

				//oof<<idx<<' '<<ii<<' '<<' '<<idx*(paras->num_all)+ii<<' '<<jacobians[0][0]<<std::endl;
				////oof<<jacobians[2][0]<<' '<<jacobians[2][1]<<std::endl;
				////oof<<jacobians[0][idx*(paras->num_all)+ii]<<std::endl;
				//oof<<(R*v).dot(n) * residuals[idx]<<std::endl;
				jacobians[0][idx*(paras->num_all)+ii] = (R*v).dot(n) * paras->nearestWeight[i];// * residuals[idx];//s
				//oof<<"i91: " << i<<' '<<idx<<std::endl;
				ii++;
				/*
				jacobians[0][idx*(paras->num_all)+ii] = (-v(2)*n(1)+v(1)*n(2))*s[0] * paras->nearestWeight[i];// * residuals[idx];  //r alpha
				//oof<<"i92: " << i<<' '<<idx<<std::endl;
				//oof<<"all: " << paras->num_all<<std::endl;
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = (-v(0)*n(2)+v(2)*n(0))*s[0] * paras->nearestWeight[i];// * residuals[idx];   //r beta
				//oof<<"i93: " << i<<' '<<idx<<std::endl;
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = (-v(1)*n(0)+v(0)*n(1))*s[0] * paras->nearestWeight[i];// * residuals[idx];  //r gamma
				ii++;*/
				jacobians[0][idx*(paras->num_all)+ii] = getr1(r, v, n) * s[0] * paras->nearestWeight[i];// * residuals[idx];  //r alpha
				//oof<<"i92: " << i<<' '<<idx<<std::endl;
				//oof<<"all: " << paras->num_all<<std::endl;
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = getr2(r, v, n) * s[0] * paras->nearestWeight[i];// * residuals[idx];   //r beta
				//oof<<"i93: " << i<<' '<<idx<<std::endl;
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = getr3(r, v, n) * s[0] * paras->nearestWeight[i];// * residuals[idx];  //r gamma
				ii++;
			
				for(int j=0; j<3; j++, ii++) //t
				{//oof<<"i10: " << i <<' '<<j<<' '<<ii<<std::endl;
					jacobians[0][idx*(paras->num_all)+ii] = n(j) * paras->nearestWeight[i];
				}
			
				for(int j=0; j<paras->num_y; j++,ii++)//y
				{////oof<<"i11: " << i<<' '<<j<<' '<<ii<<std::endl;
					Eigen::VectorXd Pi(3);
					for(int k=0; k<3; k++)
					{
						Pi(k) = paras->P_eigenvectors.coeff(j, register_idx*3+k);
					}
					jacobians[0][idx*(paras->num_all)+ii] = n.dot(R*Pi*s[0]) * paras->nearestWeight[i];
				}
			
				for(int j=0; j<paras->num_z; j++, ii++)//z
				{////oof<<"i12: " << i<<' '<<j<<' '<<ii<<std::endl;
					Eigen::VectorXd Ei(3);
					//oof<<"z "<<j<<' ';
					for(int k=0; k<3; k++)
					{
						Ei(k) = paras->E_eigenvectors.coeff(j, register_idx*3+k);
					}
					//oof<<j<<' ';
					jacobians[0][idx*(paras->num_all)+ii] = 0;//n.dot(R*Ei*s[0]) * paras->nearestWeight[i];
					//oof<<j<<std::endl;
				}
			}
		}
	/*
		for(int i=0; i<paras->marker_registerIndex.size(); i++,idx++)
		{
			//oof<<"i1: " << i<<std::endl;
			int register_idx = paras->marker_registerIndex[i];//oof<<"i2: " << i<<std::endl;
			int nearest_idx = paras->marker_nearestIndex[i];//oof<<"i3: " << i<<std::endl;
			Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
			Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
			Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
			for(int j=0; j<3; j++)
			{
				//oof<<"i7: " << i<<std::endl;
				v(j) = paras->mean_neutral(register_idx*3+j);
				for(int k=0; k<paras->num_y; k++)
				{
					v(j) = v(j) + paras->P_eigenvectors(k, register_idx*3+j)*y[k];
				}
				for(int k=0; k<paras->num_z; k++)
				{
					v(j) = v(j) + paras->E_eigenvectors(k, register_idx*3+j)*z[k];
				}
				//v(j) = paras->neutral_blendshape(register_idx*3+j);
				c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
				n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			}
			Eigen::VectorXd new_v = (R*v)*s[0]+T-c;
			//oof<<"i8: " << i<<std::endl;
			residuals[idx] = new_v.dot(n);// * paras->marker_nearestWeight[i];
			
			//oof<<"i9: " << i<<' '<<idx<<std::endl;
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
			int ii=0;
			//oof<<R<<std::endl;
			//oof<<v<<std::endl;
			//oof<<n<<std::endl;
			//oof<<residuals[idx]<<std::endl;
			//oof<<(R*v).dot(n)<<std::endl<<std::endl;

			//oof<<idx<<' '<<ii<<' '<<' '<<idx*(paras->num_all)+ii<<' '<<jacobians[0][0]<<std::endl;
			////oof<<jacobians[2][0]<<' '<<jacobians[2][1]<<std::endl;
			////oof<<jacobians[0][idx*(paras->num_all)+ii]<<std::endl;
			//oof<<(R*v).dot(n) * residuals[idx]<<std::endl;
			jacobians[0][idx*(paras->num_all)+ii] = (R*v).dot(n) * paras->marker_nearestWeight[i];// * residuals[idx];//s
			//oof<<"i91: " << i<<' '<<idx<<std::endl;
			ii++;
			
			jacobians[0][idx*(paras->num_all)+ii] = (-n(0)*v(0)-n(1)*v(2)+n(2)*v(1)) * 0.5 * sin(2*r[0]) * s[0] * paras->marker_nearestWeight[i];// * residuals[idx];  //r alpha
			//oof<<"i92: " << i<<' '<<idx<<std::endl;
			//oof<<"all: " << paras->num_all<<std::endl;
			ii++;
			jacobians[0][idx*(paras->num_all)+ii] = (n(0)*v(2)-n(1)*v(1)-n(2)*v(0)) * 0.5 * sin(2*r[1]) * s[0] * paras->marker_nearestWeight[i];// * residuals[idx];   //r beta
			//oof<<"i93: " << i<<' '<<idx<<std::endl;
			ii++;
			jacobians[0][idx*(paras->num_all)+ii] = (-n(0)*v(1)+n(1)*v(0)-n(2)*v(2)) * 0.5 * sin(2*r[2]) * s[0] * paras->marker_nearestWeight[i];// * residuals[idx];  //r gamma
			ii++;
			
			for(int j=0; j<3; j++, ii++) //t
			{//oof<<"i10: " << i <<' '<<j<<' '<<ii<<std::endl;
				jacobians[0][idx*(paras->num_all)+ii] = n(j) * paras->marker_nearestWeight[i];
			}
			
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{////oof<<"i11: " << i<<' '<<j<<' '<<ii<<std::endl;
				Eigen::VectorXd Pi(3);
				for(int k=0; k<3; k++)
				{
					Pi(k) = paras->P_eigenvectors.coeff(j, register_idx*3+k);
				}
				jacobians[0][idx*(paras->num_all)+ii] = n.dot(R*Pi*s[0]) * paras->marker_nearestWeight[i];
			}
			
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{////oof<<"i12: " << i<<' '<<j<<' '<<ii<<std::endl;
				Eigen::VectorXd Ei(3);
				for(int k=0; k<3; k++)
				{
					Ei(k) = paras->E_eigenvectors.coeff(j, register_idx*3+k);
				}
				jacobians[0][idx*(paras->num_all)+ii] = n.dot(R*Ei*s[0]) * paras->marker_nearestWeight[i];
			}
			}
		}*/
		//||Dpy||
		for(int i=0; i<paras->num_y; i++,idx++)
		{
			residuals[idx] = y[i]/ paras->P_eigenvalues[i];
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				int ii=0;
				jacobians[0][idx*(paras->num_all)+ii] = 0;//s
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;  //r alpha
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;   //r beta
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;  //r gamma
				ii++;
				for(int j=0; j<3; j++, ii++) //t
				{
					jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
				for(int j=0; j<paras->num_y; j++,ii++)//y
				{
					if(j==i)jacobians[0][idx*(paras->num_all)+ii] = 1.0/paras->P_eigenvalues[j] * paras->beta1 ;
					else jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
				for(int j=0; j<paras->num_z; j++, ii++)//z
				{
					jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
			}
		}
		// ||Ez||
		for(int i=0; i<paras->num_z; i++,idx++)
		{
			residuals[idx] = paras->E_eigenvalues[i]*z[i];
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				int ii=0;
				jacobians[0][idx*(paras->num_all)+ii] = 0;//s
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;  //r alpha
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;   //r beta
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;  //r gamma
				ii++;
				for(int j=0; j<3; j++, ii++) //t
				{
					jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
				for(int j=0; j<paras->num_y; j++,ii++)//y
				{
					jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
				for(int j=0; j<paras->num_z; j++, ii++)//z
				{
					if(j==i)jacobians[0][idx*(paras->num_all)+ii] = paras->E_eigenvalues[j] * (paras->beta2) ;
					else jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
			}
		}
		//||z||
		for(int i=0; i<paras->num_z; i++,idx++)
		{
			residuals[idx] = z[i];
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				int ii=0;
				jacobians[0][idx*(paras->num_all)+ii] = 0;//s
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;  //r alpha
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;   //r beta
				ii++;
				jacobians[0][idx*(paras->num_all)+ii] = 0;  //r gamma
				ii++;
				for(int j=0; j<3; j++, ii++) //t
				{
					jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
				for(int j=0; j<paras->num_y; j++,ii++)//y
				{
					jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
				for(int j=0; j<paras->num_z; j++, ii++)//z
				{
					if(j==i)jacobians[0][idx*(paras->num_all)+ii] = 1.0 * paras->beta3;
					else jacobians[0][idx*(paras->num_all)+ii] = 0;
				}
			}
		}
		//oof<<"iter finish"<<std::endl;
		oof.close();
    return true;
  }
	
};

/*
virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
		int ii=0;
		std::ofstream oof("error_evalue.txt",std::ios::app);
		
		oof<<"Evaluate"<<std::endl;
		for(int i=0; i<paras->num_s; i++,ii++)s[i] = parameters[0][ii];
		for(int i=0; i<paras->num_r; i++,ii++)r[i] = parameters[0][ii];
		for(int i=0; i<paras->num_t; i++,ii++)t[i] = parameters[0][ii];
		for(int i=0; i<paras->num_y; i++,ii++)y[i] = parameters[0][ii];
		for(int i=0; i<paras->num_z; i++,ii++)z[i] = parameters[0][ii];
		oof<<"copy finish"<<std::endl;
		Eigen::VectorXd T(3);
		Eigen::MatrixXd R = getR(r);
		oof<<R<<std::endl;
		oof<<T<<std::endl;
		int idx=0;
		// ||A0(Rb0+t)-c0||
		//oof<<"jac: "<<jacobians[2][0]<<std::endl;
		//int cnt=0;

		//oof<<"e: " <<jacobians[0][208]<<std::endl;
		for(int i=0; i<paras->registerIndex.size(); i++,idx++)
		{
			oof<<"i1: " << i<<std::endl;
			int register_idx = paras->registerIndex[i];oof<<"i2: " << i<<std::endl;
			int nearest_idx = paras->nearestIndex[i];oof<<"i3: " << i<<std::endl;
			Eigen::VectorXd v(3);oof<<"i4: " << i<<std::endl;
			Eigen::VectorXd c(3);oof<<"i5: " << i<<std::endl;
			Eigen::VectorXd n(3);oof<<"i6: " << i<<std::endl;
			for(int j=0; j<3; j++)
			{
				oof<<"i7: " << i<<std::endl;
				v(j) = paras->neutral_blendshape(register_idx*3+j);
				c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
				n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			}
			Eigen::VectorXd new_v = (R*v)*s[0]+T;
			oof<<"i8: " << i<<std::endl;
			residuals[idx] = new_v.dot(n);
			
			oof<<"i9: " << i<<' '<<idx<<std::endl;
			int ii=0;
			oof<<R<<std::endl;
			oof<<v<<std::endl;
			oof<<n<<std::endl;
			oof<<residuals[idx]<<std::endl;
			oof<<(R*v).dot(n)<<std::endl<<std::endl;
			oof<<idx<<' '<<ii<<' '<<std::endl;
			oof<<jacobians[2][0]<<' '<<jacobians[2][1]<<std::endl;
			oof<<jacobians[0][idx*paras->num_all+ii]<<std::endl;
			oof<<(R*v).dot(n) * residuals[idx]<<std::endl;
			jacobians[0][idx*paras->num_all+ii] = (R*v).dot(n) * residuals[idx];//s
			oof<<"i91: " << i<<' '<<idx<<std::endl;
			ii++;
			jacobians[0][idx*paras->num_all+ii] = (-v(2)*c(1)+v(1)*c(2))*s[0] * residuals[idx];  //r alpha
			oof<<"i92: " << i<<' '<<idx<<std::endl;
			oof<<"all: " << paras->num_all<<std::endl;
			ii++;
			jacobians[0][idx*paras->num_all+ii] = (-v(0)*c(2)+v(2)*c(0))*s[0] * residuals[idx];   //r beta
			oof<<"i93: " << i<<' '<<idx<<std::endl;
			ii++;
			jacobians[0][idx*paras->num_all+ii] = (-v(1)*c(0)+v(0)*c(1))*s[0] * residuals[idx];  //r gamma
			ii++;
			
			for(int j=0; j<3; j++, ii++) //t
			{oof<<"i10: " << i <<' '<<j<<' '<<ii<<std::endl;
				jacobians[0][idx*paras->num_all+ii] = n(j) * residuals[idx];
			}
			
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{oof<<"i11: " << i<<' '<<j<<' '<<ii<<std::endl;
				Eigen::VectorXd Pi(3);
				for(int k=0; k<3; k++)
				{
					Pi(k) = paras->P_eigenvectors.coeff(j, i*3+k);
				}
				jacobians[0][idx*paras->num_all+ii] = n.dot(R*Pi*s[0]) * residuals[idx];
			}
			
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{oof<<"i12: " << i<<' '<<j<<' '<<ii<<std::endl;
				Eigen::VectorXd Ei(3);
				for(int k=0; k<3; k++)
				{
					Ei(k) = paras->E_eigenvectors.coeff(j, i*3+k);
				}
				jacobians[0][idx*paras->num_all+ii] = n.dot(R*Ei*s[0]) * residuals[idx];
			}
		}

		//||Dpy||
		for(int i=0; i<paras->num_y; i++,idx++)
		{
			residuals[idx] = paras->P_eigenvalues[i]*y[i];
			
			int ii=0;
			jacobians[0][idx*paras->num_all+ii] = 0;//s
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;  //r alpha
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;   //r beta
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;  //r gamma
			ii++;
			for(int j=0; j<3; j++, ii++) //t
			{
				jacobians[0][idx*paras->num_all+ii] = 0;
			}
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{
				if(j==i)jacobians[0][idx*paras->num_all+ii] = paras->P_eigenvalues[j] * residuals[idx] * paras->beta1;
				else jacobians[0][idx*paras->num_all+ii] = 0;
			}
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{
				jacobians[0][idx*paras->num_all+ii] = 0;
			}
		}
		// ||Ez||
		for(int i=0; i<paras->num_z; i++,idx++)
		{
			residuals[idx] = paras->E_eigenvalues[i]*z[i];
					
			int ii=0;
			jacobians[0][idx*paras->num_all+ii] = 0;//s
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;  //r alpha
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;   //r beta
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;  //r gamma
			ii++;
			for(int j=0; j<3; j++, ii++) //t
			{
				jacobians[0][idx*paras->num_all+ii] = 0;
			}
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{
				jacobians[0][idx*paras->num_all+ii] = 0;
			}
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{
				if(j==i)jacobians[0][idx*paras->num_all+ii] = paras->E_eigenvalues[j] * residuals[idx] * paras->beta2;
				else jacobians[0][idx*paras->num_all+ii] = 0;
			}
		}
		//||z||
		for(int i=0; i<paras->num_z; i++,idx++)
		{
			residuals[idx] = z[i];
			
			int ii=0;
			jacobians[0][idx*paras->num_all+ii] = 0;//s
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;  //r alpha
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;   //r beta
			ii++;
			jacobians[0][idx*paras->num_all+ii] = 0;  //r gamma
			ii++;
			for(int j=0; j<3; j++, ii++) //t
			{
				jacobians[0][idx*paras->num_all+ii] = 0;
			}
			for(int j=0; j<paras->num_y; j++,ii++)//y
			{
				jacobians[0][idx*paras->num_all+ii] = 0;
			}
			for(int j=0; j<paras->num_z; j++, ii++)//z
			{
				if(j==i)jacobians[0][idx*paras->num_all+ii] = 1.0 * residuals[idx] * paras->beta3;
				else jacobians[0][idx*paras->num_all+ii] = 0;
			}
		}
		oof<<"iter finish"<<std::endl;
		oof.close();
    return true;
  }
*/
#endif