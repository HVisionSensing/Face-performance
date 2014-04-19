#ifndef RIGIDLM_H_INCLUDED
#define RIGIDLM_H_INCLUDED
#include <vector>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include <TriMesh.h>
#include <fstream>
using ceres::CostFunction;
using ceres::SizedCostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class RigidLMParameters
{
public:
	
	int outIters;
	int maxInIters;
	double registerRate;
	double s;
	Eigen::VectorXd r;
	Eigen::VectorXd t;
	int num_s;
	int num_r;
	int num_t;
	int num_all;
	int start_s;
	int start_r;
	int start_t;
	
	
	std::vector<int> init_registerIndex;
	std::vector<int> init_nearestIndex;
	std::vector<double> init_nearestWeight;
	
	std::vector<int> registerIndex;
	std::vector<int> nearestIndex;
	std::vector<double> nearestWeight;

	Eigen::VectorXd expression_blendshape;
	trimesh::TriMesh* pointCloudMesh;
	
};
class RigidLMCostFunction: public CostFunction
{
 public:
	 std::ofstream of;
	 RigidLMParameters* paras;
	 double *s;
	 double *r;
	 double *t;
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
	 RigidLMCostFunction(RigidLMParameters* paras)//Eigen::VectorXd& neutral_blendshape, const std::vector<int>& registerIndex, trimesh::TriMesh* pointCloudMesh, const std::vector<int>& nearestIndex, int num_s, int num_r, int num_t, int num_y, int num_z)//s, r, t, y, z
	 {
		 of.open("error_RigidLM.txt");
		 of<<"init "<<' '<<paras->registerIndex.size()<<' '<<paras->num_all<<std::endl;
		this->paras = paras;
		set_num_residuals(paras->registerIndex.size());
		mutable_parameter_block_sizes()->push_back(paras->num_all);//paras->num_s + paras->num_r + paras->num_t + paras->num_y + paras->num_z);
		
		s = new double[paras->num_s];
		r = new double[paras->num_r];
		t = new double[paras->num_t];
		of.close();
	 }
	virtual ~RigidLMCostFunction() 
	{
		if(s!=NULL)delete[]s;
		if(r!=NULL)delete[]r;
		if(t!=NULL)delete[]t;
	}
	
	
	virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
		int ii=0;
		std::ofstream oof("error_rigidevalue.txt",std::ios::app);
		//oof<<paras->num_all<<' '<<paras->
		//oof<<"Evaluate"<<std::endl;
		for(int i=0; i<paras->num_s; i++,ii++)s[i] = parameters[0][ii];
		for(int i=0; i<paras->num_r; i++,ii++)r[i] = parameters[0][ii];
		for(int i=0; i<paras->num_t; i++,ii++)t[i] = parameters[0][ii];
		
		//oof<<"copy finish"<<std::endl;
		Eigen::VectorXd T(3);
		Eigen::MatrixXd R = getR(r);
		//oof<<R<<std::endl;
		//oof<<T<<std::endl;
		int idx=0;
		for(int i=0; i<3; i++)T(i) = t[i];

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
				v(j) = paras->expression_blendshape(register_idx*3+j);
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