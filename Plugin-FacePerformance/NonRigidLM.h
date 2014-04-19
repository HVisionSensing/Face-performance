#ifndef NONRIGIDLM_H_INCLUDED
#define NONRIGIDLM_H_INCLUDED
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

class NonRigidLMParameters
{
public:
	int allIters;
	int outIters;
	int maxInIters;
	double registerRate;
	double sigma;//image
	double alpha;
	double beta;
	double gamma;
	double tao;
	double primal_residual_sqr_norm_eps;
	double dual_residual_sqr_norm_eps;

	Eigen::VectorXd w;
	Eigen::VectorXd v;
	Eigen::VectorXd d_w;
	Eigen::VectorXd d_v;
	double mu_w;
	double mu_v;

	double s;
	Eigen::MatrixXd R;
	Eigen::VectorXd T;

	int start_x;
	int num_x;
	int num_all;

	cv::Mat pre_colorImage;
	cv::Mat colorImage;
	double * colorCoordinate;
	double * pre_colorCoordinate;
	double maxDepth; //900
	double fx;
	double fy;
	Eigen::MatrixXd _3d_to_2d_matrix;
	Eigen::MatrixXd pre_3d_to_2d_matrix;
	int num_frame;
	//double imgCol; //
	//double imgRow;
	/*
	int num_s;
	int num_r;
	int num_t;
	
	int start_s;
	int start_r;
	int start_t;*/
	
	std::vector<int> init_registerIndex;
	std::vector<int> init_nearestIndex;
	std::vector<double> init_nearestWeight;
	
	std::vector<int> registerIndex;
	std::vector<int> nearestIndex;
	std::vector<double> nearestWeight;

	std::vector<int> marker_registerIndex;
	std::vector<int> marker_nearestIndex;
	std::vector<double> marker_nearestWeight;

	double markerWeight;
	double * marker_Coordinates;
	
	std::vector<Eigen::VectorXd> I;
	std::vector<Eigen::VectorXd> J;
	std::vector<Eigen::MatrixXd> derive_J;
	std::vector<Eigen::MatrixXd> derive_f;

	Eigen::VectorXd init_expression_blendshape;
	Eigen::VectorXd expression_blendshape;
	Eigen::VectorXd neutral_blendshape;
	Eigen::MatrixXd delta_B;
	Eigen::VectorXd pre_x;
	Eigen::VectorXd ppre_x;
	trimesh::TriMesh* pointCloudMesh;
	std::string meshIndexFilename;
	std::string depthIndexFilename;
};

class NonRigidLMCostFunction: public CostFunction
{
 public:
	 std::ofstream of;
	 NonRigidLMParameters* paras;
	 double* x;
	 static Eigen::MatrixXd derive_f(double x, double y, double z, cv::Mat& colorImage, double fx, double fy, double maxDepth, double * colorCoordinate, Eigen::MatrixXd& _3d_to_2d_matrix)
	 {
		 Eigen::MatrixXd deriveMat(2, 3);
		 for(int i=0; i<2; i++)
		 {
			 for(int j=0; j<3; j++)
			 {
				 deriveMat(i,j) = _3d_to_2d_matrix(j,i);
			 }
		 }
		 return deriveMat;
		 /*
		 Eigen::MatrixXd deriveMat(2, 3);
		 deriveMat(0,1) = deriveMat(1,0) = 0.0;
		 deriveMat(0,0) = 1.0/fx/(maxDepth-z)*colorImage.cols;
		 deriveMat(0,2) = -x/fx/(maxDepth-z)/(maxDepth-z)*colorImage.cols;
		 deriveMat(1,1) = -1.0/fy/(maxDepth-z)*colorImage.rows;
		 deriveMat(1,2) = y/fy/(maxDepth-z)/(maxDepth-z)*colorImage.rows;
		 return deriveMat;*/
	 }
	 static Eigen::VectorXd f_color(double x, double y, double z, cv::Mat& colorImage, double fx, double fy, double maxDepth, double * colorCoordinate, Eigen::MatrixXd& _3d_to_2d_matrix)
	 {
		 Eigen::VectorXd color(3);//三通道
		 double texture_x = x * _3d_to_2d_matrix(0,0) + y * _3d_to_2d_matrix(1,0) + z * _3d_to_2d_matrix(2,0) + _3d_to_2d_matrix(3,0);
		 double texture_y = x * _3d_to_2d_matrix(0,1) + y * _3d_to_2d_matrix(1,1) + z * _3d_to_2d_matrix(2,1) + _3d_to_2d_matrix(3,1);
		 texture_x *= 90;
		 texture_y *= 90;
		 if(texture_x<0)texture_x = 0;
		 else if(texture_x>=colorImage.cols)texture_x = colorImage.cols-1;
		 if(texture_y<0)texture_y = 0;
		 else if(texture_y>=colorImage.rows)texture_y = colorImage.rows-1;

		 for(int i=0; i<3; i++)
		 {
			// output[ll_l].at<Vec<unsigned char,3> >(ii,jj)[kk];
			 color(i) = colorImage.at<cv::Vec<unsigned char,3> >((int)texture_y, (int)texture_x)[i];
		 }
		 return color;
		 /*Eigen::VectorXd color(3);//三通道
		 int m = (0.5-y/fy/(maxDepth-z))*colorImage.rows+0.5;
		 int n = (x/fx/(maxDepth-z)+0.5)*colorImage.cols+0.5;
		 if(m>=colorImage.rows) m = colorImage.rows-1;
		 if(m<0)m = 0;
		 if(n>=colorImage.cols) n = colorImage.cols-1;
		 if(n<0)n = 0;

		 int k = m*colorImage.cols + n;
		 n = colorCoordinate[2*k];
		 m = colorCoordinate[2*k+1];
		 for(int i=0; i<3; i++)
		 {
			// output[ll_l].at<Vec<unsigned char,3> >(ii,jj)[kk];
			 color(i) = colorImage.at<cv::Vec<unsigned char,3> >(m,n)[i];
		 }
		 return color;*/
	 }
	 static Eigen::MatrixXd derive_color(double x, double y, double z, cv::Mat& colorImage, double fx, double fy, double maxDepth, double * colorCoordinate, Eigen::MatrixXd& _3d_to_2d_matrix)
	 {
		 Eigen::MatrixXd derive_color(3, 2);//三通道
		 double texture_x = x * _3d_to_2d_matrix(0,0) + y * _3d_to_2d_matrix(1,0) + z * _3d_to_2d_matrix(2,0) + _3d_to_2d_matrix(3,0);
		 double texture_y = x * _3d_to_2d_matrix(0,1) + y * _3d_to_2d_matrix(1,1) + z * _3d_to_2d_matrix(2,1) + _3d_to_2d_matrix(3,1);
		 texture_x *= 90;
		 texture_y *= 90;
		 if(texture_x<0)texture_x = 0;
		 else if(texture_x>=colorImage.cols)texture_x = colorImage.cols-1;
		 if(texture_y<0)texture_y = 0;
		 else if(texture_y>=colorImage.rows)texture_y = colorImage.rows-1;

		 int n = texture_x;
		 int m = texture_y;
		 for(int i=0; i<3; i++)
		 {
			 if(n+1<colorImage.cols && n>=0)
				derive_color(i,0) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n+1)[i] - (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n)[i];
			 else derive_color(i,0) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n)[i] - (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n-1)[i];
			 if(m+1<colorImage.rows && m>=0)
				derive_color(i,1) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m+1, n)[i]-(double)colorImage.at<cv::Vec<unsigned char, 3> >(m,n)[i];
			 else derive_color(i,1) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n)[i] - (double)colorImage.at<cv::Vec<unsigned char, 3> >(m-1, n)[i];
		 }
		 return derive_color;
		 /*Eigen::MatrixXd derive_color(3, 2);//三通道
		 int m = (0.5-y/fy/(maxDepth-z))*colorImage.rows+0.5;
		 int n = (x/fx/(maxDepth-z)+0.5)*colorImage.cols+0.5;
		 if(m>=colorImage.rows) m = colorImage.rows-1;
		 if(m<0)m = 0;
		 if(n>=colorImage.cols) n = colorImage.cols-1;
		 if(n<0)n = 0;
		 int k = m*colorImage.cols + n;

		 n = colorCoordinate[2*k];
		 m = colorCoordinate[2*k+1];
		 for(int i=0; i<3; i++)
		 {
			 if(n+1<colorImage.rows && n>=0)
				derive_color(i,0) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n+1)[i] - (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n)[i];
			 else derive_color(i,0) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n)[i] - (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n-1)[i];
			if(m+1<colorImage.cols && m>=0)
				derive_color(i,1) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m+1, n)[i]-(double)colorImage.at<cv::Vec<unsigned char, 3> >(m,n)[i];
			else derive_color(i,1) = (double)colorImage.at<cv::Vec<unsigned char, 3> >(m, n)[i] - (double)colorImage.at<cv::Vec<unsigned char, 3> >(m-1, n)[i];
		 }
		 return derive_color;*/
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
	 NonRigidLMCostFunction(NonRigidLMParameters* paras)//Eigen::VectorXd& neutral_blendshape, const std::vector<int>& registerIndex, trimesh::TriMesh* pointCloudMesh, const std::vector<int>& nearestIndex, int num_s, int num_r, int num_t, int num_y, int num_z)//s, r, t, y, z
	 {
		 of.open("error_NonRigidLM.txt");
		 of<<"init "<<' '<<paras->registerIndex.size()<<' '<<paras->num_all<<std::endl;
		this->paras = paras;
		set_num_residuals(paras->registerIndex.size() + paras->registerIndex.size()*3 + (paras->num_x)*3);
		mutable_parameter_block_sizes()->push_back(paras->num_all);//paras->num_s + paras->num_r + paras->num_t + paras->num_y + paras->num_z);
		
		x = new double[paras->num_x];
		of.close();
	 }
	virtual ~NonRigidLMCostFunction() 
	{
		if(x!=NULL)delete[]x;
	}
	
	
	virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const {
		std::ofstream oof("error_nonrigidevalue.txt");//,std::ios::app);
		//oof<<paras->num_all<<' '<<paras->
		//oof<<"Evaluate"<<std::endl;
		for(int i=0; i<paras->num_x; i++)x[i] = parameters[0][i];

		int idx=0;
		for(int i=0; i<paras->registerIndex.size(); i++,idx++)
		{
			//oof<<"i1: " << i<<std::endl;
			int register_idx = paras->registerIndex[i]*3;//oof<<"i2: " << i<<std::endl;
			int nearest_idx = paras->nearestIndex[i];//oof<<"i3: " << i<<std::endl;
			Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
			Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
			Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
			for(int j=0; j<3; j++)
			{
				//oof<<"i7: " << i<<std::endl;
				v(j) = paras->expression_blendshape(register_idx+j);
				c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
				n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			}
			Eigen::VectorXd new_v = (paras->R * v) * paras->s + paras->T - c;
			//oof<<"i8: " << i<<std::endl;
			residuals[idx] = new_v.dot(n);// * paras->nearestWeight[i];
			
			//oof<<"i9: " << i<<' '<<idx<<std::endl;
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				for(int j=0; j<paras->num_x; j++)
				{
					Eigen::VectorXd b(3);
					for(int k=0; k<3; k++)
					{
						b(k) = paras->delta_B.coeff(register_idx+k, j);
					}
					jacobians[0][idx*(paras->num_all)+j] = (paras->R * b).dot(n) * paras->s; 
				}
			}
		}
		//Image
		for(int i=0; i<paras->registerIndex.size(); i++)
		{
			oof<<"i1: " << i<<std::endl;
			int register_idx = paras->registerIndex[i]*3;//oof<<"i2: " << i<<std::endl;
			int nearest_idx = paras->nearestIndex[i];//oof<<"i3: " << i<<std::endl;
			Eigen::VectorXd init_v(3);
			Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
			Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
			Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
			for(int j=0; j<3; j++)
			{
				oof<<"i4: " << j<<std::endl;
				init_v(j) = paras->init_expression_blendshape(register_idx+j);oof<<"i4: " << j<<std::endl;
				v(j) = paras->expression_blendshape(register_idx+j);oof<<"i4: " << j<<std::endl;
				c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];oof<<"i4: " << j<<std::endl;
				n(j) = paras->pointCloudMesh->normals[nearest_idx][j];oof<<"i4: " << j<<std::endl;
			}
			init_v = (paras->R*init_v)*paras->s + paras->T;
			oof<<"i5: " << i<<std::endl;
			v = (paras->R*v)*paras->s + paras->T;
			oof<<"i6: " << i<<std::endl;
			Eigen::VectorXd init_color = f_color(init_v(0), init_v(1), init_v(2), paras->pre_colorImage, paras->fx, paras->fy, paras->maxDepth, paras->pre_colorCoordinate, paras->pre_3d_to_2d_matrix);
			oof<<"i7: " << i<<std::endl;
			Eigen::VectorXd color = f_color(v(0), v(1), v(2), paras->colorImage, paras->fx, paras->fy, paras->maxDepth, paras->colorCoordinate, paras->_3d_to_2d_matrix);
			oof<<"i8: " << i<<std::endl;
			Eigen::MatrixXd derive_mat = derive_color(v(0), v(1), v(2), paras->colorImage, paras->fx, paras->fy, paras->maxDepth, paras->colorCoordinate, paras->_3d_to_2d_matrix)
									   * derive_f(v(0), v(1), v(2), paras->colorImage, paras->fx, paras->fy, paras->maxDepth, paras->colorCoordinate, paras->_3d_to_2d_matrix);
			oof<<"i9: " << i<<std::endl;
			Eigen::VectorXd new_v(3);
			for(int j=0; j<3; j++) 
			{
				oof<<"i10: " << i<<std::endl;
				new_v(j) = paras->neutral_blendshape(register_idx+j);
				for(int k=0; k<paras->num_x; k++)
				{
					new_v(j) += paras->delta_B(register_idx+j, k)*x[k];
				}
			}
			new_v = paras->R * new_v * paras->s + paras->T;oof<<"i11: " << i<<std::endl;
			oof<<"derive_mat: "<<derive_mat<<std::endl;
			oof<<"init_color: "<<init_color.transpose()<<std::endl;
			oof<<"color: "<<color.transpose()<<std::endl;
			oof<<"new_v: "<<new_v.transpose()<<std::endl;
			oof<<"v: "<<v.transpose()<<std::endl;
			Eigen::VectorXd d_color = init_color - color - derive_mat*(new_v-v);oof<<"i12: " << i<<std::endl;
			for(int h=0; h<3; h++,idx++)
			{
				residuals[idx] = d_color(h);
				if(jacobians!=NULL && jacobians[0] != NULL)
				{
					for(int j=0; j<paras->num_x; j++)
					{
						Eigen::VectorXd b(3);
						for(int k=0; k<3; k++)
						{
							b(k) = paras->delta_B.coeff(register_idx+k, j);
						}
						Eigen::VectorXd vec = -derive_mat*paras->R*paras->s*b;
						jacobians[0][idx*(paras->num_all)+j] = vec(h) / (paras->sigma*paras->sigma);//(paras->R * b).dot(n) * paras->s; 
					}
				}
			}
		}
		//smooth
		for(int i=0; i<paras->num_x; i++,idx++)
		{
			residuals[idx] = paras->ppre_x[i] - paras->pre_x[i]*2 + x[i];
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				for(int j=0; j<paras->num_x; j++)
				{
					if(i==j)jacobians[0][idx*(paras->num_all)+j] = paras->alpha;
					else jacobians[0][idx*(paras->num_all)+j] = 0.0;
				}
			}
		}
		for(int i=0; i<paras->num_x; i++, idx++)
		{
			residuals[idx] = paras->w[i] - x[i] + paras->d_w[i];
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				for(int j=0; j<paras->num_x; j++)
				{
					if(i==j)jacobians[0][idx*(paras->num_all)+j] = -paras->mu_w/2;
					else jacobians[0][idx*(paras->num_all)+j] = 0.0;
				}
			}
		}
		for(int i=0; i<paras->num_x; i++, idx++)
		{
			residuals[idx] = paras->v[i] - x[i] + paras->d_v[i];
			if(jacobians!=NULL && jacobians[0] != NULL)
			{
				for(int j=0; j<paras->num_x; j++)
				{
					if(i==j)jacobians[0][idx*(paras->num_all)+j] = -paras->mu_v/2;
					else jacobians[0][idx*(paras->num_all)+j] = 0.0;
				}
			}
		}
		
		oof.close();
		return true;
  }
	
};

#endif