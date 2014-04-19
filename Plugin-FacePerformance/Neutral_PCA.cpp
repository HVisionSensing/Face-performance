#include "Neutral_PCA.h"
Neutral_PCA::Neutral_PCA()
{
	FILE * of_mean = fopen("neutral_mean_face.txt","w");
	FILE * of_eigenvalue = fopen("neutral_eigenvalue.txt","w");
	FILE * of_eigenvector = fopen("neutarl_eigenvector.txt", "w");

	
	numNeutralMeshes = 150;
	numFeature_pca = 50;
	TriMesh* neutral_mesh = new TriMesh[numNeutralMeshes];
	for(int i=1; i<=numNeutralMeshes; i++)
	{
		QString prePath = "E:/download/FaceWareHouse/FaceWarehouse_Data_0_triangle/Tester_";
		QString sufPath = "/Blendshape/shape_0.obj";
		QString fileName = prePath.append(QString::number(i)).append(sufPath);
		if(OpenMesh::IO::read_mesh(neutral_mesh[i-1], fileName.toStdString()))
		{
			//fprintf(of_mean,"1");fflush(of_mean);
		}
		else 
		{
			fprintf(of_mean, "0");
			fprintf(of_mean, "%s\n", fileName.toStdString().c_str());
			fflush(of_mean);
		}
	}
	int numFeature_org = neutral_mesh[0].n_vertices();
	cv::Mat neutral_meanMat(1, numFeature_org*3, CV_64FC1);  //平均脸
	cv::Mat neutral_orgMat(numNeutralMeshes, numFeature_org*3, CV_64FC1);
	cv::Mat neutral_pcaMat(numFeature_pca, numFeature_org*3, CV_64FC1);
	fprintf(of_mean, "before");fflush(of_mean);
	Mesh2Mat(neutral_orgMat, neutral_mesh);
	fprintf(of_mean, "after");fflush(of_mean);

	//计算平均中性脸
	for(int i=0; i<neutral_meanMat.cols; i++)
	{
		neutral_meanMat.at<double>(1,i) = 0;
	}
	for(int i=0; i<neutral_orgMat.rows; i++)
	{
		double * p = neutral_orgMat.ptr<double>(i);
		for(int j=0; j<neutral_orgMat.cols; j++)
		{
			neutral_meanMat.at<double>(1,j) += p[j];
		}
	}
	for(int i=0; i<neutral_meanMat.cols; i++)
	{
		neutral_meanMat.at<double>(1,i) /= neutral_orgMat.rows;
		fprintf(of_mean, "%lf ", neutral_meanMat.at<double>(1,i));fflush(of_mean);
	}

	fprintf(of_mean, "\n%d %d \n", neutral_orgMat.rows, neutral_orgMat.cols);fflush(of_mean);
	neutral_pca = cv::PCA(neutral_orgMat, cv::Mat(), CV_PCA_DATA_AS_ROW, numFeature_pca);
	fprintf(of_mean, "\n%d %d \n", neutral_orgMat.rows, neutral_orgMat.cols);fflush(of_mean);
	neutral_pcaMat = neutral_pca.eigenvectors;
	//neutral_pca.mean;

	mean_neutral = Eigen::VectorXd(numFeature_org*3);
	for(int i=0; i<neutral_meanMat.cols; i++)
	{
		mean_neutral(i) = neutral_meanMat.at<double>(1,i);
	}
	P_eigenvalues = Eigen::VectorXd(numFeature_pca);
	P_eigenvectors = Eigen::MatrixXd(numFeature_pca, numFeature_org*3);
	for(int i=0; i<numFeature_pca; i++)
	{
		P_eigenvalues(i) = neutral_pca.eigenvalues.at<double>(i);
		double * p = neutral_pca.eigenvectors.ptr<double>(i);
		fprintf(of_eigenvalue, "%lf ", P_eigenvalues(i));fflush(of_eigenvalue);
		for(int j=0; j<neutral_pcaMat.cols; j++)
		{
			P_eigenvectors(i,j) = p[j];//neutral_pcaMat.at<double>(i,j);
			fprintf(of_eigenvector, "%lf ", p[j]);fflush(of_eigenvector);
		}
		fprintf(of_eigenvector, "\n");
	}
	fprintf(of_eigenvalue, "\n");
	fclose(of_mean);
	fclose(of_eigenvalue);
	fclose(of_eigenvector);
	delete [] neutral_mesh;
}
void Neutral_PCA::Mesh2Mat(cv::Mat& mat, TriMesh* mesh)
{
	for(int i=0; i<mat.rows; i++)
	{
		TriMesh::VertexIter v_it, v_end = mesh[i].vertices_end();
		for(v_it = mesh[i].vertices_begin(); v_it != v_end; ++v_it)
		{
			TriMesh::Point point = mesh[i].point(v_it);
			int idx = v_it->idx();
			for(int j=0; j<3; j++)
			{
				mat.at<double>(i,idx*3+j) = point[j];
			}
		}
	}
	FILE * of = fopen("neutral_mat.txt","w");
	for(int i=0; i<mat.rows; i++)
	{
		for(int j=0; j<mat.cols; j++)
			fprintf(of, "%lf ", mat.at<double>(i,j));
		fprintf(of, "\n");
	}
	fflush(of);
	fclose(of);
}

Neutral_PCA::~Neutral_PCA()
{
}
