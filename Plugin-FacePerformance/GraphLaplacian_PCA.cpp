#include "GraphLaplacian_PCA.h"
GraphLaplacian_PCA::GraphLaplacian_PCA()
{
	FILE * of_eigenvector = fopen("gl_eigenvector.txt", "w");
	FILE * of_eigenvalue = fopen("gl_eigenvalue.txt", "w");
	FILE * of  = fopen("error_laplacian.txt","w");
	fprintf(of, "0\n");
	fflush(of);
	std::string fileName = "E:/download/FaceWareHouse/FaceWarehouse_Data_0_triangle/Tester_1/Blendshape/shape_0.obj";
	fprintf(of, "1\n");
	fflush(of);
	TriMesh mesh;
	fprintf(of, "before readmesh\n");fflush(of);
	OpenMesh::IO::read_mesh(mesh, fileName);
	fprintf(of, "after readmesh\n");fflush(of);
	int num_vertices = mesh.n_vertices();
	int numFeature_org = mesh.n_vertices();
	int numFeature_pca = 50;
	fprintf(of, "before gl_orgMat");fflush(of);
	cv::Mat gl_orgMat = cv::Mat::zeros(num_vertices*1, num_vertices*1, CV_32FC1);
	fprintf(of, "after gl_orgMat");fflush(of);
	cv::Mat gl_pcaMat(num_vertices*1, num_vertices*1, CV_32FC1);
	fprintf(of, "before gl2mat");fflush(of);
	Mesh2GLMat(gl_orgMat, mesh);
	fprintf(of, "after gl2mat");fflush(of);
	gl_pca = cv::PCA(gl_orgMat, cv::Mat(), CV_PCA_DATA_AS_ROW, numFeature_org*1); //保留全部特征向量
	fprintf(of, "after pca");fflush(of);
	gl_pcaMat = gl_pca.eigenvectors;
	//gl_pca.project(gl_orgMat, gl_pcaMat);

	E_eigenvectors =Eigen::MatrixXd(numFeature_pca*1, num_vertices*1);
	E_eigenvalues = Eigen::VectorXd(numFeature_pca*1);
	//for(int i=gl_pcaMat.rows-numFeature_pca*1,ii=0; i<gl_pcaMat.rows; i++,ii++)
	for(int i=0,ii=0; i<numFeature_pca; i++,ii++)
	{
		E_eigenvalues(ii) = gl_pca.eigenvalues.at<float>(i);
		fprintf(of_eigenvalue, "%lf ", E_eigenvalues(i));fflush(of_eigenvalue);
		float * p = gl_pca.eigenvectors.ptr<float>(i);
		for(int j=0; j<gl_pcaMat.cols; j++)
		{
			E_eigenvectors(ii,j) = p[j];//gl_pcaMat.at<double>(i,j);
			fprintf(of_eigenvector, "%f ", p[j]);fflush(of_eigenvector);
		}
		fprintf(of_eigenvector, "\n");
	}
	fprintf(of_eigenvalue, "\n");
	fclose(of);
	fclose(of_eigenvalue);
	fclose(of_eigenvector);
	/*CvMat *Vector1;  
    CvMat *AvgVector;  
    CvMat *EigenValue_Row;  
    CvMat *EigenVector;  
  
    Vector1=cvCreateMat(10,2,CV_32FC1);  
    cvSetData(Vector1,Coordinates,Vector1->step);  
    AvgVector=cvCreateMat(1,2,CV_32FC1);  
    EigenValue_Row=cvCreateMat(2,1,CV_32FC1);  
    EigenVector=cvCreateMat(2,2,CV_32FC1); */
  
    //cvCalcPCA(Vector1,AvgVector,EigenValue_Row,EigenVector,CV_PCA_DATA_AS_ROW);  
	//for(int i=0; i<gl_
}
GraphLaplacian_PCA::~GraphLaplacian_PCA()
{
}
void GraphLaplacian_PCA::Mesh2GLMat(cv::Mat& mat, TriMesh& mesh)
{
	TriMesh::VertexIter v_it, v_end = mesh.vertices_end();
	for(v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		TriMesh::VertexOHalfedgeIter voh(mesh, v_it);
		float valence = mesh.valence(v_it);
		int row_idx = v_it->idx();
		for(int i=0; i<1; i++)
		{
			mat.at<float>(row_idx*1+i, row_idx*1+i) = valence;
		}
		for( ; voh; ++voh)
		{
			for(int i=0; i<1; i++)
			{
				mat.at<float>(row_idx*1+i, (mesh.to_vertex_handle(voh)).idx()) = -1.0;
			}
		}
	}
	FILE * of = fopen("mat.txt","w");
	//fprintf(of,"0\n");
	for(int i=0; i<mat.rows; i++)
	{
		for(int j=0; j<mat.cols; j++)
			fprintf(of,"%f ",mat.at<float>(i,j));
		fprintf(of,"\n");
	}
	//fprintf(of,"1\n");
	fflush(of);
	fclose(of);
}