#include "TrackingMesh.h"
/*
TrackingMesh::TrackingMesh(trimesh::TriMesh* pointCloudMesh)
{
	of.open("error_TrackingMesh.txt");
	beta1 = 0.5;
	beta2 = 0.1;
	beta3 = 0.001;
	of<<"before neutral_pca"<<std::endl;
	//neutral_pca = Neutral_PCA();
	of<<"after neutral_pca"<<std::endl;
	//gl_pca = GraphLaplacian_PCA();
	of<<"after gl_pca"<<std::endl;
	pointCloudMesh->write("pvM.obj");////////////////////////
	this->pointCloudMesh = pointCloudMesh;
	pclCloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);
	trimesh2pcl(pointCloudMesh, pclCloud);
	pclCloudKdtree.setInputCloud(pclCloud);
	of<<"before load1"<<std::endl;
	readRegisterIndexFromFile(registerIndex, "non_rigid_selection_index.txt");of<<"before load2"<<std::endl;
	loadEigenVector(mean_neutral, "neutral_meanface.txt");of<<"before load3"<<std::endl;// = neutral_pca.mean_neutral;
	loadEigenVector(P_eigenvalues, "neutral_eigenvalue.txt");of<<"before load4"<<std::endl;// = neutral_pca.P_eigenvalues;
	loadEigenMatrix(P_eigenvectors, "neutral_eigenvector.txt");of<<"before load5"<<std::endl;// = neutral_pca.P_eigenvectors;
	loadEigenVector_gl(E_eigenvalues, "gl_eigenvalue.txt");of<<"before load6"<<std::endl;// = gl_pca.E_eigenvalues;
	loadEigenMatrix_gl(E_eigenvectors, "gl_eigenvector.txt");of<<"before load7"<<std::endl;// = gl_pca.E_eigenvectors;
	of<<"after load"<<std::endl;

	treeCloud = new trimesh::KDtree(pointCloudMesh->vertices);of<<"after load1"<<std::endl;
	pointCloudMesh->need_normals();of<<"after load2"<<std::endl;
	pointCloudMesh->need_bsphere();of<<"after load3"<<std::endl;

	//Î´ÖªÊý R(3), t(3), y(50), z(150)
	numR = 3;
	numt = 3;
	numy = P_eigenvalues.size();of<<"P_eigenvalues size: "<<numy<<std::endl;
	numz = E_eigenvalues.size();of<<"E_eigenvalues size: "<<numz<<std::endl;
	numX = numy + numz + numR + numt;
	startR = 0;
	startt = startR + numR;
	starty = startt + numt;
	startz = starty + numy;
	
	x = new double[numX];
	
	neutral_blendshape.resize(mean_neutral.size());
	numz = E_eigenvalues.size();
	of<<"before initRtyz: "<<numy<<std::endl;
	initRtyz();
	update_all();
	trimesh::TriMesh * cur_mesh = new trimesh::TriMesh;
	blendshape2trimesh(cur_mesh);
	cur_mesh->write("cur_mesh.obj");
	delete cur_mesh;
	for(int i=0; i<6; i++)of<<x[i]<<std::endl;
	of<<"Construct finish"<<std::endl;
}

TrackingMesh::~TrackingMesh()
{
	if(treeCloud != NULL)
		delete treeCloud;
	if(x != NULL)
		delete [] x;
}
void TrackingMesh::readRegisterIndexFromFile(std::vector<int>& index, std::string filename)
{
	index.clear();

	ifstream in(filename);
	int idx;
	while(in>>idx)
	{
		index.push_back(idx);
	}
	in.close();
}
void TrackingMesh::initRtyz()
{
	std::ifstream fmesh("MeshIndex.txt");
	std::ifstream fdepth("DepthIndex.txt");
	std::vector<int> m_DepthIndex;
	std::vector<int> m_MeshIndex;
	for(int i=0; i<7; i++)
	{
		int m_d, m_m;
		fdepth>>m_d;
		fmesh>>m_m;
		if(m_d != -1)
		{
			m_MeshIndex.push_back(m_m);
			m_DepthIndex.push_back(m_d);
		}
		of<<m_d<<' '<<m_m<<std::endl;
	}
	fmesh.close();fdepth.close();

	trimesh::point min_val, max_val, mean_val(0,0,0);
	for(int i=0; i<3; i++)
	{
		min_val[i] = 1.0e6f;
		max_val[i] = -1.0e6f;
	}
	for(int ii=0; ii<m_MeshIndex.size(); ii++)
	{
		int i = m_MeshIndex[ii];
		of<<"mean_neutral: ";
		for(int j=0; j<3; j++)
		{
			of<<mean_neutral(i*3+j)<<' ';
			if(min_val[j]>mean_neutral(i*3+j))
			{
				min_val[j] = mean_neutral(i*3+j);
			}
			if(max_val[j]<mean_neutral(i*3+j))
			{
				max_val[j] = mean_neutral(i*3+j);
			}
			mean_val[j] += mean_neutral(i*3+j);
		}
		of<<std::endl;
	}
	for(int i=0; i<3; i++) mean_val[i] /= m_MeshIndex.size();

	trimesh::point min_val_pointcloud(1.0e6f, 1.0e6f, 1.0e6f), max_val_pointcloud(-1.0e6f, -1.0e6f, -1.0e6f), mean_val_pointcloud(0,0,0);
	//mean_val_pointcloud = pointCloudMesh->bsphere.center;
	for(int i=0; i<3; i++)
	{
		min_val_pointcloud[i] = 1.0e6f;
		max_val_pointcloud[i] = -1.0e6f;
	}
	for(int ii=0; ii<m_DepthIndex.size(); ii++)
	{
		int i = m_DepthIndex[ii];
		of<<"pointCloud: "<<pointCloudMesh->vertices[i]<<std::endl;
		for(int j=0; j<3; j++)
		{
			if(min_val_pointcloud[j]>pointCloudMesh->vertices[i][j])
			{
				min_val_pointcloud[j]=pointCloudMesh->vertices[i][j];
			}
			if(max_val_pointcloud[j]<pointCloudMesh->vertices[i][j])
			{
				max_val_pointcloud[j]=pointCloudMesh->vertices[i][j];
			}
			mean_val_pointcloud[j] += pointCloudMesh->vertices[i][j];
		}
	}
	for(int i=0; i<3; i++) mean_val_pointcloud[i] /= m_DepthIndex.size();

	
	for(int i=0; i<3; i++)
	{
		of<<i<<": "<<max_val_pointcloud[i]<<' '<<min_val_pointcloud[i]<<' '<<max_val[i]<<' '<<min_val[i]<<' '
			<<mean_val_pointcloud[i]<<' '<<mean_val[i]<<std::endl;
		x[i] = (max_val_pointcloud[i]-min_val_pointcloud[i])/(max_val[i]-min_val[i]);//R
		if(i==2)x[i] = x[1];////////////////////¶Ô½ÇÕó
		x[i+3] = mean_val_pointcloud[i]-mean_val[i]*x[i];//t
	}
	for(int i=6; i<numX; i++)x[i] = 0;//y z
	of<<"init finish"<<std::endl;

}
void TrackingMesh::function_f(double * x, double* f_sum)
{
	Eigen::VectorXd f_val;
	double sum = 0;
	std::vector<trimesh::point>& pointcloud_normals = pointCloudMesh->normals;
	std::vector<trimesh::point>& pointcloud_vertices = pointCloudMesh->vertices;

	f_val.resize(registerIndex.size() + P_eigenvalues.size() + E_eigenvalues.size()*2);
	int idx = 0;
	
	for(int i=0; i<registerIndex.size(); i++,idx++)
	{
		f_val(idx) = 0;
		int register_index = registerIndex[i];
		int nearest_index = nearestIndex[i];
		for(int j=0; j<3; j++)
		{
			//(Rv+t-c)n
			f_val(idx) += (x[startR+j]*neutral_blendshape(register_index*3+j)+x[startt+j]-pointcloud_vertices[nearest_index][j]) * pointcloud_normals[nearest_index][j];
		}
		f_val(idx) *= f_val(idx);
		sum += f_val(idx);
	}
	for(int i=0; i<numy; i++,idx++)
	{
		f_val(idx) = x[starty+i]*P_eigenvalues(i);
		f_val(idx) *= f_val(idx);
		f_val(idx) *= beta1;
		sum += f_val(idx);
	}
	for(int i=0; i<numz; i++, idx++)
	{
		f_val(idx) = x[startz+i]*E_eigenvalues(i);
		f_val(idx) *= f_val(idx);
		f_val(idx) *= beta2;
		sum += f_val(idx);
	}
	for(int i=0; i<numz; i++, idx++)
	{
		f_val(idx) = x[startz+i];
		f_val(idx) *= f_val(idx);
		f_val(idx) *= beta3;
		sum += f_val(idx);
	}
	*f_sum = sum;
	//double a = f_val.dot(f_val);
}
void TrackingMesh::function_df(double *x, double *g)
{
	Eigen::VectorXd g_val;
	std::vector<trimesh::point>& pointcloud_normals = pointCloudMesh->normals;
	std::vector<trimesh::point>& pointcloud_vertices = pointCloudMesh->vertices;
	g_val.resize(numX);
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<registerIndex.size(); j++)
		{	
			int register_index = registerIndex[j];
			int nearest_index = nearestIndex[j];
			double val = 0;
			for(int k=0; k<3; k++)
			{
				val += (x[startR+k]*neutral_blendshape(register_index*3+k)+x[startt+k]-pointcloud_vertices[nearest_index][k])*pointcloud_normals[nearest_index][k];
			}
			g_val(i) = 2 * val * (neutral_blendshape(register_index*3+i)*pointcloud_normals[nearest_index][i]);//R
			g_val(i+3) = 2 * val * pointcloud_normals[nearest_index][i];//t
		}
	}

	Eigen::VectorXd dy(numy);
	Eigen::VectorXd dz(numz);

	for(int i=0; i<numy; i++)dy(i) = 0;
	for(int i=0; i<numz; i++)dz(i) = 0;
	
	for(int i=0; i<registerIndex.size(); i++)
	{
		double val = 0;
		int register_index = registerIndex[i];
		int nearest_index = nearestIndex[i];
		for(int j=0; j<3; j++)
		{
			val += pointcloud_normals[nearest_index][j]*(x[startR+j]*neutral_blendshape(register_index*3+j)+x[startt+j]-pointcloud_vertices[nearest_index][j]);
		}
		for(int j=0; j<3; j++)
		{
			dy = dy + 2*val*x[startR+j]*pointcloud_normals[nearest_index][j] * P_eigenvectors.col(register_index*3+j);
			dz = dz + 2*val*x[startR+j]*pointcloud_normals[nearest_index][j] * E_eigenvectors.col(register_index*3+j);
		}
	}
	for(int i=0; i<numy; i++)
	{
		dy(i) += 2*beta1*x[starty+i]*P_eigenvalues[i]*P_eigenvalues[i];
	}
	for(int i=0; i<numz; i++)
	{
		dz(i) += 2*beta2*x[startz+i]*E_eigenvalues[i]*E_eigenvalues[i];
		dz(i) += 2*beta3*x[startz+i];
	}
	for(int i=0; i<numy; i++)
	{
		g_val(starty+i) = dy.coeff(i);
	}
	for(int i=0; i<numz; i++)
	{
		g_val(startz+i) = dz.coeff(i);
	}
	for(int i=0; i<numX; i++)
	{
		g[i] = g_val.coeff(i);
	}
	for(int i=0; i<numX; i++)of<<g[i]<<' ';of<<std::endl;
}

void TrackingMesh::update_all()
{
	update_blendshape();
	find_nearest();

}
void TrackingMesh::update_blendshape()
{
	neutral_blendshape = mean_neutral;
	//of<<"neutral_blend: "<<neutral_blendshape.size()<<' '<<neutral_blendshape.cols()<<' '<<neutral_blendshape.rows()<<endl;
	for(int i=0; i<numX; i++)of<<x[i]<<' ';of<<std::endl;
	for(int i=0; i<numy; i++)
	{
		for(int j=0; j<neutral_blendshape.size(); j++)
			neutral_blendshape(j) = neutral_blendshape(j) + P_eigenvectors(i,j)*x[starty+i];
	}
	for(int i=0; i<numz; i++)
	{
		for(int j=0; j<neutral_blendshape.size(); j++)
			neutral_blendshape(j) = neutral_blendshape(j) + E_eigenvectors(i,j)*x[startz+i];
	}

}
void TrackingMesh::find_nearest()
{
	nearestIndex.clear();
	
	for(int i=0; i<registerIndex.size(); i++)
	{
		int index = registerIndex[i];
		pcl::PointXYZ searchPoint;
		searchPoint.x = neutral_blendshape(index*3)*x[startR+0]+x[startt+0];
		searchPoint.y = neutral_blendshape(index*3+1)*x[startR+1]+x[startt+1];
		searchPoint.z = neutral_blendshape(index*3+2)*x[startR+2]+x[startt+2];
		
		vector<int> indexs(1,0);
		vector<float> pointNKNDistance (1,0);
		pclCloudKdtree.nearestKSearch (searchPoint, 1, indexs, pointNKNDistance);
		double d;
		nearestIndex.push_back(indexs[0]);
	}
}

void TrackingMesh::trimesh2pcl(const trimesh::TriMesh *mesh,pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud)
{
	pointcloud->width = mesh->vertices.size();
	pointcloud->height = 1;
	pointcloud->is_dense = false;
	pointcloud->resize(pointcloud->width*pointcloud->height);
	for(int i=0;i<pointcloud->points.size();i++)
	{
		pointcloud->points[i].x = mesh->vertices[i][0];
		pointcloud->points[i].y = mesh->vertices[i][1];
		pointcloud->points[i].z = mesh->vertices[i][2];
	}
}
void TrackingMesh::blendshape2trimesh(trimesh::TriMesh * tmesh)
{
	tmesh->clear();
	int j=0;
	for(int i=0; i<neutral_blendshape.size(); i+=3)
	{
		trimesh::point p;
		
			for(int k=0; k<3; k++)
			{
				p[k] = neutral_blendshape(i+k)*x[startR+k]+x[startt+k];
			}
		//}
		tmesh->vertices.push_back(p);
	}
}
void TrackingMesh::loadEigenMatrix(Eigen::MatrixXd& mat, string filename)
{
	ifstream in(filename);
	int m,n;
	in>>m>>n;
	mat.resize(m,n);
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<n; j++)
		{
			double val;
			in>>val;
			mat(i,j) = val;
		}
	}
	in.close();
}
void TrackingMesh::loadEigenVector(Eigen::VectorXd& vec, string filename)
{
	ifstream in(filename);
	int m;
	in>>m;
	vec.resize(m);
	for(int i=0; i<m; i++)
	{
		double val;
		in>>val;
		vec(i) = val;
	}
}
void TrackingMesh::loadEigenMatrix_gl(Eigen::MatrixXd& mat, string filename)
{
	ifstream in(filename);
	int m,n;
	in>>m>>n;
	mat.resize(m*3,n*3);
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<n; j++)
		{
			double val;
			in>>val;
			for(int k=0; k<3; k++)
			{
				for(int g=0; g<3; g++)
				{
					if(k==g)
						mat(i*3+k, j*3+g) = val;
					else mat(i*3+k, j*3+g) = 0.0;
				}
			}
			
		}
	}
	in.close();
}
void TrackingMesh::loadEigenVector_gl(Eigen::VectorXd& vec, string filename)
{
	ifstream in(filename);
	int m;
	in>>m;
	vec.resize(m*3);
	for(int i=0; i<m; i++)
	{
		double val;
		in>>val;
		for(int j=0; j<3; j++)
			vec(i*3+j) = val;
	}
}*/