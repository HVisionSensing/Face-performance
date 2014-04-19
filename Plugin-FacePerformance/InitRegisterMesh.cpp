#include "InitRegisterMesh.h"
InitRegisterMesh::InitRegisterMesh(InitLMParameters* paras)
{
	of.open("error_InitRegisterMesh.txt");
	of<<"before neutral_pca"<<std::endl;
	of<<"after neutral_pca"<<std::endl;
	of<<"after gl_pca"<<std::endl;
	//pointCloudMesh->write("pvM.obj");////////////////////////
	this->paras = paras;
	pclCloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);
	trimesh2pcl(paras->pointCloudMesh, pclCloud);
	pclCloudKdtree.setInputCloud(pclCloud);
	of<<"before load1"<<std::endl;
	readRegisterIndexFromFile(paras->init_registerIndex, "non_rigid_selection_index.txt");of<<"before load2"<<std::endl;
	loadEigenVector(paras->mean_neutral, "neutral_meanface.txt");of<<"before load3"<<std::endl;// = neutral_pca.mean_neutral;
	loadEigenVector(paras->P_eigenvalues, "neutral_eigenvalue.txt");of<<"before load4"<<std::endl;// = neutral_pca.P_eigenvalues;
	loadEigenMatrix(paras->P_eigenvectors, "neutral_eigenvector.txt");of<<"before load5"<<std::endl;// = neutral_pca.P_eigenvectors;
	loadEigenVector_gl(paras->E_eigenvalues, "gl_eigenvalue.txt");of<<"before load6"<<std::endl;// = gl_pca.E_eigenvalues;
	loadEigenMatrix_gl(paras->E_eigenvectors, "gl_eigenvector.txt");of<<"before load7"<<std::endl;// = gl_pca.E_eigenvectors;
	of<<"after load"<<std::endl;

	treeCloud = new trimesh::KDtree(paras->pointCloudMesh->vertices);of<<"after load1"<<std::endl;
	paras->pointCloudMesh->need_normals();of<<"after load2"<<std::endl;
	paras->pointCloudMesh->need_bsphere();of<<"after load3"<<std::endl;

	//未知数s(1), r(3), t(3), y(50), z(150)
	//paras->beta1 = 0.5;
	//paras->beta2 = 0.01;
	//paras->beta3 = 0.0001;
	paras->num_s = 1;
	paras->num_r = 3;
	paras->num_t = 3;
	paras->num_y = paras->P_eigenvalues.size();
	paras->num_z = paras->E_eigenvalues.size();
	paras->num_all = paras->num_s + paras->num_r + paras->num_t + paras->num_y + paras->num_z;
	x = new double[paras->num_all];
	paras->start_s = 0;
	paras->start_r = paras->start_s + paras->num_s;
	paras->start_t = paras->start_r + paras->num_r;
	paras->start_y = paras->start_t + paras->num_t;
	paras->start_z = paras->start_y + paras->num_y;
	
	paras->neutral_blendshape.resize(paras->mean_neutral.size());
	of<<"before initRtyz: "<<paras->num_y<<std::endl;
	init_srtyz(paras->meshIndexFilename, paras->depthIndexFilename);
	update_all();
	get_3D_to_2D_Matrix(paras->_3d_to_2d_matrix);//3d到2d的映射矩阵
	of<<"init update all finish!"<<std::endl;
	trimesh::TriMesh * cur_mesh = new trimesh::TriMesh;
	blendshape2trimesh(cur_mesh);
	cur_mesh->write("cur_mesh.obj");
	delete cur_mesh;
	//for(int i=0; i<6; i++)of<<x[i]<<std::endl;
	of<<"Construct finish"<<std::endl;
}

InitRegisterMesh::~InitRegisterMesh()
{
	if(treeCloud != NULL)
		delete treeCloud;
	if(x != NULL)
		delete [] x;
}
void InitRegisterMesh::readRegisterIndexFromFile(std::vector<int>& index, std::string filename)
{
	index.clear();
	/*
	FILE * f = fopen(filename.c_str(),"r");
	while(feof(f)!=EOF)
	{
		int idx;
		fscanf(f, "%d;", &idx);
		index.push_back(idx);
	}*/
	ifstream in(filename);
	int idx;
	while(in>>idx)
	{
		index.push_back(idx);
	}
	in.close();
}
void InitRegisterMesh::init_srtyz(std::string& meshIndexFilename, std::string& depthIndexFilename)
{
	//读取配置参数
	std::ifstream if_paras("init-register-parameters.txt");
	string str;
	if_paras>>str>>paras->maxInIters;
	if_paras>>str>>paras->outIters;
	if_paras>>str>>paras->markerWeight;
	if_paras>>str>>paras->registerRate;
	if_paras>>str>>paras->beta1;
	if_paras>>str>>paras->beta2;
	if_paras>>str>>paras->beta3;
	if_paras.close();

	std::ifstream fmesh(meshIndexFilename);
	std::ifstream fdepth(depthIndexFilename);
	paras->marker_registerIndex.clear();
	paras->marker_nearestIndex.clear();
	paras->marker_nearestWeight.clear();

	std::vector<int> m_DepthIndex;
	std::vector<int> m_MeshIndex;
	int m_d, m_m;
	while(fdepth>>m_d)
	{
		//fdepth>>m_d;
		fmesh>>m_m;
		if(m_d != -1)
		{
			m_MeshIndex.push_back(m_m);
			m_DepthIndex.push_back(m_d);

			paras->marker_registerIndex.push_back(m_m);
			paras->marker_nearestIndex.push_back(m_d);
			paras->marker_nearestWeight.push_back(paras->markerWeight);/////marker权重为10
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
			of<<paras->mean_neutral(i*3+j)<<' ';
			if(min_val[j] > paras->mean_neutral(i*3+j))
			{
				min_val[j] = paras->mean_neutral(i*3+j);
			}
			if(max_val[j] < paras->mean_neutral(i*3+j))
			{
				max_val[j] = paras->mean_neutral(i*3+j);
			}
			mean_val[j] += paras->mean_neutral(i*3+j);
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
		of<<"pointCloud: "<<paras->pointCloudMesh->vertices[i]<<std::endl;
		for(int j=0; j<3; j++)
		{
			if(min_val_pointcloud[j]>paras->pointCloudMesh->vertices[i][j])
			{
				min_val_pointcloud[j]=paras->pointCloudMesh->vertices[i][j];
			}
			if(max_val_pointcloud[j]<paras->pointCloudMesh->vertices[i][j])
			{
				max_val_pointcloud[j]=paras->pointCloudMesh->vertices[i][j];
			}
			mean_val_pointcloud[j] += paras->pointCloudMesh->vertices[i][j];
		}
	}
	for(int i=0; i<3; i++) mean_val_pointcloud[i] /= m_DepthIndex.size();

	double size[3];
	double mean_size = 0.0;
	for(int i=0; i<3; i++)
	{
		of<<i<<": "<<max_val_pointcloud[i]<<' '<<min_val_pointcloud[i]<<' '<<max_val[i]<<' '<<min_val[i]<<' '
			<<mean_val_pointcloud[i]<<' '<<mean_val[i]<<std::endl;
		size[i] = (max_val_pointcloud[i]-min_val_pointcloud[i])/(max_val[i]-min_val[i]);//R
		mean_size += size[i];
		//if(i==2)x[i] = x[1];////////////////////对角阵
		//x[i+3] = mean_val_pointcloud[i]-mean_val[i]*x[i];//t
	}
	mean_size /= 3;
	for(int i=0; i<paras->num_all; i++)x[i] = 0;//y z
	x[0] = mean_size;
	for(int i=0; i<3; i++)
	{
		//x[paras->start_r+i] = 1.0e-8;
		x[paras->start_t+i] = mean_val_pointcloud[i]-mean_val[i]*x[0];
	}
	of<<"init finish"<<std::endl;

}
/*
void InitRegisterMesh::function_f(double * x, double* f_sum)
{
	Eigen::VectorXd f_val;
	double sum = 0;
	std::vector<trimesh::point>& pointcloud_normals = paras->pointCloudMesh->normals;
	std::vector<trimesh::point>& pointcloud_vertices = paras->pointCloudMesh->vertices;

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
void InitRegisterMesh::function_df(double *x, double *g)
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
*/

void InitRegisterMesh::update_all()
{
	of<<"update_start"<<std::endl;
	update_blendshape();
	of<<"update_blendshape finish"<<std::endl;
	find_nearest();
	of<<"find nearest finish"<<std::endl;
	/*
	for(int i=0; i<registerIndex.size(); i++)
	{
		trimesh::point p;
	}*/
}
void InitRegisterMesh::update_blendshape()
{
	paras->neutral_blendshape = paras->mean_neutral;
	//of<<"neutral_blend: "<<neutral_blendshape.size()<<' '<<neutral_blendshape.cols()<<' '<<neutral_blendshape.rows()<<endl;
	for(int i=0; i<paras->num_all; i++)of<<x[i]<<' ';of<<std::endl;
	for(int i=0; i<paras->num_y; i++)
	{
		for(int j=0; j<paras->neutral_blendshape.size(); j++)
		{
			paras->neutral_blendshape(j) += paras->P_eigenvectors.coeff(i,j)*x[paras->start_y+i];
		}
			/*
		for(int j=0; j<paras->registerIndex.size(); j++)
		{
			int idx = paras->registerIndex[j]*3;
			for(int k=0; k<3; k++)
				paras->neutral_blendshape(idx+k) += paras->P_eigenvectors.coeff(i,idx+k)*x[paras->start_y+i];
		}*/
	}
	for(int i=0; i<paras->num_z; i++)
	{
		for(int j=0; j<paras->neutral_blendshape.size(); j++)
		{
			paras->neutral_blendshape(j) += paras->E_eigenvectors.coeff(i,j)*x[paras->start_z+i];
		}
			/*
		for(int j=0; j<paras->registerIndex.size(); j++)
		{
			int idx = paras->registerIndex[j]*3;
			for(int k=0; k<3; k++)
				paras->neutral_blendshape(idx+k) += paras->E_eigenvectors.coeff(i,idx+k)*x[paras->start_z+i];
		}*/
	}
	/*
	for(int i=0; i<neutral_blendshape.size(); i+=3)
	{
		for(int j=0; j<3; j++)
		{
			neutral_blendshape(i+j) = neutral_blendshape(i+j)*x[startR+j]+x[startt+j]; 
		}
	}*/
}
void InitRegisterMesh::find_nearest()
{
	paras->registerIndex.clear();
	paras->nearestIndex.clear();
	paras->nearestWeight.clear();
	paras->init_nearestIndex.clear();
	paras->init_nearestWeight.clear();
	/*for(int i=0; i<registerIndex.size(); i++)
	{
		int index = registerIndex[i];
		trimesh::point p;
		for(int j=0; j<3; j++)
		{
			p[j] = neutral_blendshape(index*3+j);
		}
		double d;
		treeCloud->closest_to_pt_distance_and_index(p, d, index);
		nearestIndex.push_back(index);
	}*/
	Eigen::MatrixXd R = InitLMCostFunction::getR(&x[paras->start_r]);
	of<<R<<std::endl;
	double max_nearest_distance = 0;
	of<<"nearest_index:"<<std::endl;
	for(int i=0; i<paras->init_registerIndex.size(); i++)
	{
		int index = paras->init_registerIndex[i];
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			v(j) = paras->neutral_blendshape(index*3+j);
		}
		//of<<"i: "<<std::endl;
		v = R*v*x[paras->start_s];
		pcl::PointXYZ searchPoint;
		searchPoint.x = v.coeff(0)+x[paras->start_t+0];
		searchPoint.y = v.coeff(1)+x[paras->start_t+1];
		searchPoint.z = v.coeff(2)+x[paras->start_t+2];
		
		vector<int> indexs(1,0);
		vector<float> pointNKNDistance (1,0);
		pclCloudKdtree.nearestKSearch (searchPoint, 1, indexs, pointNKNDistance);
		double d;
		of<<indexs[0]<<' '<<pointNKNDistance[0]<<' ';
		paras->init_nearestIndex.push_back(indexs[0]);
		paras->init_nearestWeight.push_back(pointNKNDistance[0]); //权重
		if(max_nearest_distance<((double)(pointNKNDistance[0])))
		{
			max_nearest_distance = double(pointNKNDistance[0]);
		}
	}
	of<<std::endl<<"end nearest_index: "<<max_nearest_distance<<std::endl;
	for(int i=0; i<paras->init_nearestWeight.size(); i++)
	{
		paras->init_nearestWeight[i] = 1.0 - paras->init_nearestWeight[i]/max_nearest_distance;
		of<<paras->init_nearestWeight[i]<<' ';
	}

	std::vector<double> sort_weight = paras->init_nearestWeight;
	sort(sort_weight.begin(), sort_weight.end());

	double top = sort_weight[sort_weight.size()*(100-paras->registerRate)/100.0];
	
	//只提取最近匹配点的前30%
	for(int i=0; i<paras->init_nearestWeight.size(); i++)
	{
		//if(paras->nearestWeight[i] <top10)paras->nearestWeight[i] = 0.0;//1.0 - paras->nearestWeight[i]/max_nearest_distance;
		if(paras->init_nearestWeight[i] >= top)
		{
			paras->registerIndex.push_back(paras->init_registerIndex[i]);
			paras->nearestIndex.push_back(paras->init_nearestIndex[i]);
			paras->nearestWeight.push_back(paras->init_nearestWeight[i]);
		}
		of<<paras->init_nearestWeight[i]<<' ';
	}
	paras->markerWeight = double(paras->registerIndex.size());
	for(int i=0; i<paras->marker_registerIndex.size(); i++)
	{
		paras->registerIndex.push_back(paras->marker_registerIndex[i]);
		paras->nearestIndex.push_back(paras->marker_nearestIndex[i]);
		paras->nearestWeight.push_back(paras->marker_nearestWeight[i]);
	}
	of<<std::endl;
	
}
void InitRegisterMesh::fit_mesh()
{
	for(int i=0; i<paras->outIters; i++)
	{
		of<<"iter: "<<i<<std::endl;
		Problem problem;
		// Set up the only cost function (also known as residual).
		CostFunction* cost_function = new InitLMCostFunction(paras);
		problem.AddResidualBlock(cost_function, NULL, x);
		of<<"after add cost_fnuction"<<std::endl;
		Solver::Options options;
		options.max_num_iterations = paras->maxInIters;
		options.linear_solver_type = ceres::DENSE_QR;
		//options.minimizer_progress_to_stdout = true;
		Solver::Summary summary;
		of<<"before solve"<<std::endl;
		Solve(options, &problem, &summary);
		of<<"after solve"<<std::endl;
		of << summary.BriefReport() << "\n";
		update_all(); //更新
		
		//of << "x : " << initial_x << " -> " << x << "\n";
	}
}
void InitRegisterMesh::trimesh2pcl(const trimesh::TriMesh *mesh,pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud)
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
void InitRegisterMesh::getinitParas(Eigen::VectorXd& vec, double& s, Eigen::VectorXd& r, Eigen::VectorXd& t, Eigen::MatrixXd& _3d_to_2d_matrix)
{
	r.resize(3);
	t.resize(3);
	vec = paras->neutral_blendshape;
	s = x[paras->start_s];
	for(int i=0; i<3; i++)
	{
		r(i) = (x[paras->start_r+i]);
		t(i) = (x[paras->start_t+i]);
	}
	_3d_to_2d_matrix = paras->_3d_to_2d_matrix;
}
void InitRegisterMesh::get_3D_to_2D_Matrix(Eigen::MatrixXd & RT)
{
	double fx = 1.2017214;
	double fy = 0.9030345;
	double maxDepth = paras->maxDepth;
	//Eigen::MatrixXd RT(4,2);
	Eigen::MatrixXd A(paras->marker_nearestIndex.size(), 4);
	Eigen::MatrixXd b(paras->marker_nearestIndex.size(), 2);
	for(int i=0; i<paras->marker_nearestIndex.size(); i++)
	{
		int nearest_idx = paras->marker_nearestIndex[i];
		trimesh::point p = paras->pointCloudMesh->vertices[nearest_idx];
		for(int j=0; j<3; j++) A(i,j) = p[j];
		A(i,3) = 1;
/*
		int n = ( p[0]/fx/(900-p[2])+0.5 ) * paras->texture.cols;
		int m = ( 0.5-p[1]/fy/(900-p[2]) ) * paras->texture.rows;
		if(n<0)n=0;
		else if(n>=paras->texture.cols)n = paras->texture.cols-1;
		if(m<0)m=0;
		else if(m>=paras->texture.rows)m = paras->texture.rows-1;
		int k = m*paras->texture.cols + n;
		int texture_x = paras->colorCoordinate[2*k];
		int texture_y = paras->colorCoordinate[2*k+1];
		b(i,0) = texture_x;
		b(i,1) = texture_y;*/
		b(i,0) = paras->marker_Coordinates[i*2];
		b(i,1) = paras->marker_Coordinates[i*2+1];
	}
	
	RT = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}

void InitRegisterMesh::blendshape2trimesh(trimesh::TriMesh * tmesh)
{
	static std::string textureName = "1";
	static double max_n;
	std::vector<double> nearestWeight;
	ofstream of_tex(textureName+"vt");
	tmesh->clear();
	int j=0;
	Eigen::MatrixXd R=InitLMCostFunction::getR(&x[paras->start_r]);
	Eigen::VectorXd T(3);
	for(int i=0; i<3; i++)
	{
		T(i) = x[paras->start_t+i];
	}

	double max_nearest_distance = -1.0e10;
	for(int i=0; i<paras->neutral_blendshape.size(); i+=3)
	{
		trimesh::point p;
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			v(j) = paras->neutral_blendshape(i+j);
		}
		v = R*v*x[paras->start_s]+T;
			for(int k=0; k<3; k++)
			{
				p[k] = v(k);//paras->neutral_blendshape(i+k)*x[startR+k]+x[startt+k];
			}
		//}
		tmesh->vertices.push_back(p);
		Eigen::Vector2d vt = InitLMCostFunction::f_vt(v(0), v(1), v(2), paras->texture, paras->colorCoordinate, paras->_3d_to_2d_matrix);
		of_tex<<"vt "<<vt(0)<<' '<<vt(1)<<std::endl;
		/*
		bool mark = false;
		for(int j=0; j<paras->init_registerIndex.size(); j++)
		{
			if(paras->init_registerIndex[j]==i/3)mark = true;
			else if(paras->init_registerIndex[j]>i/3)break;
		}
		if(mark)
		{
			
		}
		else of_tex<<"vt 0 0"<<std::endl;*/
	}
	/*
	if(textureName == "1")max_n = max_nearest_distance;
	//tmesh->need_normals();
	of_tex<<"dis: "<<max_nearest_distance<<std::endl;
	for(int i=0; i<paras->neutral_blendshape.size(); i+=3)
	{
		if(nearestWeight[i/3] == -1 || nearestWeight[i/3]>=0.07)of_tex<<"vt"<<'1'<<' '<<'1'<<std::endl;
		else of_tex<<"vt"<<' '<<nearestWeight[i/3]/max_n<<' '<<nearestWeight[i/3]/max_n<<std::endl;
	}*/
	of_tex.close();
	textureName[0] += 1;
}
/*
void InitRegisterMesh::blendshape2trimesh(trimesh::TriMesh * tmesh)
{
	static std::string textureName = "1";
	static double max_n;
	std::vector<double> nearestWeight;
	ofstream of_tex(textureName+"vt");
	tmesh->clear();
	int j=0;
	Eigen::MatrixXd R=InitLMCostFunction::getR(&x[paras->start_r]);
	Eigen::VectorXd T(3);
	for(int i=0; i<3; i++)
	{
		T(i) = x[paras->start_t+i];
	}

	double max_nearest_distance = -1.0e10;
	for(int i=0; i<paras->neutral_blendshape.size(); i+=3)
	{
		trimesh::point p;
		
		
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			v(j) = paras->neutral_blendshape(i+j);
		}
		v = R*v*x[paras->start_s]+T;
			for(int k=0; k<3; k++)
			{
				p[k] = v[k];//paras->neutral_blendshape(i+k)*x[startR+k]+x[startt+k];
			}
		//}
		tmesh->vertices.push_back(p);
		bool mark = false;
		for(int j=0; j<paras->init_registerIndex.size(); j++)
		{
			if(paras->init_registerIndex[j]==i/3)mark = true;
			else if(paras->init_registerIndex[j]>i/3)break;
		}
		if(mark)
		{
			pcl::PointXYZ searchPoint;
			searchPoint.x = p[0];
			searchPoint.y = p[1];
			searchPoint.z = p[2];
			vector<int> indexs(1,0);
			vector<float> pointNKNDistance (1,0);
			pclCloudKdtree.nearestKSearch (searchPoint, 1, indexs, pointNKNDistance);
			double d;
			of<<indexs[0]<<' '<<pointNKNDistance[0]<<' ';
			//paras->init_nearestIndex.push_back(indexs[0]);
			nearestWeight.push_back(pointNKNDistance[0]); //权重
			if(max_nearest_distance<((double)(pointNKNDistance[0])))
			{
				max_nearest_distance = double(pointNKNDistance[0]);
			}
		}
		else nearestWeight.push_back(-1.0);
	}
	if(textureName == "1")max_n = max_nearest_distance;
	//tmesh->need_normals();
	of_tex<<"dis: "<<max_nearest_distance<<std::endl;
	for(int i=0; i<paras->neutral_blendshape.size(); i+=3)
	{
		if(nearestWeight[i/3] == -1 || nearestWeight[i/3]>=0.07)of_tex<<"vt"<<'1'<<' '<<'1'<<std::endl;
		else of_tex<<"vt"<<' '<<nearestWeight[i/3]/max_n<<' '<<nearestWeight[i/3]/max_n<<std::endl;
	}
	of_tex.close();
	textureName[0] += 1;
}
*/
void InitRegisterMesh::loadEigenMatrix(Eigen::MatrixXd& mat, string filename)
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
void InitRegisterMesh::loadEigenVector(Eigen::VectorXd& vec, string filename)
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
void InitRegisterMesh::loadEigenMatrix_gl(Eigen::MatrixXd& mat, string filename)
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
void InitRegisterMesh::loadEigenVector_gl(Eigen::VectorXd& vec, string filename)
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
}