#include "RigidMesh.h"
RigidMesh::RigidMesh(RigidLMParameters* rigidParas)
{
	of.open("error_RigidMesh.txt",std::ios::app);
	paras = rigidParas;
	//paras->pointCloudMesh = pointCloudMesh;
	pclCloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);
	trimesh2pcl(paras->pointCloudMesh, pclCloud);
	pclCloudKdtree.setInputCloud(pclCloud);
	of<<"before load"<<std::endl;
	readRegisterIndexFromFile(paras->init_registerIndex, "rigid_selection_index.txt");
	of<<"after load"<<std::endl;

	treeCloud = new trimesh::KDtree(paras->pointCloudMesh->vertices);of<<"after load1"<<std::endl;
	paras->pointCloudMesh->need_normals();of<<"after load2"<<std::endl;
	paras->pointCloudMesh->need_bsphere();of<<"after load3"<<std::endl;

	//未知数s(1), r(3), t(3)
	paras->num_s = 1;
	paras->num_r = 3;
	paras->num_t = 3;
	paras->num_all = paras->num_s + paras->num_r + paras->num_t;
	x = new double[paras->num_all];
	paras->start_s = 0;
	paras->start_r = paras->start_s + paras->num_s;
	paras->start_t = paras->start_r + paras->num_r;
	x[paras->start_s] = paras->s;
	for(int i=0; i<3; i++)
	{
		x[paras->start_r+i] = paras->r(i);
		x[paras->start_t+i] = paras->t(i);
	}

	
	
	//paras->expression_blendshape = expression_blendshape;
	of<<"before read_configure "<<std::endl;
	read_configure();
	update_all();
	of<<"init update all finish!"<<std::endl;
	trimesh::TriMesh * cur_mesh = new trimesh::TriMesh;
	blendshape2trimesh(cur_mesh);
	cur_mesh->write("cur_mesh_rigid.obj");
	delete cur_mesh;
	for(int i=0; i<7; i++)of<<x[i]<<std::endl;
	of<<"Construct finish"<<std::endl;
}

RigidMesh::~RigidMesh()
{
	if(treeCloud != NULL)
		delete treeCloud;
	if(x != NULL)
		delete [] x;
}
void RigidMesh::readRegisterIndexFromFile(std::vector<int>& index, std::string filename)
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
void RigidMesh::read_configure()
{
	//读取配置参数
	std::ifstream if_paras("rigid-register-parameters.txt");
	string str;
	if_paras>>str>>paras->maxInIters;
	if_paras>>str>>paras->outIters;
	if_paras>>str>>paras->registerRate;
	if_paras.close();
	of<<"read_configure finish!"<<std::endl;
}

void RigidMesh::update_all()
{
	of<<"find_nearest_start"<<std::endl;
	find_nearest();
	of<<"find nearest finish"<<std::endl;
}

void RigidMesh::find_nearest()
{
	paras->registerIndex.clear();
	paras->nearestIndex.clear();
	paras->nearestWeight.clear();
	paras->init_nearestIndex.clear();
	paras->init_nearestWeight.clear();
	
	Eigen::MatrixXd R = RigidLMCostFunction::getR(&x[paras->start_r]);
	of<<R<<std::endl;
	double max_nearest_distance = 0;
	of<<"nearest_index:"<<paras->expression_blendshape.size()<<' '<<std::endl;
	for(int i=0; i<paras->init_registerIndex.size(); i++)
	{
		int index = paras->init_registerIndex[i];
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			v(j) = paras->expression_blendshape(index*3+j);
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
		paras->init_nearestWeight[i] = 1.0 - paras->init_nearestWeight[i]/(1.0e-20+max_nearest_distance);
		of<<paras->init_nearestWeight[i]<<' ';
	}

	std::vector<double> sort_weight = paras->init_nearestWeight;
	sort(sort_weight.begin(), sort_weight.end());

	double top = sort_weight[sort_weight.size()*(100-paras->registerRate)/100.0];
	
	//只提取最近匹配点的前%
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
	
	of<<std::endl;
	
}
void RigidMesh::fit_mesh()
{
	for(int i=0; i<paras->outIters; i++)
	{
		of<<"iter: "<<i<<std::endl;
		Problem problem;
		// Set up the only cost function (also known as residual).
		CostFunction* cost_function = new RigidLMCostFunction(paras);
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
	}
	//注册后返回
	paras->s = x[paras->start_s];
	for(int i=0; i<3; i++)
	{
		paras->r(i) = x[paras->start_r+i];
		paras->t(i) = x[paras->start_t+i];
	}
	//for(int i=0; i<paras->num_all; i++) paras->x.push_back(x[i]);
}
void RigidMesh::trimesh2pcl(const trimesh::TriMesh *mesh,pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud)
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
void RigidMesh::blendshape2trimesh(trimesh::TriMesh * tmesh)
{
	tmesh->clear();
	int j=0;
	Eigen::MatrixXd R = RigidLMCostFunction::getR(&x[paras->start_r]);
	Eigen::VectorXd T(3);
	for(int i=0; i<3; i++)
	{
		T(i) = x[paras->start_t+i];
	}
	double max_nearest_distance = -1.0e10;
	for(int i=0; i<paras->expression_blendshape.size(); i+=3)
	{
		trimesh::point p;
		
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			v(j) = paras->expression_blendshape(i+j);
		}
		v = R*v*x[paras->start_s]+T;
		for(int k=0; k<3; k++)
		{
			p[k] = v[k];
		}
		//}
		tmesh->vertices.push_back(p);
	}
}
