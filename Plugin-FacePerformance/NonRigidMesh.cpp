#include "NonRigidMesh.h"
NonRigidMesh::NonRigidMesh(NonRigidLMParameters* nonRigidparas)
{
	of.open("error_NonRigidMesh.txt",std::ios::app);
	//of.clear();
	//of<<paras->expression_blendshape.size()<<std::endl;
	//paras->pointCloudMesh = pointCloudMesh;
	paras = nonRigidparas;
	pclCloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);
	trimesh2pcl(paras->pointCloudMesh, pclCloud);
	pclCloudKdtree.setInputCloud(pclCloud);
	of<<"before load"<<std::endl;
	readRegisterIndexFromFile(paras->init_registerIndex, "nonrigid_selection_index.txt");
	of<<"init_registerIndex size: "<<paras->init_registerIndex.size()<<std::endl;
	of<<"after load"<<std::endl;

	treeCloud = new trimesh::KDtree(paras->pointCloudMesh->vertices);of<<"after load1: "<<paras->pointCloudMesh->vertices.size()<<std::endl;
	paras->pointCloudMesh->need_normals();of<<"after load2"<<std::endl;
	paras->pointCloudMesh->need_bsphere();of<<"after load3"<<std::endl;

	//未知数s(1), r(3), t(3), y(50), z(150)
	//paras = nonRigidparas;
	paras->start_x = 0;
	paras->num_x = paras->pre_x.size();
	paras->num_all = paras->num_x;
	x = paras->pre_x;  //初始x
	pre_x = x;
	init();
	
	of<<"before read_configure "<<std::endl;
	read_configure();
	update_all();
	get_3D_to_2D_Matrix(paras->_3d_to_2d_matrix);//3d到2d的投影矩阵
	paras->init_expression_blendshape = paras->expression_blendshape; //存上一frame的表情
	
	of<<"Construct finish"<<std::endl;
}

NonRigidMesh::~NonRigidMesh()
{
	if(treeCloud != NULL)
		delete treeCloud;
}
void NonRigidMesh::readRegisterIndexFromFile(std::vector<int>& index, std::string filename)
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
void NonRigidMesh::init()
{
	paras->d_w.resize(x.size());
	paras->d_w.setZero();
	paras->d_v.resize(x.size());
	paras->d_v.setZero();
	paras->w = x;
	paras->v = x;
	
}

void NonRigidMesh::read_configure()
{
	//读取配置参数
	std::ifstream if_paras("nonrigid-register-parameters.txt");
	string str;
	if_paras>>str>>paras->allIters;
	if_paras>>str>>paras->maxInIters;
	if_paras>>str>>paras->outIters;
	if_paras>>str>>paras->registerRate;
	if_paras>>str>>paras->sigma;
	if_paras>>str>>paras->alpha;
	if_paras>>str>>paras->beta;

	if_paras>>str>>paras->gamma;
	if_paras>>str>>paras->tao;
	if_paras>>str>>paras->dual_residual_sqr_norm_eps;
	if_paras>>str>>paras->primal_residual_sqr_norm_eps;
	if_paras>>str>>paras->mu_w;
	if_paras>>str>>paras->mu_v;
	if_paras>>str>>paras->markerWeight;
	if_paras.close();
	of<<"read_configure finish!"<<std::endl;

	std::ifstream fmesh(paras->meshIndexFilename);
	std::ifstream fdepth(paras->depthIndexFilename);
	paras->marker_registerIndex.clear();
	paras->marker_nearestIndex.clear();
	paras->marker_nearestWeight.clear();
	//std::vector<int> m_DepthIndex;
	//std::vector<int> m_MeshIndex;
	//for(int i=0; i<paras->markerSize; i++)
	{
		int m_d, m_m;

		while(fdepth>>m_d)
		{
			fmesh>>m_m;
			if(m_d != -1)
			{
				//m_MeshIndex.push_back(m_m);
				//m_DepthIndex.push_back(m_d);
				paras->marker_registerIndex.push_back(m_m);
				paras->marker_nearestIndex.push_back(m_d);
				paras->marker_nearestWeight.push_back(paras->markerWeight);/////marker权重为1500
			}
			of<<m_d<<' '<<m_m<<std::endl;
		}
	}
	fmesh.close();fdepth.close();
}

void NonRigidMesh::update_all()
{
	update_expression();
	of<<"find_nearest_start"<<std::endl;
	find_nearest();
	of<<"find nearest finish"<<std::endl;
}
void NonRigidMesh::update_expression()
{
	paras->expression_blendshape = paras->neutral_blendshape;
	for(int i=0; i<paras->num_x; i++)
	{
		paras->expression_blendshape += paras->delta_B.col(i) * x[i];
	}
}
void NonRigidMesh::find_nearest()
{
	paras->registerIndex.clear();
	paras->nearestIndex.clear();
	paras->nearestWeight.clear();
	paras->init_nearestIndex.clear();
	paras->init_nearestWeight.clear();
	paras->I.clear();
	paras->J.clear();
	paras->derive_f.clear();
	paras->derive_J.clear();

	double max_nearest_distance = 0;
	//of<<"nearest_index:"<<std::endl;
	for(int i=0; i<paras->init_registerIndex.size(); i++)
	{
		int index = paras->init_registerIndex[i];
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			//of<<"j: "<<j<<' '<<' '<<(paras->expression_blendshape).size()<<std::endl;
			v(j) = (paras->expression_blendshape)(index*3+j);
		}
		//of<<"i: "<<std::endl;
		v = paras->R*v*paras->s;
		pcl::PointXYZ searchPoint;
		searchPoint.x = v.coeff(0)+paras->T(0);
		searchPoint.y = v.coeff(1)+paras->T(1);
		searchPoint.z = v.coeff(2)+paras->T(2);
		
		vector<int> indexs(1,0);
		vector<float> pointNKNDistance (1,0);
		pclCloudKdtree.nearestKSearch (searchPoint, 1, indexs, pointNKNDistance);
		double d;
		//of<<indexs[0]<<' '<<pointNKNDistance[0]<<' ';
		paras->init_nearestIndex.push_back(indexs[0]);
		paras->init_nearestWeight.push_back(pointNKNDistance[0]); //权重
		if(max_nearest_distance<((double)(pointNKNDistance[0])))
		{
			max_nearest_distance = double(pointNKNDistance[0]);
		}
	}
	//of<<std::endl<<"end nearest_index: "<<max_nearest_distance<<std::endl;
	for(int i=0; i<paras->init_nearestWeight.size(); i++)
	{
		paras->init_nearestWeight[i] = 1.0 - paras->init_nearestWeight[i]/(1.0e-20+max_nearest_distance);
		//of<<paras->init_nearestWeight[i]<<' ';
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

			int register_idx = paras->init_registerIndex[i]*3;
			int nearest_idx = paras->init_nearestIndex[i];
			//paras->I.push_back(NonRigidLMCostFunction::f_color( paras
		}
		//of<<paras->init_nearestWeight[i]<<' ';
	}
	//of<<std::endl;
	/*
	for(int i=0; i<paras->marker_registerIndex.size(); i++)
	{
		paras->registerIndex.push_back(paras->marker_registerIndex[i]);
		paras->nearestIndex.push_back(paras->marker_nearestIndex[i]);
		paras->nearestWeight.push_back(paras->markerWeight);
	}*/
	//of<<std::endl;
}
void NonRigidMesh::update_x()
{
	/*
	for(int i=0; i<paras->outIters; i++)
	{
		of<<"iter: "<<i<<std::endl;
		Problem problem;
		// Set up the only cost function (also known as residual).
		CostFunction* cost_function = new NonRigidLMCostFunction(paras);
		problem.AddResidualBlock(cost_function, NULL, &x(0));
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
		of <<"x :"<<std::endl;
		of<<x.transpose()<<std::endl;
		update_all(); //更新
	}*/
	of <<"before x :"<<std::endl;
	of<<x.transpose()<<std::endl;
	Eigen::MatrixXd A(paras->num_x, paras->num_x);
	Eigen::VectorXd b(paras->num_x);
	A.setZero();
	b.setZero();
	//3d
	of<<"before 3d"<<std::endl;
	for(int i=0; i<paras->registerIndex.size(); i++)//point-plane
	{
		int register_idx = paras->registerIndex[i]*3;
		int nearest_idx = paras->nearestIndex[i];//oof<<"i3: " << i<<std::endl;

		Eigen::MatrixXd cur_A(paras->num_x, paras->num_x);
		Eigen::VectorXd cur_b(paras->num_x);
		Eigen::MatrixXd cur_B(3, paras->num_x);
		Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
		Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
		Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
		Eigen::VectorXd b0(3);
		for(int j=0; j<3; j++)
		{
			//oof<<"i7: " << i<<std::endl;
			v(j) = paras->expression_blendshape(register_idx+j);
			c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
			n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			b0(j) = paras->neutral_blendshape(register_idx+j);
			cur_B.row(j) = paras->delta_B.row(register_idx+j);
		}
		cur_A = (n.transpose() * paras->s * paras->R * cur_B).transpose() * n.transpose() * paras->s * paras->R * cur_B * 2;
		cur_b = (n.transpose() * paras->s * paras->R * cur_B).transpose() * n.transpose() * (paras->s * paras->R * b0 + paras->T - c) * 2;
		A += cur_A * paras->nearestWeight[i];
		b -= cur_b * paras->nearestWeight[i];
	}
	of<<"after 3d"<<std::endl;
	
	//2d
	of<<"before 2d"<<std::endl;
	for(int i=0; i<paras->marker_registerIndex.size(); i++)//point-point
	{
		int register_idx = paras->marker_registerIndex[i]*3;
		int nearest_idx = paras->marker_nearestIndex[i];//oof<<"i3: " << i<<std::endl;

		Eigen::MatrixXd cur_A(paras->num_x, paras->num_x);
		Eigen::VectorXd cur_b(paras->num_x);
		Eigen::MatrixXd cur_B(3, paras->num_x);
		Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
		Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
		Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
		Eigen::VectorXd b0(3);
		for(int j=0; j<3; j++)
		{
			//oof<<"i7: " << i<<std::endl;
			v(j) = paras->expression_blendshape(register_idx+j);
			c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
			n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			b0(j) = paras->neutral_blendshape(register_idx+j);
			cur_B.row(j) = paras->delta_B.row(register_idx+j);
		}
		cur_A = (paras->s * paras->R * cur_B).transpose();
		cur_A *= paras->s * paras->R * cur_B * 2;
		cur_b = (paras->s * paras->R * cur_B).transpose();
		cur_b *= (paras->s * paras->R * b0 + paras->T - c) * 2;
		if(i<paras->marker_registerIndex.size()-1)
		{
			A += cur_A * paras->marker_nearestWeight[i];
			b -= cur_b * paras->marker_nearestWeight[i];
		}
		else
		{
			A += cur_A * paras->marker_nearestWeight[i] * 20;
			b -= cur_b * paras->marker_nearestWeight[i] * 20;
		}
	}
	of<<"after 2d"<<std::endl;
	/*
	Eigen::MatrixXd P3(2,3);
	Eigen::VectorXd P4(2);
	for(int i=0; i<2; i++)
	{
		for(int j=0; j<3; j++)
		{
			P3(i,j) = paras->_3d_to_2d_matrix(j,i);
		}
		P4(i) = paras->_3d_to_2d_matrix(3,i);
	}
	of<<P3<<std::endl;
	of<<P4<<std::endl;
	of<<"before 2d"<<std::endl;
	double value_2d = 0;
	for(int i=0; i<paras->marker_registerIndex.size(); i++)
	{
		int register_idx = paras->marker_registerIndex[i]*3;

		//int nearest_idx = paras->marker_nearestIndex[i];

		Eigen::MatrixXd cur_A(paras->num_x, paras->num_x);
		Eigen::VectorXd cur_b(paras->num_x);
		Eigen::MatrixXd cur_B(3, paras->num_x);
		Eigen::VectorXd v(3);//of<<"i4: " << i<<std::endl;
		//Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
		//Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
		Eigen::VectorXd b0(3);//of<<"i6: " << i<<std::endl;
		for(int j=0; j<3; j++)
		{
			v(j) = paras->expression_blendshape(register_idx+j);//of<<"i7: " << j<<std::endl;
			b0(j) = paras->neutral_blendshape(register_idx+j);//of<<"i7: " << j<<std::endl;
			cur_B.row(j) = paras->delta_B.row(register_idx+j);//of<<"i7: " << j<<std::endl;
		}
		//of<<"i8: "<<std::endl;
		Eigen::VectorXd colorCoord(2);
		colorCoord(0) = P4(0) - paras->marker_Coordinates[i*2];
		//of<<"colorCoord0: "<<colorCoord(0)<<"   ";
		colorCoord(1) = P4(1) - paras->marker_Coordinates[i*2+1];
		//of<<"colorCoord1: "<<colorCoord(1)<<std::endl;

		cur_A = 2 * (P3 * paras->s * paras->R * cur_B).transpose();//of<<"1"<<std::endl;
		cur_A *= (P3 * paras->s * paras->R * cur_B);//of<<"2"<<std::endl;
		cur_b = 2 * (P3 * paras->s * paras->R * cur_B).transpose();//of<<"3"<<std::endl;
		cur_b *= (P3 * (paras->s * paras->R * b0 + paras->T) + colorCoord);//of<<"4"<<std::endl;
		A += cur_A * paras->markerWeight;//of<<"5"<<std::endl;
		b -= cur_b * paras->markerWeight;//of<<"6"<<std::endl;

		Eigen::VectorXd val = P3*(paras->s * paras->R *(b0+cur_B*x)+paras->T)+colorCoord;
		value_2d += val.squaredNorm()*paras->markerWeight;
	}
	of<<"value_2d: "<<value_2d<<std::endl;

	of<<"after 2d"<<std::endl;
	/*of<<"A: "<<std::endl;
	of<<A<<std::endl;
	of<<"b:"<<std::endl;
	of<<b<<std::endl;*/

	//of<<"A: "<<std::endl;
	//of<<A<<std::endl;
	//of<<"b:"<<std::endl;
	//of<<b<<std::endl;
	for(int i=0; i<A.rows(); i++)
	{
		A(i,i) += (paras->mu_w + paras->mu_v + paras->alpha*2);
	}
	b += (paras->w + paras->d_w) * paras->mu_w;
	b += (paras->v + paras->d_v) * paras->mu_v;
	b += (paras->pre_x *2 - paras->ppre_x) * paras->alpha*2;
	Eigen::LDLT<Eigen::MatrixXd> ldlt;
	of<<"before ldlt"<<std::endl;
	ldlt.compute(A);
	
	x = ldlt.solve(b);
	//x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	/*of<<"A: "<<std::endl;
	of<<A<<std::endl;
	of<<"b:"<<std::endl;
	of<<b<<std::endl;*/
	of <<"after x :"<<std::endl;
	of<<x.transpose()<<std::endl;
}
void NonRigidMesh::update_w()
{
	of<<"before w "<<std::endl;
	of<<paras->w.transpose()<<std::endl;
	double norm_w = (x-paras->d_w).norm();
	if(norm_w>paras->beta/paras->mu_w)
	{
		paras->w = (1-paras->beta/paras->mu_w*norm_w)*(x-paras->d_w);
	}
	else
	{
		for(int i=0; i<paras->w.size(); i++) paras->w(i) = 0.0;
	}
	of<<"after w "<<std::endl;
	of<<paras->w.transpose()<<std::endl;
}
void NonRigidMesh::update_v()
{
	of<<"before v "<<std::endl;
	of<<paras->v.transpose()<<std::endl;
	paras->v = x-paras->d_v;
	for(int i=0; i<paras->v.size(); i++)
	{
		if(paras->v(i)>1.0) paras->v(i)=1.0;
		else if(paras->v(i)<0.0) paras->v(i)=0.0;
	}
	of<<"after v "<<std::endl;
	of<<paras->v.transpose()<<std::endl;
}
void NonRigidMesh::update_dual_variable()
{
	paras->d_w += (paras->w-x);
	paras->d_v += (paras->v-x);
}
bool NonRigidMesh::check_convergence_and_update_penalty()
{
	Eigen::VectorXd primal_residual_w = paras->w-x;
	Eigen::VectorXd primal_residual_v = paras->v-x;
	Eigen::VectorXd dual_residual_x = x-pre_x;

	double r_p = primal_residual_w.squaredNorm() + primal_residual_v.squaredNorm();
	double r_d = dual_residual_x.squaredNorm();
	of<<"r_p: "<<r_p<<' '<<paras->primal_residual_sqr_norm_eps<<std::endl;
	of<<"r_d: "<<r_d<<' '<<paras->dual_residual_sqr_norm_eps<<std::endl;
	if(r_p < paras->primal_residual_sqr_norm_eps && r_d < paras->dual_residual_sqr_norm_eps)
	{
		return true;
	}
	//update mu_w
	if( primal_residual_w.norm() > paras->gamma*sqrt(r_d) )
	{
		paras->mu_w *= paras->tao;
	}
	else if(sqrt(r_d)>primal_residual_w.norm()*paras->gamma)
	{
		paras->mu_w /= paras->tao;
	}
	//update mu_v
	if( primal_residual_v.norm() > paras->gamma*sqrt(r_d) )
	{
		paras->mu_v *= paras->tao;
	}
	else if(sqrt(r_d)>primal_residual_v.norm()*paras->gamma)
	{
		paras->mu_v /= paras->tao;
	}
	of<<"mu_w: "<<std::endl;
	of<<paras->mu_w<<std::endl;
	of<<"mu_v: "<<std::endl;
	of<<paras->mu_v<<std::endl;
	return false;
}
void NonRigidMesh::get_3D_to_2D_Matrix(Eigen::MatrixXd & RT)
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

		/*int n = ( p[0]/fx/(900-p[2])+0.5 ) * paras->colorImage.cols;
		int m = ( 0.5-p[1]/fy/(900-p[2]) ) * paras->colorImage.rows;
		if(n<0)n=0;
		else if(n>=paras->colorImage.cols)n = paras->colorImage.cols-1;
		if(m<0)m=0;
		else if(m>=paras->colorImage.rows)m = paras->colorImage.rows-1;
		int k = m*paras->colorImage.cols + n;
		double texture_x = paras->colorCoordinate[2*k];
		double texture_y = paras->colorCoordinate[2*k+1];
		b(i,0) = texture_x;
		b(i,1) = texture_y;*/
		b(i,0) = paras->marker_Coordinates[i*2];
		b(i,1) = paras->marker_Coordinates[i*2+1];
	}
	
	RT = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}
void NonRigidMesh::fit_mesh()
{
	of<<"----------------------------------------------------------"<<std::endl;
	of<<"Frame "<<paras->num_frame + 1<<" Nonrigid Start!"<<std::endl;
	of<<"----------------------------------------------------------"<<std::endl;
	of<<"f_val: "<<f_val()<<std::endl;
	for(int i=0; i<paras->allIters; i++)
	{
		pre_x = x;
		
		of<<"before update_w"<<std::endl;
		update_w();
		of<<"before update_v"<<std::endl;
		update_v();
		of<<"before update_x"<<std::endl;
		update_x();
		of<<"before update_all"<<std::endl;
		update_all();
		of<<"before update_dual"<<std::endl;
		update_dual_variable();
		of<<"iter: "<<i+1<<std::endl;
		of<<"f_val: "<<f_val()<<std::endl;
		if(check_convergence_and_update_penalty())
		{
			of<<"total iters: "<<i+1<<std::endl;
			break;
		}
	}
	for(int i=0; i<paras->num_x; i++)
	{
		if(x(i)<0)x(i)=0;
		else if(x(i)>1)x(i) = 1;
	}
	paras->ppre_x = paras->pre_x;
	paras->pre_x = x;
	of<<"----------------------------------------------------------"<<std::endl;
	of<<"Frame "<<paras->num_frame + 1<<" Nonrigid Finish!"<<std::endl;
	of<<"----------------------------------------------------------"<<std::endl;
	paras->num_frame += 1;
}

void NonRigidMesh::trimesh2pcl(const trimesh::TriMesh *mesh,pcl::PointCloud<pcl::PointXYZ>::Ptr pointcloud)
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
void NonRigidMesh::blendshape2trimesh(trimesh::TriMesh * tmesh)
{
	tmesh->clear();
	//int j=0;
	Eigen::MatrixXd R = paras->R;
	Eigen::VectorXd T = paras->T;
	//double max_nearest_distance = -1.0e10;
	for(int i=0; i<paras->expression_blendshape.size(); i+=3)
	{
		trimesh::point p;
		
		Eigen::VectorXd v(3);
		for(int j=0; j<3; j++)
		{
			v(j) = paras->expression_blendshape(i+j);
		}
		v = R*v*paras->s+T;
		for(int k=0; k<3; k++)
		{
			p[k] = v[k];
		}
		//}
		tmesh->vertices.push_back(p);
	}
}

double NonRigidMesh::f_val()
{
	//double max_nearest_distance = 0;
	double value_3d = 0.0;
	double value_2d = 0.0;
	double value_smooth = 0.0;
	double value_sparse = 0.0;
	//3d
	for(int i=0; i<paras->registerIndex.size(); i++)//point-plane
	{
		int register_idx = paras->registerIndex[i]*3;
		int nearest_idx = paras->nearestIndex[i];//oof<<"i3: " << i<<std::endl;

		Eigen::MatrixXd cur_A(paras->num_x, paras->num_x);
		Eigen::VectorXd cur_b(paras->num_x);
		Eigen::MatrixXd cur_B(3, paras->num_x);
		Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
		Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
		Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
		Eigen::VectorXd b0(3);
		for(int j=0; j<3; j++)
		{
			//oof<<"i7: " << i<<std::endl;
			v(j) = paras->expression_blendshape(register_idx+j);
			c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
			n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			b0(j) = paras->neutral_blendshape(register_idx+j);
			cur_B.row(j) = paras->delta_B.row(register_idx+j);
		}
		Eigen::VectorXd val = n.transpose() * (paras->s * paras->R * (b0 + cur_B * x) + paras->T- c);
		value_3d += paras->nearestWeight[i] * val.squaredNorm();
	}
	for(int i=0; i<paras->marker_registerIndex.size(); i++)//point-point
	{
		int register_idx = paras->marker_registerIndex[i]*3;
		int nearest_idx = paras->marker_nearestIndex[i];//oof<<"i3: " << i<<std::endl;

		Eigen::MatrixXd cur_A(paras->num_x, paras->num_x);
		Eigen::VectorXd cur_b(paras->num_x);
		Eigen::MatrixXd cur_B(3, paras->num_x);
		Eigen::VectorXd v(3);//oof<<"i4: " << i<<std::endl;
		Eigen::VectorXd c(3);//oof<<"i5: " << i<<std::endl;
		Eigen::VectorXd n(3);//oof<<"i6: " << i<<std::endl;
		Eigen::VectorXd b0(3);
		for(int j=0; j<3; j++)
		{
			//oof<<"i7: " << i<<std::endl;
			v(j) = paras->expression_blendshape(register_idx+j);
			c(j) = paras->pointCloudMesh->vertices[nearest_idx][j];
			n(j) = paras->pointCloudMesh->normals[nearest_idx][j];
			b0(j) = paras->neutral_blendshape(register_idx+j);
			cur_B.row(j) = paras->delta_B.row(register_idx+j);
		}
		Eigen::VectorXd val = (paras->s * paras->R * (b0 + cur_B * x) + paras->T- c);
		if(i<paras->marker_registerIndex.size()-1)value_2d += paras->marker_nearestWeight[i] * val.squaredNorm();
		else value_2d += paras->marker_nearestWeight[i] * val.squaredNorm()*20; //下巴点
	}
	//2d
	/*Eigen::MatrixXd P3(2,3);
	Eigen::VectorXd P4(2);
	for(int i=0; i<2; i++)
	{
		for(int j=0; j<3; j++)
		{
			P3(i,j) = paras->_3d_to_2d_matrix(j,i);
		}
		P4(i) = paras->_3d_to_2d_matrix(3,i);
	}
	for(int i=0; i<paras->marker_registerIndex.size(); i++)
	{
		int register_idx = paras->marker_registerIndex[i]*3;
		Eigen::MatrixXd cur_A(paras->num_x, paras->num_x);
		Eigen::VectorXd cur_b(paras->num_x);
		Eigen::MatrixXd cur_B(3, paras->num_x);
		Eigen::VectorXd v(3);//of<<"i4: " << i<<std::endl;
		Eigen::VectorXd b0(3);//of<<"i6: " << i<<std::endl;
		for(int j=0; j<3; j++)
		{
			v(j) = paras->expression_blendshape(register_idx+j);//of<<"i7: " << j<<std::endl;
			b0(j) = paras->neutral_blendshape(register_idx+j);//of<<"i7: " << j<<std::endl;
			cur_B.row(j) = paras->delta_B.row(register_idx+j);//of<<"i7: " << j<<std::endl;
		}
		//of<<"i8: "<<std::endl;
		Eigen::VectorXd colorCoord(2);
		colorCoord(0) = P4(0) - paras->marker_Coordinates[i*2];
		colorCoord(1) = P4(1) - paras->marker_Coordinates[i*2+1];
		Eigen::VectorXd val = P3*(paras->s * paras->R *(b0+cur_B*x)+paras->T)+colorCoord;
		value_2d += val.squaredNorm()*paras->markerWeight;
	}*/
	Eigen::VectorXd xx = x - paras->pre_x * 2 + paras->ppre_x;
	value_smooth = paras->alpha * xx.squaredNorm();

	for(int i=0; i<x.size(); i++)
	{
		value_sparse += paras->beta * abs(x(i));
	}
	double value_fit = value_2d + value_3d + value_smooth + value_sparse;
	of<<"value_2d: "<<value_2d<<std::endl;
	of<<"value_3d: "<<value_3d<<std::endl;
	of<<"value_smooth: "<<value_smooth<<std::endl;
	of<<"value_sparse: "<<value_sparse<<std::endl;
	return value_fit;
}