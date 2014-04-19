#include"InitRegister.h"
//InitRegisterMesh * initRM_global;
/*
InitRegister::InitRegister(trimesh::TriMesh* pointCloud, std::string& meshIndexFilename, std::string& depthIndexFilename)
{
	of.open("error_InitRegister.txt");
	this->pointCloud = pointCloud;
	of<<"before construct"<<std::endl;
	initRM = new InitRegisterMesh(pointCloud, meshIndexFilename, depthIndexFilename, cv::Mat());
	of<<"after construct"<<std::endl;
	//initRM_global = initRM;
}
InitRegister::~InitRegister()
{
	if(initRM != NULL)
	{
		delete initRM;
	}
}*/
/*
void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	initRM_global->of<< iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
	std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
}
void evalfunc(int N, double* x, double *prev_x, double* f, double* g)
{
	//*f = 0;
	//int numf = initRM->numX;
	//Eigen::VectorXd fval;
	//fval.setConstant(numf,0);
	initRM_global->function_f(x, f);
	//*f=fval.norm();
	//*f *= *f;
	initRM_global->function_df(x, g);
}
void InitRegister::Optimize_by_HLBFGS(int N, double *init_x, int num_iter, int M, int T)
{
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	info[4] = num_iter;
	info[6] = T;
	info[7] = 0;
	info[10] = 0;
	info[11] = 1;
	HLBFGS(N, M, init_x, evalfunc, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
}

void InitRegister::Optimize_by_LM(int N, double *init_x, int num_iter, int M, int T)
{

}
*/
/*
void InitRegister::register_mesh(int num_iter, trimesh::TriMesh* tmesh)
{
	
	//double * init_x = initRM->x;
	//for(int i=0; i<num_iter; i++)
	//{
		
		//of<<"iter: "<<i<<std::endl;
		initRM->fit_mesh();
		//initRM->update_all();
	//}
	initRM->blendshape2trimesh(tmesh);
}*/