#include "RegisterMesh.h"
/*bool isbigger(float &a,float &b)
{
	return (a<b);
}

RegisterMesh::RegisterMesh(trimesh::TriMesh* pointCloudMesh)
{
	this->pointCloudMesh = pointCloudMesh;
	templateMesh = trimesh::TriMesh::read("template.obj");
	//OpenMesh::IO::read_mesh(*templateMesh, "template.obj");

	//对模板mesh和点云mesh做处理
	pointCloudMesh->need_bsphere();
	trimesh::point& center = pointCloudMesh->bsphere.center;
	float  r = pointCloudMesh->bsphere.r*0.5;
	std::vector<trimesh::point>& ps = pointCloudMesh->vertices;
	for (size_t i=0; i != pointCloudMesh->vertices.size(); ++i)
	{
		ps[i] -= center;
		ps[i] /= r;
	}
	meshes.push_back(pointCloudMesh);////////////
	templateMesh->need_faces();
	templateMesh->need_normals();
	templateMesh->need_neighbors();
	templateMesh->need_bsphere();
	std::vector<trimesh::point>& pts = templateMesh->vertices;
	for (size_t i=0; i != templateMesh->vertices.size(); ++i)
	{
		pts[i] -= center;
		pts[i] /= r;
	}
	meshes.push_back(templateMesh);
	treeCloud = new trimesh::KDtree(pointCloudMesh->vertices);
	treeTemplate = new trimesh::KDtree(templateMesh->vertices);
	//treeCloud->closest_to_pt_distance_and_index(
}

RegisterMesh::~RegisterMesh()
{
	if(templateMesh != NULL)
		delete templateMesh;
	if(treeCloud != NULL)
		delete treeCloud;
	if(treeTemplate != NULL)
		delete treeTemplate;
}

void RegisterMesh::Transformation()
{
	FILE *f1, *f2, *f3;
	f1 = fopen("Feature.txt", "r");
	f2 = fopen("TemplateFeatureRe.txt", "r");
	f3 = fopen("1.txt", "r");

	m_kinectcloudindex.clear();
	m_templatecloudindex.clear();
	m_templateAdditionFeature.clear();

	std::vector<int> fFinalSelectPts;
	fFinalSelectPts.clear();
	int index1, index2;
	int iter = 0;
	float b, c;
	int time=0;
	while(!feof(f1))
	{
		fscanf(f1, "%d", &index1);
		fscanf(f2, "%d", &index2);
		fscanf(f3, "%f%f", &b, &c);
		if (index1 && index2)
		{
			fFinalSelectPts.push_back(index1);
			m_kinectcloudindex.push_back(index1);
			m_templatecloudindex.push_back(index2);
			trimesh::ivec2 p;
			p[0] = (int)c;
			p[1] = (int)b;
			m_imageArray.push_back(p);
		}
		else
		{
			if (fFinalSelectPts.size()<53)
			{
				fFinalSelectPts.push_back(0);
			}

			m_templateAdditionFeature.push_back(index2);
		}
		time++;
	}
	//cout<<time<<endl;
	//cout<<"fFinalSelectPts   "<<fFinalSelectPts.size()<<endl;


	//for (int i=0; i<fFinalSelectPts.size(); i++)
	{
		//cout<<fFinalSelectPts[i]<<endl;
	}
	fclose(f1);
	fclose(f2);
	fclose(f3);

	std::vector<trimesh::point> m_KscalePts;
	std::vector<trimesh::point> m_TscalePts;
	for (int i = 0; i != m_kinectcloudindex.size(); i++)
	{
		m_KscalePts.push_back(meshes[0]->vertices[m_kinectcloudindex[i]]);
		m_TscalePts.push_back(meshes[1]->vertices[m_templatecloudindex[i]]);
	}
	//iterative 
	float R[3][3], T[3];
	MySVD( m_KscalePts, m_TscalePts, R, T);
	//cout<<"here1"<<endl;
	//cout<<m_KscalePts.size()<<endl;
	std::vector<float> fDisArray;
	float fDisArrayAverage = 0;
	fDisArray.clear();
	for (int i=0; i<m_KscalePts.size(); i++)
	{
		trimesh::point p;
		p[0] = m_KscalePts[i].dot(trimesh::point(R[0][0], R[0][1], R[0][2])) + T[0];
		p[1] = m_KscalePts[i].dot(trimesh::point(R[1][0], R[1][1], R[1][2])) + T[1];
		p[2] = m_KscalePts[i].dot(trimesh::point(R[2][0], R[2][1], R[2][2])) + T[2];
		float temp = trimesh::len(p - m_TscalePts[i]);
		fDisArray.push_back(temp);
		fDisArrayAverage = fDisArrayAverage + temp;
	}

	fDisArrayAverage = fDisArrayAverage *1.0/m_KscalePts.size();
	float fVar = 0;
	for (int i =0; i<m_KscalePts.size(); i++)
	{
		fDisArray[i] -= fDisArrayAverage;
		fVar +=  (fDisArray[i] ) * (fDisArray[i] );
	}
	fVar = sqrt(fVar/(m_KscalePts.size() - 1));

	std::vector<int> nSelectPtIdxArray;
	nSelectPtIdxArray.clear();

	m_cloudDeformFeature.clear();
	m_templateDeformFeature.clear();
	for (int i=0; i<m_KscalePts.size(); i++)
	{
		if (fabs(fDisArray[i]) < 1*fVar )
		{
			nSelectPtIdxArray.push_back(i);
			m_templateDeformFeature.push_back(m_templatecloudindex[i]);
			m_cloudDeformFeature.push_back(m_kinectcloudindex[i]);
		}
		else
		{
			m_templateAdditionFeature.push_back(m_templatecloudindex[i]);
		}
	}

	m_KscalePts.clear();
	m_TscalePts.clear();
	for (int i = 0; i != m_kinectcloudindex.size(); i++)
	{
		m_KscalePts.push_back(meshes[0]->vertices[m_kinectcloudindex[i]]);
		m_TscalePts.push_back(meshes[1]->vertices[m_templatecloudindex[i]]);
	}
	std::vector<trimesh::point> m_KPtsSecond;
	std::vector<trimesh::point> m_TPtsSecond;
	for (int i=0; i<nSelectPtIdxArray.size(); i++)
	{
		m_KPtsSecond.push_back(m_KscalePts[i]);
		m_TPtsSecond.push_back(m_TscalePts[i]);
	}
	//cout<<"here2"<<endl;
	//cout<<nSelectPtIdxArray.size()<<endl;
	MySVD( m_KPtsSecond, m_TPtsSecond, R, T);

	for (int i = 0; i != meshes[0]->vertices.size(); i++)
	{
		trimesh::point p;
		p[0] = meshes[0]->vertices[i].dot(trimesh::point(R[0][0], R[0][1], R[0][2])) + T[0];
		p[1] = meshes[0]->vertices[i].dot(trimesh::point(R[1][0], R[1][1], R[1][2])) + T[1];
		p[2] = meshes[0]->vertices[i].dot(trimesh::point(R[2][0], R[2][1], R[2][2])) + T[2];
		meshes[0]->vertices[i] = p;
	}
}
void RegisterMesh::MySVD(std::vector<trimesh::point>& KPts, std::vector<trimesh::point>& TPts, float R[3][3], float T[3])
{
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			R[i][j] = 0;
		}
	}

	for (int i=0; i<3; i++)
	{
		T[i] = 0;
	}

	trimesh::point vKmean(0, 0, 0);
	trimesh::point vTmean(0, 0, 0);
	std::vector<trimesh::point> m_KscalePts;
	std::vector<trimesh::point> m_TscalePts;
	for (int i=0; i<KPts.size(); i++)
	{
		m_KscalePts.push_back(KPts[i]);
		m_TscalePts.push_back(TPts[i]);
	}
	for (int i = 0; i != m_kinectcloudindex.size(); i++)
	{
		m_KscalePts.push_back(meshes[0]->vertices[m_kinectcloudindex[i]]);
		m_TscalePts.push_back(meshes[1]->vertices[m_templatecloudindex[i]]);
	}
	for (int i=0; i<m_KscalePts.size(); i++)
	{
		vKmean = vKmean + m_KscalePts[i];
		vTmean = vTmean + m_TscalePts[i];
	}
	vKmean[0] = vKmean[0] * (1.0/m_KscalePts.size());
	vKmean[1] = vKmean[1] * (1.0/m_KscalePts.size());
	vKmean[2] = vKmean[2] * (1.0/m_KscalePts.size());
	vTmean[0] = vTmean[0] * (1.0/m_TscalePts.size());
	vTmean[1] = vTmean[1] * (1.0/m_TscalePts.size());
	vTmean[2] = vTmean[2] * (1.0/m_TscalePts.size());

	std::vector<float>	scaleArray;
	for (int i = 1; i < m_KscalePts.size(); i = i + 2)
	{
		float fKDis = sqrt((m_KscalePts[i] - m_KscalePts[i-1]).dot(m_KscalePts[i] - m_KscalePts[i-1]));
		float fTDis = sqrt((m_TscalePts[i] - m_TscalePts[i-1]).dot(m_TscalePts[i] - m_TscalePts[i-1]));
		scaleArray.push_back(fTDis*1.0/fKDis);
	}
	std::sort(scaleArray.begin(),scaleArray.end(),isbigger);
	float averageScale = 0;
	for (int i=0; i<scaleArray.size(); i++)
	{
		averageScale = averageScale + scaleArray[i];
	}
	averageScale = averageScale * 1.0 / scaleArray.size();
	averageScale = scaleArray[int(scaleArray.size()/2)];
	for (int i = 0; i != meshes[0]->vertices.size(); i++)
	{
		trimesh::point q;
		q[0] = (meshes[0]->vertices[i][0] - vKmean[0]) * averageScale + vTmean[0];
		q[1] = (meshes[0]->vertices[i][1] - vKmean[1]) * averageScale + vTmean[1];
		q[2] = (meshes[0]->vertices[i][2] - vKmean[2]) * averageScale + vTmean[2];
		meshes[0]->vertices[i] = q;
	}

	for (int i=0; i<3; i++)
	{
		vKmean[i] = 0;
		vTmean[i] = 0;
	}
	m_KscalePts.clear();
	m_TscalePts.clear();
	for (int i=0; i<KPts.size(); i++)
	{
		m_KscalePts.push_back(KPts[i]);
		m_TscalePts.push_back(TPts[i]);
	}
	for (int i=0; i<m_KscalePts.size(); i++)
	{
		vKmean = vKmean + m_KscalePts[i];
		vTmean = vTmean + m_TscalePts[i];
	}
	vKmean[0] = vKmean[0] * (1.0/m_KscalePts.size());
	vKmean[1] = vKmean[1] * (1.0/m_KscalePts.size());
	vKmean[2] = vKmean[2] * (1.0/m_KscalePts.size());
	vTmean[0] = vTmean[0] * (1.0/m_TscalePts.size());
	vTmean[1] = vTmean[1] * (1.0/m_TscalePts.size());
	vTmean[2] = vTmean[2] * (1.0/m_TscalePts.size());
	for (int i = 0; i != m_KscalePts.size(); i++)
	{
		m_KscalePts[i] = m_KscalePts[i] - vKmean;
		m_TscalePts[i] = m_TscalePts[i] - vTmean;
	}
	double H[3][3];
	for(int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			H[i][j] = 0;
		}
	}
	for (int i=0; i<m_KscalePts.size(); i++)
	{
		float tempH[3][3];
		float vKinect[3][1] =	{ m_KscalePts[i][0], m_KscalePts[i][1], m_KscalePts[i][2] };
		float vTemp[1][3]	=	{ m_TscalePts[i][0], m_TscalePts[i][1], m_TscalePts[i][2] };
		mklMulti(vKinect, 3, 1, CblasNoTrans, vTemp, 3, CblasNoTrans, tempH, CblasRowMajor);
		for(int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				H[i][j] = H[i][j] + tempH[i][j];
			}
		}
	}
	float u[3][3], w[3], vt[3][3];
	mklSvdcmp(H, 3, 3, u, w, vt);
	mklMulti(vt, 3, 3, CblasTrans, u, 3, CblasTrans, R, CblasRowMajor);
	float kinectMean[3] =	{vKmean[0], vKmean[1], vKmean[2]};
	float tempMean[3]	=	{vTmean[0], vTmean[1], vTmean[2]};
	for (int i=0; i<3; i++)
	{
		T[i] = tempMean[i] - (R[i][0] * kinectMean[0] + R[i][1] * kinectMean[1] + R[i][2] * kinectMean[2]);
	}
	m_KscalePts.clear();
	m_TscalePts.clear();
	scaleArray.clear();
}


void RegisterMesh::Generate()
{
	meshes[0]->vertices = m_deformCloud;
	//cout<<meshes[0]->vertices.size()<<"　"<<meshes[1]->vertices.size()<<endl;
	InitProject();
	CalProjectMat();
	CalColorIndex();
}

void RegisterMesh::InitProject()
{
	m_pointArray.clear();
	m_imageArray.clear();
	m_imageArrayy.clear();
	m_templatecloudindex.clear();
	FILE *f1, *f2;
	f1 = fopen("templatefeature.txt", "r");
	f2 = fopen("1.txt", "r");
	int a;
	float b, c;
	while (!feof(f2))
	{
		fscanf(f1, "%d", &a);
		fscanf(f2, "%f%f", &b, &c);
		m_templatecloudindex.push_back(a);
		trimesh::ivec2 p;
		p[0] = (int)(c);
		p[1] = (int)(b);
		m_imageArray.push_back(p);
	}
	m_templatecloudindex.push_back(35446);
	m_templatecloudindex.push_back(20079);
// 	m_templatecloudindex.push_back(12278);
// 	m_templatecloudindex.push_back(4409);
	for (int i = 0; i != pictureextrafeature.size(); i++)
	{
		trimesh::ivec2 p;
		p[0] = (int)(pictureextrafeature[i].y);
		p[1] = (int)(pictureextrafeature[i].x);
		m_imageArray.push_back(p);
	}	
	fclose(f1);
	fclose(f2);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[54]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[55]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[12]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[14]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[17]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[19]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[25]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[27]]);
// 	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[56]]);
//     m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[57]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[33]]);
	m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[39]]);
	//m_pointArray.push_back(meshes[0]->vertices[m_templatecloudindex[55]]);
	m_imageArrayy.push_back(m_imageArray[54]);
	m_imageArrayy.push_back(m_imageArray[55]);
	m_imageArrayy.push_back(m_imageArray[12]);
	m_imageArrayy.push_back(m_imageArray[14]);
	m_imageArrayy.push_back(m_imageArray[17]);
	m_imageArrayy.push_back(m_imageArray[19]);
	m_imageArrayy.push_back(m_imageArray[25]);
	m_imageArrayy.push_back(m_imageArray[27]);
// 	m_imageArrayy.push_back(m_imageArray[56]);
// 	m_imageArrayy.push_back(m_imageArray[57]);
	m_imageArrayy.push_back(m_imageArray[33]);
	m_imageArrayy.push_back(m_imageArray[39]);
	//m_imageArrayy.push_back(m_imageArray[55]);
	std::vector<int>::iterator it1, it2;
	std::vector<trimesh::ivec2>::iterator it3, it4;
	it1 = m_templatecloudindex.begin();
	it3 = m_imageArray.begin();
// 	it2 = it1 + 57;
// 	m_templatecloudindex.erase(it2);
// 	it2 = it1 + 56;
// 	m_templatecloudindex.erase(it2);
	it2 = it1 + 55;
	m_templatecloudindex.erase(it2);
	it2 = it1 + 54;
	m_templatecloudindex.erase(it2);
// 	it4 = it3 + 57;
// 	m_imageArray.erase(it4);
// 	it4 = it3 + 56;
// 	m_imageArray.erase(it4);
	it4 = it3 + 55;
	m_imageArray.erase(it4);
	it4 = it3 + 54;
	m_imageArray.erase(it4);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			m_fProjectMat[i][j] = 0;
		}
	}
}

void RegisterMesh::CalProjectMat()
{
	int nrow = m_pointArray.size();
	int ncol = 4;
	RowMat<double> rmat(nrow, ncol);
	for (int i = 0; i < nrow; i++)
	{
		rmat[i][0] = m_pointArray[i][0];
		rmat[i][1] = m_pointArray[i][1];
		rmat[i][2] = m_pointArray[i][2];
		rmat[i][3] = 1;
	}

	std::vector<double>	b(nrow, 0);
	std::vector<double>	x(ncol, 0);

	CSRMatrix<double> M(rmat);

	for (int i = 0; i < nrow; i++)
	{
		//b[i] = m_imageArray[i][0];
		b[i] = m_imageArrayy[i][0];
	}
	LeastSquareSolve(M, &b[0], &x[0]);

	for (int i = 0; i < ncol; i++)
	{
		m_fProjectMat[0][i] = x[i];
	}

	for (int i = 0; i < nrow; i++)
	{
		//b[i] = m_imageArray[i][1];
		b[i] = m_imageArrayy[i][1];
	}

	LeastSquareSolve(M, &b[0], &x[0]);

	for (int i = 0; i < ncol; i++)
	{
		m_fProjectMat[1][i] = x[i];
	}

	for (int i = 0; i < nrow; i++)
	{
		b[i] = 1;
	}
	LeastSquareSolve(M, &b[0], &x[0]);
	for (int i = 0; i < ncol; i++)
	{
		m_fProjectMat[2][i] = x[i];
	}
}

void RegisterMesh::CalColorIndex()
{
	IplImage *Img = cvLoadImage("texture.JPG", 1);
	float p[4];
	float u, v;
	trimesh::point m_color_index;
	meshes[0]->colors.clear();
	for (size_t i = 0; i != meshes[0]->vertices.size(); i++)
	{
		p[0] = meshes[0]->vertices[i][0];
		p[1] = meshes[0]->vertices[i][1];
		p[2] = meshes[0]->vertices[i][2];
		p[3] = 1;
		u = m_fProjectMat[0][0] * p[0] + m_fProjectMat[0][1] * p[1] + m_fProjectMat[0][2] * p[2] + m_fProjectMat[0][3] * p[3];
		v = m_fProjectMat[1][0] * p[0] + m_fProjectMat[1][1] * p[1] + m_fProjectMat[1][2] * p[2] + m_fProjectMat[1][3] * p[3];
		u = (int)(u + 0.5);
		v = (int)(v + 0.5);
		if (u > 0 && u < Img->height && v < Img->width && v > 0)
		{
			m_color_index[0] = cvGet2D(Img, (int)(u), (int)(v)).val[2];
			m_color_index[1] = cvGet2D(Img, (int)(u), (int)(v)).val[1];
			m_color_index[2] = cvGet2D(Img, (int)(u), (int)(v)).val[0];
			meshes[0]->colors.push_back(m_color_index);
		}
		else
		{
			m_color_index[0] = 0;
			m_color_index[1] = 0;
			m_color_index[2] = 0;
			meshes[0]->colors.push_back(m_color_index);
		}
	}
}

void RegisterMesh::ICP(const xform& xf1)
{
	int verbose = 0;
	bool do_scale	= false;
	bool do_affine	= false;
	bool bulkmode	= false;
	trimesh::TriMesh::set_verbose(verbose);

	xform  xf2;
	std::vector<float> weights1, weights2;

	trimesh::TriMesh * tmpMesh = new trimesh::TriMesh;
	trimesh::TriMesh * pmesh = meshes[1];
	const std::vector<trimesh::point>& vs = pmesh->vertices;
	const std::vector<trimesh::point>& ns = pmesh->normals;

	std::vector<trimesh::point>&  tns = tmpMesh->normals;

	int nsample = 100;
	tmpMesh->vertices.resize(vs.size()/nsample);
	tns.resize(vs.size()/nsample);

	for (size_t i=0; i<tmpMesh->vertices.size(); ++i)
	{
		tmpMesh->vertices[i] = vs[i*100];
		tns[i] = ns[i*100];
	}

	trimesh::KDtree * tmpTree = new trimesh::KDtree(tmpMesh->vertices);
	
	float err = trimesh::ICP(tmpMesh, meshes[0], xf1, xf2, tmpTree, m_treeCloud,
		weights1,weights2,verbose, do_scale, do_affine);
	if (err >= 0.0f)
	{
		err = trimesh::ICP(tmpMesh, meshes[0], xf1, xf2,  tmpTree, m_treeCloud,
			weights1,weights2,verbose, do_scale, do_affine);
	}
	if (err < 0.0f)
	{
		trimesh::TriMesh::eprintf("ICP failed\n");
		exit(1);
	}

	xform xf12 = inv(xf2) * xf1;
	WrapMesh(meshes[1], xf12);

	delete tmpMesh;
	delete tmpTree;
}
void RegisterMesh::WrapMesh(trimesh::TriMesh * pmesh, const xform& xf)
{
	std::vector<trimesh::point>& vs = pmesh->vertices;
	for (size_t i=0; i != vs.size(); ++i)
	{
		double oldx(vs[i][0]), oldy(vs[i][1]), oldz(vs[i][2]);
		vs[i][0] = xf(0,0)*oldx + xf(0,1)*oldy + xf(0,2)*oldz + xf(0,3);
		vs[i][1] = xf(1,0)*oldx + xf(1,1)*oldy + xf(1,2)*oldz + xf(1,3);
		vs[i][2] = xf(2,0)*oldx + xf(2,1)*oldy + xf(2,2)*oldz + xf(2,3);
	}
}
*/