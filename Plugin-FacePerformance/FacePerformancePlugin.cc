#include "FacePerformancePlugin.hh"
void add_face(std::string filename)
{
	std::ofstream of(filename, std::ios::app);
	std::ifstream in("face.txt");
	std::string str;
	while(in>>str)
	{
		of<<str<<' ';
		in>>str;
		of<<str<<' ';
		in>>str;
		of<<str<<' ';
		in>>str;
		of<<str<<std::endl;
	}
	in.close();
	of.close();
}
void PrintSparseMatrix(const SparseMatrixXd& matrix, const std::string & filename)
{
	std::ofstream of(filename.c_str());
	for(int i=0; i<100; i++)
	{
		for(int j=0; j<100; j++)
		{
			of<<matrix.coeff(i,j)<<' ';
		}
		of<<std::endl;
	}
	of.close();
}
void saveSparseMatrix(SparseMatrixXd& matrix, const std::string & filename)
{
	std::ofstream of(filename);
	int nonZeros = matrix.nonZeros();
	//matrix.
	double* value = matrix.valuePtr();
	int* row = matrix.innerIndexPtr();
	int* col = matrix.outerIndexPtr();
	of<<matrix.rows()<<' '<<matrix.cols()<<std::endl;
	int j=-1;
	for(int i=0; i<nonZeros; i++)
	{
		while(i>=col[j+1])j++;
		of<<row[i]<<' '<<j<<' '<<value[i]<<std::endl;
		
	}
	of.close();
}
void loadSparseMatrix(SparseMatrixXd& matrix, const std::string & filename)
{
	int rows,cols;
	int row,col;
	double val;
	std::ifstream in(filename);
	in>>rows>>cols;
	matrix.resize(rows, cols);
	std::vector<T> triplets;
	while(in>>row>>col>>val)
	{
		triplets.push_back(T(row,col,val));
	}
	matrix.setFromTriplets(triplets.begin(), triplets.end());
	in.close();
}
void load_slection_idx_to_sparsematrix(const std::string& filename, SparseMatrixXd & spm)
{
	FILE * fin = fopen(filename.c_str(),"r");
	int idx;
	std::vector<T> triplets_T;
	while(fscanf(fin,"%d;",&idx) != EOF)
	{
		for(int i=0; i<3; i++)
		{
			triplets_T.push_back(T(idx*3+i, idx*3+i, 1.0));
		}
	}
	spm.setFromTriplets(triplets_T.begin(), triplets_T.end());
}
FacePerformancePlugin::FacePerformancePlugin()
{
}

void FacePerformancePlugin::initializePlugin()
{
	of.open("error_faceperformance.txt");
	toolBox = new QWidget();
	loadButton = new QPushButton("&loadBlendshape",toolBox);
	openKinectButton = new QPushButton("&OpenKinect", toolBox);
	initRegisterButton = new QPushButton("InitRegister", toolBox);
	readFrameButton = new QPushButton("stopReadFrame", toolBox);
	trackingButton = new QPushButton("tracking", toolBox);

	layout = new QGridLayout(toolBox);
	layout->addWidget(loadButton,0,0);
	layout->addWidget(openKinectButton,1,0);
	layout->addWidget(readFrameButton,2,0);
	layout->addWidget(initRegisterButton,3,0);
	layout->addWidget(trackingButton,4,0);
	connect(loadButton, SIGNAL(clicked()), this, SLOT(loadBlendshapes()));
	connect(openKinectButton, SIGNAL(clicked()), this, SLOT(openKinect()));
	connect(readFrameButton, SIGNAL(clicked()), this, SLOT(stopReadFrame()));
	connect(initRegisterButton, SIGNAL(clicked()), this, SLOT(initRegisterMesh()));
	connect(trackingButton, SIGNAL(clicked()), this, SLOT(trackingMesh()));
	timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(ReadFrame()));
	colorImage.create(480, 640, CV_8UC3);
	depthImage.create(480, 640, CV_64FC1);
	

	emit addToolbox(tr("load"), toolBox);
	depthMesh = new trimesh::TriMesh;
	fx = 1.2017214;
	fy = 0.9030345;
	
	GTHGF.resize(num_blendshape-1);
	of<<"init finish"<<std::endl;
	
}

void FacePerformancePlugin::compute_S_face(TriMesh * mesh, OpenMesh::FPropHandleT<Eigen::MatrixXd>& S, bool isInverse)//计算每个面的变换矩阵S
{
	mesh->add_property(S, "Rt_matrix");
	TriMesh::FaceIter f_it, f_end = mesh->faces_end();
	for( f_it = mesh->faces_begin(); f_it != f_end; ++f_it )
	{
		TriMesh::Point point[4];
		TriMesh::FaceHalfedgeIter fv_it(*mesh, f_it);
		int i = 0;
		for( ; fv_it; ++fv_it )
		{
			point[i++] = mesh->point(mesh->from_vertex_handle(fv_it));
		}
		TriMesh::Point cross_dot = (point[1] - point[0]) % (point[2] - point[0]);
		point[i] = point[0] + cross_dot/ sqrt(cross_dot.norm());

		Eigen::MatrixXd cur_S(3,3);
		TriMesh::Point cur_point = point[1] - point[0];
		cur_S(0,0) = cur_point[0];
		cur_S(1,0) = cur_point[1];
		cur_S(2,0) = cur_point[2];
		cur_point = point[2] - point[0];
		cur_S(0,1) = cur_point[0];
		cur_S(1,1) = cur_point[1];
		cur_S(2,1) = cur_point[2];
		cur_point = point[3] - point[0];
		cur_S(0,2) = cur_point[0];
		cur_S(1,2) = cur_point[1];
		cur_S(2,2) = cur_point[2];
		if(isInverse)cur_S = cur_S.inverse();
		mesh->property(S, f_it) = cur_S;
	}
}

void FacePerformancePlugin::compute_H_Star(TriMesh * neutral_mesh, OpenMesh::FPropHandleT<Eigen::MatrixXd>& S_neutral, TriMesh * blendshape_mesh, OpenMesh::FPropHandleT<Eigen::MatrixXd>& S_blendshape, SparseMatrixXd & H_Star)
{
	TriMesh::FaceIter f_it_neutral, f_it_blendshape, f_end = neutral_mesh->faces_end();
	int H_row = 0;
	std::vector<T> triplets_H;
	//std::ofstream of("out.txt");
	for( f_it_neutral = neutral_mesh->faces_begin(), f_it_blendshape = blendshape_mesh->faces_begin(); f_it_neutral != f_end; ++f_it_neutral, ++f_it_blendshape)
	{
		Eigen::MatrixXd m_neutral = neutral_mesh->property(S_neutral, f_it_neutral);
		Eigen::MatrixXd m_blendshape = blendshape_mesh->property(S_blendshape, f_it_blendshape);
		Eigen::MatrixXd S_Star =  m_blendshape * m_neutral;
		for(int k=0; k<2; k++)
		{
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					triplets_H.push_back(T(H_row+i, H_row+j, S_Star(i,j)));
					//of<<S_Star(i,j)<<' '<<m_neutral(i,j)<<' '<<m_blendshape<<' ';
				}
				//of<<std::endl;
			}
			H_row += 3;
		}
		emit log(LOGERR, QString::number(H_row));
	}
	H_Star.setFromTriplets(triplets_H.begin(), triplets_H.end());
		
}
void FacePerformancePlugin::computeBlendShapeByNeutral(TriMesh* neutral_mesh, int id_mesh, Eigen::SimplicialLDLT<SparseMatrixXd>& LDLT_solver, const SparseMatrixXd& GTHGF)//SparseMatrixXd& T_Star)
{
	int num_vertices = neutral_mesh->n_vertices();
	Eigen::VectorXd vec(num_vertices*3);
	TriMesh::VertexIter v_it, v_end = neutral_mesh->vertices_end();
	for(v_it = neutral_mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		TriMesh::Point point = neutral_mesh->point(v_it);
		int idx = v_it->idx();
		for(int i=0; i<3; i++)
		{
			vec(idx*3+i) = point[i];
		}
	}
	vec = GTHGF * vec;
	vec = LDLT_solver.solve(vec);
	for(v_it = neutral_mesh->vertices_begin(); v_it != v_end; ++v_it)
	{
		TriMesh::Point point = neutral_mesh->point(v_it);
		int idx = v_it->idx();
		neutral_mesh->point(v_it) = TriMesh::Point(vec(idx*3), vec(idx*3+1), vec(idx*3+2));
	}
	OpenMesh::IO::write_mesh(*neutral_mesh, "transfer.obj");
	emit updateObject(id_mesh, UPDATE_GEOMETRY);
	emit updateView();
	//emit save(id_mesh, "result.obj");
}
void FacePerformancePlugin::computeDeltaBByNeutral(Eigen::VectorXd& neutral_mesh, Eigen::MatrixXd& delta_B)//Eigen::SimplicialLDLT<SparseMatrixXd>& LDLT_solver, const SparseMatrixXd& GTHGF)//SparseMatrixXd& T_Star)
{
	delta_B.resize(neutral_mesh.size(), num_blendshape-1);
	for(int i=0; i<delta_B.cols(); i++)
	{
		//delta_B.col(i) = GTHGF[i] * neutral_mesh;
		delta_B.col(i) = LDLT_solver.solve(GTHGF[i] * neutral_mesh)-neutral_mesh;
	}
	/*
	for(int i=0; i<delta_B.cols(); i++)
	{
		QString dir("delta_B/shape_");
		ofstream out(dir.append(QString::number(i+1)).append(".obj").toStdString());
		out<<"# OBJ"<<std::endl;
		for(int j=0; j<delta_B.rows(); j+=3)
		{
			out<<"v "<<delta_B(j, i)<<' '<<delta_B(j+1, i)<<' '<<delta_B(j+2, i)<<std::endl;
		}
		out.close();
	}
	for(int i=0; i<delta_B.cols(); i++)
	{
		//delta_B.col(i) = GTHGF[i] * neutral_mesh;
		delta_B.col(i) = delta_B.col(i);
	}*/
}
void FacePerformancePlugin::compute_T_Star()
{
	QString tstar_dir("T-Star/");
	DataType dataType = DATA_TRIANGLE_MESH;
	TriMesh * blendshape_mesh = new TriMesh;
	TriMesh * neutral_mesh = new TriMesh;
	int i=0;
	QString path = "E:/download/FaceWareHouse/FaceWarehouse_Data_0_triangle/Tester_1/Blendshape/shape_";
	QString num = QString::number(i);
	QString suffix = ".obj";
	QString fileName = path.append(num).append(suffix);
	OpenMesh::IO::read_mesh(*neutral_mesh, fileName.toStdString());
	if(neutral_mesh != NULL)
	{
		double mu = 100;
		//emit log(LOGERR, num);
		emit log(LOGERR, QString("load: ").append(QString::number(i)));
		int num_vertices = neutral_mesh->n_vertices();
		int num_faces = neutral_mesh->n_faces();
		
		SparseMatrixXd G(num_faces*3*2, num_vertices*3);
		SparseMatrixXd F(num_vertices*3, num_vertices*3);
		std::string F_selection_filename = "F-selection.txt";
		load_slection_idx_to_sparsematrix(F_selection_filename, F);//加载约束点idx

		//PrintSparseMatrix(F, "F.txt");
		TriMesh::FaceIter f_it, f_end = neutral_mesh->faces_end();
		int G_row = 0;
		std::vector<T> triplets_G;
		//std::ofstream of1("Gtest.txt");
		for(f_it = neutral_mesh->faces_begin(); f_it != f_end; ++f_it)  //构造G
		{
			TriMesh::FaceHalfedgeIter fv_it(*neutral_mesh, f_it);
			int first_idx = -1;
			for( ; fv_it; ++fv_it)
			{
				if(first_idx == -1)
				{
					first_idx = neutral_mesh->from_vertex_handle(fv_it).idx();
					//of1<<first_idx<<std::endl;
				}
				else 
				{
					int other_idx = neutral_mesh->from_vertex_handle(fv_it).idx();
					//of1<<other_idx<<std::endl;
					for(int i=0; i<3; i++)
					{
						triplets_G.push_back( T(G_row, other_idx*3+i, 1.0) );
						triplets_G.push_back( T(G_row, first_idx*3+i, -1.0) );
						G_row++;
					}
				}
			}
		}
		//of1.close();
		
		G.setFromTriplets(triplets_G.begin(), triplets_G.end());
		
		//PrintSparseMatrix(G, "G.txt");////////////////////////////////////////////////////////////
		GTGF.resize(num_vertices*3, num_vertices*3);
		GTGF = G.transpose()*G+mu*F;
		saveSparseMatrix(GTGF, QString("").append(tstar_dir).append("GTGF.txt").toStdString());//保存
		LDLT_solver.compute(GTGF);
		//PrintSparseMatrix(GTGF, "GTGF.txt");////////////////////////////////////////////////////////////
		OpenMesh::FPropHandleT<Eigen::MatrixXd> S_neutral, S_blendshape;
		compute_S_face(neutral_mesh, S_neutral, true); //计算中性脸各个面的变换矩阵
		for(int i=1; i<num_blendshape; i++) //计算blendshape各个面的变换矩阵，并计算其H*矩阵
		{
			of<<"I: "<<i<<std::endl;
			GTHGF[i-1].resize(num_vertices*3, num_vertices*3);
			SparseMatrixXd T_Star(num_vertices*3, num_vertices*3);
			SparseMatrixXd H_Star(num_faces*3*2, num_faces*3*2);
			//PluginFunctions::getMesh(id_blendshape[i], blendshape_mesh);
			QString path = "E:/download/FaceWareHouse/FaceWarehouse_Data_0_triangle/Tester_1/Blendshape/shape_";
			QString num = QString::number(i);
			QString suffix = ".obj";
			QString fileName = path.append(num).append(suffix);
			OpenMesh::IO::read_mesh(*blendshape_mesh, fileName.toStdString());

			emit log(LOGERR, QString("S_face Start: ").append(QString::number(i)));
			compute_S_face(blendshape_mesh, S_blendshape); //计算blendshape各个面的变换矩阵
			emit log(LOGERR, QString("H_star Start: ").append(QString::number(i)));
			compute_H_Star(neutral_mesh, S_neutral, blendshape_mesh, S_blendshape, H_Star);//计算H*
			//PrintSparseMatrix(H_Star, "H.txt");
			GTHGF[i-1] = G.transpose() * H_Star * G + mu * F;
			of<< QString("").append(tstar_dir).append("GTHGF_").append(QString::number(i-1)).append(".txt").toStdString()<<std::endl;
			saveSparseMatrix(GTHGF[i-1], QString("").append(tstar_dir).append("GTHGF_").append(QString::number(i-1)).append(".txt").toStdString());//保存
			of<< QString("").append(tstar_dir).append("GTHGF_").append(QString::number(i-1)).append(".txt").toStdString()<<std::endl;
			//PrintSparseMatrix(GTHGF, "GTHGF.txt");////////////////////////////////////////////////////////////
			
			//std::ofstream of("out2.txt");
			////////////////////////////////////////////
			//TriMesh* tmesh = new TriMesh;
			//OpenMesh::IO::read_mesh(*tmesh,"E:/download/FaceWareHouse/FaceWarehouse_Data_0_triangle/Tester_2/Blendshape/shape_0.obj");
			//////////////////////////////////////////////
			//computeBlendShapeByNeutral(tmesh, id_blendshape[0], LDLT_solver, GTHGF[i-1]);//T_Star);
			blendshape_mesh->remove_property(S_blendshape);
		}
		neutral_mesh->remove_property(S_neutral);
		emit log(LOGERR, QString("G norm:").append(QString::number(G.norm())));
	}
	else 
	{
		log(LOGERR, "Is not mesh");
	}
	delete neutral_mesh;
	delete blendshape_mesh;
}

void FacePerformancePlugin::load_T_Star()
{
	QString tstar_dir("T-Star/");
	loadSparseMatrix(GTGF, QString("").append(tstar_dir).append("GTGF.txt").toStdString());
	of<<QString("").append(tstar_dir).append("GTGF.txt").toStdString()<<std::endl;
	LDLT_solver.compute(GTGF);
	for(int i=1; i<num_blendshape; i++)
	{
		of<<QString("").append(tstar_dir).append("GTHGF_").append(QString::number(i-1)).append(".txt").toStdString()<<std::endl;
		loadSparseMatrix(GTHGF[i-1], QString("").append(tstar_dir).append("GTHGF_").append(QString::number(i-1)).append(".txt").toStdString());
	}
}

void FacePerformancePlugin::loadBlendshapes()
{
	DataType dataType = DATA_TRIANGLE_MESH;
	for(int i=0; i<num_blendshape; i++)
	{
		QString path = "E:/download/FaceWareHouse/FaceWarehouse_Data_0_triangle/Tester_1/Blendshape/shape_";
		QString num = QString::number(i);
		QString suffix = ".obj";
		QString fileName = path.append(num).append(suffix);
		emit load(fileName, dataType, id_blendshape[i]);
		if(id_blendshape[i]==-1)
		{
			emit log(LOGERR, num);
		}
		else
		{
			emit log(LOGERR, QString("load: ").append(QString::number(i)));
		}
	}
	compute_T_Star();
}
void FacePerformancePlugin::openKinect()
{
	//ui.OpenKinect->setEnabled(false);
	//ui.DepthCollection->setEnabled(true);
	//m_DepthCollection = 0;
	//m_depthiter = 0;
	m_colorRGBX = new BYTE[640 * 480 * 4];
	m_colorCoordinates = new LONG[640 * 480 * 2];
	
	NuiCreateSensorByIndex(0, &m_pNuiSensor);
	DWORD NuiFlag = NUI_INITIALIZE_FLAG_USES_COLOR | NUI_INITIALIZE_FLAG_USES_DEPTH;
	m_pNuiSensor->NuiInitialize(NuiFlag);
	h1 = CreateEvent(NULL, TRUE, FALSE, NULL);
	h2 = NULL;
	h3 = CreateEvent(NULL, TRUE, FALSE, NULL);
	h4 = NULL;
	m_pNuiSensor->NuiImageStreamOpen(NUI_IMAGE_TYPE_COLOR, NUI_IMAGE_RESOLUTION_640x480, 0, 2, h1, &h2);
	m_pNuiSensor->NuiImageStreamOpen(NUI_IMAGE_TYPE_DEPTH, NUI_IMAGE_RESOLUTION_640x480, 0, 2, h3, &h4);
	m_pNuiSensor->NuiImageStreamSetImageFrameFlags(h4,NUI_IMAGE_STREAM_FLAG_ENABLE_NEAR_MODE|NUI_IMAGE_STREAM_FLAG_DISTINCT_OVERFLOW_DEPTH_VALUES);
	isDepthCollection = true;
	cv::namedWindow("colorImage", CV_WINDOW_AUTOSIZE);
	cv::namedWindow("depthImage", CV_WINDOW_AUTOSIZE);
	//ReadFrame();
	//timer->start(1000);
	isStopReadFrame = false;
	featureExtraction = new FeatureExtraction(640, 480);
	timer->start(34);
	/*FeatureExtraction featureExtraction;
	cv::Mat colorImage;
	colorImage.create(480, 640, CV_8UC3); 
	cv::Mat depthImage;
	depthImage.create(480, 480, CV_8UC3); 

	HANDLE colorEvent = CreateEvent( NULL, TRUE, FALSE, NULL ); 
    HANDLE depthEvent = CreateEvent( NULL, TRUE, FALSE, NULL ); 

 
    HANDLE colorStreamHandle = NULL; 
    HANDLE depthStreamHandle = NULL; 
	
	 HRESULT hr;
	 int iSensorCount=0;
	 hr=NuiGetSensorCount(&iSensorCount);

	 if(FAILED(hr))
	 {
		 emit log(LOGERR, "fail!");
		 return;
	 }
	 for(int i=0;i<iSensorCount;i++)
	 {
		 hr=NuiCreateSensorByIndex(i,&pNuiSensor);

		 if(FAILED(hr))
		  {
			 emit log(LOGERR, "fail!");
			 return;
	      }
	 

	     hr=pNuiSensor->NuiStatus();
	     if(hr==S_OK){break;}
		 pNuiSensor->Release();
     }
     
	//HRESULT hr = NuiInitialize(NUI_INITIALIZE_FLAG_USES_COLOR | NUI_INITIALIZE_FLAG_USES_DEPTH_AND_PLAYER_INDEX);   
	 hr = NuiInitialize(NUI_INITIALIZE_FLAG_USES_COLOR | NUI_INITIALIZE_FLAG_USES_DEPTH_AND_PLAYER_INDEX); 
    if( hr != S_OK )   
    {   
		emit log(LOGERR, "NuiInitialize failed");   
        return;   
    } 
 
    hr = NuiImageStreamOpen(NUI_IMAGE_TYPE_COLOR, NUI_IMAGE_RESOLUTION_640x480, NULL, 4, colorEvent, &colorStreamHandle); 
    if( hr != S_OK )   
    {   
        emit log(LOGERR,"Open the color Stream failed"); 
        NuiShutdown(); 
        return;   
    } 
    hr = NuiImageStreamOpen(NUI_IMAGE_TYPE_DEPTH_AND_PLAYER_INDEX, NUI_IMAGE_RESOLUTION_640x480, NULL, 2, depthEvent, &depthStreamHandle); 
    if( hr != S_OK )   
    {   
        emit log(LOGERR, "Open the depth Stream failed"); 
        NuiShutdown(); 
        return;   
    } 

	cv::namedWindow("colorImage", CV_WINDOW_AUTOSIZE);
	cv::namedWindow("depthImage", CV_WINDOW_AUTOSIZE);
	
    while (1) 
    { 
        if(WaitForSingleObject(colorEvent, 0)==0) 
		{
            getColorImage(colorEvent, colorStreamHandle, colorImage); 
			featureExtraction.GetFacePoint(colorImage);
		}
        if(WaitForSingleObject(depthEvent, 0)==0) 
            getDepthImage(depthEvent, depthStreamHandle, depthImage); 
           
        imshow("colorImage", colorImage); 
        imshow("depthImage", depthImage); 
		
        if(cvWaitKey(1)==27) 
            break; 
    }
    NuiShutdown(); */
}
void FacePerformancePlugin::getColorImage(HANDLE &colorEvent, HANDLE &colorStreamHandle, cv::Mat &colorImage)
{
	const NUI_IMAGE_FRAME *colorFrame = NULL; 
 
    NuiImageStreamGetNextFrame(colorStreamHandle, 0, &colorFrame); 
    INuiFrameTexture *pTexture = colorFrame->pFrameTexture;   
 
    NUI_LOCKED_RECT LockedRect; 
    pTexture->LockRect(0, &LockedRect, NULL, 0);   
 
    if( LockedRect.Pitch != 0 ) 
    { 
		for (int i=0; i<colorImage.rows; i++) 
        {
			uchar *ptr = colorImage.ptr<uchar>(i); 
            uchar *pBuffer = (uchar*)(LockedRect.pBits) + i * LockedRect.Pitch;
            for (int j=0; j<colorImage.cols; j++) 
            { 
                ptr[3*j] = pBuffer[4*j];  
                ptr[3*j+1] = pBuffer[4*j+1]; 
                ptr[3*j+2] = pBuffer[4*j+2]; 
            } 
		} 
    } 
    else 
    { 
        emit log(LOGERR, "Get Color Image Failed!"); 
    }

	pTexture->UnlockRect(0); 
    NuiImageStreamReleaseFrame(colorStreamHandle, colorFrame );
}
void FacePerformancePlugin::getDepthImage(HANDLE &depthEvent, HANDLE &depthStreamHandle, cv::Mat &depthImage)
{
	 const NUI_IMAGE_FRAME *depthFrame = NULL; 
 
    NuiImageStreamGetNextFrame(depthStreamHandle, 0, &depthFrame); 
    INuiFrameTexture *pTexture = depthFrame->pFrameTexture;   
 
    NUI_LOCKED_RECT LockedRect; 
    pTexture->LockRect(0, &LockedRect, NULL, 0);   
 
    RGBQUAD q; 

    if( LockedRect.Pitch != 0 ) 
    { 
        for (int i=0; i<depthImage.rows; i++) 
        { 
            uchar *ptr = depthImage.ptr<uchar>(i); 
            uchar *pBuffer = (uchar*)(LockedRect.pBits) + i * LockedRect.Pitch;
			USHORT *pBufferRun = (USHORT*) pBuffer; 
			
            for (int j=0; j<depthImage.cols; j++) 
            { 
                int player = pBufferRun[j]&7; 
                int data = (pBufferRun[j]&0xfff8) >> 3; 
                 
                uchar imageData = 255-(uchar)(256*data/0x0fff); 
                q.rgbBlue = q.rgbGreen = q.rgbRed = 0; 
 
                switch(player) 
                { 
					case 0:   
						q.rgbRed = imageData / 2;   
						q.rgbBlue = imageData / 2;   
						q.rgbGreen = imageData / 2;   
						break;   
					case 1:    
						q.rgbRed = imageData;   
						break;   
					case 2:   
						q.rgbGreen = imageData;   
						break;   
					case 3:   
						q.rgbRed = imageData / 4;   
						q.rgbGreen = q.rgbRed*4;  
						q.rgbBlue = q.rgbRed*4;  
						break;   
					case 4:   
						q.rgbBlue = imageData / 4;  
						q.rgbRed = q.rgbBlue*4;   
						q.rgbGreen = q.rgbBlue*4;   
						break;   
					case 5:   
						q.rgbGreen = imageData / 4;  
						q.rgbRed = q.rgbGreen*4;   
						q.rgbBlue = q.rgbGreen*4;   
						break;   
					case 6:   
						q.rgbRed = imageData / 2;   
						q.rgbGreen = imageData / 2;    
						q.rgbBlue = q.rgbGreen*2;   
						break;   
					case 7:   
						q.rgbRed = 255 - ( imageData / 2 );   
						q.rgbGreen = 255 - ( imageData / 2 );   
						q.rgbBlue = 255 - ( imageData / 2 ); 
				} 	 
				ptr[3*j] = q.rgbBlue; 
				ptr[3*j+1] = q.rgbGreen; 
				ptr[3*j+2] = q.rgbRed; 
            } 
        } 
    } 
    else 
    { 
        emit log(LOGERR, "Get Depth Image Failed!"); 
    } 
	
	pTexture->UnlockRect(0);
    NuiImageStreamReleaseFrame(depthStreamHandle, depthFrame);  
}
void FacePerformancePlugin::ReadFrame()
{
	static int cnt = 0;
	QString data_dir("kinect-data/");
	std::vector<double> marker_Coordinates; //标记点的坐标

	if(isStopReadFrame)return;
	NUI_IMAGE_FRAME *pColorImageFrame;
	NUI_IMAGE_FRAME *pDepthImageFrame;
	pColorImageFrame = new NUI_IMAGE_FRAME[1];
	pDepthImageFrame = new NUI_IMAGE_FRAME[1];
	if (!WaitForSingleObject(h1, INFINITE))
	{
		m_pNuiSensor->NuiImageStreamGetNextFrame(h2,0,&pColorImageFrame[0]);
		INuiFrameTexture * pTexture = pColorImageFrame[0].pFrameTexture;
		NUI_LOCKED_RECT LockedRect;
		pTexture->LockRect(0,&LockedRect,NULL,0);
		if (LockedRect.Pitch != 0)
		{
			if (isDepthCollection)
			{
				BYTE* pBuffer = (BYTE*)LockedRect.pBits;
				memcpy(m_colorRGBX, LockedRect.pBits, LockedRect.size);
				//QImage colorimg = QImage(pBuffer, 640, 480, QImage::Format_RGB32);
				for (int i=0; i<colorImage.rows; i++)
				{
					uchar *ptr = colorImage.ptr<uchar>(i); 
					uchar *pBuffer = (uchar*)(LockedRect.pBits) + i * LockedRect.Pitch;
					for (int j=0; j<colorImage.cols; j++) 
					{ 
						ptr[3*j] = pBuffer[4*j];  
						ptr[3*j+1] = pBuffer[4*j+1]; 
						ptr[3*j+2] = pBuffer[4*j+2]; 
					}
				}
				//if(cnt<5)
				{
					cv::imwrite(QString("kinect-data/colorImage").append(QString::number(cnt)).append(".jpg").toStdString(), colorImage);
					featureExtraction->GetFacePoint(colorImage);
				}
				
				//cv::write(colorImage, "Texture.jpg");
				//featureExtraction.feature;
			}
		}
		m_pNuiSensor->NuiImageStreamReleaseFrame(h2, &pColorImageFrame[0]);
	}
	
	if (!WaitForSingleObject(h3, INFINITE))
	{
		m_pNuiSensor->NuiImageStreamGetNextFrame(h4, 0, &pDepthImageFrame[0]);
		INuiFrameTexture *pTexture = pDepthImageFrame[0].pFrameTexture;
		NUI_LOCKED_RECT LockedRect;
		pTexture->LockRect(0, &LockedRect, NULL, 0);
		if (LockedRect.Pitch != 0)
		{
			if (isDepthCollection)
			{
				BYTE* pBuffer = (BYTE*)LockedRect.pBits;
				cv::Mat DepthImage(480, 640, CV_16UC1, pBuffer);
				//if (m_depthiter == 10)
				{
					USHORT* depthD16;
					depthD16 = new USHORT[640*480];
					memcpy(depthD16, LockedRect.pBits, LockedRect.size);
					//获取三维点云对应在color空间的坐标
					//m_pNuiSensor->NuiGetCoordinateMapper(INuiCoordinateMapper **pMapping);
					m_pNuiSensor->NuiImageGetColorPixelCoordinateFrameFromDepthPixelFrameAtResolution(
						NUI_IMAGE_RESOLUTION_640x480,
						NUI_IMAGE_RESOLUTION_640x480,
						640*480,
						depthD16,
						640*480*2,
						m_colorCoordinates);
					delete [] depthD16;
				}
				USHORT* depth_line = (USHORT*)DepthImage.data;
				UINT stride = DepthImage.step1();
				USHORT max_depth = 900;
				USHORT min_depth = 0;
				int k = 0;
				
				for (DWORD m = 0; m < 480; m++)
				{
					for (DWORD n = 0; n < 640; n++)
					{
						int real_depth = (depth_line[n] >> 3);
						k = m * 640 + n;
						if (real_depth < max_depth && real_depth > min_depth)
						{
							//m_depth[k] = real_depth;
							depthImage.at<double>(m,n) = real_depth;
						}
						else depthImage.at<double>(m,n) = -1;
					}
					depth_line += stride;
				}
				depthMesh->clear();
				m_DepthIndex.clear();////
				
				//if(cnt<5)
				{
					int i=0;
					for(i=0; i<featureExtraction->feature.size(); i++)
					{
						//if (featureExtraction.feature[i].x && featureExtraction.feature[i].y)
						bool isFind = false;
						int m_data = 0;
						int x = (int)(featureExtraction->feature[i].x + 0.5);
						int y = (int)(featureExtraction->feature[i].y + 0.5);
						for (int m = 0; m < 480; m++)//添加trimesh的vertices中
						{
							for (int n = 0; n < 640; n++)
							{
								int k = 640 * m + n;
								LONG colorInDepthX = m_colorCoordinates[k * 2];
								LONG colorInDepthY = m_colorCoordinates[k * 2 + 1];
								if (colorInDepthX >= 0 && colorInDepthX < 640 && colorInDepthY >= 0 && colorInDepthY < 480 && depthImage.at<double>(m,n) != -1)
								{	
									if (colorInDepthX == x && colorInDepthY == y)
									{
										m_DepthIndex.push_back(m_data);//[i] = m_DataIndex;
										//featureExtraction->feature[i].x = 0;
										//featureExtraction->feature[i].y = 0;
										marker_Coordinates.push_back(colorInDepthX*1.0/90.0);
										marker_Coordinates.push_back(colorInDepthY*1.0/90.0);
										isFind = true;
										break;
									}
									m_data ++;
								}
							}
							if(isFind)break;
						}
						if(!isFind)
						{
							m_DepthIndex.push_back(-1);
							//marker_Coordinates.push_back(-1);
							//marker_Coordinates.push_back(-1);
						}
					}
					//if(i>=featureExtraction.feature.size())m_DepthIndex.push_back(-1);
				}
				int m_data = 0;

				for (int m = 0; m < 480; m++)//添加trimesh的vertices中
				{
					for (int n = 0; n < 640; n++)
					{
						int k = 640 * m + n;
						LONG colorInDepthX = m_colorCoordinates[k * 2];
						LONG colorInDepthY = m_colorCoordinates[k * 2 + 1];
						if (colorInDepthX >= 0 && colorInDepthX < 640 && colorInDepthY >= 0 && colorInDepthY < 480 && depthImage.at<double>(m,n) != -1)
						{
							//
							trimesh::point p,q;
							p[0] = (n/640.0 - 0.5) * (depthImage.at<double>(m,n))* fx;// 
							p[1] = (0.5 - m/480.0) * (depthImage.at<double>(m,n))* fy;// 
							p[2] = 900-depthImage.at<double>(m,n);
							p[0] /= 90.0;
							p[1] /= 90.0;
							p[2] /= 90.0;
							depthMesh->vertices.push_back(p);
							
							m_data ++;
						}
					}
				}
				//depthMesh->write("depthMesh.obj");
				depthMesh->write(QString("kinect-data/depthMesh").append(QString::number(cnt)).append(".obj").toStdString());
				//if(cnt<5)
				{
					FILE *f;
					f = fopen(QString("kinect-data/DepthIndex").append(QString::number(cnt)).append(".txt").toStdString().c_str(), "w");
					for (int i = 0; i < m_DepthIndex.size(); i++)
					{
						fprintf(f, "%d\n", m_DepthIndex[i]);
					}
					fclose(f);
				}
				//计算3d_to_2d的投影矩阵
				Eigen::MatrixXd A(depthMesh->vertices.size(), 4);
				Eigen::MatrixXd b(depthMesh->vertices.size(), 2);
				m_data = 0;
				for (int m = 0; m < 480; m++)//添加trimesh的vertices中
				{
					for (int n = 0; n < 640; n++)
					{
						int k = 640 * m + n;
						LONG colorInDepthX = m_colorCoordinates[k * 2];
						LONG colorInDepthY = m_colorCoordinates[k * 2 + 1];
						if (colorInDepthX >= 0 && colorInDepthX < 640 && colorInDepthY >= 0 && colorInDepthY < 480 && depthImage.at<double>(m,n) != -1)
						{
							for(int j=0; j<3; j++) A(m_data, j) = depthMesh->vertices[m_data][j];
							A(m_data, 3) = 1;
							b(m_data, 0) = colorInDepthX*1.0/90.0;
							b(m_data, 1) = colorInDepthY*1.0/90.0;
							//
							//trimesh::point p,q;
							//p[0] = (n/640.0 - 0.5) * (depthImage.at<double>(m,n))* fx;// 
							//p[1] = (0.5 - m/480.0) * (depthImage.at<double>(m,n))* fy;// 
							//p[2] = 900-depthImage.at<double>(m,n);
							//depthMesh->vertices.push_back(p);
							m_data ++;
							
						}
					}
				}
				Eigen::MatrixXd RT = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
				ofstream of_3d_2_2d(QString("kinect-data/2d_2_3d_").append(QString::number(cnt)).append(".txt").toStdString());
				of_3d_2_2d<<RT<<std::endl;
				of_3d_2_2d.close();
			}
		}
		m_pNuiSensor->NuiImageStreamReleaseFrame(h4, &pDepthImageFrame[0]);
	}
	//delete [] m_colorCoordinates;
	//delete [] m_colorRGBX;
	//delete [] m_Cout;
	//delete [] m_Sum;
	imshow("colorImage", colorImage);
    imshow("depthImage", depthImage);
	
	cv::imwrite(QString("kinect-data/colorImage_landmark").append(QString::number(cnt)).append(".jpg").toStdString(), colorImage);
	cv::imwrite(QString("kinect-data/depthImage").append(QString::number(cnt)).append(".jpg").toStdString(), depthImage);

	ofstream of_marker_coord(QString("kinect-data/marker_Coordinate").append(QString::number(cnt)).append(".txt").toStdString());
	for(int i=0; i<marker_Coordinates.size(); i++)
	{
		of_marker_coord<<marker_Coordinates[i]<<' ';
	}
	of_marker_coord<<std::endl;
	of_marker_coord.close();
	//ofstream of_coord(QString("kinect-data/colorCoordinate").append(QString::number(cnt)).append(".txt").toStdString());
	//int num_mn = 480*640*2;
	//for(int k=0 ; k<num_mn; k++) of_coord<<m_colorCoordinates[k]<<' ';of_coord<<std::endl;
	/*
	for (int m = 0; m < 480; m++)//colorCoordinate
	{
		for (int n = 0; n < 640; n++)
		{
			int k = 640 * m + n;
			LONG colorInDepthX = m_colorCoordinates[k * 2];
			LONG colorInDepthY = m_colorCoordinates[k * 2 + 1];
			of_coord<<colorInDepthX<<' '<<colorInDepthY<<std::endl;
		}
	}*/
	//of_coord.close();
	cnt++;
	delete [] pColorImageFrame;
	delete [] pDepthImageFrame;
}
void FacePerformancePlugin::initRegisterMesh()
{
	/*of<<"in"<<std::endl;
	trimesh::TriMesh* depthMesh = trimesh::TriMesh::read("depthMesh.obj");
	trimesh::TriMesh* finalMesh = new trimesh::TriMesh;
	
	
	of<<"initRegister initR"<<std::endl;///////////////
	string meshIndexFilename = "MeshIndex.txt";
	string depthIndexFilename = "DepthIndex.txt";
	InitRegisterMesh initR(depthMesh, meshIndexFilename, depthIndexFilename);
	initR.fit_mesh();
	initR.blendshape2trimesh(finalMesh);
	//std::ifstream fmesh("MeshIndex.txt");
	//std::ifstream fdepth("DepthIndex.txt");

	//of<<"register start"<<std::endl;///////////////////
	//initR.register_mesh(1, finalMesh);

	//of<<"register_end"<<std::endl;////////////////////
	finalMesh->write("finalMesh.obj");

	delete depthMesh;
	delete finalMesh;
	//InitRegisterMesh initRM(depthMesh);*/
}
void FacePerformancePlugin::stopReadFrame()
{
	isStopReadFrame = true;
}
void FacePerformancePlugin::trackingMesh()
{
	of<<"before load T*"<<std::endl;
	//compute_T_Star();//计算T*
	load_T_Star(); //加载T_star, 已经预处理好，只需从文件读出
	of<<"computer_T_Star finish"<<std::endl;
	QString kinect_dir = "kinect-data/";
	string meshIndexFilename = "MeshIndex.txt";
	Eigen::VectorXd register_mesh;
	Eigen::MatrixXd delta_B;
	//double s;
	//Eigen::VectorXd r;
	//Eigen::VectorXd t;
	RigidLMParameters rigidParas;
	NonRigidLMParameters nonrigidParas;
	InitLMParameters initParas;
	nonrigidParas.fx = fx;
	nonrigidParas.fy = fy;
	initParas.maxDepth = nonrigidParas.maxDepth = 90.0;

	double* colorCoordinate = new double[480*640*2];
	nonrigidParas.colorCoordinate = new double[480*640*2];
	nonrigidParas.pre_colorCoordinate = new double[480*640*2];

	double * marker_Coordinates = new double[87*2];
	nonrigidParas.marker_Coordinates = new double[87*2];
	
	//int colorCoordinate[480*640*2];
	//for(int i=4; i<29; i++)
	for(int i=0; i<=60; i++)
	{
		of<<i<<std::endl;
		QString depthIndexFilename = QString("").append(kinect_dir).append("DepthIndex").append(QString::number(i)).append(".txt");
		QString depthMeshFilename = QString("").append(kinect_dir).append("depthMesh").append(QString::number(i)).append(".obj");
		
		if(i==0) //初始注册
		{
			//dianyun
			/*ifstream if_coord(QString("").append(kinect_dir).append("colorCoordinate").append(QString::number(i)).append(".txt").toStdString());
			int num_coord = 480*640*2;
			for(int j=0; j<num_coord; j++)if_coord>>colorCoordinate[j];
			if_coord.close();*/
			/*paras.pointCloudMesh = pointCloudMesh;
	paras.texture = texture;
	paras.colorCoordinate = colorCoordinate;
	paras.marker_Coordinates = marker_Coordinates;*/
			
			ifstream if_3d_2_2d(QString("").append(kinect_dir).append("2d_2_3d_").append(QString::number(i)).append(".txt").toStdString());
			initParas._3d_to_2d_matrix.resize(4,2);
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<2; k++)
				{
					if_3d_2_2d>>initParas._3d_to_2d_matrix(j,k);
				}
			}
			if_3d_2_2d.close();

			ifstream if_coord(QString("").append(kinect_dir).append("marker_Coordinate").append(QString::number(i)).append(".txt").toStdString());
			double index;
			int j=0;
			while(if_coord>>index)
			{
				marker_Coordinates[j++] = index;
			}
			if_coord.close();
			//of<<"read coord finish"<<std::endl;

			cv::Mat initTexture = cv::imread( QString("").append(kinect_dir).append("colorImage").append(QString::number(i)).append(".jpg").toStdString(), CV_LOAD_IMAGE_COLOR);
			trimesh::TriMesh * depthMesh = trimesh::TriMesh::read(depthMeshFilename.toStdString());of<<"read coord finish1"<<std::endl;
			initParas.pointCloudMesh = depthMesh;
			initParas.meshIndexFilename = meshIndexFilename;
			initParas.depthIndexFilename = depthIndexFilename.toStdString();
			initParas.texture = initTexture;
			initParas.colorCoordinate = colorCoordinate;
			initParas.marker_Coordinates = marker_Coordinates;
			InitRegisterMesh initR(&initParas);//depthMesh, meshIndexFilename, depthIndexFilename.toStdString(), initTexture, colorCoordinate, marker_Coordinates);of<<"read coord finish2"<<std::endl;
			initR.fit_mesh();
			initR.getinitParas(nonrigidParas.expression_blendshape, rigidParas.s, rigidParas.r, rigidParas.t, nonrigidParas.pre_3d_to_2d_matrix);
			//获取delta_B
			computeDeltaBByNeutral(nonrigidParas.expression_blendshape, nonrigidParas.delta_B); //register_mesh为初始化后的中性脸
			trimesh::TriMesh * finalMesh = new trimesh::TriMesh;
			initR.blendshape2trimesh(finalMesh);
			finalMesh->write("finalMesh.obj");
			delete finalMesh;
			of<<"init RegisterMesh finish"<<std::endl;

			nonrigidParas.pre_x.resize(num_blendshape-1);
			nonrigidParas.ppre_x.resize(num_blendshape-1);
			nonrigidParas.pre_x.setZero();
			nonrigidParas.ppre_x.setZero();

			nonrigidParas.neutral_blendshape = nonrigidParas.expression_blendshape;
			nonrigidParas.pre_colorImage = cv::imread( QString("").append(kinect_dir).append("colorImage").append(QString::number(i)).append(".jpg").toStdString(), CV_LOAD_IMAGE_COLOR);
			//cv::imwrite(QString("colorImage").append(QString::number(i+1)).append(".jpg").toStdString(), nonrigidParas.pre_colorImage);
			nonrigidParas.pre_colorCoordinate = colorCoordinate;
			nonrigidParas.meshIndexFilename = "MeshIndex.txt";
			nonrigidParas.num_frame = 0;
			
		}
		else //realtime tracking
		{
			
			ifstream if_3d_2_2d(QString("").append(kinect_dir).append("2d_2_3d_").append(QString::number(i)).append(".txt").toStdString());
			nonrigidParas._3d_to_2d_matrix.resize(4,2);
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<2; k++)
				{
					if_3d_2_2d>>nonrigidParas._3d_to_2d_matrix(j,k);
				}
			}
			if_3d_2_2d.close();

			trimesh::TriMesh * depthMesh = trimesh::TriMesh::read(depthMeshFilename.toStdString());

			//刚性对齐
			rigidParas.pointCloudMesh = depthMesh;
			rigidParas.expression_blendshape = nonrigidParas.expression_blendshape; //复制
			of<<"before rigid"<<std::endl;
			of<<"pointCloudSize: "<<rigidParas.pointCloudMesh->vertices.size()<<std::endl;
			of<<"s:\n"<<rigidParas.s<<std::endl;
			for(int j=0; j<3; j++)of<<rigidParas.r(j)<<' ';of<<std::endl;
			for(int j=0; j<3; j++)of<<rigidParas.t(j)<<' ';of<<std::endl;

			RigidMesh rigidMesh(&rigidParas);
			rigidMesh.fit_mesh();
			//trimesh::TriMesh * finalMesh = new trimesh::TriMesh;
			//rigidMesh.blendshape2trimesh(finalMesh);
			//finalMesh->write(QString("finalMesh_nonRigid_").append(QString::number(i)).append(".obj").toStdString().c_str());
			//delete finalMesh;
			of<<"after rigid"<<std::endl;
			of<<"s:\n"<<rigidParas.s<<std::endl;
			for(int j=0; j<3; j++)of<<rigidParas.r(j)<<' ';of<<std::endl;
			for(int j=0; j<3; j++)of<<rigidParas.t(j)<<' ';of<<std::endl;

			
			//非刚性注册
			nonrigidParas.depthIndexFilename = depthIndexFilename.toStdString();
			//nonrigidParas.expression_blendshape = rigidParas.expression_blendshape;
			/*ifstream if_coord(QString("").append(kinect_dir).append("colorCoordinate").append(QString::number(i)).append(".txt").toStdString());
			int num_coord = 480*640*2;
			for(int j=0; j<num_coord; j++)if_coord>>nonrigidParas.colorCoordinate[j];
			if_coord.close();*/
			ifstream if_coord(QString("").append(kinect_dir).append("marker_Coordinate").append(QString::number(i)).append(".txt").toStdString());
			double index;
			int j=0;
			while(if_coord>>index)
			{
				nonrigidParas.marker_Coordinates[j++] = index;
			}
			if_coord.close();

			nonrigidParas.colorImage = cv::imread( QString("").append(kinect_dir).append("colorImage").append(QString::number(i)).append(".jpg").toStdString(), CV_LOAD_IMAGE_COLOR);
			cv::imwrite("tex.jpg", nonrigidParas.colorImage);///测试读写。
			nonrigidParas.pointCloudMesh = depthMesh;
			nonrigidParas.R = RigidLMCostFunction::getR(&rigidParas.r(0));
			nonrigidParas.T = rigidParas.t;
			nonrigidParas.s = rigidParas.s;
			of<<"before nonrigid1"<<std::endl;
			NonRigidMesh nonrigidMesh(&nonrigidParas);
			of<<"before nonrigid"<<std::endl;
			nonrigidMesh.fit_mesh();
			of<<"after nonrigid"<<std::endl;
			
			trimesh::TriMesh * finalMesh = new trimesh::TriMesh;
			nonrigidMesh.blendshape2trimesh(finalMesh);
			finalMesh->write(QString("finalMesh_nonRigid_").append(QString::number(i)).append(".obj").toStdString().c_str());
			add_face(QString("finalMesh_nonRigid_").append(QString::number(i)).append(".obj").toStdString());
			delete finalMesh;
			nonrigidParas.pre_colorImage = nonrigidParas.colorImage;
			nonrigidParas.pre_colorCoordinate = nonrigidParas.colorCoordinate;
			nonrigidParas.pre_3d_to_2d_matrix = nonrigidParas._3d_to_2d_matrix;
		}
	}
	//delete [] colorCoordinate;
	//delete [] nonrigidParas.colorCoordinate;
	//delete [] nonrigidParas.pre_colorCoordinate;
}
Q_EXPORT_PLUGIN2( facePerformancePlugin , FacePerformancePlugin );