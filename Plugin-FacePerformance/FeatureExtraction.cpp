#include "FeatureExtraction.h"

FeatureExtraction::FeatureExtraction(int nImgWidth, int nImgHeight)
{
	/*
	face_cascade = new cv::CascadeClassifier;
	if(!face_cascade->load("haarcascade_frontalface_alt.xml"))
	{
		printf("Cascadeclassifier load false!\n");
		return;
	}
	if(fit_asm.Read("my68-1d.amf") == false)
	{
		printf("ASM特征点模版载入失败,退出!\n");
		return ;
	}*/
	//of.open("error_featureExtraction.txt");
	//of = fopen("error_featureExtraction.txt", "w");
	std::ifstream in("key_point.txt");
	int index;
	while(in>>index) {key_point.push_back(index-1);}// fprintf(of, "%d\n", index);fflush(of);}
	in.close();
	/*key_point.push_back(20-1);
	key_point.push_back(18-1);
	key_point.push_back(28-1);
	key_point.push_back(30-1);
	key_point.push_back(49-1);
	key_point.push_back(52-1);
	key_point.push_back(63-1);
	key_point.push_back(67-1);
	key_point.push_back(58-1);
	key_point.push_back(55-1);
	key_point.push_back(78-1);*/

	this->nImgHeight = nImgHeight;
	this->nImgWidth = nImgWidth;
	this->nImgSize = nImgHeight * nImgWidth;
	bReset = true;

	char pcMdlPath[512];
    sprintf(pcMdlPath, "./aam-Data/model_map_2d_data.bin");
    FILE* fp = fopen(pcMdlPath, "rb");
    fseek(fp, 0, SEEK_SET);
    fseek(fp, 0, SEEK_END);
    long len_bytes = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    unsigned char* pLoadData = (unsigned char*)malloc(len_bytes);
    fread(pLoadData, sizeof(unsigned char), len_bytes, fp);
    if (!EiInitialize(nImgWidth, nImgHeight, pLoadData, len_bytes) ||
            !EiInitialize_NewReso(nImgWidth, nImgHeight))
    {
		//fprintf(of,"initialize failed");fflush(of);
		//of<<"initialize failed"<<std::endl;
        printf("initialize failed\r\n");
		free(pLoadData);
		//if (pCapture)
			//cvReleaseCapture(&pCapture);
        exit(-1);
    }
    free(pLoadData);
	//fprintf(of, "initialize successed\n");fflush(of);
	//of<<"initialize successed"<<std::endl;
}

FeatureExtraction::~FeatureExtraction()
{
	//delete  face_cascade;
	EiDestroy_NewReso();
    EiDestroy();
}

void FeatureExtraction::GetFacePoint(cv::Mat& img)
{
	feature.clear();//fprintf(of, "1\n");fflush(of);//of<<1<<std::endl;
	IplImage *outImg = cvCreateImage(cvSize(nImgWidth, nImgHeight), 8, 3);//fprintf(of, "2\n");fflush(of);//of<<2<<std::endl;
	IplImage pImg = img;//fprintf(of, "3\n");fflush(of);//of<<3<<std::endl;
	IplImage* pGrayImg = cvCreateImage(cvGetSize(&pImg), 8, 1);//fprintf(of, "4\n");fflush(of);//of<<4<<std::endl;
	cvCvtColor(&pImg, pGrayImg, CV_BGR2GRAY);//fprintf(of, "5\n");fflush(of);//of<<5<<std::endl;
	//cvSaveImage("kinect-data/gray.jpg", pGrayImg);
	AAM_OUTPUT_STRUCT aam_ret = {0};//fprintf(of, "6\n");fflush(of);//of<<6<<std::endl;
	aam_ret.nstate = -1;//fprintf(of, "7\n");fflush(of);//of<<7<<std::endl;
	aam_ret.nstate = EiGetExpression((unsigned char*)pGrayImg->imageData, &aam_ret, bReset);//fprintf(of, "8\n");fflush(of);//of<<8<<std::endl;
	if (aam_ret.nstate == 1)
	{
		//of<<"aam_ret_nstate"<<std::endl;
		for (int i = 0; i < aam_ret.n2DNum; ++i)
		{
			cvDrawCircle(&pImg, cvPoint(int(aam_ret.pKeyPoint2DOut[i].x + 0.5f), int(aam_ret.pKeyPoint2DOut[i].y + 0.5f)), 1, cvScalar(0, 255, 0), -1);
		}
		for(int i=0; i<key_point.size(); i++)
		{
			int index = key_point[i];
			feature.push_back(cvPoint(cvRound(aam_ret.pKeyPoint2DOut[index].x + 0.5f),cvRound(aam_ret.pKeyPoint2DOut[index].y + 0.5f)));
			cvDrawCircle(&pImg, cvPoint(int(aam_ret.pKeyPoint2DOut[index].x + 0.5f), int(aam_ret.pKeyPoint2DOut[index].y + 0.5f)), 1, cvScalar(255, 0, 0), -1);
		}
		bReset = false;
	}
	else
		bReset = true;
	//fprintf(of, "9\n");fflush(of);
	//of<<9<<std::endl;
	cvReleaseImage(&pGrayImg);
	//fprintf(of, "getFacePoint finish\n");fflush(of);
    //of<<"getFacePoint finish"<<std::endl;
    
	/*feature.clear();
	cv::Mat pic_RGB= img, pic_GRAY;
	//pic_RGB ;
	vector<cv::Rect> face_vec;
	cvtColor(pic_RGB,pic_GRAY,CV_RGB2GRAY);
	equalizeHist(pic_GRAY,pic_GRAY);
	face_cascade->detectMultiScale(pic_GRAY,face_vec,1.1,2,0|CV_HAAR_SCALE_IMAGE,cv::Size(5,5));
	if(!face_vec.size())
	{
		printf("检测不到人脸\n");
		return;
	}
	else
	{
		printf("检测到%d张人脸\n",face_vec.size());
	}

	for(int i=0;i<(int)face_vec.size();i++)
	{
		cv::Point p1((int)face_vec[i].x,(int)face_vec[i].y);
		cv::Point p2((int)(face_vec[i].x+face_vec[i].width),(int)(face_vec[i].y+face_vec[i].height*1.1));
		cv::rectangle(img, p1, p2, cv::Scalar( 255, 0,0 ), 1, 8, 0 );
		//cvRectangle( img, p1,p2, cv::Scalar( 255, 0,0 ),1,8, 0 );
	}

	asm_shape shape, detshape;
	detshape.Resize(2);

	for(int i=(int)face_vec.size() - 1; i< (int)face_vec.size(); i++)
	{
		detshape[0].x=(float)face_vec[i].x;
		detshape[0].y=(float)face_vec[i].y;
		detshape[1].x=(float)(face_vec[i].x+face_vec[i].width);
		detshape[1].y=(float)(face_vec[i].y+face_vec[i].height);

		InitShapeFromDetBox(shape, detshape, fit_asm.GetMappingDetShape(), fit_asm.GetMeanFaceWidth());
		IplImage ipl_img = img;
		bool b=fit_asm.ASMSeqSearch(shape, &ipl_img , 0, true, 30);

		for(int j=0; j<68; j++)
		{
			if(j==27 || j==29 || j==34 || j==32 || j==67 || j==48 || j==54) //主要特征点|| j==66 
			{
				cv::circle( img,cvPoint(cvRound(shape[j].x),cvRound(shape[j].y)), 2, cvScalar(255,0,0),1, 8, 0 );
				feature.push_back(cvPoint(cvRound(shape[j].x),cvRound(shape[j].y)));
			}
			else cv::circle( img,cvPoint(cvRound(shape[j].x),cvRound(shape[j].y)), 2, cvScalar(0,255,0),1, 8, 0 );
			//cvCircle( img,cvPoint(cvRound(shape[j].x),cvRound(shape[j].y)), 2, cvScalar(255,0,0),1, 8, 0 );
			
		}
	}*/
}