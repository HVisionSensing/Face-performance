#ifndef FEATUREEXTRACTION_H_INCLUDED
#define FEATUREEXTRACTION_H_INCLUDED
#include <opencv\cv.h>
#include <opencv2/highgui/highgui.hpp>
#include "asmfitting.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "FaceTrackingDll.h"

class FeatureExtraction
{
public:
	std::vector< cv::Point2d > feature;
	//asmfitting fit_asm;
	std::vector<int> key_point;
	int nImgWidth;
    int nImgHeight;
    int nImgSize;
    bool bReset;

	//FILE* of;
public:
	FeatureExtraction(){};
	FeatureExtraction(int nImgWidth, int nImgHeight);
	~FeatureExtraction();
	void GetFacePoint(cv::Mat& img);
	//void PushPoint(cv::Point2d p);
};
#endif