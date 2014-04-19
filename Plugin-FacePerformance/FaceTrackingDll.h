// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the FACETRACKINGDLL_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// FACETRACKINGDLL_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifndef _AAM_FACE_TRACKING_DLL_OUTPUT
#define _AAM_FACE_TRACKING_DLL_OUTPUT

#define FACETRACKINGDLL_EXPORTS

#ifdef FACETRACKINGDLL_EXPORTS
#ifdef WIN32
#define FACETRACKINGDLL_API __declspec(dllexport)
#else
#define FACETRACKINGDLL_API
#endif
#else
#define FACETRACKINGDLL_API
#endif

#include "aam_interface_struct.h"
#define AVATAR_AAM
// Initialize
// nW and nH is the width and height of the input image
FACETRACKINGDLL_API bool EiInitialize(int nW,
                                      int nH,
                                      unsigned char* pLoadData,
                                      int len_bytes);
FACETRACKINGDLL_API bool EiInitialize_NewReso(int nW, int nH);
// Destroy
FACETRACKINGDLL_API bool EiDestroy(void);
FACETRACKINGDLL_API bool EiDestroy_NewReso(void);

// Check whether or not the tracking has been initialized
FACETRACKINGDLL_API bool EiIsEnable(void);

// Tracking and get parameters of the model
// return the tracking state
// pBuffer  - (i) image data, GRAY
// pPara    - (o) model parameters
// bNormalPose  - (i) indicates whether or not reinitialize the tracking, leave it as false always
FACETRACKINGDLL_API int EiGetExpression(unsigned char* pBuffer,
                                        AAM_OUTPUT_STRUCT *pAAM_Output,
                                        bool bNormalPose);
#endif