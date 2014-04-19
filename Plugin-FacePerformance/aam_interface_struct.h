#ifndef AAM_INTERFACE_STRUCT_HEADER_201301051458
#define AAM_INTERFACE_STRUCT_HEADER_201301051458
typedef struct AAM_POINT_3D
{
	float x;
	float y;
	float z;
}AAM_POINT_3D;



typedef struct AAM_POINT_2D
{
	float x;
	float y;
}AAM_POINT_2D;


struct SFFIPara
{
	float fAngleY;		// out plane left/right	(rad)
	float fAngleX;		// out plane up/down (rad)
	float fAngleZ;		// in plane (rad)
	float fScale;			
	float fTranslateX;		
	float fTranslateY;		
	float fShapePara[15];	
	float fAnimPara[11];	
} ;

typedef struct AAM_OUTPUT_STRUCT
{
	AAM_POINT_2D *pKeyPoint2DOut;
	AAM_POINT_2D *pKeyPoint3DProjTo2D;
	SFFIPara *pPara;

	int nViewType;
	int nstate;
	int n3DPNum;
	int n2DNum;
} AAM_OUTPUT_STRUCT;

//坐标轴参数
typedef struct AAM_POINT_AXIS
{
	AAM_POINT_2D Center;
	AAM_POINT_2D X;
	AAM_POINT_2D Y;
	AAM_POINT_2D Z;
} AAM_POINT_AXIS;

//角度

typedef struct AAM_POINT_RAD
{
	float xr;
	float yr;
	float zr;
} AAM_POINT_RAD;
#endif
