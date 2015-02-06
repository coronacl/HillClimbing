#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <list>
#include <vector>
#include <algorithm>
#include <new>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

#include <windows.h> 

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#ifdef _DEBUG
//Debugモードの場合
#pragma comment(lib,"opencv_core249d.lib")
#pragma comment(lib,"opencv_imgproc249d.lib")
#pragma comment(lib,"opencv_highgui249d.lib")
//#pragma comment(lib,"opencv_objdetect249d.lib")
//#pragma comment(lib,"opencv_contrib249d.lib")
//#pragma comment(lib,"opencv_features2d249d.lib")
//#pragma comment(lib,"opencv_flann249d.lib")
//#pragma comment(lib,"opencv_gpu249d.lib")
//#pragma comment(lib,"opencv_haartraining_engined.lib")
//#pragma comment(lib,"opencv_legacy249d.lib")
//#pragma comment(lib,"opencv_ts249d.lib")
//#pragma comment(lib,"opencv_video249d.lib")
#else
//Releaseモードの場合
#pragma comment(lib,"opencv_core249.lib")
#pragma comment(lib,"opencv_imgproc249.lib")
#pragma comment(lib,"opencv_highgui249.lib")
//#pragma comment(lib,"opencv_objdetect249.lib")
//#pragma comment(lib,"opencv_contrib249.lib")
//#pragma comment(lib,"opencv_features2d249.lib")
//#pragma comment(lib,"opencv_flann249.lib")
//#pragma comment(lib,"opencv_gpu249.lib")
//#pragma comment(lib,"opencv_haartraining_engined.lib")
//#pragma comment(lib,"opencv_legacy249.lib")
//#pragma comment(lib,"opencv_ts249.lib")
//#pragma comment(lib,"opencv_video249.lib")
#endif

using namespace cv;

//マクロ定義
#define M_PI 3.141592653589793238463
#define Delete(p) delete[] p; p=NULL;
#define PData(p,x,y) p.data[ (y)*(p.step) + (x)*(p.elemSize())]
#define FAR_DIS (134217728)

//IplImage用のアクセス
#define PIXEL(p,x,y)    ((unsigned char *)p->imageData)[p->widthStep*(y)+(x)]
#define PIXEL_R(p,x,y)  ((unsigned char *)p->imageData)[p->widthStep*(y)+(x)*3+2]
#define PIXEL_G(p,x,y)  ((unsigned char *)p->imageData)[p->widthStep*(y)+(x)*3+1]
#define PIXEL_B(p,x,y)  ((unsigned char *)p->imageData)[p->widthStep*(y)+(x)*3]
#define IPL(p,x,y)    ((int *)p.imageData)[p.widthStep*(y)+(x)]

//Mat用のアクセス
#define PIX(p,x,y)	(uchar )p.data[y*p.cols + x]
#define PIXd(p,x,y)  (double)p.ptr<double>(y)[x]
#define PIXi(p,x,y)  (int)p.ptr<int>(y)[x]
#define PIXs(p,x,y)  p.ptr<short>(y)[x]
#define PIXci(p,x,y)  (int)p.ptr<uchar>(y)[x]
#define PIX_R(p,x,y)	(uchar)p.data[y * p.step + x * p.elemSize() + 2]
#define PIX_G(p,x,y)	(uchar)p.data[y * p.step + x * p.elemSize() + 1]
#define PIX_B(p,x,y)	(uchar)p.data[y * p.step + x * p.elemSize() + 0]

////基本処理用
//#include "Pixel2i.h"
#include "method.h"
//#include "Extract.h"
//
//////DICOM関連
//#include "Dicom/PixelStk.h"
//#include "Dicom/RawData.h"
//#include "Dicom/BinaryFileManage.h"
//#include "Dicom/DicomImage.h"
//#include "ManageDicom.h"
//#include "ManageCT.h"
//#include "ManageMR.h"
//#include "SystemValue.h"
////
//////処理用
////#include "equlization.h"
////#include "MatchingFast.h"
////#include "TransPicture.h"
//#include "TemplateMatching.h"
//
////GUI用
//#include "MyForm.h"

using namespace System;
using namespace System::IO;
using namespace System::Windows::Forms;
