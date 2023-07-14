// OpencvToC.cpp: 定义应用程序的入口点。
//

#include "OpencvToC.h"
#include "sod.h"
#include "core/imgpro.h"
#include "core/boundingRect.h"
#include "core/cvAlloc.h"

using namespace std;

int main()
{
	cout << "Start CMake." << endl;

	const char* zInput = "1.bmp";
	const char* zInput2 = "2.bmp";
	const char* zOut = "out1.png";
	cout << "Load Image." << endl;
	sod_img imgIn = sod_img_load_from_file(zInput, SOD_IMG_GRAYSCALE);
	sod_img imgIn2 = sod_img_load_from_file(zInput2, SOD_IMG_GRAYSCALE);
	if (imgIn.data == 0) {
		/* Invalid path, unsupported format, memory failure, etc. */
		puts("Cannot load input image..exiting");
		return 0;
	}
	CvMat mat1;
	mat1.ptr = imgIn.charData;
	mat1.Empty = 1;
	mat1.cols = imgIn.h;
	mat1.rows = imgIn.w;
	mat1.type = CV_8U;
	mat1.channels = imgIn.c;
	mat1.step = mat1.cols * mat1.channels * sizeof(unsigned char);
	CvSize cs;
	cs.height = imgIn.h/2;
	cs.width = imgIn.w/2;
	CvMat mat12 = createMat(cs.height,cs.width, CV_8U);
	
	mat12.ptr = imgIn.charData;
	mat12.Empty = 1;
	mat12.cols = imgIn.h;
	mat12.rows = imgIn.w;
	mat12.type = CV_8U;
	mat12.channels = imgIn.c;
	mat12.step = mat1.cols * mat1.channels * sizeof(unsigned char);
	
	

	//CvMat prev_gray = mat1;

	//CvPoint2D32f flow;

	//CvMat pts_x = createMat(mat1.rows,mat1.cols,CV_32FC1);
	//zerosMat(&pts_x);
	//CvMat pts_y = createMat(mat1.rows, mat1.cols, CV_32FC1);
	//zerosMat(&pts_y);
	//CvMat mag;
	//CvMat ang, ang2;
	//CvMat hsv = createMat(mat1.rows, mat1.cols, CV_8UC3);
	//zerosMat(&hsv);

	//UMatVector* mv;


	// 光流
	//UMatVector va = { 0 };
	CvMat va = createMat(1,1, CV_8U);
	//cvmat_vector_init(&va, 10);
	//calcOpticalFlowFarneback(mat1, mat12, &va, 0.5, 1, 15, 3, 5, 1.2, 0);
	cvcalcOpticalFlowFarneback(mat1, mat12, &va,0.25, 1, 15, 3, 5, 1.2, 0);

	CvMat pts_x = createMat(1,1,CV_32F);
	CvMat pts_y = createMat(1, 1, CV_32F);
	pts_x.fl[0] = 3;
	pts_y.fl[0] = 4;
	CvMat mag = createMat(1, 1, CV_32F);
	CvMat ang = createMat(1, 1, CV_32F);
	cvCartToPolar(&pts_x,&pts_y,&mag,&ang,true);
	CvMat ang2 = createMat(ang.rows, ang.cols, ang.type);
	cvMatDev(&ang,&ang2,2);
	cvMatMul(&mag, 100);

	CvMat magdst = createMat(mag.rows,mag.cols,CV_8U);
	
	cvconvertScaleAbs(&mag, &magdst, 1, 0);
	freeMat(&mag);
	CvMat binary = createMat(mat1.rows, mat1.cols, CV_8U);
	cvThreshold(&mat1,&binary,20,255,THRESH_BINARY);

	sod_img_blob_save_as_bmp(zOut, binary.ptr, binary.cols, binary.rows, binary.channels);


	// 开运算
	sod_img s = sod_make_image(binary.cols, binary.rows, binary.channels);
	
	int i;
	for ( i = 0; i < s.w*s.h*s.c; i++)
	{
		s.data[i] = (float)s.data[i] / 255.;
	}
	s = sod_erode_image(s,3);
	s = sod_dilate_image(s,3);



	//cartToPolar

	//convertScaleAbs
	//convertScaleAbs();
	//cvThreshold
	//findbox
	// 

	//sod_img_blob_save_as_bmp(zOut, mat122->data.ptr,mat122->cols,mat122->rows,mat122->channels);
	sod_img_save_as_png(imgIn, zOut);
	//cout << "Save Image Success." << endl;
	///* Cleanup */
	//sod_free_image(imgIn);

	return 0;
}
