#ifndef IMGPRO_H
#define IMGPRO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include "cvMat.h"
#include "tracking.h"
#include <stdbool.h>
#include "types_c.h"

#define PI 3.141592654
#  define CV_SSE 1
#  define CV_SSE2 1
    typedef enum InterpolationType {
        INTER_AREA,
        INTER_LINEAR
    } InterpolationType;
    typedef struct {
        float* data;
        int size;
    } AutoBuffer;
	//int resize(CvMat* src, CvMat** dst, CvSize dsize);
	int calcOpticalFlowFarneback(CvMat _prev0,CvMat _next0, UMatVector *_flow0,double pyr_scale, int levels, int winsize,
		int iterations, int poly_n, double poly_sigma, int flags);

    static void FarnebackPolyExp(const CvMat* src, CvMat* dst, int n, double sigma);
    void initAutoBuffer(AutoBuffer* buffer, int size);
    void releaseAutoBuffer(AutoBuffer* buffer);
    void FarnebackPrepareGaussian(int n, double sigma, float* g, float* xg, float* xxg,
        double* ig11, double* ig03, double* ig33, double* ig55);
    void FarnebackUpdateMatrices(const CvMat* _R0,const CvMat* _R1,const CvMat* _flow,const CvMat* matM,int _y0,int _y1);
    void FarnebackUpdateFlow_GaussianBlur(const CvMat* _R0,const CvMat* _R1,CvMat* _flow,CvMat* matM,int block_size,bool update_matrices);
    void FarnebackUpdateFlow_Blur(const CvMat* _R0, const CvMat* _R1, CvMat* _flow, CvMat* matM, int block_size, bool update_matrices);
    void cvcalcOpticalFlowFarneback(CvMat _prev0, CvMat _next0, CvMat* _flow0, double pyr_scale, int levels, int winsize, int iterations, int poly_n, double poly_sigma, int flags);
    int cvRound(double value);
    int cvFloor(double value);
    float* alignPtr(float* ptr, int n);
    unsigned char clipPixelValue(int value);
    void resize(const CvMat* src, CvMat* dst, CvSize dsize, double inv_scale_x, double inv_scale_y, InterpolationType interpolation);
    void resize_cv(const CvMat* src, CvMat* dst, CvSize dsize, double inv_scale_x, double inv_scale_y, InterpolationType interpolation);
    void gaussian(CvMat* src, CvMat* dst, int size, float sigma);
    void get_gau_kernel(float** kernel, int size, float sigma);

    float cvmean(const CvMat image);
    
    double  cvnorm(const CvMat arr1, const CvMat arr2,
        int norm_type);
    void convertScaleAbs(
        unsigned char* src,
        unsigned char* dst,
        double alpha,  //alpha 为乘数因子
        double beta,  //beta为偏移量
        int width, 
        int height );
    void cvconvertScaleAbs(CvMat* _src, CvMat* _dst, double alpha, double beta);
    void cvThreshold(CvMat* _src,CvMat* _dst,double thresh,double maxval,int type);

    void cvCartToPolar(const CvMat* xarr, const CvMat* yarr,
        CvMat* magarr, CvMat* anglearr,
        int angle_in_degrees);
   

   

#ifdef __cplusplus
}
#endif /* end of __cplusplus */

#endif /* IMGPRO_H */
