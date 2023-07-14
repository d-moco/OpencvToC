#include "cvMat.h"
#include "cvAlloc.h"
#include <string.h>
#include "..\sod.h"

CvMat createMat(int rows, int cols, int type)
{
    CvMat _mat;
    _mat.rows = rows;
    _mat.cols = cols;
    _mat.type = type;
    _mat.Empty = 1;// 初始为空
    _mat.channels = CV_MAT_CN(type);
    
    int size = rows * cols * CV_MAT_CN(type);

    switch (CV_MAT_DEPTH(type)) {
    case CV_8U:
    case CV_8S:
        _mat.ptr = (unsigned char*)CvCalloc(size, sizeof(unsigned char));
        _mat.step = cols * _mat.channels;
        break;
    case CV_16U:
    case CV_16S:
        _mat.s = (short*)CvCalloc(size, sizeof(short));
        _mat.step = cols * _mat.channels;
        break;
    case CV_32S:
        _mat.i = (int*)CvCalloc(size, sizeof(int));
        _mat.step = cols * _mat.channels;
        break;
    case CV_32F:
        _mat.fl = (float*)CvCalloc(size, sizeof(float));
        _mat.step = cols * _mat.channels;
        break;
    case CV_64F:
        _mat.db = (double*)CvCalloc(size, sizeof(double));
        _mat.step = cols * _mat.channels;
        break;
    case CV_16F:
        _mat.fl = (float*)CvCalloc(size, sizeof(float));
        _mat.step = cols * _mat.channels;
        break;
    }
    return _mat;
}

void createMat_s(CvMat* _mat, int rows, int cols, int type)
{
    if (_mat->Empty) 
    {
        freeMat(_mat);
    }
    _mat->rows = rows;
    _mat->cols = cols;
    _mat->type = type;
    _mat->Empty = 1;
    _mat->channels = CV_MAT_CN(type);
    _mat->step = cols * _mat->channels ;
    int size = rows * cols * CV_MAT_CN(type);

    switch (CV_MAT_DEPTH(type)) {
    case CV_8U:
    case CV_8S:
        _mat->ptr = (unsigned char*)CvCalloc(size, sizeof(unsigned char));
        break;
    case CV_16U:
    case CV_16S:
        _mat->s = (short*)CvCalloc(size, sizeof(short));
        break;
    case CV_32S:
        _mat->i = (int*)CvCalloc(size, sizeof(int));
        break;
    case CV_32F:
        _mat->fl = (float*)CvCalloc(size, sizeof(float));
        break;
    case CV_64F:
        _mat->db = (double*)CvCalloc(size, sizeof(double));
        break;
    case CV_16F:
        _mat->fl = (float*)CvCalloc(size, sizeof(float));
        break;
    }
    
}

void zerosMat(CvMat* data)
{
    int i = 0;
    int type = data->type;
    for (size_t i = 0; i < data->cols*data->rows; i++)
    {
        switch (CV_MAT_DEPTH(type)) {
        case CV_8U:
        case CV_8S:
            
            data->ptr[i] = (unsigned char)0;
            break;
        case CV_16U:
        case CV_16S:
            data->s[i] = (short)0;
            break;
        case CV_32S:
            data->i[i] = (int)0;
            break;
        case CV_32F:
            data->fl[i] = (float)0.0f;
            break;
        case CV_64F:
            data->db[i] = (double)0.0;
            break;
        case CV_16F:
            data->fl[i] = (float)0.0;
            break;
        }
    }
    
}

void freeMat(CvMat* _mat)
{
    if (_mat == (CvMat*)NULL) {
        return;
    }
    _mat->Empty = 0;
    switch (CV_MAT_DEPTH(_mat->type)) {
    case CV_8U:
    case CV_8S:
       
       CvFree(_mat->ptr);
       _mat->ptr = NULL;
        break;
    case CV_16U:
    case CV_16S:
        CvFree(_mat->s);
        break;
    case CV_32S:
        CvFree(_mat->i);
        break;
    case CV_32F:
        CvFree(_mat->fl);
        _mat->fl = 0;
        break;
    case CV_64F:
        CvFree(_mat->db);
        break;
    case CV_16F:
        CvFree(_mat->fl);
        break;
    }
   // CvFree(_mat);
}

CvMat getCvMat(CvMat* mat)
{
    CvMat cm = createMat(mat->rows,mat->cols,mat->type);
    cm.channels = mat->channels;
    cm.cols = mat->cols;
    cm.step = mat->step;
    cm.refcount = cm.refcount;
    switch (CV_MAT_DEPTH(mat->type)) {
    case CV_8U:
    case CV_8S:
        cm.ptr = mat->ptr;
        break;
    case CV_16U:
    case CV_16S:
        cm.s = mat->s;
        break;
    case CV_32S:
        cm.i = mat->i;
        break;
    case CV_32F:
        cm.fl = mat->fl;
        break;
    case CV_64F:
        cm.db = mat->db;
        break;
    case CV_16F:
        cm.fl = mat->fl;
        break;
    }
    

    return cm;
}



void cvmat_vector_init(UMatVector* vec, size_t capacity)
{
    int i = 0;
    vec->size = 0;
    vec->capacity = capacity;
    vec->data = (CvMat*)CvCalloc(capacity , sizeof(CvMat));   
}

void cvmat_vector_push_back(UMatVector* vec, CvMat value)
{
    if (vec->size >= vec->capacity) {
        vec->capacity *= 2;


        vec->data = (CvMat*)CvRealloc(vec->data,vec->capacity * sizeof(CvMat));
    }
    memcpy(&(vec->data[vec->size++]), &value, sizeof(CvMat));
   
}

void cvmat_vector_free(UMatVector* vec)
{
    CvFree(vec->data);
}

void saveImage(CvMat* data, const char* zOut)
{
    if (data->type == CV_8U)
    {
        sod_img_blob_save_as_bmp(zOut, data->ptr, data->cols, data->rows, data->channels);
    }
    else if (data->type == CV_32F) 
    {
        unsigned char* data1 = 0;
        int i, k;
        data1 = calloc(data->cols * data->rows * data->channels, sizeof(unsigned char));
        for (k = 0; k < data->channels; ++k) {
            for (i = 0; i < data->cols * data->rows; ++i) {
                data1[i * data->channels + k] = (unsigned char)(255 * data->fl[i + k * data->cols * data->rows]);
            }
        }
        sod_img_blob_save_as_bmp(zOut, data1, data->cols, data->rows, data->channels);
        free(data1);
    }
   

}

CvSize cvGetMatSize(const CvMat* mat)
{
    CvSize cs;
    cs.height = mat->rows;
    cs.width = mat->cols;
    return cs;
}

double cvGetArea(CvSize cs)
{
    return cs.width*cs.height;
}

void cvMatDev(const CvMat* src, CvMat* dst,float dev)
{
    int i;
    int type = dst->type;
    for ( i = 0; i < src->rows*src->cols*src->channels; i++)
    {
        switch (CV_MAT_DEPTH(type)) {
        case CV_8U:
        case CV_8S:
            src->ptr[i] /= dev;
            break;
        case CV_16U:
        case CV_16S:
            src->s[i] /= dev;
            break;
        case CV_32S:
            src->i[i] /= dev;
            break;
        case CV_32F:
            src->fl[i] /= dev;
            break;
        case CV_64F:
            src->db[i] /= dev;
            break;
        case CV_16F:
            src->fl[i] /= dev;
            break;
        }
        
    }
}


void cvMatMul(const CvMat* src, double mul)
{
    int i;
    int type = src->type;
    for (i = 0; i < src->rows * src->cols * src->channels; i++)
    {
        switch (CV_MAT_DEPTH(type)) {
        case CV_8U:
        case CV_8S:
            src->ptr[i] *= (unsigned char)mul;
            break;
        case CV_16U:
        case CV_16S:
            src->s[i] *= ( short)mul;
            break;
        case CV_32S:
            src->i[i] *= (int)mul;
            break;
        case CV_32F:
            src->fl[i] *= (float)mul;
            break;
        case CV_64F:
            src->db[i] *= (double)mul;
            break;
        case CV_16F:
            src->fl[i] *= (float)mul;
            break;
        }
    }
}


void cvcopyTo(const CvMat* src, CvMat* dst) 
{
    if (dst->Empty) 
    {
        freeMat(dst);
    }
    createMat_s(dst,src->rows,src->cols,src->type);
    dst->channels = src->channels;
    
    int i = 0;
    for ( i = 0; i < src->rows*src->cols*src->channels; i++)
    {
        switch (CV_MAT_DEPTH(src->type)) {
        case CV_8U:
        case CV_8S:
            dst->ptr[i] = src->ptr[i];
            dst->step = dst->cols * dst->channels;
            break;
        case CV_16U:
        case CV_16S:
            dst->s[i] = src->s[i];
            dst->step = dst->cols * dst->channels;
            break;
        case CV_32S:
            dst->i[i] = src->i[i];
            dst->step = dst->cols * dst->channels;
            break;
        case CV_32F:
            dst->fl[i] = src->fl[i];
            dst->step = dst->cols * dst->channels;
            break;
        case CV_64F:
            dst->db[i] = src->db[i];
            dst->step = dst->cols * dst->channels;
            break;
        case CV_16F:
            dst->fl[i] = src->fl[i];
            dst->step = dst->cols * dst->channels;
            break;
        }
    }
}
// 目前仅支持32 8u
void cvconvertTo(const CvMat* src,CvMat* dst, int type)
{
    if (dst->Empty)
    {
        freeMat(dst);
    }
    createMat_s(dst, src->rows, src->cols, type);
    dst->channels = src->channels;

    int i = 0;
    for ( i = 0; i < src->channels*src->cols*src->rows; i++)
    {
       
        switch (CV_MAT_DEPTH(src->type)) {
        case CV_8U:
        case CV_8S:
            dst->step = dst->cols * dst->channels;
            if (dst->type == src->type) 
            {
                dst->ptr[i] = src->ptr[i];
               
            }
            else if (dst->type == CV_32F) 
            {
                dst->fl[i] = (float)src->ptr[i]/255.;
            }
            
            break;
       
        case CV_32F:
            dst->step = dst->cols * dst->channels;
            if (dst->type == src->type)
            {
                dst->fl[i] = src->fl[i];

            }
            else if (dst->type == CV_8S || dst->type == CV_8U)
            {
                dst->ptr[i] = (unsigned char)(src->fl[i] * 255.);
            }
            break;
      
        }
    }

    
}


void cvSetIdentity(CvMat* src)
{

    
    int i;
    int type = src->type;

    if (type == CV_32FC1) 
    {
        int minSize = (src->rows < src->cols) ? src->rows : src->cols;
        for (i = 0; i < minSize; i++) {
            src->fl[i * src->cols + i] = 1.0;
        }

    }
    else if (type == CV_64FC1) 
    {
        int minSize = (src->rows < src->cols) ? src->rows : src->cols;
        for (int i = 0; i < minSize; i++) {
            src->db[i * src->cols + i] = 1.0;
        }

    }
    else
    {
        zerosMat(src);
        
    }
}

double cvInvert(CvMat* src,CvMat* dst, int method)
{
    
    

    int result = false;
    int type = src->type;

    if (!(type == CV_32F || type == CV_64F)) 
    {
        printf_s("只能是32F 和 64F");
       
        zerosMat(dst);
        return result;
    }

    size_t esz = CV_ELEM_SIZE(type);
    int m = src->rows, n = src->cols;

    CvMat _dst = createMat(n,n,type);
    CvMat dst = getCvMat(&_dst);

    if (n <= 3) 
    {
        unsigned char* srcdata = src->db;
        unsigned char* dstdata = src->db;
        size_t srcstep = src->step;
        size_t dststep = dst->step;

        if (n == 2) 
        {
            if (type == CV_32FC1)
            {
                double d = det2(Sf);
                if (d != 0.)
                {
                    result = true;
                    d = 1. / d;


                    
                    double t0, t1;
                    t0 = Sf(0, 0) * d;
                    t1 = Sf(1, 1) * d;
                    Df(1, 1) = (float)t0;
                    Df(0, 0) = (float)t1;
                    t0 = -Sf(0, 1) * d;
                    t1 = -Sf(1, 0) * d;
                    Df(0, 1) = (float)t0;
                    Df(1, 0) = (float)t1;
                    

                }
                
            }
            else
            {
                double d = det2(Sd);
                if (d != 0.)
                {
                    result = true;
                    d = 1. / d;


                    double t0, t1;
                    t0 = Sd(0, 0) * d;
                    t1 = Sd(1, 1) * d;
                    Dd(1, 1) = t0;
                    Dd(0, 0) = t1;
                    t0 = -Sd(0, 1) * d;
                    t1 = -Sd(1, 0) * d;
                    Dd(0, 1) = t0;
                    Dd(1, 0) = t1;

                }
            }
        }
        else if (n == 3)
        {
            if (type == CV_32FC1)
            {
                double d = det3(Sf);

                if (d != 0.)
                {
                    double t[12];

                    result = true;
                    d = 1. / d;
                    t[0] = (((double)Sf(1, 1) * Sf(2, 2) - (double)Sf(1, 2) * Sf(2, 1)) * d);
                    t[1] = (((double)Sf(0, 2) * Sf(2, 1) - (double)Sf(0, 1) * Sf(2, 2)) * d);
                    t[2] = (((double)Sf(0, 1) * Sf(1, 2) - (double)Sf(0, 2) * Sf(1, 1)) * d);

                    t[3] = (((double)Sf(1, 2) * Sf(2, 0) - (double)Sf(1, 0) * Sf(2, 2)) * d);
                    t[4] = (((double)Sf(0, 0) * Sf(2, 2) - (double)Sf(0, 2) * Sf(2, 0)) * d);
                    t[5] = (((double)Sf(0, 2) * Sf(1, 0) - (double)Sf(0, 0) * Sf(1, 2)) * d);

                    t[6] = (((double)Sf(1, 0) * Sf(2, 1) - (double)Sf(1, 1) * Sf(2, 0)) * d);
                    t[7] = (((double)Sf(0, 1) * Sf(2, 0) - (double)Sf(0, 0) * Sf(2, 1)) * d);
                    t[8] = (((double)Sf(0, 0) * Sf(1, 1) - (double)Sf(0, 1) * Sf(1, 0)) * d);

                    Df(0, 0) = (float)t[0]; Df(0, 1) = (float)t[1]; Df(0, 2) = (float)t[2];
                    Df(1, 0) = (float)t[3]; Df(1, 1) = (float)t[4]; Df(1, 2) = (float)t[5];
                    Df(2, 0) = (float)t[6]; Df(2, 1) = (float)t[7]; Df(2, 2) = (float)t[8];
                }
            }
            else
            {
                double d = det3(Sd);
                if (d != 0.)
                {
                    result = true;
                    d = 1. / d;
                    double t[9];

                    t[0] = (Sd(1, 1) * Sd(2, 2) - Sd(1, 2) * Sd(2, 1)) * d;
                    t[1] = (Sd(0, 2) * Sd(2, 1) - Sd(0, 1) * Sd(2, 2)) * d;
                    t[2] = (Sd(0, 1) * Sd(1, 2) - Sd(0, 2) * Sd(1, 1)) * d;

                    t[3] = (Sd(1, 2) * Sd(2, 0) - Sd(1, 0) * Sd(2, 2)) * d;
                    t[4] = (Sd(0, 0) * Sd(2, 2) - Sd(0, 2) * Sd(2, 0)) * d;
                    t[5] = (Sd(0, 2) * Sd(1, 0) - Sd(0, 0) * Sd(1, 2)) * d;

                    t[6] = (Sd(1, 0) * Sd(2, 1) - Sd(1, 1) * Sd(2, 0)) * d;
                    t[7] = (Sd(0, 1) * Sd(2, 0) - Sd(0, 0) * Sd(2, 1)) * d;
                    t[8] = (Sd(0, 0) * Sd(1, 1) - Sd(0, 1) * Sd(1, 0)) * d;

                    Dd(0, 0) = t[0]; Dd(0, 1) = t[1]; Dd(0, 2) = t[2];
                    Dd(1, 0) = t[3]; Dd(1, 1) = t[4]; Dd(1, 2) = t[5];
                    Dd(2, 0) = t[6]; Dd(2, 1) = t[7]; Dd(2, 2) = t[8];
                }
            }
        }
        else
        {
           

            if (type == CV_32FC1)
            {
                double d = Sf(0, 0);
                if (d != 0.)
                {
                    result = true;
                    Df(0, 0) = (float)(1. / d);
                }
            }
            else
            {
                double d = Sd(0, 0);
                if (d != 0.)
                {
                    result = true;
                    Dd(0, 0) = 1. / d;
                }
            }
        }
        if (!result) 
        {
            zerosMat(dst);
            return result;
        }
           
        return result;
    }

    int elem_size = CV_ELEM_SIZE(type);
    CvMat src1 = {0};
    cvcopyTo(src,&src1);
    
    cvSetIdentity(&dst);

  /*  if (method == DECOMP_LU && type == CV_32F)
        result = LU((float*)src1.fl, src1.step, n, (float*)dst->fl, dst->step, n) != 0;
    else if (method == DECOMP_LU && type == CV_64F)
        result = LU((double*)src1.db, src1.step, n, (double*)dst->db, dst->step, n) != 0;
    else if (method == DECOMP_CHOLESKY && type == CV_32F)
        result = Cholesky((float*)src1.fl, src1.step, n, (float*)dst->fl, dst->step, n);
    else
        result = Cholesky((double*)src1.d, src1.step, n, (double*)dst.data, dst.step, n);*/



}



float cvGetMatValue(CvMat* src, int row, int col)
{
    int type = src->type;
    int width = src->cols;
    int height = src->rows;
    float value = 0;
    switch (CV_MAT_DEPTH(type)) {
    case CV_8U:
    case CV_8S:
        value = (float)src->ptr[row * width + col];
        break;
    case CV_16U:
    case CV_16S:
        value = (float)src->s[row * width + col] ;
        break;
    case CV_32S:
        value = (float)src->i[row * width + col] ;
        break;
    case CV_32F:
        value = (float)src->fl[row * width + col] ;
        break;
    case CV_64F:
        value = (float)src->db[row * width + col] ;
        break;
    case CV_16F:
        value = (float)src->fl[row * width + col] ;
        break;
    }
    return value;
}



void cvAppendMatAdd(CvMat* src, int row, int col,float value)
{
    int type = src->type;
    int width = src->cols;
    int height = src->rows;

    switch (CV_MAT_DEPTH(type)) {
    case CV_8U:
    case CV_8S:
        src->ptr[row*width + col] += (unsigned char)value;
        break;
    case CV_16U:
    case CV_16S:
        src->s[row * width + col] += (short)value;
        break;
    case CV_32S:
        src->i[row * width + col] += (int)value;
        break;
    case CV_32F:
        src->fl[row * width + col] += (float)value;
        break;
    case CV_64F:
        src->db[row * width + col] += (double)value;
        break;
    case CV_16F:
        src->fl[row * width + col] += (float)value;
        break;
    }
    
}
