#include "imgpro.h"
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "cvAlloc.h"
#include <float.h>

int cvRound(double value) {
    return (int)(value + 0.5);
}

int cvFloor(double value)
{
    return (int)floor(value);
}

unsigned char clipPixelValue(int value) {
    if (value < 0) {
        return 0;
    }
    else if (value > 255) {
        return 255;
    }
    else {
        return (unsigned char)value;
    }
}

void resize(const CvMat* src, CvMat* dst, CvSize dsize, double inv_scale_x, double inv_scale_y, InterpolationType interpolation)
{
    int src_rows = src->rows;
    int src_cols = src->cols;
    int dst_rows = dsize.height;
    int dst_cols = dsize.width;

    int src_step = src->step;
    int dst_step = dst->step;

    int src_channels = src->channels;
    int dst_channels = dst->channels;

    unsigned char* src_data = src->ptr;
    unsigned char* dst_data = dst->ptr;

    int i, j, k;

    double scale_x =1/ inv_scale_x;
    double scale_y =  1/inv_scale_y;

    if (interpolation == INTER_AREA) {
        double scale = 1.0 / (scale_x * scale_y);
        double inv_scale = 1.0 / scale;

        for (i = 0; i < dst_rows; i++) {
            float fy = (float)((i + 0.5) * scale_y - 0.5);
            int sy = (int)floor(fy);
            fy -= sy;

            sy = sy < 0 ? 0 : (sy >= src_rows ? src_rows - 1 : sy);

            unsigned char* dst_row = dst_data + i * dst_step;
            unsigned char* src_row1 = src_data + sy * src_step;
            unsigned char* src_row2 = src_data + (sy + 1) * src_step;

            for (j = 0; j < dst_cols; j++) {
                float fx = (float)((j + 0.5) * scale_x - 0.5);
                int sx = (int)floor(fx);
                fx -= sx;

                sx = sx < 0 ? 0 : (sx >= src_cols ? src_cols - 1 : sx);

                unsigned char* dst_pixel = dst_row + j * dst_channels;
                unsigned char* src_pixel1 = src_row1 + sx * src_channels;
                unsigned char* src_pixel2 = src_row1 + (sx + 1) * src_channels;
                unsigned char* src_pixel3 = src_row2 + sx * src_channels;
                unsigned char* src_pixel4 = src_row2 + (sx + 1) * src_channels;

                for (k = 0; k < dst_channels; k++) {
                    float val = (1.0f - fy) * (1.0f - fx) * src_pixel1[k]
                        + fy * (1.0f - fx) * src_pixel3[k]
                        + (1.0f - fy) * fx * src_pixel2[k]
                        + fy * fx * src_pixel4[k];

                    dst_pixel[k] = clipPixelValue((int)round(val * inv_scale));
                }
            }
        }
    }
    else if (interpolation == INTER_LINEAR) {
        for (i = 0; i < dst_rows; i++) {
            float fy = (float)((i + 0.5) * scale_y - 0.5);
            int sy = (int)floor(fy);
            fy -= sy;

            sy = sy < 0 ? 0 : (sy >= src_rows ? src_rows - 1 : sy);

            unsigned char* dst_row = dst_data + i * dst_step;
            unsigned char* src_row1 = src_data + sy * src_step;
            unsigned char* src_row2 = src_data + (sy + 1) * src_step;

            for (j = 0; j < dst_cols; j++) {
                float fx = (float)((j + 0.5) * scale_x - 0.5);
                int sx = (int)floor(fx);
                fx -= sx;

                sx = sx < 0 ? 0 : (sx >= src_cols ? src_cols - 1 : sx);

                unsigned char* dst_pixel = dst_row + j * dst_channels;
                unsigned char* src_pixel1 = src_row1 + sx * src_channels;
                unsigned char* src_pixel2 = src_row1 + (sx + 1) * src_channels;
                unsigned char* src_pixel3 = src_row2 + sx * src_channels;
                unsigned char* src_pixel4 = src_row2 + (sx + 1) * src_channels;

                for (k = 0; k < dst_channels; k++) {
                    float val = (1.0f - fy) * (1.0f - fx) * src_pixel1[k]
                        + fy * (1.0f - fx) * src_pixel3[k]
                        + (1.0f - fy) * fx * src_pixel2[k]
                        + fy * fx * src_pixel4[k];

                    dst_pixel[k] = clipPixelValue((int)round(val));
                }
            }
        }
    }
}

void resize_cv(const CvMat* src, CvMat* dst, CvSize dsize, double inv_scale_x, double inv_scale_y, InterpolationType interpolation)
{
    CvSize ssize= cvGetMatSize(src);
    if (cvGetArea(ssize) <= 0 ) 
    {
        printf("resize_cv size erro");
        return;
    }
    if ((inv_scale_x <= 0 || inv_scale_y <= 0) || cvGetArea(dsize))
    {
        printf("resize_cv area erro");
        return;
    }

    if (cvGetArea(dsize) == 0) 
    {
        dsize.height = ssize.height * inv_scale_y;
        dsize.width = ssize.width * inv_scale_x;
        
    }
    else 
    {
        inv_scale_x = (double)dsize.width / ssize.width;
        inv_scale_y = (double)dsize.height / ssize.height;
    }



}

int calcOpticalFlowFarneback(CvMat _prev0, CvMat _next0, UMatVector *_flow0, double pyr_scale, int levels, int winsize, int iterations, int poly_n, double poly_sigma, int flags)
{

    if ((5 != poly_n) && (7 != poly_n))
        return false;
    if (_next0.cols != _prev0.cols || _next0.rows != _prev0.rows)
        return false;

    int typePrev = _prev0.type;
    int typeNext = _next0.type;
    if ((1 != CV_MAT_CN(typePrev)) || (1 != CV_MAT_CN(typeNext)))
        return false;




    if (_flow0->size > 0)
    {

    }
    else
    {
        cvmat_vector_push_back(_flow0, createMat(_prev0.rows, _prev0.cols, _prev0.type));
        cvmat_vector_push_back(_flow0, createMat(_prev0.rows, _prev0.cols, _prev0.type));
    }

    CvMat prev0 = _prev0, next0 = _next0;

    const int min_size = 32;
    const CvMat* img[2] = { &prev0, &next0 };

    int i, k;
    double scale;
    CvMat prevFlow = { 0 }, flow, fimg;


    CvMat flow0 = _flow0->data[0];

    for (k = 0, scale = 1; k < levels; k++)
    {
        scale *= pyr_scale;
        if (prev0.cols * scale < min_size || prev0.rows * scale < min_size)
            break;
    }

    levels = k;

    for (k = levels; k >= 0; k--)
    {
        for (i = 0, scale = 1; i < k; i++)
            scale *= pyr_scale;

        double sigma = (1. / scale - 1) * 0.5;
        int smooth_sz = cvRound(sigma * 5) | 1;
        smooth_sz = max(smooth_sz, 3);


        int width = cvRound(prev0.cols * scale);
        int height = cvRound(prev0.rows * scale);

        if (k > 0)
            flow = createMat(height, width, CV_32FC2);

        else
            flow = flow0;

        if (!prevFlow.Empty)
        {
            if (flags & OPTFLOW_USE_INITIAL_FLOW)
            {
                CvSize size;
                size.height = height;
                size.width = width;

                resize(&flow0, &flow, size, 0, 0, INTER_AREA);

                // resize(&flow0, &flow, size);
                int s = 0;
                for (s = 0; s < flow.cols * flow.rows; i++)
                {
                    flow.fl[s] *= scale;
                }
            }
            else
            {
                zerosMat(&flow);
            }
        }
        else
        {
            CvSize size;
            size.height = height;
            size.width = width;
            resize(&prevFlow, &flow, size, 0, 0, INTER_LINEAR);
            // resize(&prevFlow, &flow, size);
            int s = 0;
            for (s = 0; s < flow.cols * flow.rows; i++)
            {
                flow.fl[s] *= 1. / pyr_scale;;
            }
        }

        CvMat R[2], I[2], I1[2], M;

        for (i = 0; i < 2; i++)
        {
            /*fimg = createMat(img[i]->rows, img[i]->cols, CV_32F);
            int st = 0;
            for (st = 0; st < img[i]->rows * img[i]->cols; st++)
            {
                fimg.fl[st] = (float)img[i]->ptr[st] / 255.0;
            }*/
            CvSize size;
            size.height = height;
            size.width = width;
            //            gaussian(&fimg, &fimg, smooth_sz, sigma);
            I[i] = createMat(height, width, CV_8U);
            CvSize ssize = cvGetMatSize(img[i]);

            double inv_scale_x = (double)size.width / ssize.width;
            double inv_scale_y = (double)size.height / ssize.height;
            resize(img[i], &I, size, inv_scale_x, inv_scale_y, INTER_LINEAR);
           
            R[i] = createMat(height, width, CV_32FC(5));
            zerosMat(&R[i]);
            I1[i] = createMat(height, width, CV_32F);
            int st = 0;
            for (st = 0; st < I[i].rows * I[i].cols; st++)
            {
                I1[i].fl[st] = (float)I[i].ptr[st] / 255.0;
            }
            FarnebackPolyExp(&I1[i], &R[i], poly_n, poly_sigma);
           
        }
        freeMat(&I1[1]);
        freeMat(&I1[0]);
        
        freeMat(&I[1]);
        freeMat(&I[0]);
    
        M = createMat(R[1].rows, R[1].cols, CV_32F);
        zerosMat(&M);
        FarnebackUpdateMatrices(&R[0], &R[1], &flow, &M, 0, flow.rows);

        for (i = 0; i < iterations; i++)
        {
            if (flags & OPTFLOW_FARNEBACK_GAUSSIAN)
                FarnebackUpdateFlow_GaussianBlur(&R[0], &R[1], &flow, &M, winsize, i < iterations - 1);
            else
                //FarnebackUpdateFlow_Blur(&R[0], &R[1], &flow, &M, winsize, i < iterations - 1);
                printf("");
        }
        //freeMat(&M);
        /*freeMat(&R[1]);
        freeMat(&R[0]);*/
        prevFlow = flow;

    }

}

void FarnebackPolyExp(const CvMat* src, CvMat* dst, int n, double sigma)
{
    int k, x, y;

    int width = src->cols;
    int height = src->rows;

    AutoBuffer kbuf;
    AutoBuffer _row;
    int kbufSize = n * 6 + 3;
    int rowSize = (width + n * 2) * 3;
    initAutoBuffer(&kbuf, kbufSize);
    initAutoBuffer(&_row, rowSize);

    float* g = kbuf.data + n;
    float* xg = g + n * 2 + 1;
    float* xxg = xg + n * 2 + 1;
    float* row = _row.data + n * 3;
    double ig11, ig03, ig33, ig55;

    FarnebackPrepareGaussian(n, sigma, g, xg, xxg, &ig11, &ig03, &ig33, &ig55);
    createMat_s(dst, height, width, CV_32FC(5));
    


    for (y = 0; y < height; y++)
    {
        float g0 = g[0], g1, g2;
        float* srow0 = (float*)(src->fl + src->step * y), * srow1 = 0;
        float* drow = (float*)(dst->fl + dst->step * y);

        // vertical part of convolution
        for (x = 0; x < width; x++)
        {
            row[x * 3] = srow0[x] * g0;
            row[x * 3 + 1] = row[x * 3 + 2] = 0.f;
        }

        for (k = 1; k <= n; k++)
        {
            g0 = g[k]; g1 = xg[k]; g2 = xxg[k];
            srow0 = (float*)(src->fl + src->step * max(y - k, 0));
            srow1 = (float*)(src->fl + src->step * min(y + k, height - 1));

            for (x = 0; x < width; x++)
            {
                float p = srow0[x] + srow1[x];
                float t0 = row[x * 3] + g0 * p;
                float t1 = row[x * 3 + 1] + g1 * (srow1[x] - srow0[x]);
                float t2 = row[x * 3 + 2] + g2 * p;

                row[x * 3] = t0;
                row[x * 3 + 1] = t1;
                row[x * 3 + 2] = t2;
            }
        }

        // horizontal part of convolution
        for (x = 0; x < n * 3; x++)
        {
            row[-1 - x] = row[2 - x];
            row[width * 3 + x] = row[width * 3 + x - 3];
        }

        for (x = 0; x < width; x++)
        {
            g0 = g[0];
            // r1 ~ 1, r2 ~ x, r3 ~ y, r4 ~ x^2, r5 ~ y^2, r6 ~ xy
            double b1 = row[x * 3] * g0, b2 = 0, b3 = row[x * 3 + 1] * g0,
                b4 = 0, b5 = row[x * 3 + 2] * g0, b6 = 0;

            for (k = 1; k <= n; k++)
            {
                double tg = row[(x + k) * 3] + row[(x - k) * 3];
                g0 = g[k];
                b1 += tg * g0;
                b4 += tg * xxg[k];
                b2 += (row[(x + k) * 3] - row[(x - k) * 3]) * xg[k];
                b3 += (row[(x + k) * 3 + 1] + row[(x - k) * 3 + 1]) * g0;
                b6 += (row[(x + k) * 3 + 1] - row[(x - k) * 3 + 1]) * xg[k];
                b5 += (row[(x + k) * 3 + 2] + row[(x - k) * 3 + 2]) * g0;
            }

            // do not store r1
            drow[x * 5 + 1] = (float)(b2 * ig11);
            drow[x * 5] = (float)(b3 * ig11);
            drow[x * 5 + 3] = (float)(b1 * ig03 + b4 * ig33);
            drow[x * 5 + 2] = (float)(b1 * ig03 + b5 * ig33);
            drow[x * 5 + 4] = (float)(b6 * ig55);
        }
    }

    row -= n * 3;

    releaseAutoBuffer(&kbuf);
    releaseAutoBuffer(&_row);
}

void initAutoBuffer(AutoBuffer* buffer, int size)
{
    buffer->data = (float*)malloc(sizeof(float) * size);
    buffer->size = size;
}

void releaseAutoBuffer(AutoBuffer* buffer)
{
    if (buffer->data)
    {
        free(buffer->data);
        buffer->data = NULL;
        buffer->size = 0;
    }
}

void FarnebackPrepareGaussian(int n, double sigma, float* g, float* xg, float* xxg, double* ig11, double* ig03, double* ig33, double* ig55)
{
    if (sigma < FLT_EPSILON)
        sigma = n * 0.3;

    double s = 0.;
    int x = 0,y=0;
    for ( x = -n; x <= n; x++)
    {
        g[x] = (float)exp(-x * x / (2 * sigma * sigma));
        s += g[x];
    }

    s = 1. / s;
    for ( x = -n; x <= n; x++)
    {
        g[x] = (float)(g[x] * s);
        xg[x] = (float)(x * g[x]);
        xxg[x] = (float)(x * x * g[x]);
    }

    CvMat G = createMat(6, 6, CV_64F);
    zerosMat(&G);

    for (y = -n; y <= n; y++) 
    {
        for (x = -n; x <= n; x++) 
        {
            cvAppendMatAdd(&G,0,0,g[y]*g[x]);
            cvAppendMatAdd(&G, 1, 1, g[y] * g[x] * x * x);
            cvAppendMatAdd(&G, 3, 3, g[y] * g[x] * x * x * x * x);
            cvAppendMatAdd(&G, 5, 5, g[y] * g[x] * x * x * y * y);
        }
    }

    cvAppendMatAdd(&G, 2, 2, cvGetMatValue(&G, 1, 1));
    cvAppendMatAdd(&G, 0, 3, cvGetMatValue(&G, 1, 1));
    cvAppendMatAdd(&G, 0, 4, cvGetMatValue(&G, 1, 1));
    cvAppendMatAdd(&G, 4, 0, cvGetMatValue(&G, 1, 1));
    cvAppendMatAdd(&G, 4, 0, cvGetMatValue(&G, 1, 1));
    cvAppendMatAdd(&G, 4, 4, cvGetMatValue(&G, 3, 3));
    cvAppendMatAdd(&G, 3, 4, cvGetMatValue(&G, 5, 5));
    cvAppendMatAdd(&G, 4, 3, cvGetMatValue(&G, 5, 5));

    // invG:
   // [ x        e  e    ]
   // [    y             ]
   // [       y          ]
   // [ e        z       ]
   // [ e           z    ]
   // [                u ]



}

void FarnebackUpdateMatrices(const CvMat* _R0, const CvMat* _R1, const CvMat* _flow, CvMat* matM, int _y0, int _y1)
{
    const int BORDER = 5;
    static const float border[5] = { 0.14f, 0.14f, 0.4472f, 0.4472f, 0.4472f };

    int x, y, width = _flow->cols, height = _flow->rows;
    const float* R1 = (float*)_R1->fl;
    size_t step1 = _R1->step / sizeof(R1[0]);

    createMat_s(matM,height,width,CV_32FC(5));
    //matM.create(height, width, CV_32FC(5));

    for (y = _y0; y < _y1; y++)
    {
        const float* flow = (float*)(_flow->fl + y * _flow->step);
        const float* R0 = (float*)(_R0->fl + y * _R0->step);
        float* M = (float*)(matM->fl + y * matM->step);

        for (x = 0; x < width; x++)
        {
            float dx = flow[x * 2], dy = flow[x * 2 + 1];
            float fx = x + dx, fy = y + dy;

#if 1
            int x1 = cvFloor(fx), y1 = cvFloor(fy);
            const float* ptr = R1 + y1 * step1 + x1 * 5;
            float r2, r3, r4, r5, r6;

            fx -= x1; fy -= y1;

            if ((unsigned)x1 < (unsigned)(width - 1) &&
                (unsigned)y1 < (unsigned)(height - 1))
            {
                float a00 = (1.f - fx) * (1.f - fy), a01 = fx * (1.f - fy),
                    a10 = (1.f - fx) * fy, a11 = fx * fy;

                r2 = a00 * ptr[0] + a01 * ptr[5] + a10 * ptr[step1] + a11 * ptr[step1 + 5];
                r3 = a00 * ptr[1] + a01 * ptr[6] + a10 * ptr[step1 + 1] + a11 * ptr[step1 + 6];
                r4 = a00 * ptr[2] + a01 * ptr[7] + a10 * ptr[step1 + 2] + a11 * ptr[step1 + 7];
                r5 = a00 * ptr[3] + a01 * ptr[8] + a10 * ptr[step1 + 3] + a11 * ptr[step1 + 8];
                r6 = a00 * ptr[4] + a01 * ptr[9] + a10 * ptr[step1 + 4] + a11 * ptr[step1 + 9];

                r4 = (R0[x * 5 + 2] + r4) * 0.5f;
                r5 = (R0[x * 5 + 3] + r5) * 0.5f;
                r6 = (R0[x * 5 + 4] + r6) * 0.25f;
            }
#else
            int x1 = cvRound(fx), y1 = cvRound(fy);
            const float* ptr = R1 + y1 * step1 + x1 * 5;
            float r2, r3, r4, r5, r6;

            if ((unsigned)x1 < (unsigned)width &&
                (unsigned)y1 < (unsigned)height)
            {
                r2 = ptr[0];
                r3 = ptr[1];
                r4 = (R0[x * 5 + 2] + ptr[2]) * 0.5f;
                r5 = (R0[x * 5 + 3] + ptr[3]) * 0.5f;
                r6 = (R0[x * 5 + 4] + ptr[4]) * 0.25f;
            }
#endif
            else
            {
                r2 = r3 = 0.f;
                r4 = R0[x * 5 + 2];
                r5 = R0[x * 5 + 3];
                r6 = R0[x * 5 + 4] * 0.5f;
            }

            r2 = (R0[x * 5] - r2) * 0.5f;
            r3 = (R0[x * 5 + 1] - r3) * 0.5f;

            r2 += r4 * dy + r6 * dx;
            r3 += r6 * dy + r5 * dx;

            if ((unsigned)(x - BORDER) >= (unsigned)(width - BORDER * 2) ||
                (unsigned)(y - BORDER) >= (unsigned)(height - BORDER * 2))
            {
                float scale = (x < BORDER ? border[x] : 1.f) *
                    (x >= width - BORDER ? border[width - x - 1] : 1.f) *
                    (y < BORDER ? border[y] : 1.f) *
                    (y >= height - BORDER ? border[height - y - 1] : 1.f);

                r2 *= scale; r3 *= scale; r4 *= scale;
                r5 *= scale; r6 *= scale;
            }

            M[x * 5] = r4 * r4 + r6 * r6; // G(1,1)
            M[x * 5 + 1] = (r4 + r5) * r6;  // G(1,2)=G(2,1)
            M[x * 5 + 2] = r5 * r5 + r6 * r6; // G(2,2)
            M[x * 5 + 3] = r4 * r2 + r6 * r3; // h(1)
            M[x * 5 + 4] = r6 * r2 + r5 * r3; // h(2)
        }
    }
}
float* alignPtr(float* ptr,int n) 
{
    (float*)(((size_t)ptr + n - 1) & -n);
}
void FarnebackUpdateFlow_GaussianBlur(const CvMat* _R0, const CvMat* _R1, CvMat* _flow, CvMat* matM, int block_size, bool update_matrices)
{
    int x, y, i, width = _flow->cols, height = _flow->rows;
    int m = block_size / 2;
    int y0 = 0, y1;
    int min_update_stripe = max((1 << 10) / width, block_size);
    double sigma = m * 0.3, s = 1;
    AutoBuffer _vsum;
    initAutoBuffer(&_vsum, (width + m * 2 + 2) * 5 + 16);
    AutoBuffer _hsum;
    initAutoBuffer(&_hsum, width * 5 + 16);
    AutoBuffer _srow;
    initAutoBuffer(&_srow, m * 2 + 1);
    AutoBuffer _kernel;
    initAutoBuffer(&_kernel, (m + 1) * 5 + 16);

    float* vsum = alignPtr((float*)_vsum.data + (m + 1) * 5, 16);
    float* hsum = alignPtr((float*)_hsum.data, 16);
    float* kernel = (float*)_kernel.data;
    const float** srow = (const float**)&_srow.data[0];
    kernel[0] = (float)s;

    for (i = 0; i <= m; i++)
    {
        float t = (float)exp(-i * i / (2 * sigma * sigma));
        kernel[i] = t;
        s += t * 2;
    }

    s = 1. / s;
    for (i = 0; i <= m; i++)
    {
        kernel[i] = (float)(kernel[i] * s);
    }

    for (y = 0; y < height; y++)
    {
        double g11, g12, g22, h1, h2;
        float* flow = (float*)(_flow->fl + _flow->step * y);
        // vertical blur
        for (i = 0; i <= m; i++)
        {
            srow[m - i] = (const float*)(matM->fl + matM->step * max(y - i, 0));
            srow[m + i] = (const float*)(matM->fl + matM->step * min(y + i, height - 1));
        }

        x = 0;

        for (; x < width * 5; x++)
        {
            float s0 = srow[m][x] * kernel[0];
            for (i = 1; i <= m; i++)
                s0 += (srow[m + i][x] + srow[m - i][x]) * kernel[i];
            vsum[x] = s0;
        }

        // update borders
        for (x = 0; x < m * 5; x++)
        {
            vsum[-1 - x] = vsum[4 - x];
            vsum[width * 5 + x] = vsum[width * 5 + x - 5];
        }

        // horizontal blur
        x = 0;

        for (; x < width * 5; x++)
        {
            float sum = vsum[x] * kernel[0];
            for (i = 1; i <= m; i++)
                sum += kernel[i] * (vsum[x - i * 5] + vsum[x + i * 5]);
            hsum[x] = sum;
        }
        for (x = 0; x < width; x++)
        {
            g11 = hsum[x * 5];
            g12 = hsum[x * 5 + 1];
            g22 = hsum[x * 5 + 2];
            h1 = hsum[x * 5 + 3];
            h2 = hsum[x * 5 + 4];

            double idet = 1. / (g11 * g22 - g12 * g12 + 1e-3);

            flow[x * 2] = (float)((g11 * h2 - g12 * h1) * idet);
            flow[x * 2 + 1] = (float)((g22 * h1 - g12 * h2) * idet);
        }

        y1 = y == height - 1 ? height : y - block_size;
        if (update_matrices && (y1 == height || y1 >= y0 + min_update_stripe))
        {
            FarnebackUpdateMatrices(_R0, _R1, _flow, matM, y0, y1);
            y0 = y1;
        }
    }

    releaseAutoBuffer(&_vsum);
    releaseAutoBuffer(&_hsum);
    releaseAutoBuffer(&_srow);
    releaseAutoBuffer(&_kernel);

}

void FarnebackUpdateFlow_Blur(const CvMat* _R0, const CvMat* _R1, CvMat* _flow, CvMat* matM, int block_size, bool update_matrices)
{
    int x, y, width = _flow->cols, height = _flow->rows;
    int m = block_size / 2;
    int y0 = 0, y1;
    int min_update_stripe = max((1 << 10) / width, block_size);
    double scale = 1. / (block_size * block_size);

    AutoBuffer _vsum = { 0 };
    int _vsumSize = (width + m * 2 + 2) * 5;
    initAutoBuffer(&_vsum, _vsumSize);

    double* vsum = (double*)_vsum.data + (m + 1) * 5;

    // init vsum
    const float* srow0 = (const float*)matM->fl;
    for (x = 0; x < width * 5; x++)
        vsum[x] = srow0[x] * (m + 2);

    for (y = 1; y < m; y++)
    {
        srow0 = (float*)(matM->fl + matM->step * min(y, height - 1));
        for (x = 0; x < width * 5; x++)
            vsum[x] += srow0[x];
    }

    // compute blur(G)*flow=blur(h)
    for (y = 0; y < height; y++)
    {
        double g11, g12, g22, h1, h2;
        float* flow = (float*)(_flow->fl + _flow->step * y);

        srow0 = (const float*)(matM->fl + matM->step * max(y - m - 1, 0));
        const float* srow1 = (const float*)(matM->fl + matM->step * min(y + m, height - 1));

        // vertical blur
        for (x = 0; x < width * 5; x++)
            vsum[x] += srow1[x] - srow0[x];

        // update borders
        for (x = 0; x < (m + 1) * 5; x++)
        {
            vsum[-1 - x] = vsum[4 - x];
            vsum[width * 5 + x] = vsum[width * 5 + x - 5];
        }

        // init g** and h*
        g11 = vsum[0] * (m + 2);
        g12 = vsum[1] * (m + 2);
        g22 = vsum[2] * (m + 2);
        h1 = vsum[3] * (m + 2);
        h2 = vsum[4] * (m + 2);

        for (x = 1; x < m; x++)
        {
            g11 += vsum[x * 5];
            g12 += vsum[x * 5 + 1];
            g22 += vsum[x * 5 + 2];
            h1 += vsum[x * 5 + 3];
            h2 += vsum[x * 5 + 4];
        }

        // horizontal blur
        for (x = 0; x < width; x++)
        {
            g11 += vsum[(x + m) * 5] - vsum[(x - m) * 5 - 5];
            g12 += vsum[(x + m) * 5 + 1] - vsum[(x - m) * 5 - 4];
            g22 += vsum[(x + m) * 5 + 2] - vsum[(x - m) * 5 - 3];
            h1 += vsum[(x + m) * 5 + 3] - vsum[(x - m) * 5 - 2];
            h2 += vsum[(x + m) * 5 + 4] - vsum[(x - m) * 5 - 1];

            double g11_ = g11 * scale;
            double g12_ = g12 * scale;
            double g22_ = g22 * scale;
            double h1_ = h1 * scale;
            double h2_ = h2 * scale;

            double idet = 1. / (g11_ * g22_ - g12_ * g12_ + 1e-3);

            flow[x * 2] = (float)((g11_ * h2_ - g12_ * h1_) * idet);
            flow[x * 2 + 1] = (float)((g22_ * h1_ - g12_ * h2_) * idet);
        }

        y1 = y == height - 1 ? height : y - block_size;
        if (update_matrices && (y1 == height || y1 >= y0 + min_update_stripe))
        {
            FarnebackUpdateMatrices(_R0, _R1, _flow, matM, y0, y1);
            y0 = y1;
        }
    }

    releaseAutoBuffer(&_vsum);
}

void cvcalcOpticalFlowFarneback(CvMat _prev0, CvMat _next0, CvMat* _flow0, double pyr_scale, int levels, int winsize, int iterations, int poly_n, double poly_sigma, int flags)
{
    CvMat prev0 = getCvMat(&_prev0);
    CvMat next0 = getCvMat(&_next0);

    const int min_size = 32;
    const CvMat* img[2] = { &prev0, &next0 };
    int i, k;
    double scale;
    CvMat prevFlow = { 0 }, flow = { 0 }, fimg = {0};
    if (_flow0->Empty) 
    {
        freeMat(_flow0);
        createMat_s(_flow0, prev0.rows, prev0.cols, prev0.type);
    }
    else 
    {
        createMat_s(_flow0, prev0.rows, prev0.cols, prev0.type);
    }
    
    CvMat flow0 = getCvMat(_flow0);
    for (k = 0, scale = 1; k < levels; k++)
    {
        scale *= pyr_scale;
        if (prev0.cols * scale < min_size || prev0.rows * scale < min_size)
            break;
    }

    levels = k;
    for (k = levels; k >= 0; k--) 
    {
        for (i = 0, scale = 1; i < k; i++)
            scale *= pyr_scale;

        double sigma = (1. / scale - 1) * 0.5;
        int smooth_sz = cvRound(sigma * 5) | 1;
        smooth_sz = max(smooth_sz, 3);

        int width = cvRound(prev0.cols * scale);
        int height = cvRound(prev0.rows * scale);

        if (k > 0) 
        {
            if (flow.Empty) 
            {
                freeMat(&flow);
            }
            createMat_s(&flow, height, width, CV_32FC2);
        }
            
        else 
        {
            cvcopyTo(&flow0,&flow);
        }
           

        if (!prevFlow.Empty)
        {
            if (flags & OPTFLOW_USE_INITIAL_FLOW)
            {
                CvSize cs;
                cs.height = height;
                cs.width = width;
                resize(&flow0, &flow, cs, 0.25, 0.25, INTER_AREA);

                cvMatMul(&flow, scale);
            }
            else
                zerosMat(&flow);
        }
        else 
        {
            CvSize cs;
            cs.height = height;
            cs.width = width;
            CvSize ssize = cvGetMatSize(img[i]);
            double inv_scale_x = (double)cs.width / ssize.width;
            double inv_scale_y = (double)cs.height / ssize.height;
            resize(&prevFlow,&flow,cs, inv_scale_x, inv_scale_y,INTER_LINEAR);
            double v =1. / pyr_scale;
            cvMatMul(&flow,v);
        }

        CvMat R[2] = { 0 }, I = { 0 }, M = {0};
        for(i=0;i<2;i++)
        {
            cvconvertTo(img[i], &fimg, CV_32F);
            gaussian(&fimg, &fimg, 3, sigma);
           // unsigned char* test = calloc(fimg.cols * fimg.rows * fimg.channels, sizeof(unsigned char));
            //int ti = 0;
            /*for ( ti = 0; ti < fimg.cols * fimg.rows * fimg.channels; ti++)
            {
                test[ti] = (unsigned char)(fimg.fl[ti] * 255);

            }
            sod_img_blob_save_as_bmp("gaussian.bmp", test, fimg.cols , fimg.rows , fimg.channels);
            free(test);*/
            CvSize cs;
            cs.height = height;
            cs.width = width;
            CvMat tmpT = { 0 };

            cvconvertTo(&fimg, &tmpT, CV_8U);
            createMat_s(&I,height,width,CV_8U);
            CvSize ssize = cvGetMatSize(img[i]);
            double inv_scale_x = (double)cs.width / ssize.width;
            double inv_scale_y = (double)cs.height / ssize.height;
            resize(&tmpT,&I,cs, inv_scale_x, inv_scale_y,INTER_LINEAR);
            cvconvertTo(&I, &tmpT, CV_32F);
           

          /*  sod_img_blob_save_as_bmp("resize.bmp", I.ptr, width ,height , I.channels);
            free(test);*/

            FarnebackPolyExp(&tmpT,&R[i],poly_n,poly_sigma);
            freeMat(&I);
            freeMat(&tmpT);

        }
        FarnebackUpdateMatrices(&R[0],&R[1],&flow,&M,0,flow.rows);
        
    }
}



void get_gau_kernel(float** kernel, int size, float sigma)
{
    if (size <= 0 || sigma == 0)
        return;

    int x, y;
    int m = size / 2;
    float sum = 0;

    //get kernel
    for (y = 0; y < size; y++)
    {
        for (x = 0; x < size; x++)
        {
            kernel[y][x] = (1 / (2 * PI * sigma * sigma)) * exp(-((x - m) * (x - m) + (y - m) * (y - m)) / (2 * sigma * sigma));
            sum += kernel[y][x];
        }
    }

    //normal
    for (y = 0; y < size; y++)
    {
        for (x = 0; x < size; x++)
        {
            kernel[y][x] /= sum;
        }
    }
}

float cvmean(const CvMat image)
{
   
    int i, j, k;
    int pixelSize = image.rows * image.cols;
    float sum = 0.0f;
    if (image.channels == 1)
    {

        // 计算通道均值
        for (k = 0; k < image.channels; k++) {
            for (i = 0; i < image.rows; i++) {
                for (j = 0; j < image.cols; j++) {
                    int index = k * pixelSize + i * image.cols + j;
                    sum += (float)image.ptr[index];
                }
            }

        }
        sum /= pixelSize;
    }
    return sum;
}

double cvnorm(const CvMat arr1, const CvMat arr2, int norm_type)
{
    int rows = arr1.rows;
    int cols = arr1.cols;

    float* data1 = arr1.fl;
    float* data2 = arr2.fl;
    double result = 0.0;

    if (norm_type == NORM_L2)
    {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {

                float diff = data1[i * cols + j] - data2[i * cols + j];
                result += diff * diff;
            }
        }

        result = sqrtf(result);

        return result;
    }
    return result;
}

void convertScaleAbs(unsigned char* src, unsigned char* dst, double alpha, double beta, int width, int height)
{
    int index = 0;
    for (int h = 0; h < height; h++)
    {
        for (int w = 0; w < width; w++)
        {
            if (src[index] * alpha + beta > -255 && src[index] * alpha + beta < 255)
            {
                src[index] = abs(src[index]);
            }
            if (alpha + beta >= 255 || src[index] * alpha + beta <= -255)
            {
                src[index] = 255;
            }
            index++;
        }
    }
    dst = src;
}

void cvconvertScaleAbs(CvMat* _src, CvMat* _dst, double alpha, double beta)
{
    int rows = _src->rows;
    int cols = _src->cols ;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
           
            float srcValue = (float)_src->fl[i * cols + j];
            
           // float srcValue = ((float*)_src->i)[i * cols + j];

            float scaledValue = alpha * srcValue + beta;
            float absValue = (float)fabsf(scaledValue);
            int value = (int)(absValue + 0.5);
            if (value > 255) value = 255;
            ((unsigned char*)_dst->ptr)[i * cols + j] = (unsigned char)value;
        }
    }
}



void cvThreshold(CvMat* _src, CvMat* _dst, double thresh, double maxval, int type)
{
    if (_src->type != 0) // 只支持 CV_8UC1 类型的图像
    {
        printf("Unsupported image type\n");
        return;
    }

    int rows = _src->rows;
    int cols = _src->step;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            unsigned char pixel = _src->ptr[i * cols + j];

            if (type == THRESH_BINARY)
            {
                if (pixel > thresh)
                    _dst->ptr[i * cols + j] = maxval;
                else
                    _dst->ptr[i * cols + j] = 0;
            }
        }
    }
}

void cvCartToPolar(const CvMat* xarr, const CvMat* yarr, CvMat* magarr, CvMat* anglearr, int angle_in_degrees)
{
    int rows = xarr->rows;
    int cols = xarr->cols;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            float x = xarr->fl[i * cols + j];
            float y = yarr->fl[i * cols + j];

            float magnitude = sqrt(x * x + y * y);
            float angle = atan2(y, x);

            if (angle_in_degrees)
                angle *= 180.0f / PI;

            magarr->fl[i * cols + j] = magnitude;
            anglearr->fl[i * cols + j] = angle;
        }
    }
}



void gaussian(CvMat* src, CvMat* dst, int size, float sigma) 
{
    if (src->cols == 0 || src->rows == 0)
        return;

    int y, x;
    int i, j;
    int m = size / 2;
    float value;

    float** kernel = (float**)CvMalloc(size * sizeof(float*));
    for (i = 0; i < size; i++)
        kernel[i] = (float*)CvMalloc(size * sizeof(float));

    get_gau_kernel(kernel, size, sigma);

    float* kernel_vec = (float*)CvMalloc(size * size * sizeof(float));

   
    int k = 0;
    for (j = 0; j < size; j++)
    {
        for (i = 0; i < size; i++)
        {
            kernel_vec[k++] = kernel[j][i];
        }
    }
    if (dst->type == CV_8S || CV_8U) 
    {
        unsigned char* src_ptr = src->ptr + m * src->cols;
        unsigned char* dst_ptr = dst->ptr + m * dst->rows;

        for (y = m; y < src->rows - m; y++)
        {
            for (x = m; x < src->cols - m; x++)
            {

                value = 0;
                k = 0;
                for (j = -m; j < m; j++)
                {
                    for (i = -m; i < m; i++)
                    {
                        unsigned char temp = src_ptr[(y + j) * src->cols + (x + i)];
                        float temp1 = kernel_vec[k++];
                        value += temp * temp1;
                    }
                }

                dst_ptr[x] = (unsigned char)(value);
            }

            dst_ptr += dst->rows;
        }

        free(kernel_vec);
        for (i = 0; i < size; i++)
            free(kernel[i]);
        free(kernel);
    }
    else if (dst->type == CV_32F) 
    {
        float* src_ptr = src->fl + m * src->cols;
        float* dst_ptr = dst->fl + m * dst->rows;

        for (y = m; y < src->rows - m; y++)
        {
            for (x = m; x < src->cols - m; x++)
            {

                value = 0;
                k = 0;
                for (j = -m; j < m; j++)
                {
                    for (i = -m; i < m; i++)
                    {
                        float temp = src_ptr[(y + j) * src->cols + (x + i)];
                        float temp1 = kernel_vec[k++];
                        value += temp * temp1;
                    }
                }

                dst_ptr[x] = (float)(value);
            }

            dst_ptr += dst->rows;
        }

        free(kernel_vec);
        for (i = 0; i < size; i++)
            free(kernel[i]);
        free(kernel);
    }
   
}



