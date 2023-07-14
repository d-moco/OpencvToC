#include "boundingRect.h"
#include "cvAlloc.h"

int rect_Filter(int xmin, int ymin, int xmax, int ymax)
{
    if (xmin < xmax && ymin < ymax)
    {
        int width = xmax - xmin;
        int height = ymax - ymin;
        if (width < 10 || height < 10 || width > 200 || height > 200)
        {
            return 0;
        }

        return 1;
    }

    return 0;
}

box_rect* findbox(unsigned char* mask, int width, int height, int areaMax, int boxMax, int* ret)
{
    short* area = NULL;
    area = (short*)CvMalloc(width * height * sizeof(*area));
    memset(area, 0, width * height * sizeof(*area));

    //  ���������
    int classCount = 0;
    int pixindex = -1;
    int h = 0, w = 0;
    for (h = 0; h < height; h++)
    {
        for (w = 0; w < width; w++)
        {
            pixindex++;
            if (mask[pixindex] == 255)//ֻ�жϰ�ɫ����
            {
                if (h == 0)//��һ��
                {
                    if (pixindex % width == 0)//��һ��
                    {
                        if (mask[pixindex] == 255)
                        {
                            classCount++;
                            area[pixindex] = classCount;//�µ�������
                        }
                    }
                    else//���ǵ�һ��
                    {
                        if (mask[pixindex - 1] == 255)
                        {
                            area[pixindex] = area[pixindex - 1];//����������ص����
                        }
                        else
                        {
                            classCount++;
                            area[pixindex] = classCount;//�µ�������
                        }

                    }

                }
                else//���ǵ�һ��
                {
                    if (pixindex % width == 0)//��һ��
                    {
                        if (mask[pixindex - width] == 255)
                        {
                            area[pixindex] = area[pixindex - width];//�����������ص����
                        }
                        else
                        {
                            classCount++;
                            area[pixindex] = classCount;//�µ�������
                        }

                    }
                    else//���ǵ�һ��
                    {
                        if (mask[pixindex - width] == 255)
                        {
                            area[pixindex] = area[pixindex - width];//�����������ص����
                        }
                        else if (mask[pixindex - 1] == 255)
                        {
                            area[pixindex] = area[pixindex - 1];//����������ص����
                        }
                        else
                        {
                            classCount++;
                            area[pixindex] = classCount;//�µ�������
                        }

                    }

                }
            }
        }
    }

    //�ں�����������
    while (1)
    {
        pixindex = -1;
        int change = 0;
        for (h = 0; h < height; h++)
        {
            for (w = 0; w < width; w++)
            {
                pixindex++;
                if (area[pixindex] != 0)//ֻ�жϱ�ǹ���
                {
                    if (h == 0)//��һ��
                    {
                        if (w != 0)//���ǵ�һ��
                        {
                            if (area[pixindex - 1] != 0)
                            {
                                if (area[pixindex] > area[pixindex - 1])
                                {
                                    area[pixindex] = area[pixindex - 1];
                                    change = 1;
                                }
                                else if (area[pixindex] < area[pixindex - 1])
                                {
                                    area[pixindex - 1] = area[pixindex];
                                    change = 1;
                                }

                            }
                        }

                    }
                    else//���ǵ�һ��
                    {
                        if (w == 0)//��һ��
                        {
                            if (area[pixindex - width] != 0)
                            {
                                if (area[pixindex] > area[pixindex - width])
                                {
                                    area[pixindex] = area[pixindex - width];
                                    change = 1;
                                }
                                else if (area[pixindex] < area[pixindex - width])
                                {
                                    area[pixindex - width] = area[pixindex];
                                    change = 1;
                                }

                            }

                        }
                        else//���ǵ�һ��
                        {
                            if (area[pixindex - width] != 0)//��
                            {
                                if (area[pixindex] > area[pixindex - width])
                                {
                                    area[pixindex] = area[pixindex - width];
                                    change = 1;
                                }
                                else if (area[pixindex] < area[pixindex - width])
                                {
                                    area[pixindex - width] = area[pixindex];
                                    change = 1;
                                }

                            }

                            if (area[pixindex - 1] != 0)//��
                            {
                                if (area[pixindex] > area[pixindex - 1])
                                {
                                    area[pixindex] = area[pixindex - 1];
                                    change = 1;
                                }
                                else if (area[pixindex] < area[pixindex - 1])
                                {
                                    area[pixindex - 1] = area[pixindex];
                                    change = 1;
                                }

                            }

                        }

                    }
                }

            }
        }

        if (change == 0)
        {
            break;
        }

    }

    //  Ϊÿ�����򴴽�����
    rect** mybox = (rect**)CvMalloc(areaMax * sizeof(rect*));
    int classVal = 0;
    int num = 0, i = 0;
    for (; i < width * height; i++)
    {
        if (num >= areaMax)
        {
            break;
        }

        if (area[i] > classVal && area[i] != 0)
        {
            classVal = area[i];
            mybox[num] = (rect*)CvMalloc(sizeof(rect));
            mybox[num]->xmin = 999;
            mybox[num]->xmax = 0;
            mybox[num]->ymin = 999;
            mybox[num]->ymax = 0;
            mybox[num]->classNum = area[i];
            num++;
        }
    }

    //  ������ֵ����Сֵ
    pixindex = -1;
    for (h = 0; h < height; h++)
    {
        for (w = 0; w < width; w++)
        {
            pixindex++;
            int n = 0;
            for (; n < num; n++)
            {
                if (mybox[n]->classNum == area[pixindex])
                {
                    if (w < mybox[n]->xmin)
                    {
                        mybox[n]->xmin = w;
                    }
                    if (w > mybox[n]->xmax)
                    {
                        mybox[n]->xmax = w;
                    }
                    if (h < mybox[n]->ymin)
                    {
                        mybox[n]->ymin = h;
                    }
                    if (h > mybox[n]->ymax)
                    {
                        mybox[n]->ymax = h;
                    }
                }

            }
        }
    }

    //�ҵ����кϸ�ľ���
    rect* boxes = 0;
    boxes = (rect*)CvMalloc(sizeof(rect) * boxMax);
    memset(boxes, 0, sizeof(rect) * boxMax);

    int boxCount = 0;
    for (i = 0; i < num; i++)
    {
        int ret = rect_Filter(mybox[i]->xmin, mybox[i]->ymin, mybox[i]->xmax, mybox[i]->ymax);
        if (ret == 1)
        {
            boxes[boxCount].xmin = mybox[i]->xmin;
            boxes[boxCount].ymin = mybox[i]->ymin;
            boxes[boxCount].xmax = mybox[i]->xmax;
            boxes[boxCount].ymax = mybox[i]->ymax;
            boxes[boxCount].classNum = mybox[i]->classNum;

            boxCount++;

            if (boxCount >= boxMax)
            {
                break;
            }

        }

    }

    box_rect* boxRect = 0;
    boxRect = (box_rect*)CvMalloc(sizeof(box_rect) * boxCount);
    memset(boxRect, 0, sizeof(box_rect) * boxCount);
    for (i = 0; i < boxCount; i++)
    {
        int width = boxes[i].xmax - boxes[i].xmin;
        int height = boxes[i].ymax - boxes[i].ymin;

        boxRect[i].xmin = boxes[i].xmin;
        boxRect[i].ymin = boxes[i].ymin;
        boxRect[i].width = width;
        boxRect[i].height = height;
    }

    *ret = boxCount;

    CvFree(boxes);
    CvFree(area);
    for (i = 0; i < num; i++)
    {
        if (mybox[i] != NULL)
        {
            CvFree(mybox[i]);
        }
    }
    if (mybox != NULL)
    {
        CvFree(mybox);
    }

    return boxRect;
}