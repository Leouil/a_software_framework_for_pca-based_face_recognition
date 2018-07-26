void CLBP::CreatLBP(IplImage *src,IplImage *dst)
{
    int iTemp[8] = {0};
    CvScalar s;
    
    IplImage *ImgTemp = cvCreateImage(cvGetSize(src), IPL_DEPTH_8U, 1);
    uchar *data = (uchar*)src->imageData;
    int iStep = src-> widthStep;
    
    for (int i=1;i<src->height-1;i++)
    for(int j=1;j<src->width-1;j++){
        
        int sum=0;
        if(data[(i-1)*iStep+j-1]>data[i*iStep+j])
        iTemp[0]=1;
        else
        iTemp[0]=0;
        if(data[i*iStep+(j-1)]>data[i*iStep+j])
        iTemp[1]=1;
        else
        iTemp[1]=0;
        if(data[(i+1)*iStep+(j-1)]>data[i*iStep+j])
        iTemp[2]=1;
        else
        iTemp[2]=0;
        if (data[(i+1)*iStep+j]>data[i*iStep+j])
        iTemp[3]=1;
        else
        iTemp[3]=0;
        if (data[(i+1)*iStep+(j+1)]>data[i*iStep+j])
        iTemp[4]=1;
        else
        iTemp[4]=0;
        if(data[i*iStep+(j+1)]>data[i*iStep+j])
        iTemp[5]=1;
        else
        iTemp[5]=0;
        if(data[(i-1)*iStep+(j+1)]>data[i*iStep+j])
        iTemp[6]=1;
        else
        iTemp[6]=0;
        if(data[(i-1)*iStep+j]>data[i*iStep+j])
        iTemp[7]=1;
        else
        iTemp[7]=0;
        s.val[0] = (iTemp[0]*1+iTemp[1]*2+iTemp[2]*4+iTemp[3]*8+iTemp[4]*16+iTemp[5]*32+iTemp[6]*64+iTemp[7]*128);
        
        cvSet2D(dst,i,j,s);
        
    }
    
}    // end of CLBP::CreatLBP()
