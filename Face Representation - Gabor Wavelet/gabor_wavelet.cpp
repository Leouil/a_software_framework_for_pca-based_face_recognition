Mat WDT( const Mat &_src, const string _wname, const int _level )const
{
    int reValue = THID_ERR_NONE;
    Mat src = Mat_<float>(_src);
    Mat dst = Mat::zeros( src.rows, src.cols, src.type() );
    int N = src.rows;
    int D = src.cols;
    
    Mat lowFilter;
    Mat highFilter;
    wavelet( _wname, lowFilter, highFilter );
    
    int t=1;
    int row = N;
    int col = D;
    
    while( t<=_level )
    {
        for( int i=0; i<row; i++ )
        {
            Mat oneRow = Mat::zeros( 1,col, src.type() );
            for ( int j=0; j<col; j++ )
            {
                oneRow.at<float>(0,j) = src.at<float>(i,j);
            }
            oneRow = waveletDecompose( oneRow, lowFilter, highFilter );
            for ( int j=0; j<col; j++ )
            {
                dst.at<float>(i,j) = oneRow.at<float>(0,j);
            }
        }
        
#if 0
        //normalize( dst, dst, 0, 255, NORM_MINMAX );
        IplImage dstImg1 = IplImage(dst);
        cvSaveImage( "dst.jpg", &dstImg1 );
#endif
        for ( int j=0; j<col; j++ )
        {
            Mat oneCol = Mat::zeros( row, 1, src.type() );
            for ( int i=0; i<row; i++ )
            {
                oneCol.at<float>(i,0) = dst.at<float>(i,j);
            }
            oneCol = ( waveletDecompose( oneCol.t(), lowFilter, highFilter ) ).t();
            
            for ( int i=0; i<row; i++ )
            {
                dst.at<float>(i,j) = oneCol.at<float>(i,0);
            }
        }
        
#if 0
        //normalize( dst, dst, 0, 255, NORM_MINMAX );
        IplImage dstImg2 = IplImage(dst);
        cvSaveImage( "dst.jpg", &dstImg2 );
#endif
        
        row /= 2;
        col /=2;
        t++;
        src = dst;
    }
    
    return dst;
}
Mat IWDT( const Mat &_src, const string _wname, const int _level )const
{
    int reValue = THID_ERR_NONE;
    Mat src = Mat_<float>(_src);
    Mat dst = Mat::zeros( src.rows, src.cols, src.type() );
    int N = src.rows;
    int D = src.cols;
    
    Mat lowFilter;
    Mat highFilter;
    wavelet( _wname, lowFilter, highFilter );
    
    int t=1;
    int row = N/std::pow( 2., _level-1);
    int col = D/std::pow(2., _level-1);
    
    while ( row<=N && col<=D )
    {
        for ( int j=0; j<col; j++ )
        {
            Mat oneCol = Mat::zeros( row, 1, src.type() );
            for ( int i=0; i<row; i++ )
            {
                oneCol.at<float>(i,0) = src.at<float>(i,j);
            }
            oneCol = ( waveletReconstruct( oneCol.t(), lowFilter, highFilter ) ).t();
            
            for ( int i=0; i<row; i++ )
            {
                dst.at<float>(i,j) = oneCol.at<float>(i,0);
            }
        }
        
#if 0
        //normalize( dst, dst, 0, 255, NORM_MINMAX );
        IplImage dstImg2 = IplImage(dst);
        cvSaveImage( "dst.jpg", &dstImg2 );
#endif
        for( int i=0; i<row; i++ )
        {
            Mat oneRow = Mat::zeros( 1,col, src.type() );
            for ( int j=0; j<col; j++ )
            {
                oneRow.at<float>(0,j) = dst.at<float>(i,j);
            }
            oneRow = waveletReconstruct( oneRow, lowFilter, highFilter );
            for ( int j=0; j<col; j++ )
            {
                dst.at<float>(i,j) = oneRow.at<float>(0,j);
            }
        }
        
#if 0
        //normalize( dst, dst, 0, 255, NORM_MINMAX );
        IplImage dstImg1 = IplImage(dst);
        cvSaveImage( "dst.jpg", &dstImg1 );
#endif
        
        row *= 2;
        col *= 2;
        src = dst;
    }
    
    return dst;
}
void wavelet( const string _wname, Mat &_lowFilter, Mat &_highFilter )const
{
    if ( _wname=="haar" || _wname=="db1" )
    {
        int N = 2;
        _lowFilter = Mat::zeros( 1, N, CV_32F );
        _highFilter = Mat::zeros( 1, N, CV_32F );
        
        _lowFilter.at<float>(0, 0) = 1/sqrtf(N);
        _lowFilter.at<float>(0, 1) = 1/sqrtf(N);
        
        _highFilter.at<float>(0, 0) = -1/sqrtf(N);
        _highFilter.at<float>(0, 1) = 1/sqrtf(N);
    }
    if ( _wname =="sym2" )
    {
        int N = 4;
        float h[] = {-0.483, 0.836, -0.224, -0.129 };
        float l[] = {-0.129, 0.224,    0.837, 0.483 };
        
        _lowFilter = Mat::zeros( 1, N, CV_32F );
        _highFilter = Mat::zeros( 1, N, CV_32F );
        
        for ( int i=0; i<N; i++ )
        {
            _lowFilter.at<float>(0, i) = l[i];
            _highFilter.at<float>(0, i) = h[i];
        }
        
    }
}
Mat waveletDecompose( const Mat &_src, const Mat &_lowFilter, const Mat &_highFilter )const
{
    assert( _src.rows==1 && _lowFilter.rows==1 && _highFilter.rows==1 );
    assert( _src.cols>=_lowFilter.cols && _src.cols>=_highFilter.cols );
    Mat &src = Mat_<float>(_src);
    
    int D = src.cols;
    
    Mat &lowFilter = Mat_<float>(_lowFilter);
    Mat &highFilter = Mat_<float>(_highFilter);
    
    Mat dst1 = Mat::zeros( 1, D, src.type() );
    Mat dst2 = Mat::zeros( 1, D, src.type()  );
    
    filter2D( src, dst1, -1, lowFilter );
    filter2D( src, dst2, -1, highFilter );
    
    Mat downDst1 = Mat::zeros( 1, D/2, src.type() );
    Mat downDst2 = Mat::zeros( 1, D/2, src.type() );
    
    resize( dst1, downDst1, downDst1.size() );
    resize( dst2, downDst2, downDst2.size() );
    
    for ( int i=0; i<D/2; i++ )
    {
        src.at<float>(0, i) = downDst1.at<float>( 0, i );
        src.at<float>(0, i+D/2) = downDst2.at<float>( 0, i );
    }
    
    return src;
}
Mat waveletReconstruct( const Mat &_src, const Mat &_lowFilter, const Mat &_highFilter )const
{
    assert( _src.rows==1 && _lowFilter.rows==1 && _highFilter.rows==1 );
    assert( _src.cols>=_lowFilter.cols && _src.cols>=_highFilter.cols );
    Mat &src = Mat_<float>(_src);
    
    int D = src.cols;
    
    Mat &lowFilter = Mat_<float>(_lowFilter);
    Mat &highFilter = Mat_<float>(_highFilter);
    
    Mat Up1 = Mat::zeros( 1, D, src.type() );
    Mat Up2 = Mat::zeros( 1, D, src.type() );
    
    Mat roi1( src, Rect(0, 0, D/2, 1) );
    Mat roi2( src, Rect(D/2, 0, D/2, 1) );
    resize( roi1, Up1, Up1.size(), 0, 0, INTER_CUBIC );
    resize( roi2, Up2, Up2.size(), 0, 0, INTER_CUBIC );
    
    Mat dst1 = Mat::zeros( 1, D, src.type() );
    Mat dst2= Mat::zeros( 1, D, src.type() );
    filter2D( Up1, dst1, -1, lowFilter );
    filter2D( Up2, dst2, -1, highFilter );
    
    dst1 = dst1 + dst2;
    
    return dst1;
}
