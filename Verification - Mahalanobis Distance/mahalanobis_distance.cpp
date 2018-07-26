double CPPPCA::CalMahDistance(Mat matSrc, Mat matTest)
{
	Mat matCovar, matMean;
	CvMat cvmatSrc, cvmatTest, cvmatCovar;
	cvmatSrc = matSrc;
	cvmatTest = matTest;

	calcCovarMatrix(&matSrc, 1, matCovar, matMean, CV_COVAR_NORMAL);
	cvmatCovar = matCovar;

	double dbMahDistance = cvMahalonobis(&cvmatSrc, &cvmatTest, &cvmatCovar);

	return dbMahDistance;

}	// end of CPPPCA::CalMahDistance()
