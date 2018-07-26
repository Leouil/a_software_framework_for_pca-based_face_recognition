void CPPPCA::InitPCA(String strDirName, Mat matSrc, vector<Mat> vmatSrc)
{
    double dbnumber_principal_compent = 0.95;
    Mat pcaImg, projectImg;
    
    PCA pca(matSrc, Mat(), CV_PCA_DATA_AS_COL, dbnumber_principal_compent);
    Mat EigenVectors = pca.eigenvectors;
    
    pcaImg = normalize(pca.eigenvectors.row(0)).reshape(1, vmatSrc[0].rows);
    //    pcaFace = pca.eigenvectors.row(0).reshape(1, src[0].rows);
    imwrite((("\\") + strDirName + ("\\Pca.jpg")), pcaImg);
    
    Mat EigenValues = pca.eigenvalues;
    
    Mat dst;
    dst = pca.project(matSrc.col(0));
    projectImg = normalize(pca.backProject(dst).col(0)).reshape(1, vmatSrc[0].rows);
    imwrite((("") + strDirName + ("\\Project.jpg")), projectImg);
    
    m_pcaTrain = pca;
    
}    // end of CPPPCA::InitPCA()

void CPPPCA::GetEigenValues(String strDirName, Mat matSrc, vector<Mat> vmatSrc, int iImgNum)
{
    string strInt;
    Mat dst, projectImg;
    const string strEigenValue = "D:\\Summer Project\\EigenValue\\FaceEigenValues.xml";
    FileStorage fs(strEigenValue, FileStorage::WRITE);
    
    for(int i = 0; i < iImgNum; i++){
        strInt = inttostring(i);
        dst = m_pcaTrain.project(matSrc.col(i));
        
        //        projectImg = normalize(m_pcaTrain.backProject(dst).col(0)).reshape(1, vmatSrc[0].rows);
        //        imwrite((("D:\\Summer Project\\") + strDirName + ("\\Project") + strInt + (".jpg")), projectImg);
        fs << "eigenvalue" + strInt << dst;
    }
    fs.release();
    
}    // end of CPPPCA::GetEigenValues()
