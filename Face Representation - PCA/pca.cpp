CPCAAlgoritm::CPCAAlgoritm(int iImgWidth, int iImgHeight, int iNumSamples, PBYTE pbyImgSamples)
: CLinear2DArray(iImgWidth * iImgHeight, iNumSamples, sizeof(float)),
m_iImgWidth(iImgWidth), m_iImgHeight(iImgHeight),
m_iImgSize(iImgWidth * iImgHeight), m_iNumSamples(iNumSamples),
m_pbyImgSamples(pbyImgSamples), m_pPCAModel(NULL)
{
    ASSERT(m_iImgSize > m_iNumSamples);
}
BOOL CPCAAlgoritm::GetPCAModel(CPCAModel *pPCAModel)
{
    float fPrecision = (float)0.001, *pfAllEigenVectors, *pfTempCovarianceArray;
    int iIterationTime = 100;
    DWORD dwFullSize = m_iNumSamples * m_iNumSamples;
    
    //allocate memory for the model
    pPCAModel->SetDimensions(m_iImgWidth, m_iImgHeight, m_iNumSamples);
    
    //allocate temp memory
    pfTempCovarianceArray = new float[dwFullSize];
    pfAllEigenVectors = new float[dwFullSize];
    
    //find the mean vector and the covariance array
    CalculateCovariance(pPCAModel->MeanImgVector(), pPCAModel->CovarianceArray());
    
    //get all eigen vectors from the covariance array
    memcpy(pfTempCovarianceArray, pPCAModel->CovarianceArray(), sizeof(float) * dwFullSize);    //reserve the original covariance array
    ::MatrixEigenVectors(pfTempCovarianceArray, m_iNumSamples, pfAllEigenVectors, fPrecision, iIterationTime);
    
    //get low-dimensional eigen values and vectors
    int iNumEigenVectors = CalcalateNewSD(pfTempCovarianceArray, pfAllEigenVectors, pPCAModel);
    
    //    pPCAModel->SetNumEigens(iNumEigenVectors);
    
    delete []pfAllEigenVectors;
    delete []pfTempCovarianceArray;
    
    return TRUE;
    
}    //end of CPCAAlgoritm::DoPCA()
void CPCAAlgoritm::CalculateCovariance(PFLOAT pfMeanImgVector, PFLOAT pfCovarianceArray)
{
    int i, j, k;
    float dTmp;
    
    // pfMeanImgVector[], row by row
    for (j = 0; j < m_iImgSize; j++){
        dTmp = 0.0;
        for (i = 0; i < m_iNumSamples; i++)
        dTmp += (float)m_pbyImgSamples[m_pRows[i] + j];
        pfMeanImgVector[j] = dTmp / m_iNumSamples;
    }
    
    // m_iNumSamples * m_iNumSamples
    for (i = 0; i < m_iNumSamples; i++){    //row by row
        for (j = 0; j <= i; j ++){        //colum by colum
            dTmp = 0.0;
            for (k = 0; k < m_iImgSize; k++)
            dTmp += ((float)m_pbyImgSamples[m_pRows[i] + k] - pfMeanImgVector[k]) * ((float)m_pbyImgSamples[m_pRows[j] + k] - pfMeanImgVector[k]);
            //dTmp += (float)(m_pbyImgSamples[i][k]) * (float)(m_pbyImgSamples[j][k]);
            pfCovarianceArray[m_pCols[i] + j] = dTmp / m_iImgSize;
        }
    }
    
    //pfCovarianceArray
    for (i = 0; i < m_iNumSamples; i++){
        for (j = i + 1; j < m_iNumSamples; j++)
        pfCovarianceArray[m_pCols[i] + j] = pfCovarianceArray[m_pCols[j] + i];
    }
    
}    //end of CPCAAlgoritm::CalculateCovariance()
int CPCAAlgoritm::CalcalateNewSD(PFLOAT pfCovarianceArray, PFLOAT pfAllEigenVectors, CPCAModel *pPCAModel)
{
    DWORD    dwFullDataSize = m_iNumSamples * m_iNumSamples;
    float    dTmp, pfEigenValuesSum, dValve;
    int        i, j, k, iNumEigenVectors;
    
    PFLOAT pfAllEigenValues  = new float[m_iNumSamples * 2];
    PFLOAT pfResultTmp         = new float[dwFullDataSize];
    PFLOAT pfEigenValueArray = new float[dwFullDataSize];
    PFLOAT pfTmpEigenVectors = new float[dwFullDataSize];
    
    for (i = 0; i < (int)dwFullDataSize; i++)    pfEigenValueArray[i] = pfTmpEigenVectors[i] = 0.0;
    
    for (i = 0; i < m_iNumSamples; i++){
        pfAllEigenValues[i + i] = pfCovarianceArray[m_pCols[i] + i];    //value
        pfAllEigenValues[i + i + 1] = float(i);                                //index
    }
    
    ::BubbleSortFloat(pfAllEigenValues, m_iNumSamples, m_iNumSamples);
    
    pfEigenValuesSum = 0.0;
    for (i = 0; i < 2 * m_iNumSamples; i += 2)    pfEigenValuesSum += pfAllEigenValues[i];
    
    iNumEigenVectors = 0;
    dTmp = 0.0;
    dValve = (float)0.95;
    for (i = 0; i < m_iNumSamples; i++){
        dTmp += pfAllEigenValues[i + i] / pfEigenValuesSum;
        if (dTmp >= dValve){
            iNumEigenVectors = i + 1;
            break;
        }
    }
    
    //set number of eigens and allocate memory
    pPCAModel->SetNumEigens(iNumEigenVectors);
    
    //
    for (i = 0; i < iNumEigenVectors; i++)    pPCAModel->EigenValues()[i] = pfAllEigenValues[i << 1];
    
    // (m_iNumSamples * iNumEigenVectors)
    for (i = 0; i < iNumEigenVectors; i++){
        pfEigenValueArray[m_pCols[i] + i] = float(1.0 / sqrt (pfAllEigenValues[i << 1]));
        k = (int)pfAllEigenValues[i + i + 1];
        for (j = 0; j < m_iNumSamples; j++)
        pfTmpEigenVectors[m_pCols[j] + i] = pfAllEigenVectors[m_pCols[j] + k];
    }
    
    // (iNumEigenVectors * m_iNumSamples)
    for (i = 0; i < m_iNumSamples; i++){
        for (j = 0; j <= i - 1; j++){
            dTmp = pfTmpEigenVectors[m_pCols[i] + j];
            pfTmpEigenVectors[m_pCols[i] + j] = pfTmpEigenVectors[m_pCols[j] + i];
            pfTmpEigenVectors[m_pCols[j] + i] = dTmp;
        }
    }
    
    //    pfEigenValueArray(m_iNumSamples * iNumEigenVectors) x
    //                pfTmpEigenVectors(iNumEigenVectors * m_iNumSamples)   ==> (m_iNumSamples * m_iNumSamples)
    for (i = 0; i < m_iNumSamples; i++){
        for (j = 0; j < m_iNumSamples; j++){
            dTmp = 0.0;
            for (k = 0; k < m_iNumSamples; k++)
            dTmp += pfEigenValueArray[m_pCols[i] + k] * pfTmpEigenVectors[m_pCols[k] + j];
            pfResultTmp[m_pCols[i] + j] = dTmp;
        }
    }
    
    //pEigenImgVectors_iImgSize * iNumEigenVectors
    for (i = 0; i < iNumEigenVectors; i++){
        for (j = 0; j <= m_iImgSize - 1; j++){
            dTmp = 0.0;
            for (k = 0; k < m_iNumSamples; k++)
            dTmp += pfResultTmp[m_pCols[i] + k] * m_pbyImgSamples[m_pRows[k] + j];
            pPCAModel->EigenImgVectors()[m_pRows[i] + j] = dTmp;// + pfMeanImgVector[j];
            //            pEigenImgVectors[i][j] = dTmp / sqrt(m_iNumSamples);
        }
    }
    
    //calculate the square sum for future use
    pPCAModel->CalcEigenSQRSum();
    
    delete []pfAllEigenValues;
    delete []pfResultTmp;
    delete []pfEigenValueArray;
    delete []pfTmpEigenVectors;
    
    return iNumEigenVectors;
    
}    //end of CPCAAlgoritm::CalcalateNewSD()

void CPCAAlgoritm::GetPCACoeff(PBYTE pImg, PFLOAT pCoeff)
{
    int i, j;
    float fTemp;
    
    for (i = 0; i < m_pPCAModel->NumEigens(); i++){
        fTemp = 0.0;
        for (j = 0; j < m_iImgSize; j++)
        //            fTemp += pEigenImgVectors[i][j] * data[j];
        fTemp += (m_pPCAModel->EigenImgVector(m_pRows[i] + j)  / m_pPCAModel->EigenSQRSum(i))
        * ((float)pImg[j] - m_pPCAModel->MeanImgVector(j));
        pCoeff[i] = fTemp;
    }
    
}    //end of CPCAAlgoritm::GetPCACoeff()

void CPCAAlgoritm::ReconstructImage(PBYTE pImg, PFLOAT pCoeff)
{
    int i, j;
    float fTemp;
    
    for (i = 0; i < m_iImgSize; i++){
        fTemp = 0.0;
        for (j = 0; j < m_pPCAModel->NumEigens(); j++)
        fTemp += float(m_pPCAModel->EigenImgVector(m_pRows[j] + i) * pCoeff[j] / m_pPCAModel->EigenSQRSum(j));
        pImg[i] = (BYTE)(fTemp + m_pPCAModel->MeanImgVector(i) + 0.5);
    }
    
    
}    //end of CPCAAlgoritm::ReconstructImage()
