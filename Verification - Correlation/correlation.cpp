double CPPPCA::GetCorelation(uchar *v1, uchar *v2, int n)
{
    double dSignal = 0, dNoise = 0;
    int i;
    double e1 = 0, e2 = 0, e12 = 0, c1 = 0, c2 = 0;
    
    //calculate the expectations of V1, V2 and V1 x V2
    for (i = 0; i < n; i++){
        e1  += v1[i];
        e2  += v2[i];
        e12 += v1[i] * v2[i];
    }
    
    e1  = e1  / n;
    e2  = e2  / n;
    e12 = e12 / n;
    
    //calulate the variances of R1 and R2
    for (i = 0; i < n; i++){
        c1 += (e1 - v1[i]) * (e1 - v1[i]);
        c2 += (e2 - v2[i]) * (e2 - v2[i]);
    }
    
    c1 = sqrt(c1 / n);
    c2 = sqrt(c2 / n);
    
    //calulate the correlation
    return IsZero(c1) || IsZero(c2)? 0.0 : fabs((e12 - e1 * e2) / (c1 * c2));
    
}    //end of CPPPCA::GetCorelation()
