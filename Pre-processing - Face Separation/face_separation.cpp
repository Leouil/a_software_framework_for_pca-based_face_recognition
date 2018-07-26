Mat warp_dst = GetAffinedMat(tmpImgCopy2, frame_gray);

IplImage *imgLbpSrc = (&(IplImage)warp_dst);
IplImage *imgLbpDst = cvCreateImage(cvGetSize(imgLbpSrc), IPL_DEPTH_8U,1);;

m_lbpInst.CreatLBP(imgLbpSrc, imgLbpDst);

Mat lbp_dst(imgLbpDst);

Rect roi(51, 19, 80, 89);
Mat matRoi = lbp_dst(roi);
m_vmatLbpFace.push_back(matRoi);

Rect eyesRoi(46, 22, 84, 38);
Mat matEyes = lbp_dst(eyesRoi);
m_vmatLbpEyes.push_back(matEyes);

Rect mouthRoi(69 ,102, 38, 22);
Mat matMouth = lbp_dst(mouthRoi);
m_vmatLbpMouth.push_back(matMouth);

Rect NoseRoi(71, 57, 39, 37);
Mat matNose = lbp_dst(NoseRoi);
m_vmatLbpNose.push_back(matNose);
