void PCA::load_data(const char* data, char sep){
    
    // Read data
    unsigned int row = 0;
    ifstream file(data);
    if(file.is_open()){
        string line,token;
        while(getline(file, line)){
            stringstream tmp(line);
            unsigned int col = 0;
            while(getline(tmp, token, sep)){
                if(X.rows() < row+1){
                    X.conservativeResize(row+1,X.cols());
                }
                if(X.cols() < col+1){
                    X.conservativeResize(X.rows(),col+1);
                }
                X(row,col) = atof(token.c_str());
                col++;
            }
            row++;
        }
        file.close();
        Xcentered.resize(X.rows(),X.cols());
    }else{
        cout << "Failed to read file " << data << endl;
    }
    
}

double PCA::kernel(const VectorXd& a, const VectorXd& b){
    
    /*
     Kernels
     1 = RBF
     2 = Polynomial
     TODO - add some of these these:
     http://crsouza.blogspot.co.uk/2010/03/kernel-functions-for-machine-learning.html
     */
    switch(kernel_type){
        case 2  :
        return(pow(a.dot(b)+constant,order));
        default :
        return(exp(-gamma*((a-b).squaredNorm())));
    }
    
}

void PCA::run_kpca(){
    
    // Fill kernel matrix
    K.resize(X.rows(),X.rows());
    for(unsigned int i = 0; i < X.rows(); i++){
        for(unsigned int j = i; j < X.rows(); j++){
            K(i,j) = K(j,i) = kernel(X.row(i),X.row(j));
            //printf("k(%i,%i) = %f\n",i,j,K(i,j));
        }
    }
    //cout << endl << K << endl;
    
    EigenSolver<MatrixXd> edecomp(K);
    eigenvalues = edecomp.eigenvalues().real();
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
    vector<pair<double,VectorXd> > eigen_pairs;
    double c = 0.0;
    for(unsigned int i = 0; i < eigenvectors.cols(); i++){
        if(normalise){
            double norm = eigenvectors.col(i).norm();
            eigenvectors.col(i) /= norm;
        }
        eigen_pairs.push_back(make_pair(eigenvalues(i),eigenvectors.col(i)));
    }
    // http://stackoverflow.com/questions/5122804/sorting-with-lambda
    sort(eigen_pairs.begin(),eigen_pairs.end(), [](const pair<double,VectorXd> a, const pair<double,VectorXd> b) -> bool {return (a.first > b.first);} );
    for(unsigned int i = 0; i < eigen_pairs.size(); i++){
        eigenvalues(i) = eigen_pairs[i].first;
        c += eigenvalues(i);
        cumulative(i) = c;
        eigenvectors.col(i) = eigen_pairs[i].second;
    }
    transformed.resize(X.rows(),components);
    
    for(unsigned int i = 0; i < X.rows(); i++){
        for(unsigned int j = 0; j < components; j++){
            for (int k = 0; k < K.rows(); k++){
                transformed(i,j) += K(i,k) * eigenvectors(k,j);
            }
        }
    }
    
    /*
     cout << "Input data:" << endl << X << endl << endl;
     cout << "Centered data:"<< endl << Xcentered << endl << endl;
     cout << "Centered kernel matrix:" << endl << Kcentered << endl << endl;
     cout << "Eigenvalues:" << endl << eigenvalues << endl << endl;
     cout << "Eigenvectors:" << endl << eigenvectors << endl << endl;
     */
    cout << "Sorted eigenvalues:" << endl;
    for(unsigned int i = 0; i < eigenvalues.rows(); i++){
        if(eigenvalues(i) > 0){
            cout << "PC " << i+1 << ": Eigenvalue: " << eigenvalues(i);
            printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",eigenvalues(i)/eigenvalues.sum(),cumulative(i)/eigenvalues.sum());
        }
    }
    cout << endl;
    //cout << "Sorted eigenvectors:" << endl << eigenvectors << endl << endl;
    //cout << "Transformed data:" << endl << transformed << endl << endl;
}

void PCA::run_pca(){
    
    Xcentered = X.rowwise() - X.colwise().mean();
    C = (Xcentered.adjoint() * Xcentered) / double(X.rows());
    EigenSolver<MatrixXd> edecomp(C);
    eigenvalues = edecomp.eigenvalues().real();
    eigenvectors = edecomp.eigenvectors().real();
    cumulative.resize(eigenvalues.rows());
    vector<pair<double,VectorXd> > eigen_pairs;
    double c = 0.0;
    for(unsigned int i = 0; i < eigenvectors.cols(); i++){
        if(normalise){
            double norm = eigenvectors.col(i).norm();
            eigenvectors.col(i) /= norm;
        }
        eigen_pairs.push_back(make_pair(eigenvalues(i),eigenvectors.col(i)));
    }
    // http://stackoverflow.com/questions/5122804/sorting-with-lambda
    sort(eigen_pairs.begin(),eigen_pairs.end(), [](const pair<double,VectorXd> a, const pair<double,VectorXd> b) -> bool {return (a.first > b.first);} );
    for(unsigned int i = 0; i < eigen_pairs.size(); i++){
        eigenvalues(i) = eigen_pairs[i].first;
        c += eigenvalues(i);
        cumulative(i) = c;
        eigenvectors.col(i) = eigen_pairs[i].second;
    }
    transformed = Xcentered * eigenvectors;
    
}

void PCA::print(){
    
    cout << "Input data:" << endl << X << endl << endl;
    cout << "Centered data:"<< endl << Xcentered << endl << endl;
    cout << "Covariance matrix:" << endl << C << endl << endl;
    cout << "Eigenvalues:" << endl << eigenvalues << endl << endl;
    cout << "Eigenvectors:" << endl << eigenvectors << endl << endl;
    cout << "Sorted eigenvalues:" << endl;
    for(unsigned int i = 0; i < eigenvalues.rows(); i++){
        if(eigenvalues(i) > 0){
            cout << "PC " << i+1 << ": Eigenvalue: " << eigenvalues(i);
            printf("\t(%3.3f of variance, cumulative =  %3.3f)\n",eigenvalues(i)/eigenvalues.sum(),cumulative(i)/eigenvalues.sum());
        }
    }
    cout << endl;
    cout << "Sorted eigenvectors:" << endl << eigenvectors << endl << endl;
    cout << "Transformed data:" << endl << X * eigenvectors << endl << endl;
    //cout << "Transformed centred data:" << endl << transformed << endl << endl;
    
}

void PCA::write_transformed(string file){
    
    ofstream outfile(file);
    for(unsigned int i = 0; i < transformed.rows(); i++){
        for(unsigned int j = 0; j < transformed.cols(); j++){
            outfile << transformed(i,j);
            if(j != transformed.cols()-1) outfile << ",";
        }
        outfile << endl;
    }
    outfile.close();
    cout << "Written file " << file << endl;
    
}

void PCA::write_eigenvectors(string file){
    
    ofstream outfile(file);
    for(unsigned int i = 0; i < eigenvectors.rows(); i++){
        for(unsigned int j = 0; j < eigenvectors.cols(); j++){
            outfile << eigenvectors(i,j);
            if(j != eigenvectors.cols()-1) outfile << ",";
        }
        outfile << endl;
    }
    outfile.close();
    cout << "Written file " << file << endl;
    
}
