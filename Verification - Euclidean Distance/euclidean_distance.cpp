void Distance::distance(double a[],double b[]){
    square=0;
    for(int i=0;i<n;i++){
        square+=(a[i]-b[i])*(a[i]-b[i]);
    }
    Result=sqrt(square/n);
}
