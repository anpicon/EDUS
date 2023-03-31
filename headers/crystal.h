template <typename T>
void VecProd(T b,T a1,T a2)
{
    b[0]=(a1[1]*a2[2]-a1[2]*a2[1]);
    b[1]=(a1[2]*a2[0]-a1[0]*a2[2]);
    b[2]=(a1[0]*a2[1]-a1[1]*a2[0]);
}


template <typename T, typename U>
double DotProd(T a1,U a2)
{
    double value=0.;
    value=a1[0]*a2[0]+a1[1]*a2[1]+a1[2]*a2[2];
    return value;
}



complex<double> DDipW(vec5x& Dx,vec5x& Dy,vec5x& Dz, vec1d& ue1, int ikkx, int ikky, int ikkz, int i, int j, double EF1)
{
    return (Dx[ikkx][ikky][ikkz][i][j]*ue1[0]+Dy[ikkx][ikky][ikkz][i][j]*ue1[1]+Dz[ikkx][ikky][ikkz][i][j]*ue1[2])*EF1;
}

complex<double> DDip(vec5d& Dx,vec5d& Dy,vec5d& Dz, vec1d& ue1, int ikkx, int ikky, int ikkz, int i, int j, double EF1)
{
    return (Dx[ikkx][ikky][ikkz][i][j]*ue1[0]+Dy[ikkx][ikky][ikkz][i][j]*ue1[1]+Dz[ikkx][ikky][ikkz][i][j]*ue1[2])*EF1;
}
