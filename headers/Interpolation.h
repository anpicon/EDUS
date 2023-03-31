/*
 *  Interpolation.h
 *  
 *  Created by Antonio Picon on 5/04/19.
 *  Copyright 2019 Universidad Autónoma de Madrid
 *
 */

//############################# FOR COMPLEX DATA #############################
template <typename T>
class FittedData {
    
public:
    string SizeData;
    vec1d _spacing;
    vec1d _xmin;
    int NX;
    int NY;
    int NZ;
    multivec2D<T> coeff1D;
    multivec3D<T> coeff2D;
    multivec4D<T> coeff3D;
    FittedData(multivec3D<T>&,vec1d&,vec1d&);
    FittedData(){};
    void construct(multivec3D<T>&,vec1d&,vec1d&);
    T operator()(double x, double y, double z);
    string get_size();
};//FittedData



template <typename T>
FittedData<T>::FittedData(multivec3D<T>& Data,vec1d& spacing,vec1d& xmin)
{
    _spacing.resize(3);
    _xmin.resize(3);   
    for (int i=0; i<3; i++) {
        _spacing[i]=spacing[i];
        _xmin[i]=xmin[i];
    }
    NX=Data.n1();
    NY=Data.n2();
    NZ=Data.n3();
    string iPrint="none"; //Control printing for debugging purposes
    
    if(NX==1 && NY==1) {
        if (iPrint=="all") printf("Fitting 1D data...\n");
        SizeData="1Ddata";
        //coeff1D.resize(NZ,vec1x(4,0.));
        coeff1D.resize(NZ,4); coeff1D.fill(0.);
        if(_spacing[2] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
    }
    else if(NX==1) {
        if (iPrint=="all") printf("Fitting 2D data...\n");
        //coeff2D.resize(NY,vec2x(NZ,vec1x(16,0.)));
        coeff2D.resize(NY,NZ,16); coeff2D.fill(0.);
        SizeData="2Ddata";
        if(_spacing[1] <= 0 && _spacing[2] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
    }
    else {
        if (iPrint=="all") printf("Fitting 3D data...\n");
        //coeff3D.resize(NX+3,vec3x(NY,vec2x(NZ,vec1x(16,0.))));
        coeff3D.resize(NX+3,NY,NZ,16); coeff3D.fill(0.);
        SizeData="3Ddata";
        if(_spacing[0] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
        if(_spacing[1] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
        if(_spacing[2] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
    }
    
    //multivector to impose the boundary conditions in the boundary points calculations
    //vec3x _data(NX+3,vec2x(NY+3,vec1x(NZ+3,0.)));
    multivec3D<T> _data(NX+3,NY+3,NZ+3); _data.fill(0.);
    for (int ix=0; ix<NX+3; ix++)
    {
        int ixx=(NX > 1 ? ix-1 : ix);
        ixx=(ixx >= 0 ? ixx % NX : NX + ixx % NX);
        for (int iy=0; iy<NY+3; iy++)
        {
            int iyy=(NY > 1 ? iy-1 : iy);
            iyy=(iyy >= 0 ? iyy % NY : NY + iyy % NY);
            for (int iz=0; iz<NZ+3; iz++)
            {
                int izz=iz-1;
                izz=(izz >= 0 ? izz % NZ : NZ + izz % NZ);
                //cout << "ix " << ix << " , ixx " << ixx << " , iy " << iy << " , iyy " << iyy << " , iz " << iz << " , izz " << izz << endl;
                _data[ix][iy][iz]=Data[ixx][iyy][izz];
            }
        }
    }
    
    if (SizeData=="1Ddata") {
        for (int ix=1; ix<NX+1; ix++)
        {
            for (int iy=1; iy<NY+1; iy++)
            {
                for (int iz=1; iz<NZ+1; iz++)
                {
                    coeff1D[iz-1][0]=-0.5*_data[ix][iy][iz-1]+1.5*_data[ix][iy][iz]-1.5*_data[ix][iy][iz+1]+0.5*_data[ix][iy][iz+2];
                    coeff1D[iz-1][1]=_data[ix][iy][iz-1]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy][iz+1]-0.5*_data[ix][iy][iz+2];
                    coeff1D[iz-1][2]=-0.5*_data[ix][iy][iz-1]+0.5*_data[ix][iy][iz+1];
                    coeff1D[iz-1][3]=_data[ix][iy][iz];
                }
            }
        }
        if (iPrint=="all") printf("Done fitting 1D data...\n");
    }//-----1Ddata
    
    else if (SizeData=="2Ddata") {
        for (int ix=1; ix<NX+1; ix++)
        {
            for (int iy=1; iy<NY+1; iy++)
            {
                for (int iz=1; iz<NZ+1; iz++)
                {
                    coeff2D[iy-1][iz-1][0]=_data[ix][iy][iz]; //a00=p[1][1]
                    coeff2D[iy-1][iz-1][1]=-0.5*_data[ix][iy][iz-1]+0.5*_data[ix][iy][iz+1]; //a01=-.5*p[1][0] + .5*p[1][2]
                    coeff2D[iy-1][iz-1][2]=_data[ix][iy][iz-1]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy][iz+1]-0.5*_data[ix][iy][iz+2]; //a02=p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3]
                    coeff2D[iy-1][iz-1][3]=-0.5*_data[ix][iy][iz-1]+1.5*_data[ix][iy][iz]-1.5*_data[ix][iy][iz+1]+0.5*_data[ix][iy][iz+2]; //a03=-.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3]
                    coeff2D[iy-1][iz-1][4]=-0.5*_data[ix][iy-1][iz]+0.5*_data[ix][iy+1][iz]; //a10 = -.5*p[0][1] + .5*p[2][1];
                    coeff2D[iy-1][iz-1][5]=0.25*_data[ix][iy-1][iz-1]-0.25*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy+1][iz-1]+0.25*_data[ix][iy+1][iz+1]; //a11 = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
                    coeff2D[iy-1][iz-1][6]=-0.5*_data[ix][iy-1][iz-1]+1.25*_data[ix][iy-1][iz]-_data[ix][iy-1][iz+1]+0.25*_data[ix][iy-1][iz+2]+0.5*_data[ix][iy+1][iz-1]-1.25*_data[ix][iy+1][iz]+_data[ix][iy+1][iz+1]-0.25*_data[ix][iy+1][iz+2]; //a12 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3]
                    coeff2D[iy-1][iz-1][7]=0.25*_data[ix][iy-1][iz-1]-0.75*_data[ix][iy-1][iz]+0.75*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy-1][iz+2]-0.25*_data[ix][iy+1][iz-1]+0.75*_data[ix][iy+1][iz]-0.75*_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+1][iz+2]; //a13 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
                    coeff2D[iy-1][iz-1][8]=_data[ix][iy-1][iz]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy+1][iz]-0.5*_data[ix][iy+2][iz]; //a20 = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
                    coeff2D[iy-1][iz-1][9]=-0.5*_data[ix][iy-1][iz-1]+0.5*_data[ix][iy-1][iz+1]+1.25*_data[ix][iy][iz-1]-1.25*_data[ix][iy][iz+1]-_data[ix][iy+1][iz-1]+_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+2][iz-1]-0.25*_data[ix][iy+2][iz+1]; //a21=-.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
                    coeff2D[iy-1][iz-1][10]= _data[ix][iy-1][iz-1] - 2.5*_data[ix][iy-1][iz] + 2.*_data[ix][iy-1][iz+1] - .5*_data[ix][iy-1][iz+2] - 2.5*_data[ix][iy][iz-1] + 6.25*_data[ix][iy][iz] - 5.*_data[ix][iy][iz+1] + 1.25*_data[ix][iy][iz+2] + 2.*_data[ix][iy+1][iz-1] - 5.*_data[ix][iy+1][iz] + 4.*_data[ix][iy+1][iz+1] - _data[ix][iy+1][iz+2] - .5*_data[ix][iy+2][iz-1] + 1.25*_data[ix][iy+2][iz] - _data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a22 = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
                    coeff2D[iy-1][iz-1][11]= -.5*_data[ix][iy-1][iz-1] + 1.5*_data[ix][iy-1][iz] - 1.5*_data[ix][iy-1][iz+1] + .5*_data[ix][iy-1][iz+2] + 1.25*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.75*_data[ix][iy][iz+1] - 1.25*_data[ix][iy][iz+2] - _data[ix][iy+1][iz-1] + 3.*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + _data[ix][iy+1][iz+2] + .25*_data[ix][iy+2][iz-1] - .75*_data[ix][iy+2][iz] + .75*_data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a23 = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
                    coeff2D[iy-1][iz-1][12]= -.5*_data[ix][iy-1][iz] + 1.5*_data[ix][iy][iz] - 1.5*_data[ix][iy+1][iz] + .5*_data[ix][iy+2][iz];//a30 = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
                    coeff2D[iy-1][iz-1][13]= .25*_data[ix][iy-1][iz-1] - .25*_data[ix][iy-1][iz+1] - .75*_data[ix][iy][iz-1] + .75*_data[ix][iy][iz+1] + .75*_data[ix][iy+1][iz-1] - .75*_data[ix][iy+1][iz+1] - .25*_data[ix][iy+2][iz-1] + .25*_data[ix][iy+2][iz+1]; //a31 = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
                    coeff2D[iy-1][iz-1][14]= -.5*_data[ix][iy-1][iz-1] + 1.25*_data[ix][iy-1][iz] - _data[ix][iy-1][iz+1] + .25*_data[ix][iy-1][iz+2] + 1.5*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.*_data[ix][iy][iz+1] - .75*_data[ix][iy][iz+2] - 1.5*_data[ix][iy+1][iz-1] + 3.75*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + .75*_data[ix][iy+1][iz+2] + .5*_data[ix][iy+2][iz-1] - 1.25*_data[ix][iy+2][iz] + _data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a32 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
                    coeff2D[iy-1][iz-1][15]= .25*_data[ix][iy-1][iz-1] - .75*_data[ix][iy-1][iz] + .75*_data[ix][iy-1][iz+1] - .25*_data[ix][iy-1][iz+2] - .75*_data[ix][iy][iz-1] + 2.25*_data[ix][iy][iz] - 2.25*_data[ix][iy][iz+1] + .75*_data[ix][iy][iz+2] + .75*_data[ix][iy+1][iz-1] - 2.25*_data[ix][iy+1][iz] + 2.25*_data[ix][iy+1][iz+1] - .75*_data[ix][iy+1][iz+2] - .25*_data[ix][iy+2][iz-1] + .75*_data[ix][iy+2][iz] - .75*_data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a33 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];
                
                }
            }
        }
        if (iPrint=="all") printf("Done fitting 2D data...\n");
    }//-----2Ddata
    
    else if (SizeData=="3Ddata") {
        for (int ix=0; ix<NX+3; ix++)
        {
            for (int iy=1; iy<NY+1; iy++)
            {
                for (int iz=1; iz<NZ+1; iz++)
                {
                    coeff3D[ix][iy-1][iz-1][0]=_data[ix][iy][iz]; //a00=p[1][1]
                    coeff3D[ix][iy-1][iz-1][1]=-0.5*_data[ix][iy][iz-1]+0.5*_data[ix][iy][iz+1]; //a01=-.5*p[1][0] + .5*p[1][2]
                    coeff3D[ix][iy-1][iz-1][2]=_data[ix][iy][iz-1]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy][iz+1]-0.5*_data[ix][iy][iz+2]; //a02=p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3]
                    coeff3D[ix][iy-1][iz-1][3]=-0.5*_data[ix][iy][iz-1]+1.5*_data[ix][iy][iz]-1.5*_data[ix][iy][iz+1]+0.5*_data[ix][iy][iz+2]; //a03=-.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3]
                    coeff3D[ix][iy-1][iz-1][4]=-0.5*_data[ix][iy-1][iz]+0.5*_data[ix][iy+1][iz]; //a10 = -.5*p[0][1] + .5*p[2][1];
                    coeff3D[ix][iy-1][iz-1][5]=0.25*_data[ix][iy-1][iz-1]-0.25*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy+1][iz-1]+0.25*_data[ix][iy+1][iz+1]; //a11 = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
                    coeff3D[ix][iy-1][iz-1][6]=-0.5*_data[ix][iy-1][iz-1]+1.25*_data[ix][iy-1][iz]-_data[ix][iy-1][iz+1]+0.25*_data[ix][iy-1][iz+2]+0.5*_data[ix][iy+1][iz-1]-1.25*_data[ix][iy+1][iz]+_data[ix][iy+1][iz+1]-0.25*_data[ix][iy+1][iz+2]; //a12 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3]
                    coeff3D[ix][iy-1][iz-1][7]=0.25*_data[ix][iy-1][iz-1]-0.75*_data[ix][iy-1][iz]+0.75*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy-1][iz+2]-0.25*_data[ix][iy+1][iz-1]+0.75*_data[ix][iy+1][iz]-0.75*_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+1][iz+2]; //a13 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
                    coeff3D[ix][iy-1][iz-1][8]=_data[ix][iy-1][iz]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy+1][iz]-0.5*_data[ix][iy+2][iz]; //a20 = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
                    coeff3D[ix][iy-1][iz-1][9]=-0.5*_data[ix][iy-1][iz-1]+0.5*_data[ix][iy-1][iz+1]+1.25*_data[ix][iy][iz-1]-1.25*_data[ix][iy][iz+1]-_data[ix][iy+1][iz-1]+_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+2][iz-1]-0.25*_data[ix][iy+2][iz+1]; //a21=-.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
                    coeff3D[ix][iy-1][iz-1][10]= _data[ix][iy-1][iz-1] - 2.5*_data[ix][iy-1][iz] + 2.*_data[ix][iy-1][iz+1] - .5*_data[ix][iy-1][iz+2] - 2.5*_data[ix][iy][iz-1] + 6.25*_data[ix][iy][iz] - 5.*_data[ix][iy][iz+1] + 1.25*_data[ix][iy][iz+2] + 2.*_data[ix][iy+1][iz-1] - 5.*_data[ix][iy+1][iz] + 4.*_data[ix][iy+1][iz+1] - _data[ix][iy+1][iz+2] - .5*_data[ix][iy+2][iz-1] + 1.25*_data[ix][iy+2][iz] - _data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a22 = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
                    coeff3D[ix][iy-1][iz-1][11]= -.5*_data[ix][iy-1][iz-1] + 1.5*_data[ix][iy-1][iz] - 1.5*_data[ix][iy-1][iz+1] + .5*_data[ix][iy-1][iz+2] + 1.25*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.75*_data[ix][iy][iz+1] - 1.25*_data[ix][iy][iz+2] - _data[ix][iy+1][iz-1] + 3.*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + _data[ix][iy+1][iz+2] + .25*_data[ix][iy+2][iz-1] - .75*_data[ix][iy+2][iz] + .75*_data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a23 = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
                    coeff3D[ix][iy-1][iz-1][12]= -.5*_data[ix][iy-1][iz] + 1.5*_data[ix][iy][iz] - 1.5*_data[ix][iy+1][iz] + .5*_data[ix][iy+2][iz];//a30 = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
                    coeff3D[ix][iy-1][iz-1][13]= .25*_data[ix][iy-1][iz-1] - .25*_data[ix][iy-1][iz+1] - .75*_data[ix][iy][iz-1] + .75*_data[ix][iy][iz+1] + .75*_data[ix][iy+1][iz-1] - .75*_data[ix][iy+1][iz+1] - .25*_data[ix][iy+2][iz-1] + .25*_data[ix][iy+2][iz+1]; //a31 = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
                    coeff3D[ix][iy-1][iz-1][14]= -.5*_data[ix][iy-1][iz-1] + 1.25*_data[ix][iy-1][iz] - _data[ix][iy-1][iz+1] + .25*_data[ix][iy-1][iz+2] + 1.5*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.*_data[ix][iy][iz+1] - .75*_data[ix][iy][iz+2] - 1.5*_data[ix][iy+1][iz-1] + 3.75*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + .75*_data[ix][iy+1][iz+2] + .5*_data[ix][iy+2][iz-1] - 1.25*_data[ix][iy+2][iz] + _data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a32 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
                    coeff3D[ix][iy-1][iz-1][15]= .25*_data[ix][iy-1][iz-1] - .75*_data[ix][iy-1][iz] + .75*_data[ix][iy-1][iz+1] - .25*_data[ix][iy-1][iz+2] - .75*_data[ix][iy][iz-1] + 2.25*_data[ix][iy][iz] - 2.25*_data[ix][iy][iz+1] + .75*_data[ix][iy][iz+2] + .75*_data[ix][iy+1][iz-1] - 2.25*_data[ix][iy+1][iz] + 2.25*_data[ix][iy+1][iz+1] - .75*_data[ix][iy+1][iz+2] - .25*_data[ix][iy+2][iz-1] + .75*_data[ix][iy+2][iz] - .75*_data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a33 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];
                }
            }
        }
        if (iPrint=="all") printf("Done fitting 3D data...\n");

    }//-----3Ddata
    
} //----End   FittedData::FittedData

template <typename T>
T FittedData<T>::operator()(double x, double y, double z)
{
    // Map x,y,z to a point dx,dy,dz in the cube [0,NX) x [0,NY) x [0,NZ)
    double dx(std::fmod((x-_xmin[0])/_spacing[0],NX)), dy(std::fmod((y-_xmin[1])/_spacing[1],NY)), dz(std::fmod((z-_xmin[2])/_spacing[2],NZ));
    if(dx < 0) dx += NX;
    if(dy < 0) dy += NY;
    if(dz < 0) dz += NZ;
    //discrete integer values giving the locatio of a small cube in the [0,NX) x [0,NY) x [0,NZ) grid
    int xi = (int)std::floor(dx);
    int yi = (int)std::floor(dy);
    int zi = (int)std::floor(dz);
    //double values within a particular small cube at (xi,yi,zi)
    dx -= xi;
    dy -= yi;
    dz -= zi;
    
    /*
    cout.precision(13);
    cout << z-_xmin[2] << " " << (z-_xmin[2])/_spacing[2]  << " " << std::fmod((z-_xmin[2])/_spacing[2],NZ) << endl;
    cout << dz << " " << std::floor(dz) << " " << zi << endl;
    */
    if (SizeData=="3Ddata") {
        double y2 = dy * dy;
        double y3 = y2 * dy;
        double z2 = dz * dz;
        double z3 = z2 * dz;
        xi++;
        
        T p0=(coeff3D[xi-1][yi][zi][0] + coeff3D[xi-1][yi][zi][1] * dz + coeff3D[xi-1][yi][zi][2] * z2 + coeff3D[xi-1][yi][zi][3] * z3) +
        (coeff3D[xi-1][yi][zi][4] + coeff3D[xi-1][yi][zi][5] * dz + coeff3D[xi-1][yi][zi][6] * z2 + coeff3D[xi-1][yi][zi][7] * z3) * dy +
        (coeff3D[xi-1][yi][zi][8] + coeff3D[xi-1][yi][zi][9] * dz + coeff3D[xi-1][yi][zi][10] * z2 + coeff3D[xi-1][yi][zi][11] * z3) * y2 +
        (coeff3D[xi-1][yi][zi][12] + coeff3D[xi-1][yi][zi][13] * dz + coeff3D[xi-1][yi][zi][14] * z2 + coeff3D[xi-1][yi][zi][15] * z3) * y3;
        
        T p1=(coeff3D[xi][yi][zi][0] + coeff3D[xi][yi][zi][1] * dz + coeff3D[xi][yi][zi][2] * z2 + coeff3D[xi][yi][zi][3] * z3) +
        (coeff3D[xi][yi][zi][4] + coeff3D[xi][yi][zi][5] * dz + coeff3D[xi][yi][zi][6] * z2 + coeff3D[xi][yi][zi][7] * z3) * dy +
        (coeff3D[xi][yi][zi][8] + coeff3D[xi][yi][zi][9] * dz + coeff3D[xi][yi][zi][10] * z2 + coeff3D[xi][yi][zi][11] * z3) * y2 +
        (coeff3D[xi][yi][zi][12] + coeff3D[xi][yi][zi][13] * dz + coeff3D[xi][yi][zi][14] * z2 + coeff3D[xi][yi][zi][15] * z3) * y3;
        
        T p2=(coeff3D[xi+1][yi][zi][0] + coeff3D[xi+1][yi][zi][1] * dz + coeff3D[xi+1][yi][zi][2] * z2 + coeff3D[xi+1][yi][zi][3] * z3) +
        (coeff3D[xi+1][yi][zi][4] + coeff3D[xi+1][yi][zi][5] * dz + coeff3D[xi+1][yi][zi][6] * z2 + coeff3D[xi+1][yi][zi][7] * z3) * dy +
        (coeff3D[xi+1][yi][zi][8] + coeff3D[xi+1][yi][zi][9] * dz + coeff3D[xi+1][yi][zi][10] * z2 + coeff3D[xi+1][yi][zi][11] * z3) * y2 +
        (coeff3D[xi+1][yi][zi][12] + coeff3D[xi+1][yi][zi][13] * dz + coeff3D[xi+1][yi][zi][14] * z2 + coeff3D[xi+1][yi][zi][15] * z3) * y3;
        
        T p3=(coeff3D[xi+2][yi][zi][0] + coeff3D[xi+2][yi][zi][1] * dz + coeff3D[xi+2][yi][zi][2] * z2 + coeff3D[xi+2][yi][zi][3] * z3) +
        (coeff3D[xi+2][yi][zi][4] + coeff3D[xi+2][yi][zi][5] * dz + coeff3D[xi+2][yi][zi][6] * z2 + coeff3D[xi+2][yi][zi][7] * z3) * dy +
        (coeff3D[xi+2][yi][zi][8] + coeff3D[xi+2][yi][zi][9] * dz + coeff3D[xi+2][yi][zi][10] * z2 + coeff3D[xi+2][yi][zi][11] * z3) * y2 +
        (coeff3D[xi+2][yi][zi][12] + coeff3D[xi+2][yi][zi][13] * dz + coeff3D[xi+2][yi][zi][14] * z2 + coeff3D[xi+2][yi][zi][15] * z3) * y3;
        
        return p1 + 0.5 * dx*(p2 - p0 + dx*(2.0*p0 - 5.0*p1 + 4.0*p2 - p3 + dx*(3.0*(p1 - p2) + p3 - p0)));
    }
    if (SizeData=="2Ddata") {
        double y2 = dy * dy;
        double y3 = y2 * dy;
        double z2 = dz * dz;
        double z3 = z2 * dz;
        
        return (coeff2D[yi][zi][0] + coeff2D[yi][zi][1] * dz + coeff2D[yi][zi][2] * z2 + coeff2D[yi][zi][3] * z3) +
        (coeff2D[yi][zi][4] + coeff2D[yi][zi][5] * dz + coeff2D[yi][zi][6] * z2 + coeff2D[yi][zi][7] * z3) * dy +
        (coeff2D[yi][zi][8] + coeff2D[yi][zi][9] * dz + coeff2D[yi][zi][10] * z2 + coeff2D[yi][zi][11] * z3) * y2 +
        (coeff2D[yi][zi][12] + coeff2D[yi][zi][13] * dz + coeff2D[yi][zi][14] * z2 + coeff2D[yi][zi][15] * z3) * y3;
    }
    else if (SizeData=="1Ddata") {
        return (coeff1D[zi][0]*dz*dz*dz+coeff1D[zi][1]*dz*dz+coeff1D[zi][2]*dz+coeff1D[zi][3]);
    }
    else return 0.;
} //----End    FittedData::operator()

template <typename T>
string FittedData<T>::get_size()
{
    return SizeData;
}



template <typename T>
void FittedData<T>::construct(multivec3D<T>& Data,vec1d& spacing,vec1d& xmin)
{
    _spacing.resize(3);
    _xmin.resize(3);   
    for (int i=0; i<3; i++) {
        _spacing[i]=spacing[i];
        _xmin[i]=xmin[i];
    }
    NX=Data.n1();
    NY=Data.n2();
    NZ=Data.n3();
    string iPrint="none"; //Control printing for debugging purposes
    
    if(NX==1 && NY==1) {
        if (iPrint=="all") printf("Fitting 1D data...\n");
        SizeData="1Ddata";
        //coeff1D.resize(NZ,vec1x(4,0.));
        coeff1D.resize(NZ,4); coeff1D.fill(0.);
        if(_spacing[2] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
    }
    else if(NX==1) {
        if (iPrint=="all") printf("Fitting 2D data...\n");
        //coeff2D.resize(NY,vec2x(NZ,vec1x(16,0.)));
        coeff2D.resize(NY,NZ,16); coeff2D.fill(0.);
        SizeData="2Ddata";
        if(_spacing[1] <= 0 && _spacing[2] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
    }
    else {
        if (iPrint=="all") printf("Fitting 3D data...\n");
        //coeff3D.resize(NX+3,vec3x(NY,vec2x(NZ,vec1x(16,0.))));
        coeff3D.resize(NX+3,NY,NZ,16); coeff3D.fill(0.);
        SizeData="3Ddata";
        if(_spacing[0] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
        if(_spacing[1] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
        if(_spacing[2] <= 0) {
            printf("Bad datacube grid spacing.\n");
            exit(1);
        }
    }
    
    //multivector to impose the boundary conditions in the boundary points calculations
    //vec3x _data(NX+3,vec2x(NY+3,vec1x(NZ+3,0.)));
    multivec3D<T> _data(NX+3,NY+3,NZ+3); _data.fill(0.);
    for (int ix=0; ix<NX+3; ix++)
    {
        int ixx=(NX > 1 ? ix-1 : ix);
        ixx=(ixx >= 0 ? ixx % NX : NX + ixx % NX);
        for (int iy=0; iy<NY+3; iy++)
        {
            int iyy=(NY > 1 ? iy-1 : iy);
            iyy=(iyy >= 0 ? iyy % NY : NY + iyy % NY);
            for (int iz=0; iz<NZ+3; iz++)
            {
                int izz=iz-1;
                izz=(izz >= 0 ? izz % NZ : NZ + izz % NZ);
                //cout << "ix " << ix << " , ixx " << ixx << " , iy " << iy << " , iyy " << iyy << " , iz " << iz << " , izz " << izz << endl;
                _data[ix][iy][iz]=Data[ixx][iyy][izz];
            }
        }
    }
    
    if (SizeData=="1Ddata") {
        for (int ix=1; ix<NX+1; ix++)
        {
            for (int iy=1; iy<NY+1; iy++)
            {
                for (int iz=1; iz<NZ+1; iz++)
                {
                    coeff1D[iz-1][0]=-0.5*_data[ix][iy][iz-1]+1.5*_data[ix][iy][iz]-1.5*_data[ix][iy][iz+1]+0.5*_data[ix][iy][iz+2];
                    coeff1D[iz-1][1]=_data[ix][iy][iz-1]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy][iz+1]-0.5*_data[ix][iy][iz+2];
                    coeff1D[iz-1][2]=-0.5*_data[ix][iy][iz-1]+0.5*_data[ix][iy][iz+1];
                    coeff1D[iz-1][3]=_data[ix][iy][iz];
                }
            }
        }
        if (iPrint=="all") printf("Done fitting 1D data...\n");
    }//-----1Ddata
    
    else if (SizeData=="2Ddata") {
        for (int ix=1; ix<NX+1; ix++)
        {
            for (int iy=1; iy<NY+1; iy++)
            {
                for (int iz=1; iz<NZ+1; iz++)
                {
                    coeff2D[iy-1][iz-1][0]=_data[ix][iy][iz]; //a00=p[1][1]
                    coeff2D[iy-1][iz-1][1]=-0.5*_data[ix][iy][iz-1]+0.5*_data[ix][iy][iz+1]; //a01=-.5*p[1][0] + .5*p[1][2]
                    coeff2D[iy-1][iz-1][2]=_data[ix][iy][iz-1]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy][iz+1]-0.5*_data[ix][iy][iz+2]; //a02=p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3]
                    coeff2D[iy-1][iz-1][3]=-0.5*_data[ix][iy][iz-1]+1.5*_data[ix][iy][iz]-1.5*_data[ix][iy][iz+1]+0.5*_data[ix][iy][iz+2]; //a03=-.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3]
                    coeff2D[iy-1][iz-1][4]=-0.5*_data[ix][iy-1][iz]+0.5*_data[ix][iy+1][iz]; //a10 = -.5*p[0][1] + .5*p[2][1];
                    coeff2D[iy-1][iz-1][5]=0.25*_data[ix][iy-1][iz-1]-0.25*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy+1][iz-1]+0.25*_data[ix][iy+1][iz+1]; //a11 = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
                    coeff2D[iy-1][iz-1][6]=-0.5*_data[ix][iy-1][iz-1]+1.25*_data[ix][iy-1][iz]-_data[ix][iy-1][iz+1]+0.25*_data[ix][iy-1][iz+2]+0.5*_data[ix][iy+1][iz-1]-1.25*_data[ix][iy+1][iz]+_data[ix][iy+1][iz+1]-0.25*_data[ix][iy+1][iz+2]; //a12 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3]
                    coeff2D[iy-1][iz-1][7]=0.25*_data[ix][iy-1][iz-1]-0.75*_data[ix][iy-1][iz]+0.75*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy-1][iz+2]-0.25*_data[ix][iy+1][iz-1]+0.75*_data[ix][iy+1][iz]-0.75*_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+1][iz+2]; //a13 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
                    coeff2D[iy-1][iz-1][8]=_data[ix][iy-1][iz]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy+1][iz]-0.5*_data[ix][iy+2][iz]; //a20 = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
                    coeff2D[iy-1][iz-1][9]=-0.5*_data[ix][iy-1][iz-1]+0.5*_data[ix][iy-1][iz+1]+1.25*_data[ix][iy][iz-1]-1.25*_data[ix][iy][iz+1]-_data[ix][iy+1][iz-1]+_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+2][iz-1]-0.25*_data[ix][iy+2][iz+1]; //a21=-.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
                    coeff2D[iy-1][iz-1][10]= _data[ix][iy-1][iz-1] - 2.5*_data[ix][iy-1][iz] + 2.*_data[ix][iy-1][iz+1] - .5*_data[ix][iy-1][iz+2] - 2.5*_data[ix][iy][iz-1] + 6.25*_data[ix][iy][iz] - 5.*_data[ix][iy][iz+1] + 1.25*_data[ix][iy][iz+2] + 2.*_data[ix][iy+1][iz-1] - 5.*_data[ix][iy+1][iz] + 4.*_data[ix][iy+1][iz+1] - _data[ix][iy+1][iz+2] - .5*_data[ix][iy+2][iz-1] + 1.25*_data[ix][iy+2][iz] - _data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a22 = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
                    coeff2D[iy-1][iz-1][11]= -.5*_data[ix][iy-1][iz-1] + 1.5*_data[ix][iy-1][iz] - 1.5*_data[ix][iy-1][iz+1] + .5*_data[ix][iy-1][iz+2] + 1.25*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.75*_data[ix][iy][iz+1] - 1.25*_data[ix][iy][iz+2] - _data[ix][iy+1][iz-1] + 3.*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + _data[ix][iy+1][iz+2] + .25*_data[ix][iy+2][iz-1] - .75*_data[ix][iy+2][iz] + .75*_data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a23 = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
                    coeff2D[iy-1][iz-1][12]= -.5*_data[ix][iy-1][iz] + 1.5*_data[ix][iy][iz] - 1.5*_data[ix][iy+1][iz] + .5*_data[ix][iy+2][iz];//a30 = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
                    coeff2D[iy-1][iz-1][13]= .25*_data[ix][iy-1][iz-1] - .25*_data[ix][iy-1][iz+1] - .75*_data[ix][iy][iz-1] + .75*_data[ix][iy][iz+1] + .75*_data[ix][iy+1][iz-1] - .75*_data[ix][iy+1][iz+1] - .25*_data[ix][iy+2][iz-1] + .25*_data[ix][iy+2][iz+1]; //a31 = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
                    coeff2D[iy-1][iz-1][14]= -.5*_data[ix][iy-1][iz-1] + 1.25*_data[ix][iy-1][iz] - _data[ix][iy-1][iz+1] + .25*_data[ix][iy-1][iz+2] + 1.5*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.*_data[ix][iy][iz+1] - .75*_data[ix][iy][iz+2] - 1.5*_data[ix][iy+1][iz-1] + 3.75*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + .75*_data[ix][iy+1][iz+2] + .5*_data[ix][iy+2][iz-1] - 1.25*_data[ix][iy+2][iz] + _data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a32 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
                    coeff2D[iy-1][iz-1][15]= .25*_data[ix][iy-1][iz-1] - .75*_data[ix][iy-1][iz] + .75*_data[ix][iy-1][iz+1] - .25*_data[ix][iy-1][iz+2] - .75*_data[ix][iy][iz-1] + 2.25*_data[ix][iy][iz] - 2.25*_data[ix][iy][iz+1] + .75*_data[ix][iy][iz+2] + .75*_data[ix][iy+1][iz-1] - 2.25*_data[ix][iy+1][iz] + 2.25*_data[ix][iy+1][iz+1] - .75*_data[ix][iy+1][iz+2] - .25*_data[ix][iy+2][iz-1] + .75*_data[ix][iy+2][iz] - .75*_data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a33 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];
                }
            }
        }
        if (iPrint=="all") printf("Done fitting 2D data...\n");
    }//-----2Ddata
    
    else if (SizeData=="3Ddata") {
        for (int ix=0; ix<NX+3; ix++)
        {
            for (int iy=1; iy<NY+1; iy++)
            {
                for (int iz=1; iz<NZ+1; iz++)
                {
                    coeff3D[ix][iy-1][iz-1][0]=_data[ix][iy][iz]; //a00=p[1][1]
                    coeff3D[ix][iy-1][iz-1][1]=-0.5*_data[ix][iy][iz-1]+0.5*_data[ix][iy][iz+1]; //a01=-.5*p[1][0] + .5*p[1][2]
                    coeff3D[ix][iy-1][iz-1][2]=_data[ix][iy][iz-1]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy][iz+1]-0.5*_data[ix][iy][iz+2]; //a02=p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3]
                    coeff3D[ix][iy-1][iz-1][3]=-0.5*_data[ix][iy][iz-1]+1.5*_data[ix][iy][iz]-1.5*_data[ix][iy][iz+1]+0.5*_data[ix][iy][iz+2]; //a03=-.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3]
                    coeff3D[ix][iy-1][iz-1][4]=-0.5*_data[ix][iy-1][iz]+0.5*_data[ix][iy+1][iz]; //a10 = -.5*p[0][1] + .5*p[2][1];
                    coeff3D[ix][iy-1][iz-1][5]=0.25*_data[ix][iy-1][iz-1]-0.25*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy+1][iz-1]+0.25*_data[ix][iy+1][iz+1]; //a11 = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
                    coeff3D[ix][iy-1][iz-1][6]=-0.5*_data[ix][iy-1][iz-1]+1.25*_data[ix][iy-1][iz]-_data[ix][iy-1][iz+1]+0.25*_data[ix][iy-1][iz+2]+0.5*_data[ix][iy+1][iz-1]-1.25*_data[ix][iy+1][iz]+_data[ix][iy+1][iz+1]-0.25*_data[ix][iy+1][iz+2]; //a12 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3]
                    coeff3D[ix][iy-1][iz-1][7]=0.25*_data[ix][iy-1][iz-1]-0.75*_data[ix][iy-1][iz]+0.75*_data[ix][iy-1][iz+1]-0.25*_data[ix][iy-1][iz+2]-0.25*_data[ix][iy+1][iz-1]+0.75*_data[ix][iy+1][iz]-0.75*_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+1][iz+2]; //a13 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
                    coeff3D[ix][iy-1][iz-1][8]=_data[ix][iy-1][iz]-2.5*_data[ix][iy][iz]+2.*_data[ix][iy+1][iz]-0.5*_data[ix][iy+2][iz]; //a20 = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
                    coeff3D[ix][iy-1][iz-1][9]=-0.5*_data[ix][iy-1][iz-1]+0.5*_data[ix][iy-1][iz+1]+1.25*_data[ix][iy][iz-1]-1.25*_data[ix][iy][iz+1]-_data[ix][iy+1][iz-1]+_data[ix][iy+1][iz+1]+0.25*_data[ix][iy+2][iz-1]-0.25*_data[ix][iy+2][iz+1]; //a21=-.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
                    coeff3D[ix][iy-1][iz-1][10]= _data[ix][iy-1][iz-1] - 2.5*_data[ix][iy-1][iz] + 2.*_data[ix][iy-1][iz+1] - .5*_data[ix][iy-1][iz+2] - 2.5*_data[ix][iy][iz-1] + 6.25*_data[ix][iy][iz] - 5.*_data[ix][iy][iz+1] + 1.25*_data[ix][iy][iz+2] + 2.*_data[ix][iy+1][iz-1] - 5.*_data[ix][iy+1][iz] + 4.*_data[ix][iy+1][iz+1] - _data[ix][iy+1][iz+2] - .5*_data[ix][iy+2][iz-1] + 1.25*_data[ix][iy+2][iz] - _data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a22 = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
                    coeff3D[ix][iy-1][iz-1][11]= -.5*_data[ix][iy-1][iz-1] + 1.5*_data[ix][iy-1][iz] - 1.5*_data[ix][iy-1][iz+1] + .5*_data[ix][iy-1][iz+2] + 1.25*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.75*_data[ix][iy][iz+1] - 1.25*_data[ix][iy][iz+2] - _data[ix][iy+1][iz-1] + 3.*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + _data[ix][iy+1][iz+2] + .25*_data[ix][iy+2][iz-1] - .75*_data[ix][iy+2][iz] + .75*_data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a23 = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
                    coeff3D[ix][iy-1][iz-1][12]= -.5*_data[ix][iy-1][iz] + 1.5*_data[ix][iy][iz] - 1.5*_data[ix][iy+1][iz] + .5*_data[ix][iy+2][iz];//a30 = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
                    coeff3D[ix][iy-1][iz-1][13]= .25*_data[ix][iy-1][iz-1] - .25*_data[ix][iy-1][iz+1] - .75*_data[ix][iy][iz-1] + .75*_data[ix][iy][iz+1] + .75*_data[ix][iy+1][iz-1] - .75*_data[ix][iy+1][iz+1] - .25*_data[ix][iy+2][iz-1] + .25*_data[ix][iy+2][iz+1]; //a31 = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
                    coeff3D[ix][iy-1][iz-1][14]= -.5*_data[ix][iy-1][iz-1] + 1.25*_data[ix][iy-1][iz] - _data[ix][iy-1][iz+1] + .25*_data[ix][iy-1][iz+2] + 1.5*_data[ix][iy][iz-1] - 3.75*_data[ix][iy][iz] + 3.*_data[ix][iy][iz+1] - .75*_data[ix][iy][iz+2] - 1.5*_data[ix][iy+1][iz-1] + 3.75*_data[ix][iy+1][iz] - 3.*_data[ix][iy+1][iz+1] + .75*_data[ix][iy+1][iz+2] + .5*_data[ix][iy+2][iz-1] - 1.25*_data[ix][iy+2][iz] + _data[ix][iy+2][iz+1] - .25*_data[ix][iy+2][iz+2]; //a32 = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
                    coeff3D[ix][iy-1][iz-1][15]= .25*_data[ix][iy-1][iz-1] - .75*_data[ix][iy-1][iz] + .75*_data[ix][iy-1][iz+1] - .25*_data[ix][iy-1][iz+2] - .75*_data[ix][iy][iz-1] + 2.25*_data[ix][iy][iz] - 2.25*_data[ix][iy][iz+1] + .75*_data[ix][iy][iz+2] + .75*_data[ix][iy+1][iz-1] - 2.25*_data[ix][iy+1][iz] + 2.25*_data[ix][iy+1][iz+1] - .75*_data[ix][iy+1][iz+2] - .25*_data[ix][iy+2][iz-1] + .75*_data[ix][iy+2][iz] - .75*_data[ix][iy+2][iz+1] + .25*_data[ix][iy+2][iz+2]; //a33 = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];
                }
            }
        }
        if (iPrint=="all") printf("Done fitting 3D data...\n");

    }//-----3Ddata
    
}