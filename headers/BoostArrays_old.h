#include <memory>

template <typename T>
class multivec1D {
    int _n1;
     unique_ptr<T[]> data;
    
    public:
    multivec1D(){ }
    multivec1D(int n1): _n1(n1) {
        data= make_unique<T[]>(_n1);
    }
    void resize(int n1) {
        _n1=n1;
        data= make_unique<T[]>(_n1);
    }
    int n1() const { return _n1; }

    void fill(T value){
        for(int i=0; i<_n1; i++)
        {
            data[i] = value;
        }
    }
    //T operator[](int i) { return data[i]; }
    T &operator[](int i) {return data[i]; }
};

template <typename T>
class multivec2D {
    int _n1;
    int _n2;
     unique_ptr<T[]> data;
    
    public:
    multivec2D(){}
    multivec2D(int n1, int n2): _n1(n1), _n2(n2) {
        data= make_unique<T[]>(_n1*_n2);
    }
    void resize(int n1, int n2) {
        _n1=n1;
        _n2=n2;
        data= make_unique<T[]>(_n1*_n2);
    }
    void fill(T value){
        for(int i=0; i<_n1; i++)
        {
            for(int j=0; j<_n2; j++)
                data[j+i*_n2] = value;
        }
    }
        

     unique_ptr<T[]> at(){return data;}
    int n1() const { return _n1; }
    int n2() const { return _n2; }
    T *operator[](int i) { return i * _n2 + data.get(); }
    T &operator()(int i, int j) {
      return data[i * _n2 + j];
    }
};

template <typename T>
class multivec3D {
    int _n1;
    int _n2;
    int _n3;
     unique_ptr<T[]> data1D;
     unique_ptr<T**[]> minidata;
    
    public:
    multivec3D() {} 
    multivec3D(int n1, int n2, int n3): _n1(n1), _n2(n2), _n3(n3) {
        data1D= make_unique<T[]>(_n1*_n2*_n3);
        //data1d.reset(new T [_n1*_n2*_n3];
        minidata= make_unique<T**[]>(_n1);
       minidata.reset(new T**[_n1]);
	 for (int i=0; i<_n1; i++) {
            minidata[i] = new T* [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = i*_n2*_n3 + j*_n3 + data1D.get();
            }
        }
    }
    void resize(int n1, int n2, int n3) {
        _n1=n1;
        _n2=n2;
        _n3=n3;
        data1D= make_unique<T[]>(_n1*_n2*_n3);
        minidata= make_unique<T**[]>(_n1);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T* [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = i*_n2*_n3 + j*_n3 + data1D.get();
            }
        }
    }


    void fill(T value){
        for(int i=0; i<_n1; i++)
        {
            for(int j=0; j<_n2; j++)
                for(int k=0; k<_n3; k++)
                    data1D[k+j*_n3+i*_n2*_n3] = value;
        }
    }
    ~multivec3D(void){
        if(minidata.get() != nullptr) for (int i=0; i<_n1; i++) {
            delete [] minidata[i];
        }
    }
    int n1() const { return _n1; }
    int n2() const { return _n2; }
    int n3() const { return _n3; }
    T **operator[](int i) { return minidata[i]; }
    T &operator()(int i, int j,int k) {
        return data1D[i * _n2*_n3 + j*_n3 + k];
    }
};

template <typename T>
class multivec4D {
    int _n1;
    int _n2;
    int _n3;
    int _n4;
     unique_ptr<T[]> data1D;
     unique_ptr<T***[]> minidata;
    
    public:
    multivec4D() { }
    multivec4D(int n1, int n2, int n3, int n4): _n1(n1), _n2(n2), _n3(n3), _n4(n4) {
        data1D= make_unique<T[]>(_n1*_n2*_n3*_n4);
        //data1D.reset(new T [_n1*_n2*_n3*_n4]);
	minidata= make_unique<T***[]>(_n1);
	minidata.reset(new T***[_n1]);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T** [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = new T* [_n3];
                for (int k=0; k<_n3; k++) {
                    minidata[i][j][k] = i*_n2*_n3*_n4 + j*_n3*_n4 + k*_n4 + data1D.get();
                }
            }
        }
    }
    void resize(int n1, int n2, int n3, int n4) {
        _n1=n1;
        _n2=n2;
        _n3=n3;
        _n4=n4;
        data1D= make_unique<T[]>(_n1*_n2*_n3*_n4);
        minidata= make_unique<T***[]>(_n1);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T** [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = new T* [_n3];
                for (int k=0; k<_n3; k++) {
                    minidata[i][j][k] = i*_n2*_n3*_n4 + j*_n3*_n4 + k*_n4 + data1D.get();
                }
            }
        }
    }

    void fill(T value){
        for(int i=0; i<_n1; i++)
        {
            for(int j=0; j<_n2; j++)
            {
                for(int k=0; k<_n3; k++)
                {
                    for(int l=0; l<_n4; l++)
                    data1D[l+k*_n4+j*_n3*_n4+i*_n2*_n3*_n4] = value;
                }
            }
        }
    }

    ~multivec4D(void){
       if(minidata.get() != nullptr) for (int i=0; i<_n1; i++) {
            for (int j=0; j<_n2; j++) {
                delete [] minidata[i][j];
            }
            delete [] minidata[i];
        }
    }
    int n1() const { return _n1; }
    int n2() const { return _n2; }
    int n3() const { return _n3; }
    int n4() const { return _n4; }
    T ***operator[](int i) { return minidata[i]; }
};

template <typename T>
class multivec5D {
    int _n1;
    int _n2;
    int _n3;
    int _n4;
    int _n5;
     unique_ptr<T[]> data1D;
     unique_ptr<T****[]> minidata;
    
    public:
    multivec5D() {}
    multivec5D(int n1, int n2, int n3, int n4, int n5): _n1(n1), _n2(n2), _n3(n3), _n4(n4), _n5(n5) {
        data1D= make_unique<T[]>(_n1*_n2*_n3*_n4*_n5);
        minidata= make_unique<T****[]>(_n1);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T*** [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = new T** [_n3];
                for (int k=0; k<_n3; k++) {
                    minidata[i][j][k] = new T* [_n4];
                    for (int l=0; l<_n4; l++) {
                        minidata[i][j][k][l] = i*_n2*_n3*_n4*_n5 + j*_n3*_n4*_n5 + k*_n4*_n5 + l*_n5 + data1D.get();
                    }
                }
            }
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5) {
        _n1=n1;
        _n2=n2;
        _n3=n3;
        _n4=n4;
        _n5=n5;
        data1D= make_unique<T[]>(_n1*_n2*_n3*_n4*_n5);
        minidata= make_unique<T****[]>(_n1);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T*** [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = new T** [_n3];
                for (int k=0; k<_n3; k++) {
                    minidata[i][j][k] = new T* [_n4];
                    for (int l=0; l<_n4; l++) {
                        minidata[i][j][k][l] = i*_n2*_n3*_n4*_n5 + j*_n3*_n4*_n5 + k*_n4*_n5 + l*_n5 + data1D.get();
                    }
                }
            }
        }
    }

    void fill(T value){
        for(int i=0; i<_n1; i++)
            for(int j=0; j<_n2; j++)
                for(int k=0; k<_n3; k++)
                    for(int l=0; l<_n4; l++)
                        for(int m=0; m<_n5; m++)
                            data1D[m+l*_n5+k*_n4*_n5+j*_n3*_n4*_n5+i*_n2*_n3*_n4*_n5] = value;
        }

    ~multivec5D(void){
	if(minidata.get()!=nullptr)
	for (int i=0; i<_n1; i++) {
            for (int j=0; j<_n2; j++) {
                for (int k=0; k<_n3; k++) {
                    delete [] minidata[i][j][k];
                }
                delete [] minidata[i][j];
            }
            delete [] minidata[i];
        }
    }
    int n1() const { return _n1; }
    int n2() const { return _n2; }
    int n3() const { return _n3; }
    int n4() const { return _n4; }
    int n5() const { return _n5; }
    T ****operator[](int i) { return minidata[i]; }
};

template <typename T>
class multivec6D {
    int _n1;
    int _n2;
    int _n3;
    int _n4;
    int _n5;
    int _n6;
     unique_ptr<T[]> data1D;
     unique_ptr<T*****[]> minidata;
    
    public:
    multivec6D() {}
    multivec6D(int n1, int n2, int n3, int n4, int n5, int n6): _n1(n1), _n2(n2), _n3(n3), _n4(n4), _n5(n5), _n6(n6) {
        data1D= make_unique<T[]>(_n1*_n2*_n3*_n4*_n5*_n6);
        minidata= make_unique<T*****[]>(_n1);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T**** [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = new T*** [_n3];
                for (int k=0; k<_n3; k++) {
                    minidata[i][j][k] = new T** [_n4];
                    for (int l=0; l<_n4; l++) {
                        minidata[i][j][k][l] = new T* [_n5];
                        for (int m=0; m<_n5; m++) {
                            minidata[i][j][k][l][m] = i*_n2*_n3*_n4*_n5*_n6 + j*_n3*_n4*_n5*_n6 + k*_n4*_n5*_n6 + l*_n5*_n6 + m*_n6 + data1D.get();
                        }
                    }
                }
            }
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6) {
        _n1=n1;
        _n2=n2;
        _n3=n3;
        _n4=n4;
        _n5=n5;
        _n6=n6;
        data1D= make_unique<T[]>(_n1*_n2*_n3*_n4*_n5*_n6);
        minidata= make_unique<T*****[]>(_n1);
        for (int i=0; i<_n1; i++) {
            minidata[i] = new T**** [_n2];
            for (int j=0; j<_n2; j++) {
                minidata[i][j] = new T*** [_n3];
                for (int k=0; k<_n3; k++) {
                    minidata[i][j][k] = new T** [_n4];
                    for (int l=0; l<_n4; l++) {
                        minidata[i][j][k][l] = new T* [_n5];
                        for (int m=0; m<_n5; m++) {
                            minidata[i][j][k][l][m] = i*_n2*_n3*_n4*_n5*_n6 + j*_n3*_n4*_n5*_n6 + k*_n4*_n5*_n6 + l*_n5*_n6 + m*_n6 + data1D.get();
                        }
                    }
                }
            }
        }
    }
    void fill(T value){
        for(int i=0; i<_n1; i++)
            for(int j=0; j<_n2; j++)
                for(int k=0; k<_n3; k++)
                    for(int l=0; l<_n4; l++)
                        for(int m=0; m<_n5; m++)
                            for(int n=0; n<_n6; n++)
                             data1D[n+m*_n6+l*_n5*_n6+k*_n4*_n5*_n6+j*_n3*_n4*_n5*_n6+i*_n2*_n3*_n4*_n5*_n6] = value;
        }

    ~multivec6D(void){
	if(minidata.get()!=nullptr)
	for (int i=0; i<_n1; i++) {
            for (int j=0; j<_n2; j++) {
                for (int k=0; k<_n3; k++) {
                    for (int l=0; l<_n4; l++) {
                        delete [] minidata[i][j][k][l];
                    }
                    delete [] minidata[i][j][k];
                }
                delete [] minidata[i][j];
            }
            delete [] minidata[i];
        }
    }
    int n1() const { return _n1; }
    int n2() const { return _n2; }
    int n3() const { return _n3; }
    int n4() const { return _n4; }
    int n5() const { return _n5; }
    int n6() const { return _n6; }
    T *****operator[](int i) { return minidata[i]; }
};
