#ifndef COORDINATE
#define COORDINATE


double det(vec2d& mat)//works only if mat dim is 3x3 (SARRUS)
{
	double d = mat[0][0]*mat[1][1]*mat[2][2] +
			   mat[1][0]*mat[2][1]*mat[0][2] +
			   mat[2][0]*mat[0][1]*mat[1][2] -
			   mat[0][2]*mat[1][1]*mat[2][0] -
			   mat[1][2]*mat[2][1]*mat[0][0] -
			   mat[2][2]*mat[0][1]*mat[1][0];
	return d;	
}


class Coordinate 
{
	public:
		Coordinate() {cart.resize(3); cart.fill(0.); crys.resize(3); crys.fill(0.); };
		vec1d cart;
		vec1d crys;
				
		static void set_crys_to_cart(vec2d& matrix);
		static void set_cart_to_crys(vec2d& matrix);
		//Coordinate():  _M(), _MI() {  	for(int i=0; i<3; i++) {for(int j=0; j<3; j++) _M[i][j]=0.;}    for(int i=0; i<3; i++) {for(int j=0; j<3; j++) _MI[i][j]=0.; }}
		void setcart(double cx, double cy, double cz);
		void setcrys(double c1, double c2, double c3);

		Coordinate operator+(Coordinate& k2);
		Coordinate operator-(Coordinate& k2);

		static vec2d& matrix()
		{
			return _M;
		};
		static vec2d& invmatrix()
		{
			return _MI;
		}
		void getmat()
		{
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
					cout << _M[i][j] << "         ";
				cout << endl;
			}
		};
		void getmatinv()
		{
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
					cout << _MI[i][j] << "         ";
				cout << endl;
			}
		};
		static double getJ() {return Jacobian;}




        double dot(Coordinate&);
		double norm();
	private:
		static vec2d _M;//(3, vec1d(3,0.)); 
		static vec2d _MI;//(3, vec1d(3,0.));
		static double Jacobian;

};




vec2d Coordinate::_M;
vec2d Coordinate::_MI;
double Coordinate::Jacobian;

void Coordinate::set_crys_to_cart(vec2d& matrix)
{	
	//mat MM(3,3);
	_M.resize(3,3); _M.fill(0.); _MI.resize(3,3); _MI.fill(0.);
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			_M[i][j] = matrix[j][i];
			//MM(i,j) = _M[i][j];
			_MI[i][j] = matrix[i][j];
		}
	}
	//mat invMM(3,3);
	//invMM=inv(MM);
	//LU factorization
	std::array<int,3> IPIV ={0,1,2};
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, 
                   &(_MI[0][0]), 3, &IPIV[0]);  

	//inverse from LU
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, &(_MI[0][0]),
                    3, &IPIV[0]);
	//for(int i=0; i<3; i++)
	//{
	//	for(int j=0; j<3; j++)
	//	{	
	//		_MI[i][j] = invMM(i,j);
	//	}
	//}	

	//Jacobian=abs(det(MM));
	Jacobian=abs(det(_M));

}


void Coordinate::set_cart_to_crys(vec2d& matrix)
{
_M.resize(3,3); _M.fill(0.); _MI.resize(3,3); _MI.fill(0.);
	//mat MM(3, 3);
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)	
		{
			_MI[i][j] = matrix[i][j];
			_M[i][j] = matrix[i][j];
			//MM(i,j) = _MI[i][j];
		}
	}
	//LU factorization
	std::array<int,3> IPIV ={0,1,2};
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, 
                   &(_M[0][0]), 3, &IPIV[0]);  

	//inverse from LU
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, &(_M[0][0]),
                    3, &IPIV[0]);

	//mat invMM(3,3);
	//invMM=inv(MM);
	//for(int i=0; i<3; i++)
	//{
	//	for(int j=0; j<3; j++)
	//		_M[i][j] = invMM(i,j);
	//}	
	Jacobian=abs(det(_M));

}

void Coordinate::setcart(double cx, double cy, double cz)
{
	cart[0] = cx;
	cart[1] = cy;
	cart[2] = cz;

	crys[0] = _MI[0][0]*cart[0] + _MI[0][1]*cart[1] + _MI[0][2]*cart[2];
	crys[1] = _MI[1][0]*cart[0] + _MI[1][1]*cart[1] + _MI[1][2]*cart[2];
	crys[2] = _MI[2][0]*cart[0] + _MI[2][1]*cart[1] + _MI[2][2]*cart[2];
	
}


void Coordinate::setcrys(double c1, double c2, double c3)
{
	crys[0] = c1;
	crys[1] = c2;
	crys[2] = c3;	

	for(int i=0; i<3; i++)
	{	
		cart[i] = 0.;
		for (int j=0; j<3; j++){
			cart[i] +=_M[i][j]*crys[j];
			// cout << setprecision(17) << " _Mij " << _M[i][j] << " i: " << i << " j: " << j << endl;
		}      
	}

	
}




double Coordinate::dot(Coordinate& k2)
{
	return this->cart[0]*k2.cart[0]+this->cart[1]*k2.cart[1]+this->cart[2]*k2.cart[2];
}



double Coordinate::norm()
{
	return sqrt(dot(*this));
}


Coordinate Coordinate::operator+(Coordinate& k2)
{
	Coordinate sum;
	sum.setcart(this->cart[0]+k2.cart[0], this->cart[1]+k2.cart[1], this->cart[2]+k2.cart[2]);
	return sum;
}


Coordinate Coordinate::operator-(Coordinate& k2)
{
	Coordinate sum;
	sum.setcart(this->cart[0]-k2.cart[0], this->cart[1]-k2.cart[1], this->cart[2]-k2.cart[2]);
	return sum;
}





class Coord_R: public Coordinate
{
	private:
		static vec2d _M;//(3, vec1d(3,0.)); 
		static vec2d _MI;//(3, vec1d(3,0.));
		static double Jacobian;
};


class Coord_B: public Coordinate
{
	private:
		static vec2d _M;//(3, vec1d(3,0.)); 
		static vec2d _MI;//(3, vec1d(3,0.));
		static double Jacobian;	
};


#endif
