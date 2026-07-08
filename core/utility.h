#pragma once

//////////////// Eigen is provided via the compiler include path (see setup.py / third_party/eigen)
#include <Eigen/Eigen>	// use Eigen library
#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "trifaces.h"
#include "tri_tri_intersect_moeller.h"
#define PI 3.14159265358979323846
#define LARGE 1000
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::EigenSolver;
using Eigen::MatrixXi;
/////////////////// very low-level class and struct definitions
class index_cmp 
{
	public:
	index_cmp(const std::vector<double> arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  { 
    return arr[a] < arr[b];
  }

  const std::vector<double> arr;
};

struct vertex_info{
	vertex_info(int ix, int*rr, int* rVe, int*rr1, int*rr2, int*rr3)
	:n(ix),
	r(rr),
	indx(rVe),
	r1(rr1),
	r2(rr2),
	r3(rr3)
	{}
	int n;	//the number of edges the vertex is member of
	int* r;	//the indices of faces the vertex is member of
	int* indx; // indices of other vertices this vertex is connected to
	int* r1;
	int* r2;
	int* r3;
	void display(){std::cout<<n<<"\n";}
};
////////////////// low level utility functions including utilities for use with Eigen library
//
// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
    return rand() / double(RAND_MAX);
}
//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRand(double a, double b)
{
    return (b-a)*unifRand() + a;
}
// Generate a random integer between 1 and a given value.
// param n the largest value 
// return a uniform random value in [1,...,n]
long unifRand(long n)
{
    
    if (n < 0) n = -n;
    if (n==0) return 0;
    /* There is a slight error in that this code can produce a return value of n+1
    **
    **  return long(unifRand()*n) + 1;
    */
    //Fixed code
    long guard = (long) (unifRand() * n) +1;
    return (guard > n)? n : guard;
}
/******************************************************************************/
/* randn()
 * 
 * Normally (Gaussian) distributed random numbers, using the Box-Muller 
 * transformation.  This transformation takes two uniformly distributed deviates
 * within the unit circle, and transforms them into two independently 
 * distributed normal deviates.  Utilizes the internal rand() function; this can
 * easily be changed to use a better and faster RNG.
 * 
 * The parameters passed to the function are the mean and standard deviation of 
 * the desired distribution.  The default values used, when no arguments are 
 * passed, are 0 and 1 - the standard normal distribution.
 * 
 * 
 * Two functions are provided:
 * 
 * The first uses the so-called polar version of the B-M transformation, using
 * multiple calls to a uniform RNG to ensure the initial deviates are within the
 * unit circle.  This avoids making any costly trigonometric function calls.
 * 
 * The second makes only a single set of calls to the RNG, and calculates a 
 * position within the unit circle with two trigonometric function calls.
 * 
 * The polar version is generally superior in terms of speed; however, on some
 * systems, the optimization of the math libraries may result in better 
 * performance of the second.  Try it out on the target system to see which
 * works best for you.  On my test machine (Athlon 3800+), the non-trig version
 * runs at about 3x10^6 calls/s; while the trig version runs at about
 * 1.8x10^6 calls/s (-O2 optimization).
 * 
 * 
 * Example calls:
 * randn_notrig();	//returns normal deviate with mean=0.0, std. deviation=1.0
 * randn_notrig(5.2,3.0);	//returns deviate with mean=5.2, std. deviation=3.0
 * 
 * 
 * Dependencies - requires <cmath> for the sqrt(), sin(), and cos() calls, and a
 * #defined value for PI.
 */

/******************************************************************************/
//	"Polar" version without trigonometric calls
double randn_notrig(double mu=0.0, double sigma=1.0) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double polar, rsquared, var1, var2;
	
	//	If no deviate has been stored, the polar Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose pairs of uniformly distributed deviates, discarding those 
		//	that don't fall within the unit circle
		do {
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);
		
		//	calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);
		
		//	store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;
		
		//	return second deviate
		return var2*polar*sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}


/******************************************************************************/
//	Standard version with trigonometric calls
#define PI 3.14159265358979323846

double randn_trig(double mu=0.0, double sigma=1.0) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double dist, angle;
	
	//	If no deviate has been stored, the standard Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose a pair of uniformly distributed deviates, one for the
		//	distance and one for the angle, and perform transformations
		dist=sqrt( -2.0 * log(double(rand()) / double(RAND_MAX)) );
		angle=2.0 * PI * (double(rand()) / double(RAND_MAX));
		
		//	calculate and store first deviate and set flag
		storedDeviate=dist*cos(angle);
		deviateAvailable=true;
		
		//	calcaulate return second deviate
		return dist * sin(angle) * sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}
/////////////////////////////////////////////////////////////////
void rot_mx(double ang, int flag, MatrixXd &R)
{
	// populate R with a cosine rotation matrix (angles in radians)
	//flag determines the axis around which will be rotated
	// ang is an angle in radians
	R.resize(3,3);
	double ca = std::cos(ang);
	double sa = std::sin(ang);
	if(flag == 1)//    %% rotate around z
	{
		R(0,0) = ca;
		R(0,1) = -sa;
		R(0,2) = 0.0;

		R(1,0) = sa;
		R(1,1) = ca;
		R(1,2) = 0.0;

		R(2,0) = 0.0;
		R(2,1) = 0.0;
		R(2,2) = 1.0;
	}
	if(flag == 2)//    %% rotate around y
	{
		R(0,0) = ca;
		R(0,1) = 0.0;
		R(0,2) = sa;

		R(1,0) = 0.0;
		R(1,1) = 1.0;
		R(1,2) = 0.0;

		R(2,0) = -sa;
		R(2,1) = 0.0;
		R(2,2) = ca;
	}
	if(flag == 3)//    %% rotate around x
	{
		R(0,0) = 1.0;
		R(0,1) = 0.0;
		R(0,2) = 0.0;

		R(1,0) = 0.0;
		R(1,1) = ca;
		R(1,2) = -sa;

		R(2,0) = 0.0;
		R(2,1) = sa;
		R(2,2) = ca;
	}
}

inline int kk_round(double x) { return (floor(x + 0.5)); }
inline int isin(std::vector<int> vec, int v)
{
	for(int i=0;i<vec.size();i++)
	{	
		if(vec[i]==v){return 1;}
	}
	return 0;
}
void Tokenize(const std::string& str, std::vector<std::string>& tokens,const std::string& delimiters = " ")
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
double factorial ( int n )
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
int share_vert(int f00, int f01, int f02, int f10, int f11, int f12)
{
	int count = 0;
	if(f00==f10||f00==f11||f00==f12){count++;}
	if(f01==f10||f01==f11||f01==f12){count++;}
	if(f02==f10||f02==f11||f02==f12){count++;}
	return count;
}
void disp_vector(std::vector<double> *vec){	for (int i = 0;i<int(vec->size());i++){std::cout<<double(vec->operator[](i))<<std::endl;}std::cout<<std::endl;}
void disp_vector(std::vector<size_t> *vec){	for (int i = 0;i<int(vec->size());i++){std::cout<<size_t(vec->operator[](i))<<std::endl;}std::cout<<std::endl;}
void d2eig(double *F1, MatrixXd &f1){
	for(int i=0;i<f1.rows();i++){f1(i) = F1[i];}
}
void d2eig(int *F1, MatrixXi &f1){for(int i=0;i<f1.rows();i++){f1(i,0) = F1[i];}}
void export_eigen_mx(MatrixXd &mx, const char*fn)
{
  std::ofstream file(fn);
  if (file.is_open()) {file << mx << '\n';}
}
void reshape(MatrixXd &x,unsigned int const r, unsigned int const c )	// slow:copies to a temp matrix, but does it Matlab style
{
	MatrixXd temp = x;
	x.resize(r,c);
	int counter = 0;
	int rp, cp;
	for(int cix = 0;cix<(temp.cols());cix++)
		for(int rix = 0;rix<temp.rows();rix++)
		{
			rp = int(counter%r);
			cp = int(counter%c);
			x(rp, cp) = temp(rix,cix);
			counter++;
		}
}
void reorder(std::vector<double> & unordered, std::vector<size_t> const & index_map, std::vector<double> & ordered)
{
  // copy for the reorder according to index_map, because unsorted may also be
  // sorted
  std::vector<double> copy = unordered;
  ordered.resize(index_map.size());
  for(int i = 0; i<int(index_map.size());i++)
  {
    ordered[i] = copy[index_map[i]];
  }
}

void kk_sort(std::vector<double> & unsorted,std::vector<double> & sorted,std::vector<size_t> & index_map)
{
  // Original unsorted index map
  index_map.resize(unsorted.size());
  for(size_t i=0;i<unsorted.size();i++){index_map[i] = i;}	// like Matlab's colon operator a:b
  // Sort the index map, using unsorted for comparison
  std::sort(
    index_map.begin(), 
    index_map.end(), 
    index_cmp(unsorted));

  sorted.resize(unsorted.size());
  reorder(unsorted,index_map,sorted);

}
void kk_cross(MatrixXd &u, MatrixXd &v, MatrixXd &r)
{
	//u and v are assumed to be 3-vector(s). The columns being 3.
	// u = < a , b , c > and v = < d , e , f >
	// The cross product, noted by x, of the two vectors u and v given above is another vector w given by
	// w = u x v = < a , b , c > x < d , e , f > = < x , y , z >
	// with the components x, y and z given by:
	// x = b*f - c*e , y = c*d - a*f and z = a*e - b*d 
		r.resize(u.rows(),3);
	r.col(0) = u.col(1).array()*v.col(2).array()-u.col(2).array()*v.col(1).array();	// x
	r.col(1) = u.col(2).array()*v.col(0).array()-u.col(0).array()*v.col(2).array();	// y
	r.col(2) = u.col(0).array()*v.col(1).array()-u.col(1).array()*v.col(0).array();	// z
	//y = u(:,3).*v(:,1)-u(:,1).*v(:,3);
	//z = u(:,1).*v(:,2)-v(:,1).*u(:,2);
	//r = [x y z];
}
void kk_kron(MatrixXd &a, MatrixXd &b, MatrixXd &r)
{
	//%%% a and b are assumed to be nx3 arrays
	//%%% returns an nx9 array
	r.resize(a.rows(),9);
	r.col(0) = a.col(0).array()*b.col(0).array();
	r.col(1) = a.col(1).array()*b.col(0).array();
	r.col(2) = a.col(2).array()*b.col(0).array();

	r.col(3) = a.col(0).array()*b.col(1).array();
	r.col(4) = a.col(1).array()*b.col(1).array();
	r.col(5) = a.col(2).array()*b.col(1).array();

	r.col(6) = a.col(0).array()*b.col(2).array();
	r.col(7) = a.col(1).array()*b.col(2).array();
	r.col(8) = a.col(2).array()*b.col(2).array();

}
void kk_kron_fast(MatrixXd &a, MatrixXd &b, MatrixXd &r)
{
	//%%% a and b are assumed to be nx3 arrays
	//%%% returns an nx9 array
	//r.resize(a.rows(),9);
	r.col(0) = a.col(0).array()*b.col(0).array();
	r.col(1) = a.col(1).array()*b.col(0).array();
	r.col(2) = a.col(2).array()*b.col(0).array();

	r.col(3) = a.col(0).array()*b.col(1).array();
	r.col(4) = a.col(1).array()*b.col(1).array();
	r.col(5) = a.col(2).array()*b.col(1).array();

	r.col(6) = a.col(0).array()*b.col(2).array();
	r.col(7) = a.col(1).array()*b.col(2).array();
	r.col(8) = a.col(2).array()*b.col(2).array();

}

void kk_transpose(MatrixXd &a, MatrixXd &r)
{
	//%%% a is nx9
	//%%% returns nx9
	r.resize(a.rows(), a.cols());
	r.col(0) = a.col(0);
	r.col(1) = a.col(3);
	r.col(2) = a.col(6);
	r.col(3) = a.col(1);
	r.col(4) = a.col(4);
	r.col(5) = a.col(7);
	r.col(6) = a.col(2);
	r.col(7) = a.col(5);
	r.col(8) = a.col(8);
	//r = [a(:,1) a(:,4) a(:,7) a(:,2) a(:,5) a(:,8) a(:,3) a(:,6) a(:,9)];
}
void kk_transpose_fast(MatrixXd &a, MatrixXd &r)
{
	//%%% a is nx9
	//%%% returns nx9
	//r.resize(a.rows(), a.cols());
	r.col(0) = a.col(0);
	r.col(1) = a.col(3);
	r.col(2) = a.col(6);
	r.col(3) = a.col(1);
	r.col(4) = a.col(4);
	r.col(5) = a.col(7);
	r.col(6) = a.col(2);
	r.col(7) = a.col(5);
	r.col(8) = a.col(8);
	//r = [a(:,1) a(:,4) a(:,7) a(:,2) a(:,5) a(:,8) a(:,3) a(:,6) a(:,9)];
}
void kk_mx_mult(MatrixXd &a, MatrixXd &b, MatrixXd &r)
{
	//%%% a and b are assumed to be nx9 arrays -- we are vectorizing matrix
	//%%% multiplication.
	//%%% returns an nx9 array
	r.resize(a.rows(), a.cols());
	
	r.col(0) = a.col(0).array()*b.col(0).array() + a.col(3).array() * b.col(1).array() + a.col(6).array() * b.col(2).array();
	r.col(3) = a.col(0).array()*b.col(3).array() + a.col(3).array() * b.col(4).array() + a.col(6).array() * b.col(5).array();
	r.col(6) = a.col(0).array()*b.col(6).array() + a.col(3).array() * b.col(7).array() + a.col(6).array() * b.col(8).array();

	r.col(1) = a.col(1).array()*b.col(0).array() + a.col(4).array() * b.col(1).array() + a.col(7).array() * b.col(2).array();
	r.col(4) = a.col(1).array()*b.col(3).array() + a.col(4).array() * b.col(4).array() + a.col(7).array() * b.col(5).array();
	r.col(7) = a.col(1).array()*b.col(6).array() + a.col(4).array() * b.col(7).array() + a.col(7).array() * b.col(8).array();

	r.col(2) = a.col(2).array()*b.col(0).array() + a.col(5).array() * b.col(1).array() + a.col(8).array() * b.col(2).array();
	r.col(5) = a.col(2).array()*b.col(3).array() + a.col(5).array() * b.col(4).array() + a.col(8).array() * b.col(5).array();
	r.col(8) = a.col(2).array()*b.col(6).array() + a.col(5).array() * b.col(7).array() + a.col(8).array() * b.col(8).array();
}
void kk_mx_mult_fast(MatrixXd &a, MatrixXd &b, MatrixXd &r)
{
	//%%% a and b are assumed to be nx9 arrays -- we are vectorizing matrix
	//%%% multiplication.
	//%%% returns an nx9 array
	//r.resize(a.rows(), a.cols());
	
	r.col(0) = a.col(0).array()*b.col(0).array() + a.col(3).array() * b.col(1).array() + a.col(6).array() * b.col(2).array();
	r.col(3) = a.col(0).array()*b.col(3).array() + a.col(3).array() * b.col(4).array() + a.col(6).array() * b.col(5).array();
	r.col(6) = a.col(0).array()*b.col(6).array() + a.col(3).array() * b.col(7).array() + a.col(6).array() * b.col(8).array();

	r.col(1) = a.col(1).array()*b.col(0).array() + a.col(4).array() * b.col(1).array() + a.col(7).array() * b.col(2).array();
	r.col(4) = a.col(1).array()*b.col(3).array() + a.col(4).array() * b.col(4).array() + a.col(7).array() * b.col(5).array();
	r.col(7) = a.col(1).array()*b.col(6).array() + a.col(4).array() * b.col(7).array() + a.col(7).array() * b.col(8).array();

	r.col(2) = a.col(2).array()*b.col(0).array() + a.col(5).array() * b.col(1).array() + a.col(8).array() * b.col(2).array();
	r.col(5) = a.col(2).array()*b.col(3).array() + a.col(5).array() * b.col(4).array() + a.col(8).array() * b.col(5).array();
	r.col(8) = a.col(2).array()*b.col(6).array() + a.col(5).array() * b.col(7).array() + a.col(8).array() * b.col(8).array();
}


void sort_eig_vec(Eigen::MatrixXd &mx, Eigen::MatrixXi &indx)
{
	std::vector<size_t> i;		// indices vector
	std::vector<double> a;	// unsorted vector
	std::vector<double> b;	// sorted vector
	a.resize(mx.size());
	for(int ix = 0;ix<(mx.size());ix++){a[ix] = mx(ix,0);}// read values from mx into a
	kk_sort(a,b,i); // do the sorting to get the vector of indices i
	std::sort(mx.data(), mx.data() + mx.size());	//  sort mx 
	for(int ix = 0;ix<(indx.size());ix++){indx(ix,0)= i[ix];}// read values from mx into a

}
////////////////// GAUSSQUAD ///////////////////////////////
inline static void gaussquad(const double a, const double b, MatrixXd &xx, MatrixXd &ww)
// in all steps we are reproducing the custom Matlab function gaussquad.m as closely as possible
{
	int n = xx.rows();
	MatrixXd u(n-1,1);
	for(int ix=1;ix<n;ix++){u(ix-1) = ix/double(sqrt(double(4*ix*ix-1)));}
	MatrixXd A(n,n);
	A.fill(0.0);
	for(int ix = 0;ix<n-1;ix++){
		A(ix+1,ix) = u(ix); 
		A(ix,ix+1) = u(ix);
	}
	EigenSolver<MatrixXd> eig(A); // [vec val] = eig(A)
	xx = eig.pseudoEigenvalueMatrix().diagonal(); // diag(val)
	MatrixXd vec = eig.pseudoEigenvectors(); // eigenvector matrix
	MatrixXi indx(n,1);
	sort_eig_vec(xx, indx); // xx gets sorted
	for(int ix=0;ix<n;ix++){ww(ix,0) = 2*(vec(0,indx(ix,0))) * (vec(0,indx(ix,0)));}	// make a vector v with the elements of first row in val inserted using the order in indx
	MatrixXd ones(xx.rows(), xx.cols());ones.fill(1.0);
    xx = ((b-a)/2)*xx + (a+b)/2*ones;     // Transform base points X.// Linearly transform from [-1,1] to [a,b].
    ww = ((b-a)/2)*ww;               // Adjust weigths.
}
///////// meshgrid
void meshgrid(MatrixXd &x, MatrixXd &y, MatrixXd &X, MatrixXd &Y)
// duplicates Matlab's meshgrid functionality
// needs: check that x and y are vectors (in the Matlab sense)
{
	int xdim = x.size();
	int ydim = y.size();
	if(x.cols()==1){x.transposeInPlace();}
	if(y.rows()==1){y.transposeInPlace();}
	X.resize(xdim,ydim);
	Y.resize(ydim,xdim);
	for(int rix = 0;rix<ydim;rix++){X.row(rix) = x;}
	for(int cix = 0;cix<xdim;cix++){Y.col(cix) = y;}
}

////////////////////////////////// BASIS GENERATION FUNCTIONS ///////////////////////////////
double *pm_polynomial ( int mm, int n, int m, double x[] )
{//****************************************************************************80
//
//  Purpose:
//
//    PM_POLYNOMIAL evaluates the Legendre polynomials Pm(n,m,x).
//
//  Differential equation:
//
//    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
//
//  First terms:
//
//    M = 0  ( = Legendre polynomials of first kind P(N,X) )
//
//    Pm(0,0,x) =    1
//    Pm(1,0,x) =    1 X
//    Pm(2,0,x) = (  3 X^2 -   1)/2
//    Pm(3,0,x) = (  5 X^3 -   3 X)/2
//    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
//    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
//    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
//    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
//
//    M = 1
//
//    Pm(0,1,x) =   0
//    Pm(1,1,x) =   1 * SQRT(1-X^2)
//    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
//    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
//    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
//
//    M = 2
//
//    Pm(0,2,x) =   0
//    Pm(1,2,x) =   0
//    Pm(2,2,x) =   3 * (1-X^2)
//    Pm(3,2,x) =  15 * (1-X^2) * X
//    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
//
//    M = 3
//
//    Pm(0,3,x) =   0
//    Pm(1,3,x) =   0
//    Pm(2,3,x) =   0
//    Pm(3,3,x) =  15 * (1-X^2)^1.5
//    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
//
//    M = 4
//
//    Pm(0,4,x) =   0
//    Pm(1,4,x) =   0
//    Pm(2,4,x) =   0
//    Pm(3,4,x) =   0
//    Pm(4,4,x) = 105 * (1-X^2)^2
//
//  Recursion:
//
//    if N < M:
//      Pm(N,M,x) = 0
//    if N = M:
//      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
//      all the odd integers less than or equal to N.
//    if N = M+1:
//      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
//    if M+1 < N:
//      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the maximum first index of the Legendre
//    function, which must be at least 0.
//
//    Input, int M, the second index of the Legendre function,
//    which must be at least 0, and no greater than N.
//
//    Input, double X[MM], the point at which the function is to be
//    evaluated.
//
//    Output, double PM_POLYNOMIAL[MM*(N+1)], the function values.
//

  double fact;
  int i;
  int j;
  int k;
  double *v;

  v = new double[mm*(n+1)];

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }
//
//  J = M is the first nonzero function.
//
  if ( m <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+m*mm] = 1.0;
    }

    fact = 1.0;
    for ( k = 0; k < m; k++ )
    {
      for ( i = 0; i < mm; i++ )
      {
        v[i+m*mm] = - v[i+m*mm] * fact * sqrt ( 1.0 - x[i] * x[i] );
      }
      fact = fact + 2.0;
    }
  }
//
//  J = M + 1 is the second nonzero function.
//
  if ( m + 1 <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+(m+1)*mm] = x[i] * ( double ) ( 2 * m + 1 ) * v[i+m*mm];
    }
  }
//
//  Now we use a three term recurrence.
//
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( double ) ( 2 * j     - 1 ) * x[i] * v[i+(j-1)*mm]
                  + ( double ) (   - j - m + 1 ) *        v[i+(j-2)*mm] )
                  / ( double ) (     j - m     );
    }
  }

  return v;
}

void legendre(int L, MatrixXd &ct, MatrixXd &pl)	//mimic Matlab's legendre function
{	// accepts a vector ct with elements between -1 and 1 (usually called as cos of the angle)
	// Output: pl of dimensions (Lo + 1) x (rowsxcols)
	// To do: pm_polynomial calculates other L values as well, which we throw away at the moment

	int dim = ct.cols()*ct.rows();
	if(ct.cols()>1 && ct.rows()>1) {reshape(ct,1,ct.rows()*ct.cols());}	// make sure ct and st are reshaped as row vectors
	if(ct.cols()==1){ct.transposeInPlace();}	// make sure ct is a row vector
	double* x = new double[dim];for (int i=0; i<dim; i++) { x[i] = ct(0,i);} // create and initialize to values in ct
	double* v ;//= new double[dim];for (int i=0; i<dim; i++) { v[i] = 0.0;}	// points to the array returned from
	pl.resize(L+1, ct.cols());		// make sure pl is the right size

	// loop over K values and fill in the output matrix pl
	for(int K=L;K>-1;K--)
	{
		int counter = 0;
		v = pm_polynomial (dim, L, K, x);	//call recursion formula
		//for (int i=dim; i<dim*2; i++) {std::cout<<v[i]<<std::endl;}
		for (int i=L*dim; i<dim*(L+1); i++) // we are only interested in the last section of v (see definition of output of pm_polynomial)
		{ 
			pl(K,counter) = v[i];			// fill in the values into our output matrix
			counter++;
		}
		delete []v;		// free the memory
	}
	delete []x;x = NULL;
}




void plkt(int L_max, MatrixXd &P, MatrixXd &P_T) // calculate the derivative of P w.r.t. theta and store in P_T
{

	MatrixXd tmp, t1, t2;
	double fac1, fac2;
	int ix;
	P_T.resize(P.rows(), P.cols());
	P_T.fill(0.0);
	int ia = 1;
	MatrixXd dpnm2 = P.col(1);
	for(int L = 1;L<=L_max;L++)
	{
		ia +=L;
		tmp = P.col(ia-1);
		P_T.col(ia-1) = (-1)*sqrt(double(ia-1))*P.col(ia);
//if(true){std::cout<<"(-1)*sqrt(double(ia-1))*P.col(ia)"<<std::endl<<(-1)*sqrt(double(ia-1))*P.col(ia)<<std::endl;}
		fac1 = sqrt(double(2.0*L*(L+1)));
		for(int K = 1;K<=(L-1);K++)
		{
			ix = ia+K;
			fac2 = sqrt(double((L-K)*(L+K+1)));
//if(true){std::cout<<"t1"<<std::endl<<fac1*tmp.array()<<std::endl;}
//if(true){std::cout<<"t2"<<std::endl<<(fac2*P.col(ix).array())<<std::endl;}
			t1 = (fac1*tmp.array());
			t2 = (fac2*P.col(ix).array());
			P_T.col(ix-1) = 0.5*(t1-t2);
//if(true){std::cout<<"0.5*(t1-t2)"<<std::endl<<0.5*(t1-t2)<<std::endl;}

//if(true){std::cout<<"P_T --- up to Lmax = 1"<<std::endl<<"length: "<<P_T.rows()<<std::endl<<P_T<<std::endl;}

			tmp = P.col(ix-1);
			fac1 = fac2;
		}
		P_T.col(ia+L-1) = sqrt(double(L)/2)*tmp;
//if(true){std::cout<<"P_T --- up to Lmax = 1"<<std::endl<<"length: "<<P_T.rows()<<std::endl<<P_T<<std::endl;}

	}
	P_T.col(2) = dpnm2;
//if(true){std::cout<<"P_T --- up to Lmax = 1"<<std::endl<<"length: "<<P_T.rows()<<std::endl<<P_T<<std::endl;}

}
double N_LK_bosh(int L, int K)
{
	K = std::abs(K);
	if(K>L){ return 0;}
	else{ return std::sqrt((2-int(K==0))*(2*L+1)*factorial(L-K)/factorial(L+K));}
}

void ylk_cos_sin_bosh(int L_max, MatrixXd &p, MatrixXd &t, MatrixXd &YLK, MatrixXd &PLK)
// generates a spherical harmonics basis according to the Bosh 2000 definition and normalization (as in custom Matlab code)
// Input: L_max, p and t are matrices of equal size
// Output:	YLK matrix dimensions [p.cols()xp.rows()] x (L_max + 1)^2 and holds the basis vectors
//			in the sequence 0,0   1,-1   1,0   1,1   2,-2   2,-1 ... etc.
//			PLK matrix of same dimensions as YLK, only holds the associated Legendre function values for use in derivative calculations
{
	MatrixXd ct, plk_mat, tmp, tmp2, tmp3, phi, theta;
	int counter = 0;
	int pcounter = 0;
	double NLK;
	int dim = p.rows()*p.cols();
	phi = p;
	reshape(phi,dim,1);
	theta = t;
	reshape(theta,dim,1);

	YLK.resize(dim, (L_max + 1)*(L_max + 1));		// make sure YLK has the right size
	double pdim = ((double(L_max)*double(L_max))/2) + (3*double(L_max)/2) + 1;
	PLK.resize(dim,int(pdim));  // will hold the associated Legendre function values
	ct = t.array().cos();										// precalculate cosine theta for passing to legendre
    for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
		legendre(L, ct, plk_mat);				// works just like Matlab's legendre
		for(int K = -L;K<=L;K++)						// loop over K
		{
            NLK = N_LK_bosh(L,K);

            if(K >=0)
			{
				tmp = K*phi;
				tmp2 = tmp.array().cos();
				tmp3 = plk_mat.row(K);
				PLK.col(pcounter) = NLK*(tmp3.transpose());
                YLK.col(counter) = PLK.col(pcounter).array()*tmp2.array();
				pcounter++;
			}
            else if(K<0)
			{
				tmp = std::abs(K)*phi;
				tmp2 = tmp.array().sin();
				tmp3 = plk_mat.row(std::abs(K));
				YLK.col(counter) = NLK*tmp3.transpose().array() * tmp2.array();
			}
			counter++;
		}
	}
}

void ylk_cos_sin_dphi_bosh(int L_max, MatrixXd &p, MatrixXd &t, MatrixXd &P, MatrixXd &Y_P)
// generates a spherical harmonics basis according to the Bosh 2000 definition and normalization (as in custom Matlab code)
// Input: L_max, p and t are matrices of equal size
// Output:	YLK matrix dimensions [p.cols()xp.rows()] x (L_max + 1)^2 and holds the basis vectors
//			in the sequence 0,0   1,-1   1,0   1,1   2,-2   2,-1 ... etc.
//			PLK matrix of same dimensions as YLK, only holds the associated Legendre function values for use in derivative calculations
{
	MatrixXd tmp, tmp2, tmp3, phi, theta;
	int counter = 0;
	int pindx;
	int dim = p.rows()*p.cols();
	phi = p;
	reshape(phi,dim,1);
	theta = t;
	reshape(theta,dim,1);
	Y_P.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_P has the right size
    for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
        for(int K = -L;K<=L;K++)						// loop over K
		{
            if(K==0)
			{
				Y_P.col(counter) = Y_P.col(counter).array()*0.0;
			}
			else if(K>0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = (K*phi).array().sin();
                Y_P.col(counter) = P.col(pindx).array()*(-1)*K*tmp.array();
			}
            else if(K<0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = std::abs(K)*phi;
				tmp2 = std::abs(K) * tmp.array().cos();
				tmp3 = P.col(pindx);
				Y_P.col(counter) = tmp3.array() * tmp2.array();
			}
			counter++;
		}
	}
}


void ylk_cos_sin_dphiphi_bosh(int L_max, MatrixXd &p, MatrixXd &t, MatrixXd &P, MatrixXd &Y_PP)
// generates a spherical harmonics basis according to the Bosh 2000 definition and normalization (as in custom Matlab code)
// Input: L_max, p and t are matrices of equal size
// Output:	YLK matrix dimensions [p.cols()xp.rows()] x (L_max + 1)^2 and holds the basis vectors
//			in the sequence 0,0   1,-1   1,0   1,1   2,-2   2,-1 ... etc.
//			PLK matrix of same dimensions as YLK, only holds the associated Legendre function values for use in derivative calculations
{
	MatrixXd tmp, tmp2, tmp3, phi, theta;
	int counter = 0;
	int pindx;
	int dim = p.rows()*p.cols();
	phi = p;
	reshape(phi,dim,1);
	theta = t;
	reshape(theta,dim,1);
	Y_PP.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_P has the right size
    for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
        for(int K = -L;K<=L;K++)						// loop over K
		{
            if(K==0)
			{
				Y_PP.col(counter) = Y_PP.col(counter).array()*0.0;
			}
			else if(K>0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = (K*phi).array().cos();
                Y_PP.col(counter) = P.col(pindx).array()*(-1)*K*K*tmp.array();
			}
            else if(K<0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = std::abs(K)*phi;
				tmp2 = -K*K* tmp.array().sin();
				tmp3 = P.col(pindx);
				Y_PP.col(counter) = tmp3.array() * tmp2.array();
			}
			counter++;
		}
	}
}
void ylk_cos_sin_dtheta_bosh(int L_max, MatrixXd &p, MatrixXd &t, MatrixXd &P, MatrixXd &Y_T, MatrixXd &P_T)
{
	MatrixXd tmp, tmp2, tmp3, phi, theta;
	int counter = 0;
	int pindx;
	int dim = p.rows()*p.cols();
	phi = p;
	reshape(phi,dim,1);
	theta = t;
	reshape(theta,dim,1);
	Y_T.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_T has the right size
	
	plkt(L_max, P, P_T);							// calculate derivative P_T of associated Legendre polynomials

	for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
        for(int K = -L;K<=L;K++)						// loop over K
		{
			if(K>=0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = (K*phi).array().cos();
                Y_T.col(counter) = P_T.col(pindx).array()*tmp.array();
			}
            else if(K<0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = std::abs(K)*phi;
				tmp2 = tmp.array().sin();
				tmp3 = P_T.col(pindx);
				Y_T.col(counter) = tmp3.array() * tmp2.array();
			}
			counter++;
		}
	}
	/**/
}

void ylk_cos_sin_dthetaphi_bosh(int L_max, MatrixXd &p, MatrixXd &t, MatrixXd &P, MatrixXd &P_T, MatrixXd &Y_TP)
{
	MatrixXd tmp, tmp2, tmp3, phi, theta;
	int counter = 0;
	int pindx;
	int dim = p.rows()*p.cols();
	phi = p;
	reshape(phi,dim,1);
	theta = t;
	reshape(theta,dim,1);
	Y_TP.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_T has the right size

	for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
        for(int K = -L;K<=L;K++)						// loop over K
		{
			if(K>=0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = (K*phi).array().sin();
                Y_TP.col(counter) = P_T.col(pindx).array()*-K*tmp.array();
			}
            else if(K<0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = std::abs(K)*phi;
				tmp2 = std::abs(K) * tmp.array().cos();
				tmp3 = P_T.col(pindx);
				Y_TP.col(counter) = tmp3.array() * tmp2.array();
			}
			counter++;
		}
	}
}
void ylk_cos_sin_dthetatheta_bosh(int L_max, MatrixXd &p, MatrixXd &t, MatrixXd &P_T, MatrixXd &Y_TT)
{
	MatrixXd tmp, tmp2, tmp3, phi, theta, P_TT;
	int counter = 0;
	int pindx;
	int dim = p.rows()*p.cols();
	phi = p;
	reshape(phi,dim,1);
	theta = t;
	reshape(theta,dim,1);
	Y_TT.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_T has the right size // crashes when compiled as 32 bit with L_max 48 and 10242 points --- the matrix then is 
	
	plkt(L_max, P_T, P_TT);							// calculate derivative P_T of associated Legendre polynomials

	for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
        for(int K = -L;K<=L;K++)						// loop over K
		{
			if(K>=0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = (K*phi).array().cos();
                Y_TT.col(counter) = P_TT.col(pindx).array()*tmp.array();
			}
            else if(K<0)
			{
				pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
				tmp = std::abs(K)*phi;
				tmp2 = tmp.array().sin();
				tmp3 = P_TT.col(pindx);
				Y_TT.col(counter) = tmp3.array() * tmp2.array();
			}
			counter++;
		}
	}
	/**/
}

void tri_plane_intersect(std::vector< std::vector<double> >X, const std::vector< std::vector<int> > F,const double yval, std::vector<int> &indx)
	{	
		double U0[3];
		double U1[3];
		double U2[3];
		//double fac = 100.0;
		int counter = 0;
		// make a really big triangle
		U0[0] = -LARGE; U0[1] = yval; U0[2] = -LARGE;
		U1[0] =  0.0;	U1[1] = yval; U1[2] =  LARGE;
		U2[0] =  LARGE; U2[1] = yval; U2[2] = -LARGE;
		for(int i = 0;i<F.size();i++)
		{
				if(tri_tri_intersect(&X[F[i][0]][0],&X[F[i][1]][0],&X[F[i][2]][0],U0,U1,U2)==1)
					//if(NoDivTriTriIsect(V0,V1,V2,U0,U1,U2)==1)
				{
					indx.push_back(i);
					counter++;
				}
		}
		
		indx.resize(counter);
	}
void tri_cut_plane(double **X,const int n_points, int*F0, int*F1, int*F2,const int n_faces, double *D, const std::vector< std::vector<double> > &N, const double yval, std::vector< std::vector<double> > &Xin, std::vector< std::vector<double> > &Xout,std::vector< std::vector<double> > &Xcut, std::vector< std::vector<int> > &Fc,std::vector< std::vector<int> > &Fcp)
{		// generates the faces and vertices necessary to plot a cut-through view of a shell
		// Output:
		//vector<vector<double>> Xin, Xout;			// vertex positions of outer and inner shells after cutting
		//vector<vector<int>> Fc;		// the faces of the inner and outer shell cut surface (not the cut polygons of the cut plane)
		//vector<vector<double>> Xcut;			// vertices of the cut-plane through the thickness of the shell
		//vector<vector<int>> Fcp;	// faces of the polygons corresponding to the cut-plane through the thickness of the shell
		// Note: changes X (please make a copy of X in a future revision)
		double fac = 100.0;
		// copy points into the array Xin
		Xin.resize(n_points);
		for(int i=0;i<n_points;i++)
		{
			Xin[i].resize(3);
			Xin[i][0] = round(X[i][0]*fac)/fac;
			Xin[i][1] = round(X[i][1]*fac)/fac;
			Xin[i][2] = round(X[i][2]*fac)/fac;
		}
		// copy triangles into the array Fc
		Fc.resize(n_faces);
		for(int i = 0;i<n_faces;i++)
		{	
			Fc[i].resize(3);
			Fc[i][0] = F0[i]-1;
			Fc[i][1] = F1[i]-1;
			Fc[i][2] = F2[i]-1;
		}


		std::vector<int> indx;		//
		tri_plane_intersect(Xin,Fc,yval, indx);		// find the indices of triangles intersecting the plane y = yval
		//%% of the intersecting triangles, project the vertex that is closest to the plane yval onto yval
		for(int ix = 0; ix<indx.size(); ix++)        // % loop over intersecting triangles
		{
			//%%% find the vertex closest to the plane yval
			//tri = F(int_indx(ix),:);
			double dmin = LARGE;
			int minix = 0;
			for(int dix = 0;dix<3;dix++)
			{
				if(std::abs(Xin[Fc[indx[ix]][dix]][1]-yval)<dmin)
				{
					dmin = std::abs(Xin[Fc[indx[ix]][dix]][1]-yval);
					minix = dix;
				}
			}
			//%%% set the yvalue of that vertex to yval
			Xin[Fc[indx[ix]][minix]][1] = yval;
		}

		//%%  find and then delete all triangles with vertices of values larger than yval, except those that only have one value
		// those will have their y values projected on the plane yval
		std::vector<int> del_tri;
		std::vector<int> ylix_vec;
		for(int ix = 0; ix<Fc.size();ix++)	//   % loop over the triangles
		{
			  int cum = 0;
			  int ylix = 0;
			  for(int dix = 0;dix<3;dix++)
			  {
				if(Xin[Fc[ix][dix]][1]>yval)
				{
					cum++;
					ylix = dix;
				}
			  }
		     if (cum>1)
			 {
				del_tri.push_back(ix);
			 }
			 else 
			 {
				if(cum==1)			//%%% then the point with the one yval larger in the triangle is projected onto yval
				{
					ylix_vec.push_back(Fc[ix][ylix]);		// store the index of the vertex whose y value will be projected. Do not do the projection here!!!!
				}
			 }
		}
	// project the yvalues
		for(int ix = 0;ix<ylix_vec.size();ix++)
		{
			Xin[ylix_vec[ix]][1] = yval;
		}
	
		// delete the triangles with indices del_tri;///////////////////////////////////////// Inefficiency alert ///////////////////////
		for(int i = del_tri.size()-1;i>=0;i--){Fc.erase(Fc.begin() + del_tri[i]);}
		del_tri.clear();

		//%% clean the surface of triangles that are coplanar with the y plane y = yval
		for(int ix = 0; ix<Fc.size();ix++)  // % loop over the triangles
		{
			int cum = 0;
			for(int dix = 0;dix<3;dix++)
			{
				if(Xin[Fc[ix][dix]][1]== yval){ cum++;}
			}
			if(cum == 3){ del_tri.push_back(ix);}
		}
		// delete the triangles with indices del_tri;
		for(int i = del_tri.size()-1;i>=0;i--){Fc.erase(Fc.begin() + del_tri[i]);}
		del_tri.clear();
		
		// find all indices where Xin[i][1] == yval and store them in c_ix, which we will need later
		std::vector<int> c_ix;
		for(int ix = 0;ix<Xin.size();ix++)
		{
			if(Xin[ix][1]==yval){c_ix.push_back(ix);}
		}

		// now that we have a clean smooth cut through Xin we use the normals to make Xin's final configuration and also Xout
		// X_in = X + N .* d/2;  %%% now generate the inner surface
		Xout.resize(n_points);
		for(int i=0;i<n_points;i++)
		{
			Xout[i].resize(3);
			Xout[i][0] = Xin[i][0]-N[i][0]*D[i]/2;
			Xout[i][1] = Xin[i][1]-N[i][1]*D[i]/2;
			Xout[i][2] = Xin[i][2]-N[i][2]*D[i]/2;

			Xin[i][0] = Xin[i][0]+N[i][0]*D[i]/2;
			Xin[i][1] = Xin[i][1]+N[i][1]*D[i]/2;
			Xin[i][2] = Xin[i][2]+N[i][2]*D[i]/2;
		}
		// the displacement with normals also changes the smooth appearance of the cut surface
		// therefore we project the vertices on the edge to the yval again
		for(int ix = 0;ix<c_ix.size();ix++)
		{
			Xin[c_ix[ix]][1]=yval;
			Xout[c_ix[ix]][1]=yval;
		}

		//////////////////////////////////////////////////////////////////////////////
		// now calculate the fill-in space between the two surfaces (generate Xcut and Fcp)
		///////////////////////////////////////////////////////////////////////////////


		// we have to generate a fresh array of links L the vip data will not do.
		// use this until an "edge_info.m" equivalent has been written
		std::vector< std::vector<int> > L;		// This is where we store the links for every vertex in Xin
		L.resize(n_points);
		std::vector<int>vl;	//store the indices for a particular vertex
		for(int lix = 0;lix<L.size();lix++)	// in principle we loop over the number of vertices
		{
			int v0, v1, v2;		// current face index
			for(int fix = 0;fix<Fc.size();fix++)		// loop over all faces and find the neighbors of vertex lix
			{
				v0 = Fc[fix][0];
				v1 = Fc[fix][1];
				v2 = Fc[fix][2];
				if(v0==lix|v1==lix|v2==lix)		// if this face includes our vertex
				{
					if(v0==lix)
					{	// check if the other vertices are already in the vector vl, if they are not present, then include them
						if(isin(vl,v1)==0){vl.push_back(v1);}
						if(isin(vl,v2)==0){vl.push_back(v2);}
					}
					if(v1==lix)
					{	// check if the other vertices are already in the vector vl, if they are not present, then include them
						if(isin(vl,v0)==0){vl.push_back(v0);}
						if(isin(vl,v2)==0){vl.push_back(v2);}
					}
					if(v2==lix)
					{	// check if the other vertices are already in the vector vl, if they are not present, then include them
						if(isin(vl,v1)==0){vl.push_back(v1);}
						if(isin(vl,v0)==0){vl.push_back(v0);}
					}
				}
			}
			L[lix].resize(vl.size());
			for(int i = 0;i<vl.size();i++)
			{
				L[lix][i] = vl[i];
			}
			vl.clear();
		}

		///////////////////////////////////
		Fcp.reserve(c_ix.size() * 7);
		Xcut.reserve(c_ix.size() * 28);
		std::vector<double> p1(3,0.0);
		std::vector<double> p2(3,0.0);
		std::vector<double> p3(3,0.0);
		std::vector<double> p4(3,0.0);
		std::vector<int> poly_cell(4,0);
		std::vector<int>linx;	
		int counter = 0;
		for(int i = 0;i<c_ix.size();i++)	// loop over the vertices at the edge
		{	
			linx = L[c_ix[i]];
			if(linx.size())
			{
				//std::sort(linx.begin(), linx.end());
				// find the members of a polygon making up the fill-in surface

				for(int ix = 0;ix<3;ix++)		// copy the contents of a connected edge into p1 and p2
				{
					p1[ix] = Xin[c_ix[i]][ix];
					p2[ix] = Xout[c_ix[i]][ix];
				}
				// use the freshly determined L to find the other two members
				for(int lix = 0;lix<linx.size();lix++)
				{
					if(Xin[linx[lix]][1]==yval)		//i.e. if a neighbor also is on the edge yval then let us use that neighbor
					{
						p3[0] = Xin[linx[lix]][0];
						p3[1] = Xin[linx[lix]][1];
						p3[2] = Xin[linx[lix]][2];

						p4[0] = Xout[linx[lix]][0];
						p4[1] = Xout[linx[lix]][1];
						p4[2] = Xout[linx[lix]][2];
					
				
						// add the four points to Xcut and generate entries into Fcp
						// note that this will generate duplicate points -- but we don't care if vtk doesn't! --

						poly_cell[0] = counter*4;
						poly_cell[1] = counter*4+1;
						poly_cell[2] = counter*4+3;
						poly_cell[3] = counter*4+2;
						Fcp.push_back(poly_cell);
						Xcut.push_back(p1);
						Xcut.push_back(p2);
						Xcut.push_back(p3);
						Xcut.push_back(p4);
						counter++;
					}
				}
			}
		}
		Xcut.resize((counter-1)*4);
		Fcp.resize(counter-1);


}

void over_NLK(int L_max, std::vector<double> &onlk)
{
	int counter = 0;
	for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
		for(int K = -L;K<=L;K++)						// loop over K
		{
            onlk[counter] = 1.0/N_LK_bosh(L,K);
			counter++;
		}
	}
}
void times_NLK(int L_max, std::vector<double> &tnlk)
{
	int counter = 0;
	for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
	{
		for(int K = -L;K<=L;K++)						// loop over K
		{
            tnlk[counter] = N_LK_bosh(L,K);
			counter++;
		}
	}
}


//#endif