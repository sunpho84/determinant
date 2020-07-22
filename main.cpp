#include <iostream>
#include <vector>
#include <gmpxx.h>

using namespace std;

/// Iteration to compute the submatrix determinant
template <typename T>                       //  Fundamental type
T determinantIter(const T *m,               //< Matrix in an array format
		  vector<vector<int>> &s,   //< Indirect index used to access to submatrix
		  const int n,              //< Size of the submatrix
		  const int N)              //< Total matrix suze
{
  /// Reference to current line index
  const vector<int>& l=s[N-n];
  
  if(n==1)
    return m[l[0]];
  else
    {
      /// Next line swaps
      vector<int>& nl=s[N-n+1];
      
      /// Keep separate the result of even and odd lines
      T out[2]={0.0,0.0};
      
      //prepare submatrix index
      for(int i=0;i<n-1;i++) nl[i]=l[i+1];
      
      for(int p=0;p<n;p++)
	{ 
	  /// Compute submatrix determinant
	  const T inDet=
	    determinantIter(m+N,s,n-1,N);
	  
	  out[p%2]+=m[l[p]]*inDet;
	  
	  // Put back current index
	  nl[p]=l[p];
	}
      
      // Returns the difference between even and odd line minors
      return out[0]-out[1];
    }
}

//compute the determinant of a NxN matrix through a recursive formula
template <typename T>                  //  Fundamental type
T matrixDeterminant(const T *m,        //< Matrix in an array format
		    const int n)       //< Size of the submatrix
{
  /// Index used to refer to submatrix
  vector<vector<int>> l(n);
  for(int i=0;i<n;i++)
    l[i].resize(n-i+1);
  
  /// When the submatrix is the matrix itself, the index is referring to the actual line
  for(int i=0;i<=n;i++)
    l[0][i]=i;
  
  return determinantIter(m,l,n,n);
}

int main()
{
  using T=mpf_class;
  
  int mpfPrecision;
  cout<<"Please enter precision in bits: ";
  cin>>mpfPrecision;
  mpf_set_default_prec(mpfPrecision);
  
  int n;
  cout<<"Please enter matrix size: ";
  cin>>n;
  
  T* m=new T[n*n];
  cout<<"Please enter matrix (elements separated by spaces): ";
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      cin>>m[j+n*i];
  
  cout<<matrixDeterminant(m,n)<<endl;
  
  delete[] m;
  
  return 0;
}

