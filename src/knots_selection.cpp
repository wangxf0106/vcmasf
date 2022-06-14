//  Created by xufeiwang on 21/12/19.
#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <armadillo>
using namespace std;
using namespace Rcpp;
using namespace arma;

double log_res(float m, mat xx, colvec xy, double yy)
{
  if (det(xx) < 1e-10)
    return (log(yy));
  double res = yy - as_scalar(trans(xy) * solve(xx, xy));
  return (log(res) - log(m));
}

double compute_loss(mat xx, colvec xy, double yy, colvec beta) 
{
  return (yy - 2 * as_scalar(trans(xy) * beta) + as_scalar(trans(beta) * (xx * beta)));
}

double index(int n, int i, int j)
{
  int ind = (2*n-i+1)*i/2+j-i-1;
  return (ind);
}

void M_cal(NumericVector* F, vector<int>* S, int p, double gamma)
{ 
  NumericVector M(p+1,0.0);
  for (int i=1; i<p+1; i++)
  {
    M[i] = gamma +(*F)[i-1];
  }
  
  for (int i=2;i<p+1; i++)
  {
    for (int k=1; k<i; k++)
    {
      double temp = gamma + M[k] +(*F)[index(p,k,i)];
      if (temp > M[i])
      {
        M[i] = temp;
        (*S)[i] = k;
      }
    }
  }
  
  return;
}

// [[Rcpp::export]]
NumericVector knots_selection_cpp(NumericMatrix X, NumericVector y, int m, double lam0, NumericVector Knots, NumericVector u)
{
  int n = X.nrow();
  int l = X.ncol();
  
  NumericMatrix SXX(n+1, l*l);
  NumericMatrix SXy(n+1, l);
  NumericVector Syy(n+1);
  
  for (int j=0; j<l*l; j++)
    SXX(0, j) = 0;
  for (int j=0; j<l; j++)
    SXy(0, j) = 0;
  Syy[0] = 0;
  for (int i=1; i<n+1; i++)
  {
    SXy(i, 0) = SXy(i-1, 0) + y[i-1];
    Syy[i] = Syy[i-1] + y[i-1] * y[i-1];
    for (int j=0; j<l; j++) {
      SXy(i, j) = SXy(i-1, j) + X(i-1, j) * y[i-1];
      for (int k=0; k<l; k++) 
        SXX(i, j*l+k) = SXX(i-1, j*l+k) + X(i-1, j) * X(i-1, k);
    }
  }
  
  int p;
  if (m == 0)
    p = Knots.size() + 1;
  else  
    p = n / m;
  NumericVector num(p+1, 0.0);
  if (m == 0) {
    int ind = 0;
    for (int i=0; i<n; i++) {
      if (u[i] > Knots[ind]) {
        num[ind] = i;
        ind = ind + 1;
      }
      num[p] = n;
    }
  } else {
    int r = n - p * m;
    int ind = 0 ;
    for (int i=1;i<p+1;i++) {
      if (i <= r)
        ind += m+1;
      else
        ind += m;
      num[i] = ind;
    }
  }
  
  int len = p*(p+1)/2;
  NumericVector *ls = new NumericVector(len, 0.0);
  vector<int> *S = new vector<int>(p+1, 0);
  
  NumericMatrix SXX_simple(p+1, l*l);
  NumericMatrix SXy_simple(p+1, l);
  NumericVector Syy_simple(p+1);
  
  for (int i=1;i<p+1;i++) {
    int ind = int(num[i]);
    for (int k=0; k<l*l; k++)
      SXX_simple(i,k) = SXX(ind,k);
    for (int k=0; k<l; k++)
      SXy_simple(i,k) = SXy(ind,k);
    Syy_simple[i] = Syy[ind];
  }
  
  for (int i=0; i<p; i++) {
    for (int j=i+1; j<p+1; j++){
      mat xx(l, l);
      colvec xy(l);
      for (int k=0; k<l; k++) 
        xy[k] = SXy_simple(j, k) - SXy_simple(i, k);
      for (int k1=0; k1<l; k1++) 
        for (int k2=0;k2<l; k2++)
          xx(k1, k2) = SXX_simple(j, k1*l+k2) - SXX_simple(i, k1*l+k2);
      (*ls)[index(p,i,j)] = -(num[j]-num[i])/2.0*log_res(num[j]-num[i], xx, xy, Syy_simple[j] - Syy_simple[i]);
    }
  }
  double lam = -lam0*log(n)/2;
  M_cal(ls, S, p, lam);
  int slice = 1;
  int temp = (*S)[p];
  NumericVector knots(1);
  knots[0] = n;
  while (temp != 0)
  {
    knots.push_back(num[temp]);
    temp = (*S)[temp];
    slice += 1;
  }
  delete ls;	
  delete S;
  return (knots);
}