#include <Rcpp.h>
#include <iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include <deque>
#include <sstream>

using namespace Rcpp;

typedef double floattype;
typedef unsigned int uint;
using namespace std;

// Rcpp wrapper of Daniel Lemire's O(3n) envelope algorithim. From his libimproved library.
// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.

// [[Rcpp::export]]
NumericMatrix computeEnvelope_lb(NumericVector array, int constraint) {

      //Upper bound will be first column, lower bound the second.
      NumericMatrix mxmi(array.size(),2);

      uint width = 1+  2 * constraint;
        deque<int> maxfifo, minfifo;
        maxfifo.push_back(0);
        minfifo.push_back(0);
        for(uint i = 1; i < array.size(); ++i) {
          if(i >=constraint+1) {
            //maxvalues[i-constraint-1] = array[maxfifo.front()];
            //minvalues[i-constraint-1] = array[minfifo.front()];
            mxmi(i-constraint-1,0)= array[maxfifo.front()];      //Replaced above for Rcpp
            mxmi(i-constraint-1,1)= array[minfifo.front()];
          }
          if(array[i] > array[i-1]) { //overshoot
            maxfifo.pop_back();
            while(maxfifo.size() > 0) {
              if(array[i] <= array[maxfifo.back()]) break;
              maxfifo.pop_back();
            }
          } else {
            minfifo.pop_back();
            while(minfifo.size() > 0) {
              if(array[i] >= array[minfifo.back()]) break;
              minfifo.pop_back();
            }
          }
          maxfifo.push_back(i);
          minfifo.push_back(i);
         if(i==  width+maxfifo.front()) maxfifo.pop_front();
         else if(i==  width+minfifo.front()) minfifo.pop_front();
        }
        for(uint i = array.size(); i <= array.size() + constraint; ++i) {
          	//maxvalues[i-constraint-1] = array[maxfifo.front()];
          	//minvalues[i-constraint-1] = array[minfifo.front()];
            mxmi(i-constraint-1,0)= array[maxfifo.front()];       //Replaced above for Rcpp
            mxmi(i-constraint-1,1)= array[minfifo.front()];

          	if(i-maxfifo.front() >= width) maxfifo.pop_front();
          	if(i-minfifo.front() >= width) minfifo.pop_front();
        }

        return(mxmi);
}


// [[Rcpp::export]]
double LB_Keogh_lb(NumericVector candidateArray, NumericVector upperBound, NumericVector lowerBound) {

  double bet  = 0.0;

  for(int i=0; i<candidateArray.size(); i++){

    if(candidateArray[i] > upperBound[i]) {

      bet = bet + candidateArray[i] - upperBound[i];

    } else {
      if(candidateArray[i] < lowerBound[i]){

          bet = bet + lowerBound[i] - candidateArray[i];

      }
    }

  }

  return(bet);
}


// [[Rcpp::export]]
double LB_Improved(NumericVector candidateArray, NumericVector queryArray, NumericVector ubQuery, NumericVector lbQuery, int windowSize) {

  //NumericVector candidateArrayPrime = candidateArray; //This should be a copy CAUSING PROBLEMS???????
  //NumericVector candidateArrayPrime = clone(candidateArray);
  const NumericVector  v_source(candidateArray);
  NumericVector candidateArrayPrime = Rcpp::clone(v_source);



  double bet = 0.0;      //Will eventually be lb_improved
  double lb_keogh = 0.0; //We compute this on the way

  for(int i=0; i<candidateArray.size(); i++) {

    if(candidateArray[i] > ubQuery[i]) {

      bet = bet + candidateArray[i] - ubQuery[i];
      candidateArrayPrime[i] = ubQuery[i];

    } else {
      if(candidateArray[i] < lbQuery[i]) {

        bet = bet + lbQuery[i] - candidateArray[i];
        candidateArrayPrime[i] = lbQuery[i];

      }
    }
  }
  //bet is LB_Keogh at this point.
  //Rcout << bet << endl;
  lb_keogh = bet;

  //Now do the second pass for LB_Improved, NOTE: this starts with bet == LB_Keogh
  //Compute the upper and lower bounds on the "new" array: proj of candidate array on query = candidateArrayPrime
  NumericMatrix ulb_on_proj = computeEnvelope_lb(candidateArrayPrime,windowSize); //upper and lower boundry on projection.

  for(int i=0; i<queryArray.size(); i++){
    if(queryArray[i] > ulb_on_proj(i,0) ) {
      bet = bet + queryArray[i] - ulb_on_proj(i,0);
    } else {
      if(queryArray[i] < ulb_on_proj(i,1)) {
        bet = bet + ulb_on_proj(i,1) - queryArray[i];
      }
    }
  }

  return(bet);

}
