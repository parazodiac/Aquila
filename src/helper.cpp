#include <RcppEigen.h>
#include <queue>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// taken from https://www.techiedelight.com/replace-each-element-corresponding-rank-array/
// Function to replace each element of the array by its rank in the array
void transform(std::vector<double> &vec) {
  // construct a priority queue of pairs
  std::priority_queue<std::pair <double, int>> maxHeap;

  // push all input elements with their corresponding index in the priority queue
  for (size_t i = 0; i < vec.size(); i++) {
    maxHeap.push({ vec[i], i });
  }

  // get input size
  double rank = 1.0;

  // run until max-heap is empty
  while (!maxHeap.empty())
  {
    // take next maximum element from heap and replace its value
    // in the input vector with its corresponding rank
    vec[maxHeap.top().second] = rank;
    maxHeap.pop();

    // decrement rank for next maximum element
    rank += 1.0;
  }
}

double getSpearmanRho(Eigen::SparseMatrix<double> &atac,
                      Eigen::SparseMatrix<double> &rna,
                      size_t i, size_t j) {
  size_t num_col = atac.cols();

  std::vector<double> a(num_col);
  std::vector<double> b(num_col);
  for (size_t k=0; k<num_col; k++) {
    a[k] = atac.coeff(i, k);
    b[k] = rna.coeff(j, k);
  }

  transform(a);
  transform(b);

  double squaredRankSum {0.0};
  for (size_t k=0; k<num_col; k++) {
    squaredRankSum += pow(a[k] - b[k], 2);
  }

  return (1.0 - ( (6 * squaredRankSum) / (pow(num_col, 3) - num_col) ));
}

// [[Rcpp::export]]
Eigen::MatrixXd correlationMatrix(Eigen::SparseMatrix<double> atac, Eigen::SparseMatrix<double> rna) {
  size_t numAtacRows = atac.rows();
  size_t numRNARows = rna.rows();

  Eigen::MatrixXd gamma(numAtacRows, numRNARows);
  for (size_t i=0 ;i< numAtacRows; i++) {
    std::cout<< "\r Working on row #" <<i;
    for (size_t j=0; j<numRNARows; j++) {
      gamma(i, j) = getSpearmanRho(atac, rna, i, j);
    }
  }

  return gamma;
}
