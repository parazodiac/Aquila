































#include <RcppEigen.h>
#include <queue>
#include <random>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

uint32_t key_selector(std::pair<uint32_t, std::vector<uint32_t>> pair){
  return pair.first;
};

// [[Rcpp::export]]
Eigen::MatrixXd getGamma(Eigen::MatrixXd &alpha, Eigen::MatrixXd &beta,
                         DataFrame &overlaps, size_t max_iteration ) {

  std::unordered_map<uint32_t, std::vector<uint32_t>> peaks2genes;
  std::unordered_map<uint32_t, std::vector<uint32_t>> genes2peaks;
  {
    IntegerVector olap_peaks = overlaps[0];
    IntegerVector olap_genes = overlaps[1];

    for (size_t i=0; i<olap_peaks.size(); i++) {
      uint32_t peak_idx = olap_peaks[i]-1;
      uint32_t gene_idx = olap_genes[i]-1;

      peaks2genes[peak_idx].emplace_back(gene_idx);
      genes2peaks[gene_idx].emplace_back(peak_idx);
    }
  } // extracting a map of peak to genes, genes to peaks

  size_t num_cells = alpha.cols();
  size_t num_peaks = alpha.rows();
  size_t num_genes = beta.rows();
  Eigen::MatrixXd gamma = Eigen::MatrixXd::Zero(num_peaks, num_genes);

  // preprocessing overlaps
  Rcout << "Found " << num_cells << " Cells\n"
        << "Found " << num_peaks << " Peaks\n"
        << "Found " << num_genes << " Genes\n"
        << "Found " << overlaps.nrows() << " Overlapping Regions\n"
        << std::endl << std::flush;

  std::random_device rd;
  std::mt19937 mt(rd());

  std::uniform_int_distribution<size_t> sample_cell(0, num_cells-1);
  std::uniform_int_distribution<size_t> sample_peak(0, num_peaks-1);
  std::uniform_int_distribution<size_t> sample_gene(0, num_genes-1);
  std::uniform_real_distribution<double> coin_toss(0.0, 1.0);

  std::pair<size_t, size_t> state = std::make_pair(sample_peak(mt),
                                                   sample_gene(mt));

  for(size_t it = 1; it <= max_iteration; it++) {
    if (it % 100000 == 0) {
      Rcout << "\rDone iterations: "<< it << std::flush;
    }

    // first sample a cell
    size_t given_cell = sample_cell(mt);

    // next: given gene and a cell, choose peak
    {
      size_t given_gene = state.second;

      double norm{0.0};
      std::vector<uint32_t> &possible_peaks = genes2peaks[given_gene];
      size_t num_possible_peaks = possible_peaks.size();
      if (num_possible_peaks == 1) {
        state.first = possible_peaks[0];
      } else {
        std::vector<double> counts(num_possible_peaks);
        for (size_t i=0; i<num_possible_peaks; ++i) {
          double count = alpha(possible_peaks[i], given_cell);
          counts[i] = count;
          norm += count;
        }

        if (norm == 0.0) {
          state.first = sample_peak(mt);
        } else {
          double running_sum {0.0};
          double coin_value = coin_toss(mt);
          for(size_t i=0; i<num_possible_peaks; ++i) {
            if (counts[i] == 0.0) { continue; }
            running_sum += (counts[i] / norm);
            if (running_sum > coin_value) {
              state.first = possible_peaks[i];
              break;
            }
          }
        } // end-else norm!=0
      } //end-else num_possible_peaks != 1
    } // end selecting peak

    {
      // next: given peak and a cell, choose genes
      size_t given_peak = state.first;

      double norm{0.0};
      std::vector<uint32_t> &possible_genes = peaks2genes[given_peak];
      size_t num_possible_genes = possible_genes.size();
      if (num_possible_genes == 1) {
        state.second = possible_genes[0];
      } else {
        std::vector<double> counts(num_possible_genes);
        for (size_t i=0; i<num_possible_genes; ++i) {
          double count = beta(possible_genes[i], given_cell);
          counts[i] = count;
          norm += count;
        }

        if (norm == 0.0) {
          state.second = sample_gene(mt);
        } else {
          double running_sum {0.0};
          double coin_value = coin_toss(mt);
          for(size_t i=0; i<num_possible_genes; ++i) {
            if (counts[i] == 0.0) { continue; }
            running_sum += (counts[i] / norm);
            if (running_sum > coin_value) {
              state.second = possible_genes[i];
              break;
            }
          }
        } // end-else norm!=0
      } // end-else num_possible_genes != 1
    } // end selecting gene

    gamma(state.first, state.second) += 1;
  }

  Rcout << std::endl << std::flush;
  return gamma;
}
