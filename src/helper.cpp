































#include <Rcpp.h>
#include <queue>
#include <random>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getGamma(NumericMatrix &alpha, NumericMatrix &beta,
                       DataFrame &overlaps, Int32 max_iteration ) {
  size_t num_cells = alpha.ncol();
  size_t num_peaks = alpha.nrow();
  size_t num_genes = beta.nrow();

  CharacterVector peaks = rownames(alpha);
  CharacterVector genes = rownames(beta);

  std::unordered_map<String, uint32_t> peaks_indices;
  for(auto peak: peaks) {
    uint32_t idx = peaks_indices.size();
    peaks_indices[peak] = idx;
  }

  std::unordered_map<String, uint32_t> genes_indices;
  for(auto gene: genes) {
    uint32_t idx = genes_indices.size();
    genes_indices[gene] = idx;
  }

  NumericMatrix gamma(num_peaks, num_genes);
  rownames(gamma) = peaks;
  colnames(gamma) = genes;

  // preprocessing dataframe to find overlaps
  CharacterVector olap_peaks = overlaps[0];
  CharacterVector olap_genes = overlaps[1];
  Rcout << "Found " << num_cells << " Cells\n"
        << "Found " << num_peaks << " Peaks\n"
        << "Found " << num_genes << " Genes\n"
        << "Found " << olap_peaks.size() << " Overlapping Regions\n"
        << std::endl << std::flush;

  std::unordered_map<uint32_t, std::vector<uint32_t>> peaks2genes;
  std::unordered_map<uint32_t, std::vector<uint32_t>> genes2peaks;
  for (size_t i=0; i<olap_genes.size(); i++) {
    uint32_t gene_idx = genes_indices[olap_genes[i]];
    uint32_t peak_idx = peaks_indices[olap_peaks[i]];

    peaks2genes[peak_idx].emplace_back(gene_idx);
    genes2peaks[gene_idx].emplace_back(peak_idx);
  }

  std::random_device rd;
  std::mt19937 mt(rd());

  std::uniform_int_distribution<size_t> sample_cell(0, num_cells-1);
  std::uniform_real_distribution<double> coin_toss(0.0, 1.0);

  std::pair<size_t, size_t> state = std::make_pair(0, 0);
  std::vector<double> probabilities(100000);
  for(size_t it = 1; it <= max_iteration; it++) {
    if (it % 100000 == 0) {
      Rcout << "\rDone iterations: "<< it << std::flush;
    }

    // first sample a cell
    size_t given_cell = sample_cell(mt);

    // next: given gene and a cell, choose peak
    size_t given_gene = state.first;
    auto possible_peaks = genes2peaks[given_gene];

    //for (size_t i=0; i<possible_peaks.size(); ++i) {
    //  probabilities[i] = alpha(possible_peaks[i], given_cell);
    //}

    //double norm = std::accumulate(probability_peaks.begin(),
    //                              probability_peaks.end(),
    //                              0.0);
    //std::transform(probability_peaks.begin(),
    //               probability_peaks.end(),
    //               probability_peaks.begin(),
    //               [&norm](double& c){return c/norm;});


    //// next: given peak and a cell, choose genes
    //size_t given_peak = state.second;
    //auto possible_genes = peaks2genes[given_peak];

    //std::vector<double> probability_genes(possible_genes.size());
    //for (size_t i=0; i<possible_genes.size(); ++i) {
    //  probability_genes[i] += alpha(possible_genes[i], given_cell);
    //}
    //norm = std::accumulate(probability_genes.begin(),
    //                              probability_genes.end(),
    //                              0.0);
    //std::transform(probability_genes.begin(),
    //               probability_genes.end(),
    //               probability_genes.begin(),
    //               [&norm](double& c){return c/norm;});
  }

  Rcout << std::endl << std::flush;
  return gamma;
}
