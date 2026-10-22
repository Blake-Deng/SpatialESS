#ifndef SPATIALESS_SPARSE_TYPE7_HPP
#define SPATIALESS_SPARSE_TYPE7_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace spatialess {

struct Type7Spec {
  int lower_rank;
  int upper_rank;
  double fraction;
};

inline std::array<Type7Spec, 3> trimean_specs(int total_count) {
  const std::array<double, 3> probabilities{{0.25, 0.50, 0.75}};
  std::array<Type7Spec, 3> specs;
  for (int i = 0; i < 3; ++i) {
    const double h = static_cast<double>(total_count - 1) * probabilities[i] + 1.0;
    const int lower = static_cast<int>(std::floor(h));
    specs[i] = Type7Spec{
      lower, std::min(total_count, lower + 1), h - static_cast<double>(lower)
    };
  }
  return specs;
}

inline double combine_trimean(const std::array<Type7Spec, 3>& specs,
                              int zero_count,
                              const std::array<int, 6>& selected_indices,
                              const std::array<double, 6>& selected_values,
                              int selected_count) {
  auto value_at = [&](int rank) {
    if (rank <= zero_count) return 0.0;
    const int index = rank - zero_count - 1;
    for (int i = 0; i < selected_count; ++i) {
      if (selected_indices[i] == index) return selected_values[i];
    }
    return 0.0;
  };
  std::array<double, 3> quantiles;
  for (int i = 0; i < 3; ++i) {
    const double lower = value_at(specs[i].lower_rank);
    quantiles[i] = lower + specs[i].fraction *
      (value_at(specs[i].upper_rank) - lower);
  }
  return (quantiles[0] + 2.0 * quantiles[1] + quantiles[2]) / 4.0;
}

inline int required_positive_indices(const std::array<Type7Spec, 3>& specs,
                                     int zero_count,
                                     std::array<int, 6>& indices) {
  int count = 0;
  for (const Type7Spec& spec : specs) {
    const std::array<int, 2> ranks{{spec.lower_rank, spec.upper_rank}};
    for (int rank : ranks) {
      if (rank > zero_count) indices[count++] = rank - zero_count - 1;
    }
  }
  std::sort(indices.begin(), indices.begin() + count);
  int unique_count = 0;
  for (int i = 0; i < count; ++i) {
    if (unique_count == 0 || indices[i] != indices[unique_count - 1]) {
      indices[unique_count++] = indices[i];
    }
  }
  return unique_count;
}

inline double sparse_type7_trimean_sort(std::vector<double>& positive,
                                        int total_count) {
  if (total_count <= 0 || positive.size() > static_cast<std::size_t>(total_count)) {
    return NA_REAL;
  }
  const int zero_count = total_count - static_cast<int>(positive.size());
  const std::array<Type7Spec, 3> specs = trimean_specs(total_count);
  std::array<int, 6> indices{{0, 0, 0, 0, 0, 0}};
  const int count = required_positive_indices(specs, zero_count, indices);
  if (count == 0) return 0.0;
  std::sort(positive.begin(), positive.end());
  std::array<double, 6> values{{0, 0, 0, 0, 0, 0}};
  for (int i = 0; i < count; ++i) values[i] = positive[indices[i]];
  return combine_trimean(specs, zero_count, indices, values, count);
}

inline double sparse_type7_trimean_select(std::vector<double>& positive,
                                          int total_count) {
  if (total_count <= 0 || positive.size() > static_cast<std::size_t>(total_count)) {
    return NA_REAL;
  }
  const int zero_count = total_count - static_cast<int>(positive.size());
  const std::array<Type7Spec, 3> specs = trimean_specs(total_count);
  std::array<int, 6> indices{{0, 0, 0, 0, 0, 0}};
  const int count = required_positive_indices(specs, zero_count, indices);
  if (count == 0) return 0.0;

  std::array<double, 6> values{{0, 0, 0, 0, 0, 0}};
  int first = 0;
  for (int i = 0; i < count; ++i) {
    const int target = indices[i];
    std::nth_element(
      positive.begin() + first,
      positive.begin() + target,
      positive.end()
    );
    values[i] = positive[target];
    first = target + 1;
  }
  return combine_trimean(specs, zero_count, indices, values, count);
}

}  // namespace spatialess

#endif
