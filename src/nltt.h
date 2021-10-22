
std::vector<float> create_normalized_brts(const std::vector<float>& v) {
  std::vector< float > output = v;
  std::sort(output.begin(), output.end());
  if (output.front() == 0.f) {
    // t goes from 0 to T
    float maxT = output.back();
    output.push_back(maxT);
    for (auto& i : output) {
      i = i / maxT;
    }
  } else {
    output.push_back(0.f);

    // we assume here that the brts go from -T to 0 ...
    auto min_brts = output.front();
    for (auto& i : output) {
      i = 1.f - i / min_brts;
    }
  }
  return output;
}

std::vector< float > create_normalized_lins(size_t num_lin) {
  std::vector< float > output(num_lin - 1);
  std::iota(output.begin(), output.end(), 2.f);
  auto max_num = output.back();
  output.push_back(max_num);
  for (auto& i : output) {
    i *= 1.f / max_num;
  }
  return(output);
}

int get_index(const std::vector<float>& brts,
              float tim) {

  auto index = std::lower_bound(brts.begin(), brts.end(), tim);
  if (index != brts.begin()) index--;

  return std::distance(brts.begin(), index);
}

float calc_nltt_from_data(const std::vector<float>& b1,
                          const std::vector<float>& b2,
                          const std::vector<float>& n1,
                          const std::vector<float>& n2,
                          const std::vector<float>& all_b) {

  float nltt = 0.f;
  for (size_t k = 1; k < all_b.size(); ++k) {
    float tim = all_b[k];
    auto index1 = get_index(b1, tim);
    auto index2 = get_index(b2, tim);

    auto num_lin1 = n1[index1];
    auto num_lin2 = n2[index2];
    auto diff_lin = num_lin1 - num_lin2;

    if (diff_lin < 0) diff_lin *= -1;
    float dt = all_b[k] - all_b[k-1];
    nltt += dt * diff_lin;
    // std::cerr << tim << " " << num_lin1 << " " << num_lin2 << " " << dt << "\n";
  }
  return nltt;
}

// please note that the branching times have to be from -T to 0
float calc_nltt(const std::vector<float>& v1,
                const std::vector<float>& v2) {

  std::vector< float > b_times_1 = create_normalized_brts(v1);
  std::vector< float > b_times_2 = create_normalized_brts(v2);

  std::vector< float > lin_1 = create_normalized_lins(v1.size());
  std::vector< float > lin_2 = create_normalized_lins(v2.size());

  std::vector< float > all_branching_times(b_times_1.size() + b_times_2.size());
  std::merge(b_times_1.begin(), b_times_1.end(),
             b_times_2.begin(), b_times_2.end(),
             all_branching_times.begin());
  // notice that this might introduce duplicate branching times, but this does not
  // affect nltt, as dt = 0.
  return calc_nltt_from_data(b_times_1, b_times_2, lin_1, lin_2, all_branching_times);
}
