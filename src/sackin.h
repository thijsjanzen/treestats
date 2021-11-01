#ifndef sackin_h
#define sackin_h

#include <vector>

size_t find_parent(const std::vector< std::vector< float >>& ltable,
                   int focal_id,
                   int start_index) {

  for(int i = start_index; i >= 0; i--) {
    if (static_cast<int>(ltable[i][2]) == focal_id) {
      return i;
    }
  }

  if (start_index != ltable.size()) {
   // force_output("could not find parent, retrying\n");
    return find_parent(ltable, focal_id, ltable.size());
  } else {
  //  force_output("could not find parent at all\n");
    return -1; // trigger access violation --< update to throw
  }
}


size_t calc_sackin(const std::vector< std::vector< float >>& ltable) {
  std::vector< int > s_values(ltable.size(), 0);
  s_values[0] = 1;
  s_values[1] = 1;

  // ltable:
  // 0 = branching time // not used here
  // 1 = parent
  // 2 = id
  // 3 = extinct time // not used here
  for (size_t i = 2; i < ltable.size(); ++i) {
    int parent_id = ltable[i][1];
    int parent_index = find_parent(ltable, parent_id, i);
    s_values[parent_index]++;
    s_values[i] = s_values[parent_index];
  }
  // verified with R for correct values
  // compared with apTreeShape results, based on ltable
  // 24-09-2021
  return(std::accumulate(s_values.begin(), s_values.end(), 0));
}

#endif
