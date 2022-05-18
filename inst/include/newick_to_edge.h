#ifndef NEWICK_TO_EDGE_H
#define NEWICK_TO_EDGE_H

////////////////////////////////////////////////////////////////////////////////
////////  Code below was exported from the Castor R package
////////  found here: https://cran.r-project.org/web/packages/castor/index.html
////////  Following the GNU licencse, we have extracted the relevant parts,
////////  required to convert newick to an edge table.




#include <string>
#include <vector>

long count_occurrences( const std::string    &haystack,
                        const char            needle) {    // (INPUT) if true, only occurrences not bracketed by single or double quotes are counted
  long count = 0;
  bool open_single = false, open_double = false;
  for(long i=0; i<haystack.size(); ++i){
    if(haystack[i]==needle){
      if(((!open_single) && (!open_double))){
        ++count;
      }
    } else if((!open_single) && (haystack[i]=='"')){
      open_double = !open_double;
    } else if((!open_double) && (haystack[i]=='\'')){
      open_single = !open_single;
    }
  }
  return count;
}

bool extract_next_edge( const std::string     &input,
                        long                  &pointer,            // (INPUT/OUTPUT) will move towards the left
                        std::string           &child_name) {        // (OUTPUT) child name. Will be empty ("") if not available

  long left = -1, split=-1;
  for(long i=pointer; i>=0; --i){
    if (input[i]==':'){
      split = i;
    } else if ((input[i]=='(') || (input[i]==')') || (input[i]==',')){
      left = i;
      break;
    }
  }


  if (left < 0){
    //  error = "Missing terminal character '(', ')' or ','";
    return false;
  }

  if(left == pointer){
    // no child name nor edge information available
    child_name     = "";
    return true;
  }

  if(split < 0){
    // no edge information available, interpret whole specifier as child_name
    child_name = input.substr(left + 1, pointer - left);
  } else {
    child_name = input.substr(left + 1, split - left - 1);
  }
  pointer = left;
  return true;
}

void reindex_clades(    const long                Nclades,            // (INPUT) number of clades (tips or nodes) in the tree
                        const long                Nedges,                // (INPUT) number of edges in the tree
                        const std::vector<long>    &tree_edge,            // (INPUT) 2D array of size Nedges x 2, in row-major format
                        const bool                root_convention,    // (INPUT) If true, the root of the tree (if existent) is ensured to obtain the index = Ntips
                        long                    &Ntips,                // (OUTPUT) the inferred number of tips in the tree
                        long                    &Nnodes,            // (OUTPUT) the inferred number of nodes in the tree
                        std::vector<long>        &old2new_clade){    // (OUTPUT) 1D array of size Nclades, mapping old-->new clade indices
  // determine tips & nodes
  std::vector<bool> clade_is_tip(Nclades,true);
  for(long edge=0; edge<Nedges; ++edge){
    clade_is_tip[tree_edge[edge*2+0]] = false;
  }
  Ntips = Nnodes = 0;
  for(long clade=0; clade<Nclades; ++clade){
    if(clade_is_tip[clade]) ++Ntips;
    else ++Nnodes;
  }

  // re-index clades
  old2new_clade.resize(Nclades);
  long next_tip=0, next_node=1;
  for(long clade=0; clade<Nclades; ++clade){
    if(clade_is_tip[clade]) old2new_clade[clade] = (next_tip++);
    else old2new_clade[clade] = Ntips+(next_node++);
  }

  // make sure root is indexed = Ntips, if requested
  if(root_convention){
    std::vector<bool> clade_is_root(Nclades,true);
    for(long edge=0; edge<Nedges; ++edge){
      clade_is_root[tree_edge[edge*2+1]] = false;
    }
    long root = -1, occupier=-1;
    for(long clade=0; clade<Nclades; ++clade){
      if(clade_is_root[clade]) root = clade;
      if(old2new_clade[clade]==Ntips) occupier = clade; // clade currently re-indexed to Ntips
    }
    if(root>=0){
      // tree has a root, so correct index
      long temp_root_new         = old2new_clade[root];
      old2new_clade[root]     = Ntips;
      old2new_clade[occupier] = temp_root_new;
    }
  }
}

std::vector< long > newick_to_edge(const std::string& input) {
  // estimate number of tips, nodes & edges for pre-allocation purposes
  const long estimated_Nclades = count_occurrences(input, ',') + count_occurrences(input, ')');
  const long estimated_Nedges  = estimated_Nclades - 1;


  // pre-allocate space
  std::vector<std::string> clade_names, edge_names;
  std::vector< long > tree_edge;
  clade_names.reserve(estimated_Nclades);
  tree_edge.reserve(2*estimated_Nedges);

  // prepare auxiliary data structures
  std::vector<int> clade_stack; // keep track of which node we are currently in. clade_stack[n+1] is a child of clade_stack[n]
  long pointer = input.length()-1;
  if(input[pointer]==';') --pointer;
  std::string child_name;

  // read input left<--right
  while(pointer>=0) {
    if(clade_stack.empty() && (!clade_names.empty())){
      throw std::runtime_error("Tree appears to have multiple roots");
    }

    if(!extract_next_edge(input, pointer, child_name)){
      throw std::runtime_error("Invalid child specifier to the left of position");
    }

    clade_names.push_back(child_name);
    if(!clade_stack.empty()){ // if empty, clade is root.
      tree_edge.push_back(clade_stack.back());
      tree_edge.push_back(clade_names.size()-1);
    }

    if(input[pointer]==')'){
      // moving one level deeper, into a new child
      clade_stack.push_back(clade_names.size()-1);
      --pointer;
    } else if(input[pointer]=='('){
      // finished at this level, moving up to parents
      while((pointer>=0) && ((input[pointer]=='(') || std::isspace(input[pointer]))){
        if(input[pointer]=='('){

          if(clade_stack.empty()) throw std::runtime_error("Unbalanced parentheses");

          clade_stack.pop_back();
        }
        --pointer;
      }
      if((pointer>=0) && (input[pointer]==',')) --pointer;
      else if((pointer>=0) && (input[pointer]==')')) {
        throw std::runtime_error("Unexpected opening paranthesis");
      }
    }else{
      // more clades to be extracted at this level
      --pointer;
    }
  }

  // nothing left to parse, so check if we came back to level 0
  if(!clade_stack.empty()) {
    throw std::runtime_error("Unbalanced parentheses");
  }

  const long Nclades = clade_names.size();
  const long Nedges  = tree_edge.size() / 2;
  std::vector<long> old2new_clade;
  long Ntips, Nnodes;
  reindex_clades( Nclades,
                  Nedges,
                  tree_edge,
                  true,
                  Ntips,
                  Nnodes,
                  old2new_clade);

  for(long edge=0; edge<Nedges; ++edge){
    tree_edge[2*edge+0] = old2new_clade[tree_edge[2*edge+0]];
    tree_edge[2*edge+1] = old2new_clade[tree_edge[2*edge+1]];
  }

  for(auto& i : tree_edge) i++;

  return tree_edge;
}


#endif
