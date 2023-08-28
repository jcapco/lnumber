#include <lib.h>
#include <laman_number.h>

using namespace std;

static size_t bit[] = {020000000000,010000000000,04000000000,02000000000,
                 01000000000,0400000000,0200000000,0100000000,040000000,
                 020000000,010000000,04000000,02000000,01000000,0400000,
                 0200000,0100000,040000,020000,010000,04000,02000,01000,
                 0400,0200,0100,040,020,010,04,02,01};

inline std::vector<std::vector<int>> convert_to_edgelist(mpz_ptr nptr, int n=0)
{
  using namespace std;
  vector<vector<int>> out;

  for (int i=1; i<n;++i)
    for (int j=0;j<i;++j)
      if (mpz_tstbit(nptr,idx_flat(i,j)))
      {
        int temp[]={j,i};
        out.push_back(vector<int>(temp,temp+2));
      } 
  return out;
}

size_t LIBLNUMBER_LIBRARY_INTERFACE laman_number(char* graph, size_t verts)
{
  vector<vector<int>> edge_list;
  mpz_class n(graph, 10);
  edge_list = convert_to_edgelist(n.get_mpz_t(),verts);
  return laman_number(edge_list);
}

size_t LIBLNUMBER_LIBRARY_INTERFACE laman_number_nauty(unsigned int* g, int verts)
{
  std::vector<std::vector<int>> edge_list;
  for (int i=1; i<verts;++i)
    for (int j=0;j<i;++j)
      //if (g[i] & (1 << j))
      if (g[i] & bit[j])
      {
        int temp[]={j,i};
        edge_list.push_back(vector<int>(temp,temp+2));
      }
  return laman_number(edge_list);
}

size_t LIBLNUMBER_LIBRARY_INTERFACE laman_numbern(mpz_ptr nptr, size_t verts)
{
  vector<vector<int>> edge_list;
  edge_list = convert_to_edgelist(nptr,verts);
  return laman_number(edge_list);
}

