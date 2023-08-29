#pragma once

#pragma warning( push )
#pragma warning(disable: 4018 4244) //unsigned int comparisons
#include <iostream>
#include <vector>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#ifndef _WIN32
using namespace std::tr1;
#else
using namespace boost;
#endif

int laman_number(std::vector<std::vector<int>> e1, std::vector<std::vector<int>> e2, bool first=false);

// Vertices of a graph (given by its edges)
inline std::vector<int> vertices (std::vector<std::vector<int>>& e)
{
  unordered_set<int> v;
  for (size_t i = 0; i < e.size(); i++)
  {
    v.insert(e[i][0]);
    v.insert(e[i][1]);
  }
  std::vector<int> v1(v.begin(), v.end());
  return v1;
}

/* //for debug
inline void print_edgelist(std::vector<std::vector<int>> e)
{
  using namespace std;
  for (size_t i=0; i<e.size(); ++i)
    cout << "{" << e[i][0] << "," << e[i][1] << "}, ";
  cout << endl << endl;
}
*/

inline int laman_number(std::vector<std::vector<int>>& e) //edgelist
{
  int i, j;
  // preprocessing: remove vertices of valency/degree 2.
  std::vector<int> v = vertices(e);
  int n = v.size(), ln = 1;
  std::vector<int> cnt(n);
  unordered_map<int,int> w;
  for (i = 0; i < n; i++)
    w[v[i]]=i;
  bool flag = true;
  //not important since it's done only once!
  while (flag && n > 3)
  {
    for (i = 0; i < n; i++) cnt[i] = 0;
    for (i = 0; i < e.size(); i++)
    {
      cnt[w[e[i][0]]]++;
      cnt[w[e[i][1]]]++;
    }
    flag = false;
    for (i = 0; i < n; i++)
    {
      if (cnt[i] == 2)
      {
        for (j = e.size() - 1; j >= 0; j--) if (e[j][0] == v[i] || e[j][1] == v[i]) e.erase(e.begin() + j);
        n--;
        ln *= 2;
        flag = true;
      }
    }
  }

  return ln * laman_number(e, e, true);
}


#pragma warning(pop)