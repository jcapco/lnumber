#pragma once

#pragma warning( push )
#pragma warning(disable: 4018 4244) //unsigned int comparisons
#include <iostream>
#include <vector>
#include <set>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#ifndef _WIN32
using namespace std::tr1;
#else
using namespace boost;
#endif

int laman_number(std::vector<std::vector<int> > e1, std::vector<std::vector<int> > e2, bool first=false);

//convert symmetric matrix a[i,j] to flat symmetric no diagonal upper triangular 
//nxn-matrix by getting the index of [i,j] for i>j
inline int idx_flat(int i, int j)
{
  return int(i*(i-1)/2 + j);
}

// Vertices of a graph (given by its edges)
inline std::vector<int> vertices (std::vector<std::vector<int> >& e)
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

inline int laman_number(std::vector<std::vector<int> >& e) //edgelist
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

//spherical...
template<class A> inline void powerset_loop(const std::set<A>& a, bool (*f)(std::set<A>&,const size_t, void*), void *data)
{
  std::set<A> vec;
  size_t size = a.size();
  typename std::set<A>::const_iterator it = a.begin();

  for(size_t i = 0; i < 1U<<size; ++i) //1U=1 as unsigned int
  {
    for(size_t j = 0; j < size; ++j)
      if(i & 1<<j)
      {
        std::advance(it,j);
        vec.insert(*it);
        std::advance(it, -int(j));
      }
    if (!f(vec,i,data)) break;
    vec.clear();
  }
}

struct sph_Data
{
  //N,Q
  bool cont, half;
  std::vector<size_t> q0;
  size_t sum,star;
  std::set<size_t> N;
  std::vector<std::vector<size_t> > Q; 
  
  sph_Data(size_t* NN, size_t n, std::vector<std::vector<size_t> > QQ):
    cont(false), half(false), sum(0), star(100), N(NN,NN+n),Q(QQ)
  {}

    sph_Data(std::set<size_t> NN,std::vector<std::vector<size_t> > QQ):
    cont(false), half(false), sum(0),  star(100), N(NN),Q(QQ) 
  {}

  sph_Data(): 
    cont(false), half(false), sum(0), star(100)
  {};
};

bool sph_loop_pass(std::set<size_t>& I, const size_t cntr, void* data);

size_t sph_cnt_realizations(sph_Data* data);

#pragma warning(pop)