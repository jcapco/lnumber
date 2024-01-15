#pragma warning( push )
#pragma warning(disable: 4018 4244) //unsigned int comparisons

#include <laman_number.h>
#include <omp.h>

using namespace std;

// Connected components of a graph
vector<vector<int> > components (vector<vector<int> >& e, vector<int> v)
{
  int i, j, p, a;
  vector<vector<int> > cmps(0);
  vector<int> cmp;

  while (v.size() > 0)
  {
    cmp.clear();
    cmp.push_back(v[0]);
    v.erase(v.begin());
    p = 0;

    //this is expensively done..
    while (p < cmp.size())
    {
      for (i = 0; i < e.size(); i++)
      {
        a = -1;
        if (e[i][0] == cmp[p]) a = e[i][1];
        if (e[i][1] == cmp[p]) a = e[i][0];
        if (a != -1)
        {
          for (j = 0; j < v.size(); j++) 
            if (a == v[j]) v.erase(v.begin() + j);
          for (j = 0; j < cmp.size(); j++) 
            if (cmp[j] == a) j = cmp.size();
          if (j == cmp.size()) cmp.push_back(a);
        }
      }
      p++;
    }
    cmps.push_back(cmp);
  }

  return cmps;
}

// Quotient of a graph by a subset of edges
//x encodes information of edges to remove and edges to contract
void quotient_graph(vector<vector<int> >& e2, vector<vector<int> >& e, vector<int> v, long long int x)
{

  int i, j, k, n, s;
  int nv = v.size();
  
  e2.clear();
  
  // e1 is the set of edges that will be removed
  // e2 is the set of edges that are kept
  vector<vector<int> > e1(0, vector<int>(2)); //, e2(0, vector<int>(2));
  for (i = 0; i < e.size(); i++)
  {
    if (((x >> i) & 1) == 0) e2.push_back(e[i]);
    else e1.push_back(e[i]);
  }
  k = e2.size();

  // use the connected components to create a renaming w of the vertices
  vector<vector<int> > cm = components(e1, v);
  unordered_map<int,int> w1;
  for (i = 0; i < nv; i++)
  {
    for (n = 0; n < cm.size(); n++)
    {
      s = (cm[n]).size();
      for (j = 0; j < s; j++)
      {
        if (v[i] == cm[n][j])
        {
          w1[v[i]]= n + 1; 
          j = s;
        }
      }
      if (j == s + 1) n = cm.size();
    }
  }

  // rename the vertices in the new graph
  for (i = 0; i < k; i++)
  {
    e2[i][0] = w1[e2[i][0]];
    e2[i][1] = w1[e2[i][1]];
  }

  // if we have a self-loop we return the empty graph
  for (i = 0; i < k; i++)
  {
    if (e2[i][0] == e2[i][1]) i = k;
  }
  if (i == k + 1)
    e2.clear();
  
}

// Determine multiple edges in a graph
vector<vector<int> > multiple_edges(vector<vector<int> > e)
{
  int i, j;
  vector<vector<int> > me(0);
  vector<int> m;

  for (i = 0; i < e.size() - 1; i++)
  {
    if (e[i][0] != -1)
    {
      m.clear();
      m.push_back(i);
      for (j = i + 1; j < e.size(); j++)
      {
        if ((e[j][0] == e[i][0] && e[j][1] == e[i][1] && e[j][0] != -1) ||
            (e[j][0] == e[i][1] && e[j][1] == e[i][0] && e[j][0] != -1))
        {
          m.push_back(j);
          e[j][0] = -1;
        }
      }
      if (m.size() > 1) me.push_back(m);
    }
  }
  return me;
}

// Determine triangles in a graph
vector<vector<int> > triangles(vector<vector<int> >& e)
{
  int i, j1, j2, c;
  int k = e.size();
  vector<vector<int> > tr(0, vector<int>(3));

  for (i = 0; i < k - 2; i++)
  {
    for (j1 = i + 1; j1 < k; j1++)
    {
      c = -1;
      if (e[j1][0] == e[i][0]) c = e[j1][1];
      else if (e[j1][1] == e[i][0]) c = e[j1][0];
      if (c != -1)
      {
        for (j2 = i + 1; j2 < k; j2++)
        {
          if ((e[j2][0] == e[i][1] && e[j2][1] == c) || (e[j2][1] == e[i][1] && e[j2][0] == c))
          {
            int temp[]  = {i,j1,j2};
            tr.push_back(vector<int>(temp,temp+3));
          }
        }
      }
    }
  }

  return tr;
}

// Determine quadrilaterals in a graph
vector<vector<int> > quadrilaterals(vector<vector<int> >& e)
{
  int i, j1, j2, j3, c1, c2;
  int k = e.size();
  vector<vector<int> > qu(0, vector<int>(4));

  for (i = 0; i < k - 3; i++)
  {
    for (j1 = i + 1; j1 < k; j1++)
    {
      c1 = -1;
      if (e[j1][0] == e[i][0] && e[j1][1] != e[i][1]) c1 = e[j1][1];
      else if (e[j1][1] == e[i][0] && e[j1][0] != e[i][1]) c1 = e[j1][0];
      if (c1 != -1)
      {
        for (j2 = i + 1; j2 < k; j2++)
        {
          c2 = -1;
          if (e[j2][0] == e[i][1] && e[j2][1] != e[i][0]) c2 = e[j2][1];
          else if (e[j2][1] == e[i][1] && e[j2][0] != e[i][0]) c2 = e[j2][0];
          if (c2 != -1 && c1 != c2)
          {
            for (j3 = i + 1; j3 < k; j3++)
            {
              if ((e[j3][0] == c1 && e[j3][1] == c2) || (e[j3][0] == c2 && e[j3][1] == c1))
              {
                int temp[]  = {i,j1,j3,j2};
                qu.push_back(vector<int>(temp,temp+4));
              }
            }
          }
        }
      }
    }
  }

  return qu;
}

int laman_number(vector<vector<int> > e1, vector<vector<int> > e2, bool first)
{
  //needed to distinguish globals for parallel
  vector<long long int> mask1(0), mask2(0);
  vector<vector<long long int> > keep(0, vector<long long int>(5)), excl(0, vector<long long int>(5));
  int ln=0, k = e1.size();
  vector<int> v1 = vertices(e1), v2 = vertices(e2);
  long long int one = 1;

  bool t1, t2;
  int c, i, j, sx, l1;
  long long int x, y, xm;  
  vector<vector<int> > q1, q2, s1, s2;

  
  // initial conditions (base cases)
  if (k == 1)
    return 1;
  
  if (v1.size() - (components(e1, v1)).size() + v2.size() - (components(e2, v2)).size() != k + 1)
    return 0;

  // which special edge to choose (one that is involved in the fewest number of triangles)
  // this is a heuristic that seems to work quite well, but there is no guarantee of optimality
  vector<vector<int> > tr1 = triangles(e1), tr2 = triangles(e2);
  vector<int> cnt(k, 0);
  for (i = 0; i < tr1.size(); i++) for (j = 0; j < 3; j++) cnt[tr1[i][j]]++;
  for (i = 0; i < tr2.size(); i++) for (j = 0; j < 3; j++) cnt[tr2[i][j]]++;
  c = 0;
  j = cnt[0];
  for (i = 1; i < k; i++) 
    if (cnt[i] < j)
    {
      j = cnt[i];
      c = i;
    }
  // we exchange edges c and k-1, also in the already computed triangles
  (e1[c]).swap(e1[k - 1]);
  (e2[c]).swap(e2[k - 1]);
  for (i = 0; i < tr1.size(); i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (tr1[i][j] == c) tr1[i][j] = k - 1;
      else if (tr1[i][j] == k - 1) tr1[i][j] = c;
    }
  }
  for (i = 0; i < tr2.size(); i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (tr2[i][j] == c) tr2[i][j] = k - 1;
      else if (tr2[i][j] == k - 1) tr2[i][j] = c;
    }
  }

  // subsets x that we want to keep or discard
  // keep: if (x & mask) matches one of vals, then keep x, otherwise discard x.
  // excl: if (x & mask) matches one of vals, then discard x, otherwise keep x.
  vector<long long int> vals(5);  

  // generate conditions imposed by multiple edges
  vector<vector<int> > csg;
  for (c = 0; c < 2; c++)
  {
    if (c == 0) csg = multiple_edges(e1);
    else csg = multiple_edges(e2);
    for (i = 0; i < csg.size(); i++)
    {
      xm = 0;
      for (j = 0; j < (csg[i]).size(); j++) xm |= (one << csg[i][j]);
      if (((xm >> (k - 1)) & 1) == 0)
      {
        vals[0] = 2;
        vals[1] = 0;
        vals[2] = xm;
      }
      else
      {
        xm &= ~(one << (k - 1));
        vals[0] = 1;
        if (c == 0) vals[1] = xm;
        else vals[1] = 0;
      }
      mask1.push_back(xm);
      keep.push_back(vals);
    }
  }

  // generate conditions imposed by triangles
  for (c = 0; c < 2; c++)
  {
    if (c == 0) csg = tr1;
    else csg = tr2;
    for (i = 0; i < csg.size(); i++)
    {
      xm = (one << csg[i][0]) | (one << csg[i][1]) | (one << csg[i][2]);
      if (((xm >> (k - 1)) & 1) == 0)
      {
        vals[0] = 3;
        if (c == 0)
        {
          for (j = 0; j < 3; j++) vals[j + 1] = (one << csg[i][j]);
        }
        else
        {
          for (j = 0; j < 3; j++) vals[j + 1] = xm ^ (one << csg[i][j]);
        }
      }
      else
      {
        xm &= ~(one << (k - 1));
        vals[0] = 1;
        if (c == 0) vals[1] = 0;
        else vals[1] = xm;
      }
      mask2.push_back(xm);
      excl.push_back(vals);
    }
  }

  // generate conditions imposed by quadrilaterals
  for (c = 0; c < 2; c++)
  {
    if (c == 0) csg = quadrilaterals(e1);
    else csg = quadrilaterals(e2);
    for (i = 0; i < csg.size(); i++)
    {
      xm = (one << csg[i][0]) | (one << csg[i][1]) | (one << csg[i][2]) | (one << csg[i][3]);
      if (((xm >> (k - 1)) & 1) == 0)
      {
        vals[0] = 4;
        if (c == 0)
        {
          for (j = 0; j < 4; j++) vals[j + 1] = (one << csg[i][j]);
        }
        else
        {
          for (j = 0; j < 4; j++) vals[j + 1] = xm ^ (one << csg[i][j]);
        }
      }
      else
      {
        xm &= ~(one << (k - 1));
        vals[0] = 1;
        if (c == 0) vals[1] = 0;
        else vals[1] = xm;
      }
      mask2.push_back(xm);
      excl.push_back(vals);
    }
  }

  // run through all nontrivial subsets x of the first k-1 edges
  t1 = true;
  
  long long int N=0;
  if (first) N=(one << (k - 2));
  else N=(one << (k - 1))-1;

  #pragma omp parallel for reduction(+:ln) private(i,j,x,xm, t1,y,sx,q1,q2,s1,s2,l1) 
  for (x = 1; x < N; x++)
  {
    // filter out those subsets which lead to graphs with self-loops due to triangles etc.
    // test the "keep" cases    
    for (i = 0; i < mask1.size(); i++)
    {
      xm = x & mask1[i];
      for (j = 1; j <= keep[i][0]; j++)
      {
        if (xm == keep[i][j]) j = keep[i][0] + 1;
      }
      if (j == keep[i][0] + 1) i = mask1.size();
    }
    
    if (i!=mask1.size()) continue;
    for (i = 0; i < mask2.size(); i++)
    {
      xm = x & mask2[i];
      for (j = 1; j <= excl[i][0]; j++)
      {
        if (xm == excl[i][j]) j = excl[i][0] + 1;
      }
      if (j == excl[i][0] + 2) i = mask2.size();
    }
    if (i!=mask2.size()) continue;
    
    //this is the sum of products part
    y = (~x) & ((one << (k - 1)) - 1);
    sx = 0;
    for (i = 0; i < k - 2; i++)
      if (((x >> i) & 1) == 1) sx++;
    
    t1 = true;
    if (2 * sx < k)
    {
      quotient_graph(q1,e1, v1, y);
      t1 = (q1.size() != 0);
      if (t1)
      {
        quotient_graph(q2,e2, v2, x);
        t1 = (q2.size() != 0);
      }
    }
    else
    {
      quotient_graph(q2, e2, v2, x);
      t1 = (q2.size() != 0);
      if (t1)
      {
        quotient_graph(q1, e1, v1, y);
        t1 = (q1.size() != 0);
      }
    }
    if (t1)
    {
      q1.pop_back();
      q2.pop_back();
      s1.clear();
      s2.clear();
      for (i = 0; i < k - 1; i++)
      {
        if (((x >> i) & 1) == 0) s1.push_back(e1[i]);
        if (((y >> i) & 1) == 0) s2.push_back(e2[i]);
      }
      l1 = laman_number(q1, s2);
      if (l1 != 0) ln += l1 * laman_number(q2, s1);
    }

  }

  if (first) ln *= 2;
  // extreme case (test whether the special edge is a bridge)
  quotient_graph(q1, e1, v1, (one << (k - 1)) - 1);
  quotient_graph(q2, e2, v2, (one << (k - 1)) - 1);
  t1 = (q1.size() == 0);
  t2 = (q2.size() == 0);
  if (!(t1 && t2))
  {
    e1.pop_back();
    e2.pop_back();
    l1 = laman_number(e1, e2);
    if (!t1 && !t2) l1 *= 2;
    ln += l1;
  }

  return ln;
}

#pragma warning(pop)
