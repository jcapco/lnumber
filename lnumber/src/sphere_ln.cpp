#include <iostream>
#include <algorithm>

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable: 4146 4244 4800)
#include <mpirxx.h>
#pragma warning(pop)
#else
  #ifdef __cplusplus
    #include <gmpxx.h>
  #else
    #include <gmp.h>
  #endif
#endif

#include <laman_number.h>

using namespace std;

/*
computes all the Qi i=0,1,2,3,4 in the algo, consisting a vector of 4-tuples, everything is sorted!
Note: 
- You can return array of fixed size by using static within the function!
- our I is their I -{a,b} , our J is their J -{c,d}
- I is sorted!
- elts of Q are sorted!
*/
void Qi(vector<vector<size_t> > &Q01, vector<vector<size_t> > &Q34, const size_t* q0, vector<vector<size_t> >& Q, const set<size_t> &I, size_t dJ, size_t star, bool& cont)
{
  std::vector<size_t> q_inter;
  std::vector<size_t> q_diff;

  const size_t &a = q0[0];
  const size_t &b = q0[1];

  size_t s = 0;
  cont = false;
  for (vector<vector<size_t> >::iterator q = Q.begin(); q!=Q.end(); ++q)
  {
    set_intersection(q->begin(), q->end(), I.begin(), I.end(), std::back_inserter(q_inter));    
    if (find(q->begin(),q->end(),a)!=q->end())  q_inter.push_back(a); 
    if (find(q->begin(),q->end(),b)!=q->end())  q_inter.push_back(b);
    s = q_inter.size();

    if ( s == 2 || ((s==4 || s==3) && Q34.size()> I.size())
      || ((s==0 || s==1) && Q01.size()> dJ) )
    {
      cont = true;
      break;
    }

    if (s==0) Q01.push_back(*q);
    else if (s==4) Q34.push_back(*q);
    else if (s==1)
    {
      vector<size_t> qnew = *q;
      for (size_t i = 0; i<4;++i)
        if (qnew[i]==q_inter[0])
        {
          qnew[i]=star;
          break;
        } 
      sort(qnew.begin(), qnew.end()); 
      Q01.push_back(qnew);
    }
    else if (s==3)
    {
      vector<size_t> qnew = *q;
      sort(q_inter.begin(), q_inter.end());
      set_difference(q->begin(), q->end(), q_inter.begin(), q_inter.end(), std::back_inserter(q_diff));
      for (size_t i = 0; i<4;++i)
        if (qnew[i]==q_diff[0])
        {
          qnew[i]=star;
          break;
        }
      sort(qnew.begin(), qnew.end());
      Q34.push_back(qnew);
    }
    else //s==2
    {
      cont = true;
      break;
    }

    q_inter.clear();
    q_diff.clear();
    s=0;
  }

  if (Q34.size()!= I.size() || Q01.size()!= dJ) cont = true;
}

//our I is their L, return 1???
bool loop_pass(std::set<size_t>& I, const size_t cntr, void* data)
{
  sph_Data* dta = (sph_Data*)data;
  if (dta->half && cntr > 1U<<(dta->N.size()-1) ) return false;

  vector<vector<size_t> > Q01,Q34;
  size_t dJ = dta->N.size()-I.size();
  
  Qi(Q01,Q34, &(dta->q0[0]), dta->Q, I, dJ, dta->star, dta->cont);
  if (dta->cont) return true;
  
  size_t crI = 0, crJ = 0;  
  if (I.size()+3 <= 4) crI = 1;
  if (dJ+3 <= 4) crJ = 1;

  if (!crJ) //before manipulating I
  {
    std::vector<size_t> J; 
    set_difference(dta->N.begin(), dta->N.end(), I.begin(), I.end(), std::back_inserter(J));
    J.push_back(dta->q0[2]); J.push_back(dta->q0[3]); J.push_back(dta->star);
    sph_Data dataJ(&J[0],J.size(),Q01);
    crJ = sph_cnt_realizations(&dataJ);
  }

  if (!crI) //I can now be manipulated
  { 
    I.insert(dta->q0[0]); I.insert(dta->q0[1]); I.insert(dta->star); 
    sph_Data dataI(I,Q34);
    crI = sph_cnt_realizations(&dataI);
  }


  dta->sum = dta->sum + crI*crJ;
  return true;
}

//we assume data.N has size \ge 3!
size_t sph_cnt_realizations(sph_Data* data)
{
  if (data->N.size() <= 4) return 1;
  // if =5, return 1
  // if =6, return 2

  set<size_t>::iterator it = data->N.end();
  --it;
  data->star = (*it)+1;
  data->q0 = data->Q.back();
  data->Q.pop_back(); //Q is out of scope now.. damn!
  for (size_t i=0; i<4; ++i)
    data->N.erase((data->q0)[i]); 
  //data->N.erase(&(data->q0)[0],&(data->q0)[0]+4); 
  powerset_loop<size_t>(data->N,loop_pass,data);

  return data->sum; //struct data goes out of scope, so its destructor is automatically called!
}
