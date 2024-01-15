#include <laman_number.h>

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable: 4146 4244 4800)
#include <mpirxx.h>
#pragma warning(pop)
#else
#include <cmath>
#include <gmpxx.h>
#endif

#include <fstream>
#include <sstream>

#include <omp.h>

using namespace std;

template<class T>
inline int int_ceil(T t)
{
  int i = (int)t; /* truncate */
  return i + ( i < t ); /* convert trunc to ceil */
}

inline std::vector<std::vector<int> > convert_to_edgelist(mpz_ptr nptr, int n=0)
{
  using namespace std;
  vector<vector<int> > out;

  for (int i=1; i<n;++i)
  {
    for (int j=0;j<i;++j)
    {
      if (mpz_tstbit(nptr,idx_flat(i,j)))
      {
        int temp[]={j,i};
        out.push_back(vector<int>(temp,temp+2));
      }
    }    
  }
  return out;
}

inline void sph_convert_to_data(mpz_ptr nptr, size_t n, sph_Data* data)
{  
  if (n>0)
  {
    data->N.insert(0); data->N.insert(n);
  }
  
  for (size_t i=1; i<n;++i)
  {
    data->N.insert(i);
    data->N.insert(i+n);
    for (size_t j=0;j<i;++j)
      if (mpz_tstbit(nptr,idx_flat(i,j)))
      {
        size_t temp[4]={j,i,j+n,i+n};
        (data->Q).push_back(vector<size_t>(temp,temp+4));
      }
  }
}

int usage()
{
  cout << "Usage: lnumber n   or   lnumber -s n \n  This will compute (recursively) the number of planar or spherical realizations \
of a generic laman graph encoded in n (an integer).\n";
  cout << "\nDescription of arguments:\n\
  -s: tells the program to compute the number spherical realizations, without this option the program computes the number of planar realizations.\n\
  n:  the upper half triangle (without diagonal) of the adjancy matrix of the laman graph is a sequence (reading the matrix in row-major order) of bits which \
is a binary number converted to a positive decimal integer n.\n\
  As an example consider the (unique) laman graph with four vertices, this can be encoded as n=decimal(11101_2)=61. \
To compute its number of planar realizations (laman number on the plane), we execute \"lnumber 61 \"\n\
  ";
  return 1;
}

inline void trim(std::string& s)
{
  size_t p = s.find_first_not_of(" \t");
  s.erase(0, p);

  p = s.find_last_not_of(" \t");
  if (std::string::npos != p)
    s.erase(p+1);
}

inline void read_list(std::vector<string>& lgraph, std::vector<size_t>& lverts)
{
  try
  {
    ifstream file("lnumber_list.txt");
    string str;
    getline(file, str); 
    lgraph.push_back(strtok(&str[0], " ,\t"));
    lverts.push_back(atoi(strtok(NULL, " (),\t")));

    //note: strtok is not thread-safe, but we don't need threading here
    while (std::getline(file, str))
    {
      trim(str); 
      if (str.empty()) continue;
      lgraph.push_back(strtok(&str[0], " ,\t"));
      lverts.push_back(atoi(strtok(NULL, " (),\t")));
    }
    file.close();
  }
  catch(...)
  {
    std::cout << "There was an reading lnumber_list.txt\n";
  }

}

//jcapco todo: do this with the new format!
inline size_t highest_lnumber(size_t no_verts, size_t start_index, 
  size_t end_index, size_t interval_parse, char* ifilec, char* ofilec)
{
  size_t out = 0;
  FILE *ifile=0;

  ifile = fopen(ifilec,"rb");
  if (ifile==NULL) 
  {
    cout << "Cannot open " << ifilec << endl;
    return 0;
  }

  mpz_class n;
  mpz_ptr nptr= n.get_mpz_t();
  
  size_t ln=0;
  vector<vector<int> > edge_list;

  size_t cnt = -1;
  ofstream ofile(ofilec, ios::out | ios::binary);
  while (mpz_inp_raw(nptr, ifile))
  {
    ++cnt;
    if (cnt<start_index) continue;  
    if (end_index!=0 && cnt>end_index) break;
    edge_list.clear(); ln = 0;
    edge_list = convert_to_edgelist(n.get_mpz_t(),no_verts);
    ln = laman_number(edge_list);
    ofile.write((char*)&(ln),sizeof(int32_t));
    if (ln > out) out = ln;    
    if (cnt%interval_parse==0)
    {
      cout << "Parsed index: " << cnt << endl;
      cout << "At the moment the highest laman number is: " << out << endl;
    }
  } //todo check uniqueness of highest laman number
  fclose(ifile);

  return out;
}

int main(int argc, char *argv[])
{    
  if (argc<2) return usage();
  vector<vector<int> > edge_list;
  mpz_class n(0);
  int nvertices=0, ln = 0;
  double time = 0;
  if (strcmp("h",argv[1])==0 && argc==7)
  {
    //lnumber h 13 0-100 10 laman_13.bin laman_out.bin, 10=interval, 0-100 index (0-0 = all), 13=vertex, 
    //laman_13.bin = bin file of laman graphs, laman_out.bin = bin output of laman number associated to the graph.
    size_t start_index = atoi(strtok(argv[3],"-\t"));
    size_t end_index = atoi(strtok(NULL,"-\t"));
    ln  = highest_lnumber(atoi(strtok(argv[2], " \t")),start_index,end_index,atoi(strtok(argv[4], " \t")), 
      argv[5], argv[6]);
    cout << "Highest Laman number is: " << ln << endl;
  }
  else if (argc==2 && strcmp("list",argv[1])==0)
  {
    vector<string> lgraph; vector<size_t> lverts;
    read_list(lgraph, lverts);
    for (size_t i=0; i<lgraph.size(); ++i)
    {
      ln=0; edge_list.clear();
      n.set_str(lgraph[i],10);
      edge_list = convert_to_edgelist(n.get_mpz_t(),lverts[i]);
      time = omp_get_wtime();
      ln = laman_number(edge_list);
      time = omp_get_wtime()-time;
      cout << "Laman number of " << n << ": " << ln << endl; //6180, 1.7sec in my laptop
      cout << "Elapsed time: " << time  << " seconds\n\n";
    }
  }
  else if (argc==2 || argc==3)
  { 
    if (strcmp("-s",argv[1])!=0 && argc == 3) return usage();

    if (argc==2)
      n.set_str(strtok(argv[1], " \t"),10);
    else
      n.set_str(strtok(argv[2], " \t"),10);

    nvertices = mpz_sizeinbase(n.get_mpz_t(),2)-1; //floor(log2(n))
    nvertices = int_ceil((1+std::sqrt(1+8.0f*nvertices))/2.0f);
    cout << n << " has " << nvertices << " vertices\n";

    time = omp_get_wtime();
    if (argc == 2)
    {
      edge_list = convert_to_edgelist(n.get_mpz_t(),nvertices);  
      ln = laman_number(edge_list);
    }
    else
    {
      sph_Data data;
      sph_convert_to_data(n.get_mpz_t(), nvertices, &data);
      data.half = true; //because of symmetry take half of the power set in the beginning of algo
      ln = 2*sph_cnt_realizations(&data);
    }
    time = omp_get_wtime()-time;
    
    if (argc == 2) cout << "Planar ";
    else cout << "Spherical  ";
    cout << "Laman number of " << n << ": " << ln << endl; //planar = 6180, for ca. 4.3sec. in my laptop
    cout << "Elapsed time: " << time << " seconds\n";  
    //test case 252590061719913632 12   
    //print_edgelist(edge_list);
  }

  else return usage();

  return 0;
}
