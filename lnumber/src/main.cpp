#include <laman_number.h>
#include <omp.h>

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable: 4146 4244 4800)
#include <mpirxx.h>
#pragma warning(pop)
#else
#include <gmpxx.h>
#endif

#include <fstream>
#include <sstream>

namespace Time{
#include <utils/hr_time.h>
}

using namespace std;

//convert symmetric matrix a[i,j] to flat symmetric no diagonal upper triangular 
//nxn-matrix by getting the index of [i,j] for i>j
inline int idx_flat(int i, int j)
{
  return int(i*(i-1)/2 + j);
}

inline std::vector<std::vector<int>> convert_to_edgelist(mpz_ptr nptr, int n=0)
{
  using namespace std;
  vector<vector<int>> out;

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

int usage()
{
  cout << "Usage: lnumber n v\n  This will compute (recursively) the number of (planar) realizations \
of a generic laman graph encoded in n with v vertices.\n";
  cout << "\nDescription of arguments:\n\
  n:  the upper half triangle (without diagonal) of the adjancy matrix of the laman graph is a sequence (reading the matrix in row-major order) of bits which \
is a binary number converted to a positive decimal integer n.\n\
  v: the number (positive integer) of vertices of the laman graph \n\n\
  As an example consider the (unique) laman graph with four vertices, this can be encoded as n=decimal(11101_2)=61 with v=4. \
To compute its number of realizations (laman number), we execute \"lnumber 61 4\"\n\
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
  /*
  string filename;
  sstr << "./laman_d3/laman_" << no_verts << ".bin";
  ifile = fopen(sstr.str().c_str(),"rb");
  */
  ifile = fopen(ifilec,"rb");
  if (ifile==NULL) 
  {
    //cout << "Cannot open " << sstr.str() << endl;
    cout << "Cannot open " << ifilec << endl;
    return 0;
  }
  //sstr.str("");

  mpz_class n;
  mpz_ptr nptr= n.get_mpz_t();
  
  size_t ln=0;
  vector<vector<int>> edge_list;
  
  //TODO, check for uniqueness, use a for-loop for parallelization?.. nah..
  size_t cnt = -1;
  //stringstream sstr;
  //sstr << "lnumber_d3_" << end_index << ".bin";
  //ofstream ofile(sstr.str().c_str(), ios::out | ios::binary);
  ofstream ofile(ofilec, ios::out | ios::binary);
  //sstr.str("");
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
#pragma omp parallel //8 in typical laptop, 32 for big computers
  {
    cout << "Number of Threads: " << omp_get_num_threads() << endl;
  }

  if (argc<2) return usage();
  vector<vector<int>> edge_list;
  mpz_class n(0);
  int nvertices=0, ln = 0;
  double time = 0;
  Time::CStopWatch csw;


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
  else if (argc==3)
  { 
    n.set_str(strtok(argv[1], " \t"),10);
    nvertices= atoi(strtok(argv[2], " \t"));
    
    n.set_str(strtok(argv[1], " \t"),10);
    edge_list = convert_to_edgelist(n.get_mpz_t(),nvertices);
    //print_edgelist(edge_list);    
    csw.startTimer();
    ln = laman_number(edge_list);
    csw.stopTimer();
    time += csw.getElapsedTime();
    cout << "Laman number of " << n << ": " << ln << endl; //6180, for ca. 4.3sec. in my laptop
    cout << "Elapsed time: " << time << " seconds\n";  
    //test case 252590061719913632 12   
    //print_edgelist(edge_list);
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
      csw.startTimer();
      ln = laman_number(edge_list);
      csw.stopTimer();
      time = csw.getElapsedTime();
      cout << "Laman number of " << n << ": " << ln << endl; //6180, 1.7sec in my laptop
      cout << "Elapsed time: " << time << " seconds\n\n";
    }
  }
  else return usage();

  return 0;
}
