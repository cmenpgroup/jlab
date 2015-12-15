#include "cpp_headers.hh"

const unsigned sLinesize = 1000;
vector< pair<string,string> > fAllLines;
int    fCurline;
string fCurkey;
string fCurtag;

bool parse_file(const char* filename)
{
  //cout<<"file : "<<filename<<"\t\t";
  ifstream in(filename);
  if( !in.good() )
  {
    cout<<endl<<"File \"" <<filename << "\" does not exist"<<endl<<endl;
    return false;
  }

  char data[sLinesize];
  char s[sLinesize];
  int pos=0;

  do
  {
    in.seekg(pos);
    in.getline(s,sLinesize);
    
    pos = in.tellg();     
 
    if(string(s).empty())
    {
      continue; // remove empty lines
    }
 
    istringstream lin(s);  

    string tag;
    lin>>tag;

    if(!strncmp(tag.c_str(),"//",2)) continue; // remove commented lines can be done better ...

    lin.get(data,sLinesize);
    
    fAllLines.push_back(pair<string, string>(tag, data));
  } while(in.good());
  
  if(in.eof())
  {
    return true;
  }
  else
  {
    cout<<"error"<<endl;
    return false;
  }
}

     
string get_line_data(const char* tag, const char* key)
{
  char data[sLinesize];
  bool found = false;
  for(unsigned i=0; i<fAllLines.size(); i++)
  {
    if(!fnmatch(fAllLines[i].first.c_str(), tag, 0))
    {
      istringstream in(fAllLines[i].second.c_str());
      string readkey;
      in>>readkey;

      if(readkey == key)
      {
	found=true;
	in.get(data,sLinesize);
      }
    }
  }
  if(found)
  {
    string str(data);
    int pos = str.find_last_not_of(" ");
    str.replace(pos+1, str.size(),"");
    return str;
  }
  else return string();
}


string get_next_line_data(const char* tag, const char* key)
{
  if(fCurtag != tag || fCurkey != key)
  {
    // not the same request
    fCurline = 0;
    fCurtag = tag;
    fCurkey = key;
  }
  // cout<<fCurline<<" "<<fCurtag<<" "<<fCurkey<<endl;

  char data[sLinesize];
  bool found = false;
  for(unsigned i=fCurline; i<fAllLines.size(); i++)
  {
    if( !fnmatch(fAllLines[i].first.c_str(), tag, 0) )
    {
      istringstream in(fAllLines[i].second.c_str());
      string readkey; in>>readkey;
      if(readkey == key)
      {
	found=true;
	in.get(data,sLinesize);
	fCurline=i+1;
	break;
      }
    }
  }
  if(found)
  {
    return string(data);
  }
  else
  {
    return string();
  }
}
