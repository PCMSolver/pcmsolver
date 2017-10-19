/** \file Getkw.cpp
 *
 * \date Jun 3, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 * \brief
 */

#include "Getkw.h"

#define TEST_ARRAY                                                                  \
  if (len > 1)                                                                      \
    std::cout << "Warning, invalid length of 1 for " << name << std::endl;

using namespace std;

Getkw::Getkw(const string file, bool _verbose, bool _strict)
    : verbose(_verbose), strict(_strict) {

  toplevel = 0;
  if (file.empty() != 0 || file.compare("stdin") == 0 ||
      file.compare("STDIN") == 0) {
    if (verbose) {
      cout << "Reading input from stdin " << endl;
    }
    toplevel = readSect(cin);
  } else {
    const char * fname = file.data();
    if (verbose)
      cout << "Opening input file, '" << file << "'" << endl;
    ifstream fis(fname);
    if (not fis) {
      THROW_GETKW("Open failed: " + file);
    }
    toplevel = readSect(fis);
  }
  cur = toplevel;
}

Getkw::Getkw() {
  verbose = false;
  strict = false;
  toplevel = 0;
  cur = 0;
}

Getkw::Getkw(const Getkw & kw) {
  verbose = kw.verbose;
  strict = kw.strict;
  file = kw.file;
  *toplevel = *kw.toplevel;
  cur = toplevel;
}

Getkw & Getkw::operator=(const Getkw & kw) {
  if (&kw == this) {
    return *this;
  }
  verbose = kw.verbose;
  strict = kw.strict;
  file = kw.file;
  if (toplevel == 0) {
    toplevel = new Section(*kw.toplevel);
  } else {
    *toplevel = *kw.toplevel;
  }
  cur = toplevel;
  sstack = stack<const Section *>();
  return *this;
}

Getkw::~Getkw() { delete toplevel; }

void Getkw::setStrict(bool flag) { strict = flag; }

void Getkw::setVerbose(bool flag) { verbose = flag; }

void Getkw::print() const { cout << &repr(cout) << endl; }

ostream & Getkw::repr(ostream & o) const {
  if (toplevel == 0) {
    o << "Getkw not yet initialized" << endl;
  } else {
    o << *toplevel;
  }
  return o;
}

template <class T> const T & Getkw::get(const string & path) const {
  if (cur == 0) {
    THROW_GETKW("Getkw has not been initialized!");
  }
  const Keyword<T> & key = cur->getKey<T>(path);
  return key.get();
}

template <class T> const Keyword<T> & Getkw::getKeyword(const string & path) const {
  if (cur == 0) {
    THROW_GETKW("Getkw has not been initialized!");
  }
  return cur->getKey<T>(path);
}

const Section & Getkw::getSect(const string & path) const {
  if (cur == 0) {
    THROW_GETKW("Getkw has not been initialized!");
  }
  return cur->getSect(path);
}

void Getkw::pushSection(const string & path) {
  if (cur == 0) {
    THROW_GETKW("Getkw has not been initialized!");
  }
  sstack.push(cur);
  try {
    const Section & newsec = cur->getSect(path);
    cur = &newsec;
  } catch (string err) {
    cout << err;
    if (strict)
      exit(1);
  }
}

void Getkw::popSection() {
  if (sstack.empty()) {
    cout << "Error! Getkw stack is empty!" << endl;
    if (strict)
      exit(1);
  }
  cur = sstack.top();
  sstack.pop();
}

Section * Getkw::readSect(istream & fis) {
  istringstream isi;
  string name;
  string tag;
  string dum1, dum2;
  string set;
  int nsect, nkeys, i;
  bool isSet, hasTag;

  readline(fis, isi);
  isi >> dum1 >> name >> nsect >> set;
  readline(fis, isi);
  isi >> dum1 >> tag >> dum2 >> nkeys;

  isSet = convBool(set);
  hasTag = convBool(tag);

  Section * thissect = new Section(name);
  thissect->setDefined(isSet);
  if (hasTag) {
    readline(fis, isi);
    isi >> tag;
    thissect->setTag(tag);
  } else {
    thissect->setTag(name);
  }

  for (i = 0; i < nkeys; i++) {
    readKey(thissect, fis);
  }
  for (i = 0; i < nsect; i++) {
    Section * newsect = readSect(fis);
    thissect->addSect(newsect);
  }
  return thissect;
}

bool Getkw::readKey(Section * sect, istream & fis) {
  istringstream isi;
  string type, name, set;
  int len;

  readline(fis, isi);
  isi >> type >> name >> len >> set;

  bool setf = convBool(set);
  int kind = convKind(type);

  string ss;
  vector<int> iv;
  vector<double> dv;
  vector<bool> bv;
  vector<string> sv;

  switch (kind) {
    case KeyType::Int:
      TEST_ARRAY;
      if (len == 0)
        return false;
      int is;
      readline(fis, isi);
      isi >> is;
      sect->addKey(new Keyword<int>(name, is, setf));
      break;
    case KeyType::Dbl:
      TEST_ARRAY;
      if (len == 0)
        return false;
      double ds;
      readline(fis, isi);
      isi >> ds;
      sect->addKey(new Keyword<double>(name, ds, setf));
      break;
    case KeyType::Bool:
      TEST_ARRAY;
      if (len == 0)
        return false;
      bool bs;
      readline(fis, isi);
      isi >> ss;
      bs = convBool(ss);
      sect->addKey(new Keyword<bool>(name, bs, setf));
      break;
    case KeyType::Str:
      TEST_ARRAY;
      if (len == 0)
        return false;
      getline(fis, ss);
      sect->addKey(new Keyword<string>(name, ss, setf));
      break;
    case KeyType::IntArray:
      if (len == 0)
        return false;
      int ival;
      for (int i = 0; i < len; i++) {
        readline(fis, isi);
        isi >> ival;
        iv.push_back(ival);
      }
      sect->addKey(new Keyword<vector<int> >(name, iv, setf));
      break;
    case KeyType::DblArray:
      if (len == 0)
        return false;
      double dval;
      for (int i = 0; i < len; i++) {
        readline(fis, isi);
        isi >> dval;
        dv.push_back(dval);
      }
      sect->addKey(new Keyword<vector<double> >(name, dv, setf));
      break;
    case KeyType::BoolArray:
      if (len == 0)
        return false;
      bool bval;
      for (int i = 0; i < len; i++) {
        readline(fis, isi);
        isi >> ss;
        bval = convBool(ss);
        bv.push_back(bval);
      }
      sect->addKey(new Keyword<vector<bool> >(name, bv, setf));
      break;
    case KeyType::StrArray:
    case KeyType::Data:
      if (len == 0)
        return false;
      for (int i = 0; i < len; i++) {
        getline(fis, ss);
        sv.push_back(ss);
      }
      sect->addKey(new Keyword<vector<string> >(name, sv, setf));
      break;
    default:
      THROW_GETKW("Unknown keyword type: " + name + " <> " + type);
  }
  return true;
}

void Getkw::readline(istream & fis, istringstream & isi) {
  static string buf;
  getline(fis, buf);
  isi.clear();
  isi.str(buf);
}

bool Getkw::convBool(const string & val) {
  if (val[0] == 'T' or val[0] == 't')
    return true;
  return false;
}

int Getkw::convKind(const string & typ) {
  static const string INT = "INT";
  static const string DBL = "DBL";
  static const string BOOL = "BOOL";
  static const string STR = "STR";
  static const string DATA = "DATA";
  static const string INT_ARRAY = "INT_ARRAY";
  static const string DBL_ARRAY = "DBL_ARRAY";
  static const string BOOL_ARRAY = "BOOL_ARRAY";
  static const string STR_ARRAY = "STR_ARRAY";

  if (typ.compare(INT) == 0)
    return KeyType::Int;
  if (typ.compare(DBL) == 0)
    return KeyType::Dbl;
  if (typ.compare(BOOL) == 0)
    return KeyType::Bool;
  if (typ.compare(STR) == 0)
    return KeyType::Str;
  if (typ.compare(DATA) == 0)
    return KeyType::Data;
  if (typ.compare(INT_ARRAY) == 0)
    return KeyType::IntArray;
  if (typ.compare(DBL_ARRAY) == 0)
    return KeyType::DblArray;
  if (typ.compare(BOOL_ARRAY) == 0)
    return KeyType::BoolArray;
  if (typ.compare(STR_ARRAY) == 0)
    return KeyType::StrArray;
  return -1;
}

template const int & Getkw::get<int>(const string &) const;
template const bool & Getkw::get<bool>(const string &) const;
template const double & Getkw::get<double>(const string &) const;
template const string & Getkw::get<string>(const string &) const;
template const vector<int> & Getkw::get<vector<int> >(const string &) const;
template const vector<double> & Getkw::get<vector<double> >(const string &) const;
template const vector<bool> & Getkw::get<vector<bool> >(const string &) const;
template const vector<string> & Getkw::get<vector<string> >(const string &) const;
