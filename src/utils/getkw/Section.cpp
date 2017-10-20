/** \file Section.cpp
 *
 * \date Jun 3, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 * \brief  Container for Sections and Keywords
 */

#include <iostream>

#include "GetkwError.h"
#include "Section.h"

#define IF_ANY_KEYTYPE_IS(A, T) if (A.type() == typeid(const Keyword<T> *))

#define ANY_TO_CONST_KEY_PTR(A, T) boost::any_cast<const Keyword<T> *>(A)

#define PRINT_FUNC_NAME std::cout << "@ Section::" << __func__ << std::endl;

using namespace std;

Section::Section(const string & _name, const string & _tag) : name(_name) {
  isDefd = false;
  nsect = 0;
  nkeys = 0;
  if (tag.empty()) {
    this->tag = _name;
  } else {
    this->tag = _tag;
  }
}

Section::Section(const Section & s) {
  tag = s.tag;
  nkeys = s.nkeys;
  nsect = s.nsect;

  copySects(s);
  copyKeys(s);
}

Section & Section::operator=(const Section & s) {
  if (&s == this) {
    return *this;
  }
  name = s.name;
  tag = s.tag;
  nkeys = s.nkeys;
  nsect = s.nsect;

  copySects(s);
  copyKeys(s);

  return *this;
}

Section::~Section() {
  map<string, Section *>::iterator sit;
  map<string, boost::any>::iterator kit;
  for (sit = sects.begin(); sit != sects.end(); ++sit) {
    delete sects[sit->first];
  }
  for (kit = keys.begin(); kit != keys.end(); ++kit) {

    IF_ANY_KEYTYPE_IS(kit->second, int) {
      const Keyword<int> * key = ANY_TO_CONST_KEY_PTR(kit->second, int);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, bool) {
      const Keyword<bool> * key = ANY_TO_CONST_KEY_PTR(kit->second, bool);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, double) {
      const Keyword<double> * key = ANY_TO_CONST_KEY_PTR(kit->second, double);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, string) {
      const Keyword<string> * key = ANY_TO_CONST_KEY_PTR(kit->second, string);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, vector<int>) {
      const Keyword<vector<int> > * key =
          ANY_TO_CONST_KEY_PTR(kit->second, vector<int>);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, vector<bool>) {
      const Keyword<vector<bool> > * key =
          ANY_TO_CONST_KEY_PTR(kit->second, vector<bool>);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, vector<double>) {
      const Keyword<vector<double> > * key =
          ANY_TO_CONST_KEY_PTR(kit->second, vector<double>);
      delete key;
      continue;
    }
    IF_ANY_KEYTYPE_IS(kit->second, vector<string>) {
      const Keyword<vector<string> > * key =
          ANY_TO_CONST_KEY_PTR(kit->second, vector<string>);
      delete key;
      continue;
    }
    THROW_GETKW("Error! Unknown key type!");
  }
}

void Section::copySects(const Section & s) {
  map<string, Section *>::const_iterator iter;
  for (iter = s.sects.begin(); iter != s.sects.end(); ++iter) {
    sects[iter->first] = new Section(*iter->second);
    tags[(*iter->second).tag] = sects[iter->first];
  }
}

void Section::copyKeys(const Section & s) {
  map<string, boost::any>::const_iterator iter;

  for (iter = s.keys.begin(); iter != s.keys.end(); ++iter) {
    IF_ANY_KEYTYPE_IS(iter->second, int) {
      const Keyword<int> * key = ANY_TO_CONST_KEY_PTR(iter->second, int);
      keys[iter->first] = boost::any(new const Keyword<int>(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, bool) {
      const Keyword<bool> * key = ANY_TO_CONST_KEY_PTR(iter->second, bool);
      keys[iter->first] = boost::any(new const Keyword<bool>(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, double) {
      const Keyword<double> * key = ANY_TO_CONST_KEY_PTR(iter->second, double);
      keys[iter->first] = boost::any(new const Keyword<double>(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, string) {
      const Keyword<string> * key = ANY_TO_CONST_KEY_PTR(iter->second, string);
      keys[iter->first] = boost::any(new const Keyword<string>(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<int>) {
      const Keyword<vector<int> > * key =
          ANY_TO_CONST_KEY_PTR(iter->second, vector<int>);
      keys[iter->first] = boost::any(new const Keyword<vector<int> >(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<bool>) {
      const Keyword<vector<bool> > * key =
          ANY_TO_CONST_KEY_PTR(iter->second, vector<bool>);
      keys[iter->first] = boost::any(new const Keyword<vector<bool> >(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<double>) {
      const Keyword<vector<double> > * key =
          ANY_TO_CONST_KEY_PTR(iter->second, vector<double>);
      keys[iter->first] = boost::any(new const Keyword<vector<double> >(*key));
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<string>) {
      const Keyword<vector<string> > * key =
          ANY_TO_CONST_KEY_PTR(iter->second, vector<string>);
      keys[iter->first] = boost::any(new const Keyword<vector<string> >(*key));
      continue;
    }
    THROW_GETKW("Error! Unknown key type!");
  }
}

const Section & Section::getSect(const string & pathspec) const {
  vector<string> path;
  splitPath(pathspec, path);
  const Section * sect = traversePath(path, pathspec);
  return *sect;
}

template <typename T>
const Keyword<T> & Section::getKey(const string & pathspec) const {
  vector<string> path;
  splitPath(pathspec, path);

  string name = path.back();
  const Section * sect = traversePath(path, pathspec);
  if (!sect->has_key(name)) {
    THROW_GETKW("Invalid keyword, " + pathspec);
  }
  map<string, boost::any>::const_iterator iter = sect->keys.find(name);
  return *ANY_TO_CONST_KEY_PTR(iter->second, T);
}

template <class T> const T & Section::get(const string & path) const {
  const Keyword<T> & key = this->getKey<T>(path);
  return key.get();
}

Section * Section::readSect(ifstream & /* fis */) { return 0; }

/* Sections can be multiply defined, provided they have different tags.
 * The first section section is added to the box, as well as to the tags
 * map if applicable. The following sections with the same name are added
 * to the tags map.
 */
void Section::addSect(Section * sect) {
  string name = sect->name + "<" + sect->tag + ">";
  if (has_key(name)) {
    THROW_GETKW("Section::add: Section already defined, " + name);
  }

  sects[name] = sect;
  tags[sect->tag] = sects[name];
  nsect++;
}

void Section::addSect(Section & sect) {
  string name = sect.name + "<" + sect.tag + ">";
  if (has_key(name)) {
    THROW_GETKW("Section::add: Section already defined, " + name);
  }

  sects[name] = new Section(sect);
  tags[sect.tag] = sects[name];
  nsect++;
}
//! Add an allocated key to a section
template <class T> void Section::addKey(const Keyword<T> * key) {
  const string & name = key->getName();

  if (has_key(name)) {
    THROW_GETKW("Section::add: Key already defined, " + name);
  }

  keys[name] = boost::any(key);
  nkeys++;
}

//! Copy and and add a keyword to a section
template <class T> void Section::addKey(const Keyword<T> & key) {
  string name = key.getName();

  if (has_key(name)) {
    THROW_GETKW("Section::add: Key already defined, " + name);
  }

  keys[name] = boost::any(new Keyword<T>(key));
  nkeys++;
}

const Section * Section::traversePath(vector<string> & path,
                                      const string & pathspec) const {
  string cur = path[0];

  if (path.size() == 1) {
    if (has_key(cur))
      return this;
    if (!has_sect(cur)) {
      cur = cur + "<" + cur + ">";
    }
    if (has_sect(cur)) {
      map<string, Section *>::const_iterator iter = sects.find(cur);
      return iter->second;
    }
    THROW_GETKW("traversePath: Invalid path, " + pathspec);
  }

  if (!has_sect(cur))
    cur = cur + "<" + cur + ">";
  if (!has_sect(cur)) {
    THROW_GETKW("traversePath: Invalid path, " + pathspec);
  }

  path.erase(path.begin());
  map<string, Section *>::const_iterator iter = sects.find(cur);
  return iter->second->traversePath(path, pathspec);
}

void Section::splitPath(const string & pathspec, vector<string> & path) const {
  string str = pathspec;
  string::size_type m = 0;
  while (true) {
    m = str.find('.', m);
    if (m == string::npos) {
      path.push_back(str);
      break;
    } else {
      path.push_back(str.substr(0, m));
    }
    m++;
    str = str.substr(m);
  }
}

int Section::splitTag(const string & path, string & tag) const {
  string::size_type m, n = 0;
  m = path.find('<');
  if (m == string::npos)
    return 0;
  n = path.find('>');
  if (n == string::npos or n - m - 1 < 1)
    return 0;
  tag.clear();
  tag.append(path.substr(m + 1, n - m - 1));
  return m;
}

bool Section::has_key(const string & b) const {
  if (keys.find(b) == keys.end())
    return false;
  return true;
}

bool Section::has_sect(const string & b) const {
  if (sects.find(b) == sects.end())
    return false;
  return true;
}

bool Section::has_tag(const string & b) const {
  if (tags.find(b) == tags.end())
    return false;
  return true;
}

void Section::print() const { cout << &repr(cout) << endl; }

ostream & Section::repr(ostream & o) const {
  o << endl << name;
  if (name.compare(tag) != 0) {
    o << "<" + tag + ">";
  }
  o << " {" << endl;

  map<string, boost::any>::const_iterator iter;
  for (iter = keys.begin(); iter != keys.end(); ++iter) {
    IF_ANY_KEYTYPE_IS(iter->second, int) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, int) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, bool) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, bool) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, double) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, double) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, string) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, string) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<int>) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, vector<int>) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<bool>) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, vector<bool>) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<double>) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, vector<double>) << endl;
      continue;
    }
    IF_ANY_KEYTYPE_IS(iter->second, vector<string>) {
      o << *ANY_TO_CONST_KEY_PTR(iter->second, vector<string>) << endl;
      continue;
    }
    THROW_GETKW("Unknown key type!");
  }

  map<string, Section *>::const_iterator s_it;
  for (s_it = sects.begin(); s_it != sects.end(); ++s_it) {
    o << *s_it->second << endl;
  }

  o << "}";
  return o;
}

template void Section::addKey(const Keyword<int> *);
template void Section::addKey(const Keyword<double> *);
template void Section::addKey(const Keyword<bool> *);
template void Section::addKey(const Keyword<string> *);
template void Section::addKey(const Keyword<vector<int> > *);
template void Section::addKey(const Keyword<vector<double> > *);
template void Section::addKey(const Keyword<vector<bool> > *);
template void Section::addKey(const Keyword<vector<string> > *);

template void Section::addKey(const Keyword<int> &);
template void Section::addKey(const Keyword<double> &);
template void Section::addKey(const Keyword<bool> &);
template void Section::addKey(const Keyword<string> &);
template void Section::addKey(const Keyword<vector<int> > &);
template void Section::addKey(const Keyword<vector<double> > &);
template void Section::addKey(const Keyword<vector<bool> > &);
template void Section::addKey(const Keyword<vector<string> > &);

template const Keyword<int> & Section::getKey(const string &) const;
template const Keyword<double> & Section::getKey(const string &) const;
template const Keyword<bool> & Section::getKey(const string &) const;
template const Keyword<string> & Section::getKey(const string &) const;
template const Keyword<vector<int> > & Section::getKey(const string &) const;
template const Keyword<vector<double> > & Section::getKey(const string &) const;
template const Keyword<vector<bool> > & Section::getKey(const string &) const;
template const Keyword<vector<string> > & Section::getKey(const string &) const;

template const int & Section::get<int>(const string &) const;
template const bool & Section::get<bool>(const string &) const;
template const double & Section::get<double>(const string &) const;
template const string & Section::get<string>(const string &) const;
template const vector<int> & Section::get<vector<int> >(const string &) const;
template const vector<double> & Section::get<vector<double> >(const string &) const;
template const vector<bool> & Section::get<vector<bool> >(const string &) const;
template const vector<string> & Section::get<vector<string> >(const string &) const;
