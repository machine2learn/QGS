#ifndef INCL_GZFILE_H
#define INCL_GZFILE_H

#define BOOST_IOSTREAMS_NO_LIB

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <fstream>
#include <string>
#include <memory>

class GZfile {
  bool d_gz;
  std::unique_ptr<std::ifstream> d_file;
  std::unique_ptr<boost::iostreams::filtering_streambuf<boost::iostreams::input>> d_gzf;
  std::unique_ptr<std::istream> d_handle;
  
  public:
  
  GZfile(std::string const & fname)
  :
    d_gz{fname.length() > 2 && fname.substr(fname.length() - 3) == ".gz"},
    d_file{std::unique_ptr<std::ifstream>(new std::ifstream(fname, d_gz ? std::ios::binary : std::ios::in))}, // make_unique
    d_gzf{std::unique_ptr<boost::iostreams::filtering_streambuf<boost::iostreams::input>>(new boost::iostreams::filtering_streambuf<boost::iostreams::input>())}, // make_unique
    d_handle{std::unique_ptr<std::istream>(new std::istream(d_gzf.get()))} // make_unique
  {
    if (d_gz)
      d_gzf->push(boost::iostreams::gzip_decompressor());
    d_gzf->push(*d_file);
  }

  GZfile(GZfile && tmp) :
    d_gz{tmp.d_gz},
    d_file(std::move(tmp.d_file)),
    d_gzf(std::move(tmp.d_gzf)),
    d_handle(std::move(tmp.d_handle))
  {
  }

  GZfile & operator=(GZfile && tmp) {
    std::swap(d_gz, tmp.d_gz);
    std::swap(d_file, tmp.d_file);
    std::swap(d_gzf, tmp.d_gzf);
    std::swap(d_handle, tmp.d_handle);
    return *this;
  }
  
  std::istream & handle() {
    return *d_handle;
  };
  
  operator bool() const {
    return *d_file && *d_handle;
  }
};

template <typename T>
GZfile & operator>>(GZfile & gzfile, T & t) {
  operator>>(gzfile.handle(), t);
  return gzfile;
}

#endif
