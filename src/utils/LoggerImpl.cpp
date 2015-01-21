#include "LoggerImpl.hpp"

#include <fstream>
#include <stdexcept>

namespace logging
{
    void FileLogPolicy::open_ostream(const std::string & name)
    {
        outStream_->open(name.c_str(), std::ios_base::binary | std::ios_base::out);
        if(!outStream_->is_open()) {
            throw(std::runtime_error("LOGGER: Unable to open an output stream"));
        }
    }

    void FileLogPolicy::close_ostream()
    {
        if(outStream_) {
            outStream_->close();
        }
    }

    void FileLogPolicy::write(const std::string & msg)
    {
        (*outStream_) << msg << std::endl;
    }

    FileLogPolicy::~FileLogPolicy()
    {
        if(outStream_) {
            close_ostream();
        }
    }

} // close namespace logging
