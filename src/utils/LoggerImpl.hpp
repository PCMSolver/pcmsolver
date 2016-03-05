#ifndef LOGGERIMPL_HPP
#define LOGGERIMPL_HPP

#include <fstream>
#include <memory>
#include <string>

namespace logging
{
    enum printLevel {
        timings,
        coarse,
        fine,
        everything
    };

    /*! \class ILogPolicy
     *  \brief ABC for logging policy classes
     */
    class ILogPolicy
    {
    public:
        /*! \brief Opens an output stream with the given name
         *  \param[in] name name of the stream
         */
        virtual void open_ostream(const std::string & name) = 0;
        /*! \brief Closes an output stream with the given name
         */
        virtual void close_ostream() = 0;
        /*! \brief Writes to stream
         *  \param[in] msg message to be written to stream
         */
        virtual void write(const std::string & msg) = 0;
    };

    /*! \class FileLogPolicy
     *  \brief Implementation which allows to write into a file
     */
    class FileLogPolicy : public ILogPolicy
    {
    private:
        std::unique_ptr<std::ofstream> outStream_;
    public:
        /// Constructor
        FileLogPolicy() : outStream_(new std::ofstream) {}
        /// Destructor
        virtual ~FileLogPolicy() {
            if(outStream_) {
                close_ostream();
            }
        }
        /*! \brief Opens an output stream with the given name
         *  \param[in] name name of the stream
         */
        virtual void open_ostream(const std::string & name) {
            outStream_->open(name.c_str(), std::ios_base::binary | std::ios_base::out );
            if(!outStream_->is_open()) {
                PCMSOLVER_ERROR("LOGGER: Unable to open an output stream", BOOST_CURRENT_FUNCTION);
            }
        }
        /*! \brief Closes an output stream with the given name
         */
        virtual void close_ostream() {
            if(outStream_) {
                outStream_->close();
            }
        }
        /*! \brief Writes to stream
         *  \param[in] msg message to be written to stream
         */
        virtual void write(const std::string & msg) {
            (*outStream_) << msg << std::endl;
        }
    };

} // close namespace logging

#endif // LOGGERIMPL_HPP
