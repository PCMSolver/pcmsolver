#ifndef LOGGERIMPL_HPP
#define LOGGERIMPL_HPP

#include <fstream>
#include <memory>
#include <string>

namespace logging
{
    enum severityType {
        debug = 1,
        error,
        warning
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
	 *  \param[in] name name of the stream
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
        virtual ~FileLogPolicy();
	/*! \brief Opens an output stream with the given name
	 *  \param[in] name name of the stream
	 */
        virtual void open_ostream(const std::string & name);
	/*! \brief Closes an output stream with the given name
	 *  \param[in] name name of the stream
	 */
        virtual void close_ostream();
	/*! \brief Writes to stream
	 *  \param[in] msg message to be written to stream
	 */
        virtual void write(const std::string & msg);
    };

} // close namespace logging

#endif // LOGGERIMPL_HPP
