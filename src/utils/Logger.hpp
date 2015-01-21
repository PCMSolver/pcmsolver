#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <ctime>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

#include "LoggerImpl.hpp"

namespace logging
{
    template<typename logPolicy>
    class logger
    {
    private:
        size_t logLineNumber_;
        std::stringstream logStream_;
        logPolicy * policy_;
        std::mutex writeMutex_;

        std::string getTime() {
            std::string time_str;
            time_t raw_time;

            std::time(&raw_time);
            time_str = std::ctime(&raw_time);

            //without the newline character
            return time_str.substr(0, time_str.size() - 1);
        }
        std::string getLoglineHeader() {
            std::stringstream header;

            header.str("");
            header.fill('0');
            header.width(7);
            header << logLineNumber_++ << " < " << getTime() << " - ";

            header.fill('0');
            header.width(7);
            header << clock() << " > ~ ";

            return header.str();
        }

        /*! @name Core printing functionality
         *
         *  A variadic template is used, we specify the version
         *  with a variable number of parameters/types and the version
         *  with no parameters/types.
         *  The variadic template is called recursively.
         *  The type for the first parameter is resolved and streamed
         *  to logStream_. When all the parameters have been streamed
         *  the version with no arguments is called.
         */
        /// @{
        /*! */
        void printImpl() {
            policy_->write(getLoglineHeader() + logStream_.str());
            logStream_.str("");
        }
        template<typename First, typename...Rest>
        void printImpl(First parm1, Rest...parm) {
            logStream_ << parm1;
            printImpl(parm...);
        }
        /// @}
    public:
        /*! Constructor
        *  \param[in] name name for the log file
         */
        logger(const std::string & name) : logLineNumber_(0), policy_(new logPolicy) {
            if(!policy_) {
                throw std::runtime_error("LOGGER: Unable to create the logger instance");
            }
            policy_->open_ostream(name);
        }
        /// Destructor
        ~logger() {
            if(policy_) {
                policy_->close_ostream();
                delete policy_;
            }
        }

        /// User interface for the logger class
        template<severityType severity, typename...Args>
        void print(Args...args) {
            writeMutex_.lock();
            switch(severity) {
            case severityType::debug:
                logStream_ << "<DEBUG> : \n";
                break;
            case severityType::warning:
                logStream_ << "<WARNING> : \n";
                break;
            case severityType::error:
                logStream_ << "<ERROR> : \n";
                break;
            };
            printImpl(args...);
            writeMutex_.unlock();
        }

    };
} // close namespace logging

#endif // LOGGER_HPP
