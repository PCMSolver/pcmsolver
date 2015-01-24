#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <ctime>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

#include "BuildInfo.hpp"
#include "LoggerImpl.hpp"

namespace logging
{
    /*! \brief Returns date and time
     */                                      
    std::string getTime() {
        std::string time_str;
        time_t raw_time;
                                          
        std::time(&raw_time);
        time_str = std::ctime(&raw_time);
                                          
        // Without the newline character
        return time_str;
    }

    template<typename logPolicy>
    class logger
    {
    private:
	printLevel globalPrintLevel_;
        std::stringstream logStream_;
        logPolicy * policy_;
        std::mutex writeMutex_;
    
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
            policy_->write(logStream_.str());
            logStream_.str("");
        }
        template<typename First, typename...Rest>
        void printImpl(First parm1, Rest...parm) {
	    logStream_.precision(std::numeric_limits<double>::digits10);
            logStream_ << parm1 << std::endl;
            printImpl(parm...);
        }
        /// @}
    public:
        /*! Constructor
        *  \param[in] name name for the log file
	*  
	*  The build parameters are logged first
         */
        logger(const std::string & name, printLevel print = coarse) 
		: globalPrintLevel_(print), policy_(new logPolicy) 
	{
            if(!policy_) {
                throw std::runtime_error("LOGGER: Unable to create the logger instance");
            }
            policy_->open_ostream(name);
	    // Write the logfile header
	    logStream_ << "\t\tPCMSolver execution log\n" 
		       << buildInfo() << "\n\t\tLog started : " << getTime() << std::endl;
        }
        /// Destructor
        ~logger() {
            if(policy_) {
                policy_->close_ostream();
                delete policy_;
            }
        }

	void globalPrintLevel(int print) { globalPrintLevel_ = print; }

        /// User interface for the logger class
        template<printLevel print, typename...Args>
        void print(Args...args) {
	    if (globalPrintLevel_ >= print) {
               writeMutex_.lock();   
               printImpl(args...);
               writeMutex_.unlock();
	    }
        }

    };
} // close namespace logging

#endif // LOGGER_HPP
