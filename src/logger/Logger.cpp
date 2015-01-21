#include "Logger.hpp"

#include <ctime>
#include <stdexcept>

#include "LoggerImpl.hpp"

namespace logging 
{
   template <typename logPolicy>
   logger<logPolicy>::logger(const std::string & name) 
       : logLineNumber_(0), policy_(new logPolicy) 
   {
       if(!policy_) {
           throw std::runtime_error("LOGGER: Unable to create the logger instance");
       }
       policy_->open_ostream(name);
   }

   template <typename logPolicy>
   logger<logPolicy>::~logger()
   {
            if(policy_) {
                policy_->close_ostream();
                delete policy_;
            }
   }

   template <typename logPolicy>
   std::string logger<logPolicy>::getTime()
   {
          std::string time_str;                               
	  time_t raw_time;                                    
                                                             
	  std::time(&raw_time);                                 
          time_str = std::ctime(&raw_time);                      
                                                             
          //without the newline character                     
          return time_str.substr(0, time_str.size() - 1);  
      }                                                       
   
   template <typename logPolicy>
   std::string logger<logPolicy>::getLoglineHeader() 
   {
            std::stringstream header;

            header.str("");
            header.fill('0');
            header.width(7);
            header << logLineNumber_++ << " < " <<getTime() << " - ";

            header.fill('0');
            header.width(7);
            header << clock() << " > ~ ";

            return header.str();
        }

   template <typename logPolicy>
   void logger<logPolicy>::printImpl()
   {
            policy_->write(getLoglineHeader() + logStream_.str());
            logStream_.str("");
   }
   
   /*! @name Explicit specializations
    */ 
   /// @{
   /*! Specialization for FileLogPolicy */
   template class logger<FileLogPolicy>;
   /// @}
} // close namespace logging
