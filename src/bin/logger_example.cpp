#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

int main() {
   LOG(INFO) << "Write something to the logger";
   return 0;
}
