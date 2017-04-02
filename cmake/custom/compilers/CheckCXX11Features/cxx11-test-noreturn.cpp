[[noreturn]] void do_not_return() {
  throw "error";
  // OK
}

int main() { return 0; }
