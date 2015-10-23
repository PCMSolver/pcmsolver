void do_not_throw(double aug) noexcept
{
    double d = 2.213132145525;
    d += aug;
}


int main()
{
    do_not_throw(2.0);
	return 0;
}
