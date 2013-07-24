void testsuite1();

int main( int argc, char **argv )
{
    testsuite1(); // register tests

    if        ( argc==1 ) {
        run_all();
    } else if ( argc==2 ) {
        run_one( atoi(argv[1]) );
    } else {
        printf( "too many arguments\n" );
        exit(-1);
    }
}
