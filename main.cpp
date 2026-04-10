#include <elda/linalg.hpp>
using namespace linalg;
using namespace std;
int main() {
    matrix m(3,3);
    m.arr ={{2,3,1},{4,7,7},{-2,4,5}};
    matrix l = m.lu_decomp_l();
    matrix u = m.lu_decomp_u();
    m.print();
    cout << endl;
    l.print();
    cout << endl;
    u.print();
    return 0;
}
