#include <elda/linalg.hpp>
using namespace linalg;
using namespace std;
int main() {
    matrix m (3,3);
    m.input();
    // m.print();
    // matrix q = m.qr_decomp_q();
    // matrix r = m.qr_decomp_r();
    // q.print();
    // cout << endl;
    // r.print();
    matrix eigen = m.eigenvalues();
    cout << endl;
    eigen.print();
}
