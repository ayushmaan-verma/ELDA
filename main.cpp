#include <elda/matrix.hpp>
using namespace linalg;
using namespace std;
int main() {
    matrix m (3,3);
    m.arr = {{3,4,5},{6,7,8},{8,2,3}};
    m.gaussian();
    m.print();
}
