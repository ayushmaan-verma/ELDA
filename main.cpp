#include <elda/linalg.hpp>
using namespace linalg;
using namespace std;
int main() {
    matrix m = vec4(1,0,0,1);
    matrix rot = rot_y(PI/2);
    (rot*m).print();
}
