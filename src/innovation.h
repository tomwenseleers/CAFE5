#ifndef CAFE_INNOVATION_H
#define CAFE_INNOVATION_H
#include <vector>
namespace innovation {
// Rows/columns include zero. Entries are exact unbounded-process probabilities;
// finite rows are deliberately NOT renormalized.
std::vector<double> transition(int maximum, double lambda, double nu, double time);
int run(int argc, char *const argv[]);
}
#endif
