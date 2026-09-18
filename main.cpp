#include <random>
#include <string>
#include "src/innovation.h"
#include "src/easylogging++.h"

INITIALIZE_EASYLOGGINGPP

std::random_device rd;
std::mt19937 randomizer_engine(rd()); // seeding random number engine


int cafe5(int argc, char *const argv[]);

int main(int argc, char *const argv[]) {

    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Global, el::ConfigurationType::Format, "%msg");
    defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "false");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, "WARNING: %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);

    for (int i=1; i<argc; ++i)
        if (std::string(argv[i]) == "--innovation")
            return innovation::run(argc, argv);
    return cafe5(argc, argv);
}
