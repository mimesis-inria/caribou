// Here are just several convenient functions to help user to know what contains the plugin
#include <SofaCaribou/config.h>

extern "C" {
     void        initExternalModule();
     const char* getModuleName();
     const char* getModuleVersion();
     const char* getModuleLicense();
     const char* getModuleDescription();
     const char* getModuleComponentList();
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return "Caribou";
}

const char* getModuleVersion()
{
    return "alpha 1.0";
}

const char* getModuleLicense()
{
    return "TBD";
}

const char* getModuleDescription()
{
    return "Caribou's power in SOFA Framework";
}

const char* getModuleComponentList()
{
    /// string containing the names of the classes provided by the plugin
    return "";
}