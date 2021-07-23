// Here are just several convenient functions to help user to know what contains the plugin
#include <SofaCaribou/config.h>

extern "C" {
    CARIBOU_API void        initExternalModule();
    CARIBOU_API const char* getModuleName();
    CARIBOU_API const char* getModuleVersion();
    CARIBOU_API const char* getModuleLicense();
    CARIBOU_API const char* getModuleDescription();
    CARIBOU_API const char* getModuleComponentList();
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