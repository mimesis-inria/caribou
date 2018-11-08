#include <sofa/core/ObjectFactory.h>

#ifndef WIN32
#define SOFA_EXPORT_DYNAMIC_LIBRARY
#define SOFA_IMPORT_DYNAMIC_LIBRARY
#define SOFA_CARIBOU_API
#else
#ifdef SOFA_CARIBOU_API
#define SOFA_CARIBOU_API SOFA_EXPORT_DYNAMIC_LIBRARY
#else
#define SOFA_CARIBOU_API SOFA_IMPORT_DYNAMIC_LIBRARY
#endif
#endif


// Here are just several convenient functions to help user to know what contains the plugin

extern "C" {
SOFA_CARIBOU_API void        initExternalModule();
SOFA_CARIBOU_API const char* getModuleName();
SOFA_CARIBOU_API const char* getModuleVersion();
SOFA_CARIBOU_API const char* getModuleLicense();
SOFA_CARIBOU_API const char* getModuleDescription();
SOFA_CARIBOU_API const char* getModuleComponentList();
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