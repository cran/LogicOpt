/* LINTLIBRARY */
#include "port.h"
#include "utility.h"
#include "time.h"
long util_cpu_time(void)
{
    return clock();
}
