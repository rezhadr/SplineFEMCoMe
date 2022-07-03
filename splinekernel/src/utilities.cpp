#include "utilities.hpp"

void runtime_check(bool result, const char message[])
{
    if( !result )
    {
        throw std::runtime_error{ message };
    }
}