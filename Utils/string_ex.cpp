#include "string_ex.h"

std::string str_ex_tolower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

size_t str_ex_tolower(char *s)
{
    std::string ss = s;
    std::string slow = str_ex_tolower(ss);
    size_t length = slow.length();
    strcpy(s, slow.c_str());

    return length;
}

bool str_ex_compare(char *s1, char *s2, bool case_sens)
{
    if (case_sens)
    {
        return !strcmp(s1, s2);
    }
    else
    {
        std::string ss1 = s1, ss2 = s2;
        std::string sslowl = str_ex_tolower(ss1);
        std::string sslow2 = str_ex_tolower(ss2);
        return !strcmp(sslowl.c_str(), sslow2.c_str());
    }
}
