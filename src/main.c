/*
 * @Description:
 * @Author: Zhuo Yue
 * @Date: 2021-06-02 15:57:40
 * @LastEditors: Zhuo Yue
 * @LastEditTime: 2021-07-20 18:10:13
 * @Calls:
 * @Called By:
 * @FilePath: \CDpan\src\main.c
 */

#define print(LOG, format, args...)                                                                                    \
    do {                                                                                                               \
        fprintf(stdout, format, ##args);                                                                               \
        fprintf(LOG, format, ##args);                                                                                  \
    } while (0)

#include "cmdwrapper.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    struct CmdLineParser *par = NewCmd(argc, argv);

    if (GetCmdBool(par, "quite")) {
        freopen("/dev/null", "w", stdout);
        freopen("/dev/null", "w", stderr);
    }

    // char *path_log = (char *)malloc(sizeof(char) * PATH_MAX);
    // sprintf(path_log, "%s.log", GetCmdString(par, "script"));
    // FILE *log = fopen(path_log, "w+");

    // clang-format off
    // print(log, " .----------------.  .----------------.  .----------------.  .----------------.  .-----------------.\n");
    // print(log, "| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |\n");
    // print(log, "| |     ______   | || |  ________    | || |   ______     | || |      __      | || | ____  _____  | |\n");
    // print(log, "| |   .' ___  |  | || | |_   ___ `.  | || |  |_   __ \\   | || |     /  \\     | || ||_   \\|_   _| | |\n");
    // print(log, "| |  / .'   \\_|  | || |   | |   `. \\ | || |    | |__) |  | || |    / /\\ \\    | || |  |   \\ | |   | |\n");
    // print(log, "| |  | |         | || |   | |    | | | || |    |  ___/   | || |   / ____ \\   | || |  | |\\ \\| |   | |\n");
    // print(log, "| |  \\ `.___.'\\  | || |  _| |___.' / | || |   _| |_      | || | _/ /    \\ \\_ | || | _| |_\\   |_  | |\n");
    // print(log, "| |   `._____.'  | || | |________.'  | || |  |_____|     | || ||____|  |____|| || ||_____|\\____| | |\n");
    // print(log, "| |              | || |              | || |              | || |              | || |              | |\n");
    // print(log, "| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |\n");
    // print(log, " '----------------'  '----------------'  '----------------'  '----------------'  '----------------' \n");
    // clang-format on

    setenv("CDPAN_SCRIPT", GetCmdString(par, "script"), 1);
    // setenv("CDPAN_LOG", path_log, 1);
    if (GetCmdBool(par, "quite")) setenv("CDPAN_QUITE", "1", 1);
    if (GetCmdBool(par, "debug")) setenv("CDPAN_DEBUG", "1", 1);

    char *pl_path = getenv("CDPAN_PATH");
    system(pl_path);

    return 0;
}
