/*
 * @Description:
 * @Author: Zhuo Yue
 * @Date: 2021-06-02 15:57:40
 * @LastEditors: Zhuo Yue
 * @LastEditTime: 2021-07-22 01:12:43
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
    printf(" .----------------.  .----------------.  .----------------.  .----------------.  .-----------------.\n");
    printf("| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |\n");
    printf("| |     ______   | || |  ________    | || |   ______     | || |      __      | || | ____  _____  | |\n");
    printf("| |   .' ___  |  | || | |_   ___ `.  | || |  |_   __ \\   | || |     /  \\     | || ||_   \\|_   _| | |\n");
    printf("| |  / .'   \\_|  | || |   | |   `. \\ | || |    | |__) |  | || |    / /\\ \\    | || |  |   \\ | |   | |\n");
    printf("| |  | |         | || |   | |    | | | || |    |  ___/   | || |   / ____ \\   | || |  | |\\ \\| |   | |\n");
    printf("| |  \\ `.___.'\\  | || |  _| |___.' / | || |   _| |_      | || | _/ /    \\ \\_ | || | _| |_\\   |_  | |\n");
    printf("| |   `._____.'  | || | |________.'  | || |  |_____|     | || ||____|  |____|| || ||_____|\\____| | |\n");
    printf("| |              | || |              | || |              | || |              | || |              | |\n");
    printf("| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |\n");
    printf(" '----------------'  '----------------'  '----------------'  '----------------'  '----------------' \n");
    // clang-format on

    setenv("CDPAN_SCRIPT", GetCmdString(par, "script"), 1);
    // setenv("CDPAN_LOG", path_log, 1);
    if (GetCmdBool(par, "quite")) setenv("CDPAN_QUITE", "1", 1);
    if (GetCmdBool(par, "debug")) setenv("CDPAN_DEBUG", "1", 1);

    char *pl_path = getenv("CDPAN_PATH");
    char *cmd = (char *)malloc(sizeof(char) * PATH_MAX);
    strncpy(cmd, pl_path, PATH_MAX);
    if (GetCmdBool(par, "debug")) strcat(cmd, " > /dev/null 2> /dev/null");

    return 0;
}
