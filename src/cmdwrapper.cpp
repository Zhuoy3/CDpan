/*
 * @Description:
 * @Author: Zhuo Yue
 * @Date: 2021-07-12 23:09:49
 * @LastEditors: Zhuo Yue
 * @LastEditTime: 2021-07-20 17:24:42
 * @Calls:
 * @Called By:
 * @FilePath: \CDpan\src\cmdwrapper.cpp
 */

#include "cpp/cmdline.h"
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

struct CmdLineParser {
    cmdline::parser cmd;
};

extern struct CmdLineParser *NewCmd(int argc, char *argv[]) {
    // Create a new CMDLINE class

    CmdLineParser *a = new struct CmdLineParser;

    (a->cmd).add<string>("script", 's', "script file path", true, "");
    (a->cmd).add("quite", 'q', "do not output to stdout");
    (a->cmd).add("debug", 'd', "output debug information");
    (a->cmd).add("version", 'v', "print version message");
    (a->cmd).add("help", 'h', "print help message");
    (a->cmd).set_program_name("CDpan");

    bool ok=(a->cmd).parse(argc, argv);

    if (!ok){
        cerr<<(a->cmd).error()<<endl<<(a->cmd).usage();
        exit(EXIT_FAILURE);
    }

    if (argc==1 || (a->cmd).exist("help")){
        cout<<(a->cmd).usage();
        exit(EXIT_SUCCESS);
    }

    if ((a->cmd).exist("version")){
        cout<<"CDpan v0.1.12 Jue 20 2021"<<endl;
        exit(EXIT_SUCCESS);
    }

    return a;
}

extern char *GetCmdString(struct CmdLineParser *cmd, char *par){
    return (char*)((cmd->cmd).get<string>(par)).data();
}

extern int GetCmdInt(struct CmdLineParser *cmd, char *par){
    return (cmd->cmd).get<int>(par);
}

extern int GetCmdBool(struct CmdLineParser *cmd, char *par){
    return ((cmd->cmd).exist(par))?1:0;
}

#ifdef __cplusplus
};
#endif /* __cplusplus */
