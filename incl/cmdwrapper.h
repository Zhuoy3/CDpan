/*
 * @Description:
 * @Author: Zhuo Yue
 * @Date: 2021-07-12 23:08:58
 * @LastEditors: Zhuo Yue
 * @LastEditTime: 2021-07-20 17:03:09
 * @Calls:
 * @Called By:
 * @FilePath: \CDpan\incl\cmdwrapper.h
 */

#ifndef _CMDWRAPPER_H_
#define _CMDWRAPPER_H_

struct CmdLineParser;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// Create a new CMDLINE class
extern struct CmdLineParser *NewCmd(int argc, char *argv[]);

// Get value of string argument
extern char *GetCmdString(struct CmdLineParser *cmd, char *par);

// Get value of int argument
extern int GetCmdInt(struct CmdLineParser *cmd, char *par);

// Get value of bool argument
extern int GetCmdBool(struct CmdLineParser *cmd, char *par);

#ifdef __cplusplus
};
#endif /* __cplusplus */

#endif /* _CMDWRAPPER_H_ */
