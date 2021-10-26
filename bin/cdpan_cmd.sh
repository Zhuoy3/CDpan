#!/usr/bin/env bash
###
# @Description:
# @Author: zhuoy
# @Date: 2021-10-22 17:44:24
# @LastEditors: zhuoy
# @LastEditTime: 2021-10-26 18:01:29
# @Calls:
# @Called By:
# @FilePath: \CDpan\bin\cdpan_cmd.sh
###
###
# @Description:
# @Author: zhuoy
# @Date: 2021-10-22 17:44:24
# @LastEditors: zhuoy
# @LastEditTime: 2021-10-22 18:13:25
# @Calls:
# @Called By:
# @FilePath: \CDpan\bin\cdpan_cmd.sh
###

export IFS=' '

_cdpan_cmd() {

    local cur=${COMP_WORDS[COMP_CWORD]}
    local cmd=${COMP_WORDS[COMP_CWORD - 1]}

    case $cmd in
    'cdpan')
        COMPREPLY=($(compgen -W 'filter align extract assembly mope vot soot merge location RUN-ALL RUN-DIEM' -- $cur))
        ;;
    'filter') ;&
    'align') ;&
    'extract') ;&
    'assembly') ;&
    'mope') ;&
    'vot') ;&
    'soot') ;&
    'merge') ;&
    'location') ;&
    'RUN-ALL') ;&
    'RUN-DIEM')
        COMPREPLY=($(compgen -W '-i --input -c --config -w --work_dir -l --output-level -o --output -O --output_dir --no-qc -v --version -h --help' -- $cur))
        ;;
    '-l') ;&
    '--output-level')
        COMPREPLY=($(compgen -W '0 1 2' -- $cur))
        ;;
    '--no-qc')
        COMPREPLY=($(compgen -W '-i --input -c --config -w --work_dir -l --output-level -o --output -O --output_dir --no-qc -v --version -h --help' -- $cur))
        ;;
    '-v') ;&
    '--version') ;&
    '-h') ;&
    '--help')
        COMPREPLY=()
        ;;
    '-i') ;&
    '--input') ;&
    '-c') ;&
    '--config') ;&
    '-w') ;&
    '--work_dir') ;&
    '-o') ;&
    '--output') ;&
    '-O') ;&
    '--output_dir') ;&
    '*') ;;
    esac

    return 0
}

complete -F _cdpan_cmd -o dirnames cdpan
