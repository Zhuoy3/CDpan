#!/bin/sh
###
 # @Description:
 # @Author: Zhuo Yue
 # @Date: 2021-07-03 23:57:18
 # @LastEditors: Zhuo Yue
 # @LastEditTime: 2021-07-22 01:00:59
 # @Calls:
 # @Called By:
 # @FilePath: \CDpan\bin\CDpan.sh
###

if [ -h $0 ]
then
    bin=`ls -ld $0|awk '{print $NF}'`
else
    bin=`cd $(dirname $0); pwd`'/'`basename  $0`
fi

export CDPAN_PATH=${bin}.pl
${bin}.bin "$@"
