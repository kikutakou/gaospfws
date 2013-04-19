#!/bin/bash

BINARY="gaospfws_$$"
TIMMSR="ostime_$$"

trap 'rm -f ./$BINARY ./$TIMMSR' INT


while [ -n "$1" ] ; do
if [ "${1:0:2}" == "-D" ] ; then MACRO=( ${MACRO[@]} $1 ) ; shift
elif echo "$1" | grep -E "^--thread=[0-9]+$" >/dev/null ; then THREAD="${1:9:10}" ; shift
elif echo "$1" | grep -E "^--tmlim=[0-9]+$" >/dev/null ; then MACRO=( "${MACRO[@]}" -DTMLIM="${1:8:10}" ) ; shift
elif [ "$1" == "-d" ] ; then DEBUG=set ; MACRO=( "${MACRO[@]}" -DDEBUG=0 ) ; shift
elif echo "$1" | grep -E "^[0-9]+$" >/dev/null ; then ARGS=( "${ARGS[@]}" $1 ) ; shift
else echo "ERROR : bad argment or option \"$1\"" ; exit 1 ; fi
done

[ "${#ARGS[@]}" == 3 ] || { echo "ERROR  : argments too many or less" ; exit 1 ;}
echo "${ARGS[@]}" | grep -E "^[0-9 ]*$" >/dev/null || { echo "ERROR  : non numerical argment detected" ; exit 1 ;}
[ -n "$THREAD" ] && MACRO=( "${MACRO[@]}" "-DMP=$THREAD")

#filename
PBTITLE=$(echo "${ARGS[@]}" | sed "s/ /_/g")
HEADDERFILE="$(dirname $0)/header/${PBTITLE}.h" ; mkdir -p $(dirname $0)/header
TITLE=${PBTITLE}$(echo "${MACRO[@]}" | sed "s/-D/_/g" | sed "s/ //g")
OUTPUT="$(dirname $0)/out/ecmp_${TITLE}.txt" ; mkdir -p $(dirname $0)/out
while [ -e "${OUTPUT%.txt}${PP}.txt" ] ; do PP=_$[${PP:1:5}+1] ; done ; OUTPUT="${OUTPUT%.txt}${PP}.txt"
[ -n "$DEBUG" ] && OUTPUT="/dev/stdout"
[ -n "$THREAD" ] && MACRO=( "${MACRO[@]}" "-fopenmp")

#make headder
ruby $(dirname $0)/mkheadder.rb ${ARGS[@]} "$HEADDERFILE" >/dev/null
[ "$?" == 1 ] && { echo "$0:$LINENO ($(date '+%m/%d %H:%M')): ruby error" 1>&2 ; exit 1 ;}

#compile
g++ -g -lm ${MACRO[@]} -include "$HEADDERFILE" $(dirname $0)/gaospfws.cc -o $(dirname $0)/$BINARY
[ "$?" == 1 ] && { echo "$0:$LINENO ($(date '+%m/%d %H:%M')): gcc error" 1>&2 ; exit 1 ;}

(time ./$BINARY) > "$OUTPUT" 2>./$TIMMSR && { cat ./$TIMMSR >> "$OUTPUT" ; rm -f ./$BINARY ./$TIMMSR ;}

[ -z "$DEBUG" -a -s "$OUTPUT" ] && bash $(dirname $0)/../common/resultshandling/gaospfws_onecolumn.sh $OUTPUT



