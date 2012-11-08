#!/bin/sh
if test $# -lt 1
then
  echo "Usage: $0 <cfilenamelist>"
  exit 1
fi
cat << ENDOFINCLUDE
/* 
  This file is generated. Do not edit. 
*/

#ifndef PROTODEF_H
#define PROTODEF_H

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "optdesc.h"
#include "multidef.h"
#include "mumcand.h"

ENDOFINCLUDE

cproto "$@"

cat <<EOF 

#endif /* PROTODEF_H */
EOF
