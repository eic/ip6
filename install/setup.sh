#!/bin/sh

export BEAMLINE=ip6
export BEAMLINE_PATH=/Works/Sep30/ip6/install/share/ip6
export BEAMLINE_CONFIG=ip6
export BEAMLINE_VERSION=master

## note: we will phase out the JUGGLER_* flavor of variables in the future
export JUGGLER_BEAMLINE=$BEAMLINE
export JUGGLER_BEAMLINE_CONFIG=$BEAMLINE_CONFIG
export JUGGLER_BEAMLINE_VERSION=$BEAMLINE_VERSION
export JUGGLER_BEAMLINE_PATH=$BEAMLINE_PATH

## Export beamline libraries
if [[ "$(uname -s)" = "Darwin" ]] || [[ "$OSTYPE" == "darwin"* ]]; then
	export DYLD_LIBRARY_PATH=/Works/Sep30/ip6/install/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}
else
	export LD_LIBRARY_PATH=/Works/Sep30/ip6/install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
fi
