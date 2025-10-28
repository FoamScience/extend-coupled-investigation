#!/usr/bin/env bash

## Computes metrics for optimization
## -- return NaN if metric is not available yet
## -- should run in case folder

set -e
solver=$(awk -F'[; ]' '/^solver/ {print($2)}' domainDict)

case "$1" in
pError)
    uv run --script scripts/taylor_couette_analytical.py --compare "$solver" 2>/dev/null | awk '/RMSE \(p/ {found=1; print($4)} END {if (!found) print "nan"}' || print "nan"
    ;;
UError)
    uv run --script scripts/taylor_couette_analytical.py --compare "$solver"  2>/dev/null | awk '/RMSE \(u_theta/ {found=1; print($4)} END {if (!found) print "nan"}' || print "nan"
    ;;
executionTime)
    awk '/ExecutionTime/{found=1; tm=$3} END{if (!found) print "nan"; else print tm}' "logs/log.$solver" 2>/dev/null || echo "nan"
    ;;
continuityError)
    awk -F'[, ]' '/time step continuity errors/ {
        val = ($9 < 0 ? -$9 : $9)
        if (!found || val > max) max = val
        found = 1
    }
    END {
        if (!found) print "nan"
        else print max
    }' "logs/log.$solver"  2>/dev/null|| echo "nan"
    ;;
*)
    echo "Usage, from case dir.: mestric.sh [pError|UError|executionTime|continuityError]"
    exit 1
    ;;
esac
