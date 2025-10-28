#!/usr/bin/env bash
set -e

if [ ".$WM_FORK" != ".extend" ]; then
    echo "Source Foam Extend 5.0"
    exit 1
fi

if ! command -v uv >/dev/null 2>&1; then
    echo "Install UV from https://docs.astral.sh/uv/"
fi

wmake ./solvers/MRFPUCoupledFoam

mkdir -p /tmp/trials artifacts
ln -fs "$PWD/scripts/metric.sh" /tmp/trials/metric.sh

echo "======="
echo "All set"
echo "======="
echo -e "To run the optimization: \e[1muvx foamBO --config MOO.yaml\e[0m"
echo -e "To check on the running/complete optimization: \e[1muvx foamBO --visualize --config MOO.yaml ++store.read_from=json\e[0m"
echo -e "More docs: \e[1muvx foamBO --docs\e[0m"
