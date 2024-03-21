#! /bin/bash

# run.sh
#
# Copyright (c) 2020 Carlos Braga
# This program is free software; you can redistribute it and/or modify it under
# the terms of the MIT License. See accompanying LICENSE.md or
# https://opensource.org/licenses/MIT.

# -----------------------------------------------------------------------------
# Die with message
die() {
    echo >&2 "$@"
    exit 1
}

# Run command and check exit code
run() {
    echo "$@" && "$@"
    code=$?
    [[ $code -ne 0 ]] && die "[$@] failed with error code $code"
    return 0
}

# Ask for input query
ask() {
    echo -n "$@ (y/n [n]): "
    local ans
    read ans
    [[ "$ans" != "y" ]] && return 1
    return 0
}

# -----------------------------------------------------------------------------
do_clean() {
    [[ -d build ]] && rm -rvf build
}

do_build() {
    mkdir build
    pushd build
    cmake .. && make -j16 all
    popd
}

do_execute() {
    dir="$1"
    target="$2"
    pushd "$dir"
    run "$target"
    popd
}

do_execute_mpi() {
    dir="$1"
    target="$2"
    pushd "$dir"
    run mpirun --mca shmem posix --use-hwthread-cpus -np 16 "$target"
    popd
}

[[ $# == 0 ]] && die "usage: $0 clean | build | execute"

DO_CLEAN=false
DO_BUILD=false
DO_EXECUTE=false
for arg in "$@"; do
    case $arg in
    clean)
        DO_CLEAN=true
    ;;
    build)
        DO_CLEAN=true
        DO_BUILD=true
    ;;
    execute)
        DO_CLEAN=true
        DO_BUILD=true
        DO_EXECUTE=true
    ;;
    *)
        die "unknown argument: $arg"
    ;;
    esac
done

[[ $DO_CLEAN == "true" ]] && do_clean
[[ $DO_BUILD == "true" ]] && do_build
if [[ $DO_EXECUTE == "true" ]]; then
    pushd build/algods
    do_execute arrays                               ./arrays
    do_execute hashing                              ./hashing
    do_execute linked                               ./linked
    do_execute sorting                              ./sorting
    do_execute searching                            ./searching
    do_execute unionfind                            ./unionfind
    do_execute numeric                              ./numeric
    popd
fi
