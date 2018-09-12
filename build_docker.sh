#!/usr/bin/env bash

#change the version number for each new build
docker build -t cmap/psp:latest -t cmap/psp:v0.0.1 --rm=true .